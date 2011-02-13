library(quadprog)

bedcols=c("character","integer","integer","NULL","numeric","character")
gffcols=c("character","NULL","NULL","integer","integer","numeric","character","NULL","NULL")

cross.correlate = function(counts,threshold=1,lag.max=600,val.cut=.4) {
    signal=matrix(0,ncol=2,nrow=0)
    for (chr in names(counts[[1]])) {
        I = which(counts[[2]][[chr]]<val.cut)
        for (c in counts[[1]][[chr]][I])
          signal = rbind(signal,cbind(c$minus,c$plus))
    }
    if (nrow(signal) < 2) return(list(lag=c(0),acf=c(0)))
    signal[signal<=threshold] = 0
    ccf = ccf(x=ts(signal[,1]),y=ts(signal[,2]),lag.max=lag.max,plot=FALSE)
    I = which(ccf$lag>0)
    ccf$acf[ccf$acf<0] = 0
    list(lag=ccf$lag[I],acf=ccf$acf[I])
}

kernel.types = c("geometric","gaussian")

camel.kernel = function(size,mu,lambda,len,ktype) {
    if (ktype=="geometric") {
        q = 1/lambda
### need mu>len, lambda>mu
        k = c(
          1:len*q,
          rep(len*q,mu-len),
          q*((len-1):1)+1-(1-q)^(1:(len-1)),
          1-(1-q)^len,
          (1-q)^(1:(size-len-mu))*(1-(1-q)^len))
    } else if (ktype=="gaussian") {
        k = dnorm(x=0:(size+len-2),mean=mu+len,sd=lambda)
    }
    k/sum(k)
}

fit.solution = function(prob,mu,lambda,len,ktype) {
    size = length(prob)
    kernel = camel.kernel(size=size,mu,lambda,len,ktype)
#    a1 = amp[1]
#    a2 = a1
#    if (length(amp)>1) a2 = amp[2]
    data.frame(plus=convolve(kernel,prob,type="open")[len+size:1],
               minus=convolve(kernel,rev(prob),type="open")[len+1:size])
}

make.matrix = function(data.minus,data.plus,mu,lambda,len,reg,ktype) {
    size = length(data.minus)
    kernel = camel.kernel(size=size,mu=mu,lambda=lambda,len=len,ktype)
    Kmat = matrix(0,nrow=size,ncol=size-2*len)
    for (n in 1:(size-2*len)) Kmat[n:size,n] = kernel[1:(size-n+1)]
    sdKmat=try(svd(Kmat))
    if (class(sdKmat) == "try-error") Dmat=diag(1,ncol=size-2*len,nrow=size-2*len)
    else {
        I=which(sdKmat$d^2>reg)
        Kmat=sdKmat$u[,I] %*% diag(sdKmat$d[I]) %*% t(sdKmat$v[,I])
        Kmat=sweep(Kmat,2,apply(Kmat,2,sum),"/")
        sdKmat$d[-I]=1
        K2=sdKmat$u %*% diag(sdKmat$d) %*% t(sdKmat$v)
        K2=sweep(K2,2,apply(K2,2,sum),"/")
        Dmat=t(K2) %*% K2
#        Dmat=sdKmat$v %*% diag(sdKmat$d^2) %*% t(sdKmat$v)
        Dmat=Dmat+Dmat[(size-2*len):1,(size-2*len):1]
    }
    dvec.minus = t(Kmat) %*% data.minus
    dvec.plus  = rev(t(Kmat) %*% rev(data.plus))
    dvec = dvec.plus+dvec.minus
###Constraints are: prob[1:N] = prob[Aind[2,]] >= bvec[1:N] = 0
    N = nrow(Dmat)
#    bvec = rep(0,N)
    Aind = matrix(c(rep(1,N),1:N),ncol=N,nrow=2,byrow=T)
    Amat = matrix(1,ncol=N,nrow=1)
    list(Dmat=Dmat,dvec=dvec,Amat=Amat,Aind=Aind)  #bvec=bvec,
}

solve.one = function(counts,mu,lambda,len,reg,ktype) {
    allout=list(rtn=vector(mode="list", length=length(counts)),score=0)
    names(allout$rtn)=names(counts)
    ngood=0
    for (chr in names(counts)) {
#        allout$rtn[[chr]] = vector(mode="list", length=length(counts[[chr]]))
#        score = 0
        for (n in 1:length(counts[[chr]])) {
            cnt=counts[[chr]][[n]]
            N=length(cnt$minus)
            if (N > 5000 || N < 10 || sum(c(cnt$minus,cnt$plus)^2)<2) {
                allout$rtn[[chr]][[n]] = list(prob=rep(0,N),value=1)
                next
            }
            D = make.matrix(cnt$minus,cnt$plus,mu,lambda,len,reg,ktype)
            N = nrow(D$Dmat)
            prob = solve.QP.compact(Dmat=D$Dmat,dvec=D$dvec,Amat=D$Amat,Aind=D$Aind)
                                        #,bvec=D$bvec,factorized=TRUE)
            padpr=c(rep(0,len),prob$sol,rep(0,len))
            fit = fit.solution(padpr,mu,lambda,len,ktype)
            val = sqrt(sum(c(cnt$minus-fit$minus,cnt$plus-fit$plus)^2)/sum(c(cnt$minus,cnt$plus)^2))
#            allout$score = allout$score+prob$value
            allout$score = allout$score+val
            ngood=ngood+1
            allout$rtn[[chr]][[n]] = list(prob=padpr,value=val)
        }
    }
    allout$score = allout$score/ngood
    allout
}

fit.score = function(par,counts,len,reg,ktype) {
    if (par[1] > len && par[2] > 0) {
        y = try(solve.one(counts,par[1],par[2],len,reg,ktype)$score)
        if (class(y) != "try-error") return(y)
    }
    1
}

inverse.solve = function(counts,mu=50,lambda=120,len=36,regul=1e-4,
  val.cut=.35,ktype="geometric",optimize=FALSE) {
    if (mu <= len) stop("mu must be larger than len\n")
    par = c(mu,lambda)
    if (optimize) {
        cnt.high=list()
        for (chr in names(counts[[1]])) {
            I = which(counts[[2]][[chr]]<val.cut)
            cnt.high[[chr]]=counts[[1]][[chr]][I]
        }
        par = optim(par=par,fn=fit.score,counts=cnt.high,
          len=len,reg=regul,ktype=ktype,
          gr=NULL,method='L-BFGS-B',lower=c(len+1,len+1),upper=c(1000,1000),
          control=list(maxit=100,trace=4))$par
    }
    solved = solve.one(counts[[1]],par[1],par[2],len,regul,ktype)
    return(list(sol=solved$rtn,par=list(mu=par[1],lambda=par[2],len=len)))
}

read.data = function(filename,chrlist=c(),reglist=list(),offset=FALSE,gff=FALSE) {
    print(paste("Read data from",filename))
    con=gzfile(filename)
    open(con,open="r")
    if (gff) {
        colc=gffcols
    } else {
        colc=bedcols
    }
    data=read.table(con,comment.char="#",sep="\t",colClasses=colc)
    close(con)
    if (offset) data[,2]=data[,2]+1
    dspl=split(data[,-1],data[,1])
    if (length(chrlist)) {
        dspl=dspl[intersect(chrlist,names(dspl))]
        if (length(reglist)) reglist=reglist[intersect(chrlist,names(reglist))]
    }
    if (length(reglist)) {
        for (c in names(reglist)) {
            I=c()
            for (nr in 1:nrow(reglist[[c]])) {
                I=c(I,which(dspl[[c]][,1]<reglist[[c]][nr,2]
                  & dspl[[c]][,2]>reglist[[c]][nr,1]))
            }
            dspl[[c]]=dspl[[c]][unique(I),]
        }
    }
    dspl
}

data.to.strands = function(chrlist,data,lambda=120) {
    outlist=list()
    for (chr in chrlist) {
        dm = data[[chr]][which(data[[chr]][,4]=="-"),-4]
        dp = data[[chr]][which(data[[chr]][,4]=="+"),-4]
        dpunroll = data.frame(
          pos=as.integer(unlist(apply(dp[,1:2],1,function(x){x[1]:x[2]+lambda/2}))),
          counts=as.numeric(rep(dp[,3],dp[,2]-dp[,1]+1)))
        dmunroll = data.frame(
          pos=as.integer(unlist(apply(dm[,1:2],1,function(x){x[1]:x[2]-lambda/2}))),
          counts=as.numeric(rep(dm[,3],dm[,2]-dm[,1]+1)))
        d = merge(dpunroll, dmunroll, by="pos", all=TRUE)
        d$prod = sqrt(d$counts.x*d$counts.y)
        dclean = matrix(nrow=nrow(d),ncol=3)
        dsize = 1
        clast = -1
        left = 0
        right = 0
        for (n in which(!is.na(d$prod))) {
            if (d$prod[n] == clast && d$pos[n] == right+1) {
                right = d$pos[n]
            } else {
                dclean[dsize,] = c(left,right,clast)
                dsize = dsize+1
                left = d$pos[n]
                right = d$pos[n]
                clast = d$prod[n]
            }
        }
        dclean[dsize,] = c(left,right,clast)
        outlist$plus[[chr]]=dp
        outlist$minus[[chr]]=dm
        outlist$prod[[chr]]=dclean[2:dsize,]
    }
    outlist
}

setMinCoverage = function(count,plot=FALSE) {
    w = c()
    cnt = c()
    for (chr in names(count)) {
        w = c(w,count[[chr]][,2]-count[[chr]][,1]+1)
        cnt = c(cnt,count[[chr]][,3])
    }
#    w = count[,2]-count[,1]+1
    ww = w/sum(w)
    mu = sum(cnt*ww)
    sig = sum((cnt-mu)^2*ww)
    print(paste("Mean/sd +:",mu,sqrt(sig)))
    p = 1-exp(-(10:1e4)*1e-3)
    foreground = quantile(rep(cnt, w),p)
    background = qpois(p=p,lambda=mu)
    dev = sd((foreground-background)[p<.5])
    T = background[max(which(foreground-background<2*dev))]
    if (plot) {
        plot(background,foreground-background,pch=20,main="Deviation from Poisson")
        abline(h=dev)
        abline(v=T,col='blue')
    }
    T
}

find.enriched.regions = function(pfactor,threshold,datalist) {
####### this dynamic programming approach: P. Bucher + G. Ambrosini
    penalty = -pfactor*threshold
    reg = list()
    for (chr in names(datalist)) {
        data = datalist[[chr]]
        reg[[chr]] = matrix(0,ncol=2,nrow=nrow(data))
        if (all(data[,3]<=threshold)) {
            reg[[chr]]=matrix(0,ncol=2,nrow=0)
            next
        }
        pos0 = 1
        back = rep(0,nrow=nrow(data))
        score = list(x=0,y=0)
### I should build the whole matrix first, then apply max/cumsum?
        for (n in 1:nrow(data)) {
            xs = c(score$x,score$y+penalty,score$x+2*penalty)
            increment = (data[n,3]-threshold)*(data[n,2]-data[n,1]+1)
            xs = xs + increment
            xs = xs + c(pos0-data[n,1],data[n,1]-pos0,data[n,1]-pos0-1)*threshold
            xi = which.max(xs)
            ys = c(score$x+penalty,score$y)
            ys = ys - increment + (data[n,1]-pos0)*threshold
            yi = which.max(ys)
            score$x = xs[xi]
            score$y = ys[yi]
            back[n] = xi+3*yi 
        # 4 = + -> + / + -> -
        # 5 = - -> + / + -> -
        # 6 = + -> - -> + / + -> -
        # 7 = + -> + / - -> -
        # 8 = - -> + / - -> -
        # 9 = + -> - -> + / - -> -
            pos0 = data[n,2]+1
        }
        cnts = c(0)
        rptr = nrow(reg[[chr]])
        if (score$x>score$y) {
            reg[[chr]][rptr,2] = pos0-1
        } else {
            reg[[chr]][rptr,2] = 0
        }
        for (n in nrow(data):2) {
            if (reg[[chr]][rptr,2] == 0) {
                 # in - region
                if (back[n] < 7) reg[[chr]][rptr,2] = data[n-1,2]
            } else {
                 # in + region
                if (any(back[n] == c(5,6,8,9))) {
                    reg[[chr]][rptr,1] = data[n,1]
                    cnts[1] = cnts[1]+data[n,3]*(data[n,2]-data[n,1]+1)
                    cnts = c(0,cnts)
                    rptr=rptr-1
                    if (any(back[n] == c(6,9))) reg[[chr]][rptr,2] = data[n-1,2]
                } else {
                    cnts[1] = cnts[1]+data[n,3]*(data[n,2]-data[n,1]+1)
                }
            }
        }
        reg[[chr]] = reg[[chr]][-(1:rptr),]
        if (reg[[chr]][1,2] == 0) reg[[chr]] = reg[[chr]][-1,]
        if (nrow(reg[[chr]])>0 && reg[[chr]][1,1] == 0) reg[[chr]][1,1] = data[1,1]
    }
    reg
}

merge.regions = function(regions,gap) {
    nr = nrow(regions)
    regions2 = matrix(0,nrow=nr,ncol=2) ### do it in place?
    rptr = 0
    regions = rbind(regions,c(1,1)*regions[nr,2]+gap+1)
    for (n in 1:nr+1) {
        if (regions[n,1]>regions[n-1,2]+gap) {
            if (regions[n-1,2]-regions[n-1,1]>=10) {
#                regions2 = rbind(regions2,regions[n-1,]+c(-1,1)*gap/2)
                rptr = rptr+1
                regions2[rptr,] = regions[n-1,]+c(-1,1)*gap/2
            }
        } else {
            regions[n,1] = regions[n-1,1]
        }
    }
    regions2[1,1] = max(1,regions2[1,1])
    matrix(regions2[1:rptr,],ncol=2)
}

fit.filter = function(regions,datalist,chr,
  mu=90,lambda=100,len=36,regul=1e-2,ktype="geometric",var.cut=.65) {
    nr = nrow(regions)
    if (is.null(nr) || nr==0) return(regions)
    regs = cbind(regions,rep(1,nrow(regions)))
    n1 = 0
    for (n in 1:nr) {
        allpos = regions[n,1]:regions[n,2] 
        cp = make.counts(allpos,datalist$plus[[chr]])[,2]
        cm = make.counts(allpos,datalist$minus[[chr]])[,2]
        D = make.matrix(cm,cp,mu,lambda,len,regul,ktype)
        N = nrow(D$Dmat)
        if (N > 3000 || N < 10) {next}
        prob = solve.QP.compact(Dmat=D$Dmat,dvec=D$dvec,Amat=D$Amat,Aind=D$Aind)
                                        #,bvec=D$bvec,factorized=TRUE)
        padpr=c(rep(0,len),prob$sol,rep(0,len))
        fit = fit.solution(padpr,mu,lambda,len,ktype)
        val = sqrt(sum(c(cm-fit$minus,cp-fit$plus)^2)/sum(c(cm,cp)^2))
        if (!is.na(val) && val<var.cut) {
            n1 = n1+1
            regs[n1,] = c(regions[n,1:2],val)
        } else {
            print(paste("Region [",regions[n,1],",",regions[n,2],"] rejected: val =",val))
        }
    }
    regs[1:n1,]
}

regions2counts = function(regionlist,datalist,threshold)  {
    counts = list()
    vals = list()
    for (chr in intersect(names(regionlist),names(datalist$prod))) {
        regions = regionlist[[chr]]
        nr = nrow(regions)
        if (is.null(nr) || nr==0) next
        for (n in 1:nr) {
            allpos = regions[n,1]:regions[n,2]
            cnts = cbind(
              make.counts(allpos,datalist$prod[[chr]]),
              make.counts(allpos,datalist$plus[[chr]])[,2],
              make.counts(allpos,datalist$minus[[chr]])[,2])
            colnames(cnts) = c("pos","prod.counts","plus.counts","minus.counts")
            if ((threshold < 0) ||
                (max(cnts[,2]) > threshold && all(apply(cnts[,-1],2,sum) > 0))) {
                counts[[chr]] = c(counts[[chr]],list(cnts))
                if (ncol(regions)>2) vals[[chr]] = c(vals[[chr]],regions[n,3])
                else vals[[chr]] = c(vals[[chr]],-1)
            }
        }
    }
    list(counts,vals)
}

make.counts = function(allpos,data,shift=0) {
    pmin = allpos[1]
    pmax = allpos[length(allpos)]
    I = which(data[,1]<=pmax & data[,2]>=pmin)
    if (length(I)>1) {
        unroll = data.frame(
          pos=as.integer(unlist(apply(data[I,1:2],1,function(x){x[1]:x[2]+shift}))),
          counts=as.numeric(rep(data[I,3],data[I,2]-data[I,1]+1)))
    } else if (length(I)==1) {
        unroll = data.frame(
          pos=as.integer(data[I,1]:data[I,2]+shift),
          counts=as.numeric(rep(data[I,3],data[I,2]-data[I,1]+1)))
    } else {
        unroll=data.frame(pos=allpos,counts=rep(0,length(allpos)))
    }
    cnts = merge(allpos,unroll,by.x=1,by.y="pos",all=T)
    if (any(is.na(cnts[,2])))    cnts[which(is.na(cnts[,2])),2] = 0
    if (cnts[1,1]<pmin)          cnts = cnts[-which(cnts[,1]<pmin),]
    if (cnts[nrow(cnts),1]>pmax) cnts = cnts[-which(cnts[,1]>pmax),]
    cnts
}

plot.counts = function(counts,data,cumul=TRUE) {
    if (cumul) {
        ymax = sum(data[,3])
        plot(data[,1]*1e-3,cumsum(data[,3]),type='l',ylim=c(0,ymax),
             xlab="Position [kbp]",ylab="Cumulative counts",main="Enriched regions")
        for (c in counts) 
          lines(x=range(c$pos)*1e-3,y=c(ymax,ymax)/2,col='cyan',lwd=5)
    } else {
        ymax = max(data[,3])
        plot(data[,1]*1e-3,data[,3],type='l',ylim=c(0,ymax),
             xlab="Position [kbp]",ylab="Counts",main="Enriched regions")
        for (c in counts) 
          lines(x=range(c$pos)*1e-3,y=c(0,0),col='cyan',lwd=5)
    }
}

plot.sol = function(counts,solution,par,ktype="geometric",legend=TRUE,ymax) {
    xlab=''
    if (legend) xlab="Position [kbp]"
    if (length(counts)*length(solution)*length(par) == 0) {
        plot(0,0,xlab=xlab,ylab="Counts",main='',ylim=c(0,1))
        return(0)
    }
    pos = counts$pos*1e-3
    prob = solution$prob
    if (missing(ymax)) ymax = max(c(counts$minus,counts$plus,prob))
    main = paste("Rel. fit error =",round(solution$value,digits=4))
    plot(pos,counts$minus,col='red',type='s',
         xlab=xlab,ylab="Counts",main=main,lwd=2,ylim=c(0,ymax))
    lines(pos,counts$plus,col='blue',lwd=.5,t='s')
    if (length(prob)) {
        fit.sol = fit.solution(prob,par$mu,par$lambda,par$len,ktype)
        lines(pos,prob,col='black',lwd=2)
        lines(pos,fit.sol$plus,col='cyan',lwd=2)
        lines(pos,fit.sol$minus,col='magenta',lwd=2)
        if (legend) legend("topleft",c("data +","data -","prob","fit +","fit -"),
                           lty=1,lwd=2,col=c("blue","red","black","cyan","magenta"))
    }
}

write.bed.counts = function(data,header,file,offset=TRUE) {
    bounds = t(sapply(data,function(x){range(x$pos)},simplify=TRUE))
    if (offset) bounds[,1] = bounds[,1]-1
    score = sapply(counts,function(x){(sum(x$plus)+sum(x$minus))/2},simplify=TRUE)
    name = paste("region",1:length(data),sep='')
    dbed = rbind(c(header[1],rep('',3)),cbind(header[2],bounds,name,score))
    write.table(dbed,file=file,row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}

write.bed.sol = function(sol,counts,header,file,offset=TRUE,append=FALSE,cutoff=1e-3,gff=FALSE) {
    result = matrix(nrow=0,ncol=4)
    for (n in 1:length(counts)) {
        I = which(sol[[n]]$prob>cutoff*sum(sol[[n]]$prob))
        if (length(I)<2) next
        interval = range(counts[[n]]$pos[I])
        score = sum(sol[[n]]$prob[I])
        result = rbind(result,c(interval,score,round(sol[[n]]$val,digits=4)))
    }
    if (offset) result[,1] = result[,1]-1
    name = paste("ID=",header[2],"_",1:nrow(result),";FERR=",result[,4],sep='')
    if (nrow(result)==1) {
        dbed = c(header[2],result[,1:2],name,result[,3])
    } else {
        dbed = cbind(header[2],result[,1:2],name,result[,3])
    }
    if (gff) {
        dbed = cbind(dbed[,1],"deconvolution","binding_site",dbed[,c(2,3,5,6)],".",name)
    } else {
        if (!append) dbed = rbind(c(header[1],rep('',4)),dbed)
    }
    write.table(dbed,file=file,row.names=FALSE,col.names=FALSE,quote=FALSE,
                append=append,sep="\t")
}

write.wig = function(data,header,file,offset=TRUE,append=FALSE,cutoff=1e-3) {
    if (nrow(data)<10) return(0)
#    I0 = range(which(data[,3]>cutoff))
#    I = I0[1]:I0[2]
    I=which(data[,3]>cutoff*sum(data[,3]))
    if (offset) data[,1] = data[,1]-1
    dwig = cbind(rep(as.character(header[2]),length(I)),
      as.integer(data[I,1]),as.integer(data[I,2]),as.numeric(data[I,3]))
#    I = which(dwig[,4]<=cutoff)
#    dwig[I,4] = 0
    if (!append) {dwig = rbind(c(header[1],rep('',3)),dwig)}
    write.table(dwig,file=file,row.names=FALSE,col.names=FALSE,quote=FALSE,append=append,sep="\t")
}
