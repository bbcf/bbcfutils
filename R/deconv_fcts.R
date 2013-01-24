library(quadprog)

read.data = function(file) {
    load(file)
    split(counts[,c("pos","plus","minus")],counts$name)
}

cross.correlate = function(counts,threshold=1,lag.max=600) {
    signal=matrix(0,ncol=2,nrow=0)
    for (c in counts) {
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
    for (n in names(counts)) {
        cnt=counts[[n]]
        N=length(cnt$minus)
        if (N > 5000 || N < 2*len+lambda || sum(c(cnt$minus,cnt$plus)^2)<2) {
            allout$rtn[[n]] = list(prob=rep(0,N),value=1)
            next
        }
        D = make.matrix(cnt$minus,cnt$plus,mu,lambda,len,reg,ktype)
        N = nrow(D$Dmat)
        prob = solve.QP.compact(Dmat=D$Dmat,dvec=D$dvec,Amat=D$Amat,Aind=D$Aind)
                                #,bvec=D$bvec,factorized=TRUE)
        padpr=c(rep(0,len),prob$sol,rep(0,len))
        fit = fit.solution(padpr,mu,lambda,len,ktype)
        val = sqrt(sum(c(cnt$minus-fit$minus,cnt$plus-fit$plus)^2)/sum(c(cnt$minus,cnt$plus)^2))
        #allout$score = allout$score+prob$value
        allout$score = allout$score+val
        ngood=ngood+1
        allout$rtn[[n]] = list(prob=padpr,value=val)
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

inverse.solve = function(counts,
  mu=50,lambda=120,len=36,regul=1e-4,
  ktype="geometric",optimize=FALSE) {
    if (mu <= len) stop("mu must be larger than len\n")
    par = c(mu,lambda)
    if (optimize) {
        solved = solve.one(counts,par[1],par[2],len,1e-3,ktype)
        O = order(sapply(solved$rtn,function(x)x$value))
        npeaks = min(5,length(O))
        par = optim(par=par,fn=fit.score,counts=counts[O[1:npeaks]],
          len=len,reg=1e-3,ktype=ktype,
          gr=NULL,method='L-BFGS-B',lower=c(len+1,len+1),upper=c(600,600),
          control=list(maxit=5,trace=4))$par
    }
    solved = solve.one(counts,par[1],par[2],len,regul,ktype)
    return(list(sol=solved$rtn,par=list(mu=par[1],lambda=par[2],len=len)))
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

