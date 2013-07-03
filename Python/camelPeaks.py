#!/usr/bin/env python

"""
Wrapper to run the deconvolution algorithm for ChIP-seq peak.
Need a list of peaks and a density file in input.

"""

from bbcflib.track import track
from bbcflib import genrep
from bbcflib.gfminer.common import sorted_stream
import rpy2.robjects as robjects
import sys, optparse, os

opts = (("-p", "--peaks",
         "A bed-like file with selected enriched regions where deconvolution will be performed",{}),
        ("-f", "--forward",
         "A bedgraph-like file with ChIP density on the forward strand",{}),
        ("-r", "--reverse",
         "A bedgraph-like file with ChIP density on the reverse strand",{}),
        ("-o", "--output", "Output prefix (default 'camels')",{'default':'camels'}),
        ("-c", "--chromosome", "Chromosome name",{'default':None}),
        ("-l", "--length", "Chromosome length",{'default':sys.maxint,'type':"int"}),
        ("-s", "--sizecutoff",
         "Maximum region size to consider (default 3000)",{'default': 3000,'type':"int"}),
        ("-e", "--extension",
         "Read extension used for computing the density (default 40)",{'default': 40,'type':"int"}),
        ("-x", "--mu", "Lower cutoff to fragment size (default 80)",{'default': 80,'type':"int"}),
        ("-z", "--script", "R script path (default ./)",{'default': './'}),
        ("-g", "--genome", "Genome assembly code",{'default': None}),
        ("-k", "--kernel", "Deconvolution kernel ('gaussian' or 'geometric', default 'geometric')",{'default': "geometric"})
        )

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    try:
        usage = "camelPeaks.py [OPTIONS]"
        desc = """A ChIP-seq peak deconvolution algorithm."""
        parser = optparse.OptionParser(usage=usage, description=desc)
        for opt in opts:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
        (opt, args) = parser.parse_args()
        if not(opt.peaks and os.path.exists(opt.peaks)):
            parser.print_help()
            raise Usage("Specify a valid peaks file with -p.")
        if not(opt.forward and os.path.exists(opt.forward)):
            parser.print_help()
            raise Usage("Specify a valid forward strand density file with -f.")
        if not(opt.reverse and os.path.exists(opt.reverse)):
            parser.print_help()
            raise Usage("Specify a valid reverse strand density file with -r.")
####
        if opt.chromosome and opt.length: chrmeta = {opt.chromosome: {'length': opt.length}}
        else: chrmeta = opt.genome
        peak_track = track(opt.peaks,chrmeta=chrmeta)
        chrmeta = peak_track.chrmeta
        if opt.chromosome: chrmeta = {opt.chromosome: chrmeta[opt.chromosome]}
        track_info = {'datatype': peak_track.info.get('datatype','qualitative')}
        outbed = track(opt.output+"_peaks.bed", chrmeta=chrmeta,
                             fields=["chr","start","end","name","score"])
        outwig = track(opt.output+"_deconv.bedgraph", chrmeta=chrmeta)
        outwig.open(mode='overwrite')
        topts = {'chrmeta': chrmeta, 'readonly': True}
        for chrom,cv in chrmeta.iteritems():
            peak_stream = sorted_stream(peak_track.read(selection=chrom),[chrom])
            strands = {track(opt.forward,**topts).read(chrom,fields=['start','end','score']): 'plus',
                       track(opt.reverse,**topts).read(chrom,fields=['start','end','score']): 'minus'}
            robjects.r('options(stringsAsFactors=F)')
            robjects.r('counts=data.frame()')
            for row_count,peak in enumerate(peak_stream):
                start = int(peak[peak_stream.fields.index('start')])
                end = int(peak[peak_stream.fields.index('end')])
                if end-start > opt.sizecutoff: continue
                if start < 0: start = 0
                if not(end <= cv['length']): end = cv['length']
                if 'name' in peak_stream.fields:
                    reg_name = peak[peak_stream.fields.index('name')]
                else:
                    reg_name = str(row_count+1)
                data_block = robjects.DataFrame({'pos':   robjects.IntVector(range(start+1,end+1)),
                                                 'plus':  robjects.FloatVector([0]*(end-start)),
                                                 'minus': robjects.FloatVector([0]*(end-start)),
                                                 'name':  robjects.StrVector([reg_name]*(end-start))})
                for stream,strnd in strands.iteritems():
                    for row in stream:
                        if row[0]<start: continue
                        if row[1]>end: break
                        data_block.rx2(strnd)[(row[0]-start):(row[1]-start)] = \
                            robjects.FloatVector([row[2]]*(row[1]-row[0]))
                robjects.r.assign('newblock',data_block)
                robjects.r('counts=rbind(counts,newblock)')
            robjects.r('read.length=%i' %opt.extension)
            robjects.r('chr.name="%s"' %chrom)
            robjects.r('pdf.file="%s.pdf"' %opt.output)
            robjects.r('mu=%i' %opt.mu)
            robjects.r('ktype="%s"' %opt.kernel)
            robjects.r('source("%s")' %os.path.join(opt.script,"deconv_fcts.R"))
            robjects.r("""
    counts = split(counts[,c("pos","plus","minus")],counts$name)
    pdf(file=pdf.file,title='chip-seq',paper='a4',width=8,height=11)
    par(cex=1.5,lwd=1.5)
    ccf = cross.correlate(counts,threshold=.5)
    plot(ccf$lag,ccf$acf,t='l',ylim=c(0,1),
         xlab='Lag',ylab='Cross-correlation',
         main=paste('Strand cross-correlation',chr.name))
    cut.ccf = ccf$acf
    cut.ccf[which(ccf$lag<mu)] = 0
    lambda = ccf$lag[which.max(cut.ccf)]
    sol = inverse.solve(counts,mu=mu,lambda=lambda,len=read.length,regul=1e-3,optimize=TRUE,ktype=ktype)
    col = 'red'
    lab = paste('lambda=',sol$par$lambda,sep='')
    abline(v=sol$par$lambda,col=col)
    text(sol$par$lambda,0,lab,col=col,pos=4)
    col = 'blue'
    lab = paste('mu=',sol$par$mu,sep='')
    abline(v=sol$par$mu,col=col)
    text(sol$par$mu,0.3,lab,col=col,pos=4)
    col = 'darkgreen'
    lab = paste('l=',read.length,sep='')
    abline(v=read.length,col=col)
    text(read.length,0.6,lab,col=col,pos=4)
    par(mfrow=c(4,2))
    for (n in names(counts)) {
      if (sol$sol[[n]]$value>.65) next
      plot.sol(counts[[n]],sol$sol[[n]],sol$par)
      title(sub=chr.name)
    }
    dev.off()
    bed = data.frame()
    cutoff = 1e-3
    for (n in names(counts)) {
      I = which(sol$sol[[n]]$prob>cutoff*sum(sol$sol[[n]]$prob))
      if (length(I)<2) next
      interval = range(counts[[n]]$pos[I])
      score = sum(sol$sol[[n]]$prob[I])
      name = paste('ID=',n,';FERR=',round(sol$sol[[n]]$val,digits=4),sep='')
      bed = rbind(bed,data.frame(
          start=interval[1],end=interval[2],
          name=name,score=score))
    }
    bed[,'start'] = as.integer(bed[,'start']-1)
    wig = data.frame()
    for (n in names(counts)) {
      I = which(sol$sol[[n]]$prob>cutoff*sum(sol$sol[[n]]$prob))
      wig = rbind(wig,data.frame(
      pos = as.integer(counts[[n]]$pos[I]),
      score = as.numeric(sol$sol[[n]]$prob[I])))
    }
    """)
            nrow = robjects.r("nrow(bed)")[0]
            outbed.write(((robjects.r("bed").rx2('start')[ri],
                           robjects.r("bed").rx2('end')[ri],
                           robjects.r("bed").rx2('name')[ri],
                           robjects.r("bed").rx2('score')[ri]) for ri in xrange(nrow)),
                         fields=["start","end","name","score"], chrom=chrom, mode='append')
            nrow = robjects.r("nrow(wig)")[0]
            outwig.write(((robjects.r("wig").rx2('pos')[ri]-1,
                           robjects.r("wig").rx2('pos')[ri],
                           robjects.r("wig").rx2('score')[ri]) for ri in xrange(nrow)),
                         fields=["start","end","score"], chrom=chrom, mode='append')
        outwig.close()
        print "************OUTPUT FILES**********"
        print "\n".join([opt.output+".pdf",
                         opt.output+"_peaks.bed",
                         opt.output+"_deconv.bedgraph"])
        print "************PARAMETERS**********"
        print "lambda=%f|mu=%f|len=%f" %(robjects.r("sol$par$lambda")[0],robjects.r("sol$par$mu")[0],robjects.r("read.length")[0])
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(2)

if __name__ == '__main__':
    sys.exit(main())

