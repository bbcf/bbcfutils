#!/usr/bin/env python

import optparse, os
from bbcflib import btrack as track
from bbcflib.bFlatMajor.stream import merge_scores
from bbcflib.bFlatMajor.numeric import correlation

parameters = {
'usage'   : '%prog [OPTIONS] -f fwd_file -r rev_file',
'description' : ' Merges two signal tracks after shifting them.',
'option_list' : [optparse.make_option("-f", "--forward",
                                      help="A bedgraph-like file with ChIP density on the forward strand"),
                 optparse.make_option("-r", "--reverse",
                                      help="A bedgraph-like file with ChIP density on the reverse strand"),
                 optparse.make_option("--formatf",
                                      help="Format of the forward track."),
                 optparse.make_option("--formatr",
                                      help="Format of the reverse track."),
                 optparse.make_option("-p", "--shift",
                                      help="Shift positions downstream. If <0 will compute an optimal shift using autocorrelation.", default=0, action="store", type=int),
                 optparse.make_option("-o", "--output",
                                      help="Output file", default='merged.sql'),
                 optparse.make_option("-a", "--assembly",
                                      help = "The assembly (GenRep assembly name such as 'sacCer2').")]}

parser = optparse.OptionParser(**parameters)
options, args = parser.parse_args()
if not(options.forward and os.path.exists(options.forward)):
    parser.print_help()
    raise ValueError("Specify a valid forward strand density file with -f.")
if not(options.reverse and os.path.exists(options.reverse)):
    parser.print_help()
    raise ValueError("Specify a valid reverse strand density file with -r.")


def _shift(stream,shift):
    istart = stream.fields.index('start')
    iend   = stream.fields.index('end')
    i1 = min(istart,iend)
    i2 = max(istart,iend)
    def _apply_shift(x):
        return x[:i1]+(x[i1]+shift,)+x[i1+1:i2]+(x[i2]+shift,)+x[i2+1:]
    return track.FeatureStream((_apply_shift(x) for x in stream),
                               fields=stream.fields)
fields = ['chr','start','end','score']
tfwd = track.track(options.forward,format=options.formatf,chrmeta=options.assembly)
trev = track.track(options.reverse,format=options.formatr,chrmeta=options.assembly)
if tfwd.chrmeta:
    chrmeta = tfwd.chrmeta
elif trev.chrmeta:
    chrmeta = trev.chrmeta
else:
    parser.print_help()
    raise ValueError("Specify an assembly with -a.")

shiftval = int(options.shift)
if shiftval < 0:
    slim = 300
    chrsize,chrom = sorted([(v['length'],k) for k,v in chrmeta.iteritems()],reverse=True)[0]
    xcor = correlation([tfwd.read(chrom),trev.read(chrom)],1,chrsize,limits=(-slim,slim))
    shiftval = (xcor[0].argmax()-slim-1)/2
    print "Autocorrelation shift=%i, correlation is %f." %(shiftval,xcor[0].max())

tout = track.track(options.output,fields=fields,chrmeta=chrmeta,info={'datatype': 'quantitative'})
mode = 'write'
for chrom in chrmeta.keys():
    tout.write(merge_scores([_shift(tfwd.read(selection=chrom), shiftval),
                             _shift(trev.read(selection=chrom),-shiftval)]),
               chrom=chrom,mode=mode,clip=True)
    mode = 'append'
tout.close()
trev.close()
tfwd.close()
