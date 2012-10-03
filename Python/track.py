#!/usr/bin/env python

import optparse, json, sys, os
from bbcflib import btrack
from bbcflib.bFlatMajor.stream import merge_scores
from bbcflib.bFlatMajor.numeric import correlation

functions = ["convert","read","merge"]
usage = {'all': "track.py %s [OPTIONS]"}
description = {'all': "Command-line interface to bbcflib.btrack functionalities."}
opts = {'all':
        (("-a", "--assembly", "Assembly name or id",{'default':None}),
         ("-o", "--output", "Output file (default stdout)",{'default':None}),
         ("-c", "--chrmeta", "Json-formatted chrmeta dictionary (if not an assembly)", {'default': None}))}

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

############## CONVERT ##############
f = 'convert'
usage[f] = usage['all'] %f +" file_in file_out"
description[f] = 'Converts between file types.'
opts[f] = (("-1", "--format1", "File format of first input file if extension not explicit",{'default':None}),
           ("-2", "--format2", "File format of second input file if extension not explicit",{'default':None}),
           ("-d", "--datatype", "Datatype info in output track", {}))

def convert(*args,**kw):
    if len(args) != 2:
        if len(args) == 1 and kw['output']:
            args.append(kw['output'])
        else:
            raise Usage("Please specify SOURCE and DESTINATION files")
    if kw['format1']:
        args[0] = (args[0],kw['format1'])
    if kw['format2']:
        args[1] = (args[1],kw['format2'])
    if kw['chrmeta']:
        chrmeta = json.loads(kw['chrmeta'])
    else:
        chrmeta = kw['assembly']
    if kw['datatype']:
        info = {'datatype': str(kw['datatype'])}
    else:
        info = None
    btrack.convert(*args, chrmeta=chrmeta, info=info)
    return 0

############## READ ##############
f = 'read'
usage[f] = usage['all'] %f +" file1 [file2 ...]"
description[f] = 'A generic genomic track data reader.'
opts[f] = (("-t", "--format", "File format if extension not explicit",{'default':None}),
        ("-s", "--selection", "Selection dictionary or coma-separated list of chromosome names",
         {'default':None}),
        ("-f", "--fields", "Fields in output",{'default':None}),
        ("-d", "--description", "Only print a description of the file",{'action':"store_true"}))

def read(*args,**kw):
    if len(args) < 1: raise Usage("No input file provided")
    selection = None
    if kw['selection']:
        if kw['selection'].count("{"):
            jsonargs = json.loads(kw['selection'])
            for k,v in jsonargs.iteritems():
                if isinstance(v,basestring): jsonargs[k] = str(v)
            selection = dict((str(k),v) for k,v in jsonargs.iteritems())
        else:
            selection = str(kw['selection']).split(",")
    fields = None
    if kw['fields']:
        fields = str(kw['fields']).split(",")
    if kw['output'] is None:
        output = sys.stdout
    else:
        output = open(kw['output'],'w')
    for infile in args:
        intrack = btrack.track(infile,format=kw['format'],chrmeta=kw['assembly'])
        if kw['description']:
            fileinfo = ",".join(["%s=%s" %(k,v) for k,v in intrack.info.iteritems()])
            chromlist = ",".join(sorted(intrack.chrmeta.keys()))
            fields = ",".join(intrack.fields)
            output.write(\
"""*****************************
%s (%s):
infos: %s
chromosomes: %s
fields: %s
*****************************
""" %(os.path.basename(infile), intrack.format, fileinfo, chromlist, fields))
            continue
        for x in intrack.read(selection=selection,fields=fields):
            output.write("\t".join([str(y) for y in x])+"\n")
    return 0

############## MERGE ##############
f = 'merge'
usage[f] = usage['all'] %f +" -f fwd_file -r rev_file"
description[f] = 'Merges two signal tracks after shifting each in downstream direction.'
opts[f] = (("-f", "--forward", "A bedgraph-like file with ChIP density on the forward strand", {}),
        ("-r", "--reverse", "A bedgraph-like file with ChIP density on the reverse strand", {}),
        ("-1","--formatf", "Format of the forward track.", {}),
        ("-2","--formatr", "Format of the reverse track.", {}),
        ("-p", "--shift", "Shift positions downstream. If <0 will compute an optimal shift using autocorrelation.", {'default':0, 'action':"store", 'type':int}))

def merge(*args,**kw):
    if not(kw['forward'] and os.path.exists(kw['forward'])):
        raise Usage("Specify a valid forward strand density file with -f.")
    if not(kw['reverse'] and os.path.exists(kw['reverse'])):
        raise Usage("Specify a valid reverse strand density file with -r.")

    def _shift(stream,shift):
        istart = stream.fields.index('start')
        iend   = stream.fields.index('end')
        i1 = min(istart,iend)
        i2 = max(istart,iend)
        def _apply_shift(x):
            return x[:i1]+(x[i1]+shift,)+x[i1+1:i2]+(x[i2]+shift,)+x[i2+1:]
        return btrack.FeatureStream((_apply_shift(x) for x in stream),
                                    fields=stream.fields)
    fields = ['chr','start','end','score']
    tfwd = btrack.track(kw['forward'],format=kw['formatf'],chrmeta=kw['assembly'])
    trev = btrack.track(kw['reverse'],format=kw['formatr'],chrmeta=kw['assembly'])
    if tfwd.chrmeta:
        chrmeta = tfwd.chrmeta
    elif trev.chrmeta:
        chrmeta = trev.chrmeta
    else:
        raise Usage("Specify an assembly with -a.")

    shiftval = int(kw['shift'])
    if shiftval < 0:
        slim = 300
        chrsize,chrom = sorted([(v['length'],k) for k,v in chrmeta.iteritems()],reverse=True)[0]
        xcor = correlation([tfwd.read(chrom),trev.read(chrom)],(1,chrsize),limits=(-slim,slim))
        shiftval = (xcor.argmax()-slim-1)/2
        print "Autocorrelation shift=%i, correlation is %f." %(shiftval,xcor.max())

    tout = btrack.track(kw['output'],fields=fields,chrmeta=chrmeta,info={'datatype': 'quantitative'})
    mode = 'write'
    for chrom in chrmeta.keys():
        tout.write(merge_scores([_shift(tfwd.read(chrom), shiftval),
                                 _shift(trev.read(chrom),-shiftval)]),
                   chrom=chrom,mode=mode,clip=True)
        mode = 'append'
    tout.close()
    trev.close()
    tfwd.close()
    return 0

usage['all'] = usage['all'] %str(functions)
def main(argv=None):
    parser = None
    try:
        fct = None
        if len(sys.argv) > 1: fct = sys.argv.pop(1)
        if fct not in functions:
            m = fct and "No such operation: %s, choose one of %s." %(fct,str(functions)) or ''
            fct = 'all'
            raise Usage(m)

        parser = optparse.OptionParser(usage=usage[fct],
                                       description=description[fct])
        for opt in opts['all']+opts[fct]:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
        (options, args) = parser.parse_args()
        return eval(fct)(*args,**options.__dict__)

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage[fct]
        if parser: parser.print_help()
        return 1

if __name__ == '__main__':
    sys.exit(main())

