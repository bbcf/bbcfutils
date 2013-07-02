#!/usr/bin/env python

import optparse, json, sys, os
from bbcflib import track
from bbcflib.gfminer.common import sorted_stream
from bbcflib.gfminer.stream import merge_scores
from bbcflib.gfminer.numeric import correlation

functions = ["convert","read","merge","stats","check","sort"]
usage = {'all': "track.py %s [OPTIONS]"}
description = {'all': "Command-line interface to bbcflib.track functionalities."}
opts = {'all':
        (("-a", "--assembly", "Assembly name or id. If not standard, use `-a guess` \
                             (text formats only).", {'default':None}),
         ("-o", "--output", "Output file (default stdout).", {'default':None}),
         ("-c", "--chrmeta", "(if not assembly specified) Json-formatted chrmeta dictionary, or a tab-delimited file with two columns: <chromosome name> <size in bases>.",
                            {'default': None}))}

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def _get_chrmeta(**kw):
    if kw['chrmeta']:
        if kw['chrmeta'] == "guess":
            chrmeta = kw['chrmeta']
        elif os.path.exists(kw['chrmeta']):
            chrmeta = {}
            with open(kw['chrmeta'],'r') as chr_sizes:
                for line in chr_sizes:
                    try:
                        chrname,chrsize = line.split()
                        chrmeta[chrname] = {'length':int(chrsize)}
                    except ValueError:
                        raise ValueError("Wrong line in chr.sizes file:\n%s." % line)
        else:
            chrmeta = json.loads(kw['chrmeta'])
    else:
        chrmeta = kw['assembly'] # None by default
    return chrmeta

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
    chrmeta = _get_chrmeta(**kw)
    if kw['datatype']:
        info = {'datatype': str(kw['datatype'])}
    else:
        info = None
    track.convert(*args, chrmeta=chrmeta, info=info)
    return 0

############## READ ##############
f = 'read'
usage[f] = usage['all'] %f +" file1 [file2 ...]"
description[f] = 'A generic genomic track data reader.'
opts[f] = (("-t", "--format", "File format, if extension not explicit", {'default':None}),
        ("-s", "--selection", "Selection or comma-separated list of chromosome names. Selection \
         can be a dictionary (json) or a string of the form `chr:start-end`.", {'default':None}),
        ("-f", "--fields", "Fields in output", {'default':None}),
        ("-d", "--description", "Only print a description of the file", {'action':"store_true"}))

def read(*args,**kw):
    if len(args) < 1: raise Usage("No input file provided")
    selection = None
    if kw['selection']:
        if kw['selection'].count("{"):
            jsonargs = json.loads(kw['selection'])
            for k,v in jsonargs.iteritems():
                if isinstance(v,basestring): jsonargs[k] = str(v)
            selection = dict((str(k),v) for k,v in jsonargs.iteritems())
        elif kw['selection'].count(":"):
            chr,coord = kw['selection'].split(':')
            start,end = coord.split('-')
            selection = {'chr':chr,'start':(start,end),'end':(start,end)}
        else:
            selection = str(kw['selection']).split(",")
    fields = None
    if kw['fields']:
        fields = str(kw['fields']).split(",")
    if kw['output'] is None:
        output = sys.stdout
    else:
        output = open(kw['output'],'w')
    chrmeta = _get_chrmeta(**kw)
    for infile in args:
        intrack = track.track(infile,format=kw['format'],chrmeta=chrmeta)
        if kw['description']:
            if intrack.info:
                fileinfo = ",".join(["%s=%s" %(k,v) for k,v in intrack.info.iteritems()])
            else: fileinfo = 'None'
            chromlist = ",".join(sorted(intrack.chrmeta.keys())) or "None"
            fields = ",".join(intrack.fields)
            output.write(\
"""# *****************************************
# File '%s' (%s):
# Infos: %s
# Chromosomes: %s
# Fields: %s
# *****************************************
""" %(os.path.basename(infile), intrack.format, fileinfo, chromlist, fields))
            continue
        for x in intrack.read(selection=selection,fields=fields):
            output.write("\t".join([str(y) for y in x])+"\n")
        intrack.close()
    try: output.close()
    except IOError: pass # if stdout
    return 0

############## MERGE ##############
f = 'merge'
usage[f] = usage['all'] %f +" -f fwd_file -r rev_file"
description[f] = 'Merges two signal tracks after shifting each in downstream direction.'
opts[f] = (("-f", "--forward", "A bedgraph-like file with ChIP density on the forward strand", {}),
        ("-r", "--reverse", "A bedgraph-like file with ChIP density on the reverse strand", {}),
        ("-1","--formatf", "Format of the forward track.", {}),
        ("-2","--formatr", "Format of the reverse track.", {}),
        ("-p", "--shift", "Shift positions downstream. If <0 will compute an optimal shift \
                           using autocorrelation.", {'default':0, 'action':"store", 'type':int}))

def merge(*args,**kw):
    if not(kw['forward'] and os.path.exists(kw['forward'])):
        raise Usage("Specify a valid forward strand density file with -f.")
    if not(kw['reverse'] and os.path.exists(kw['reverse'])):
        raise Usage("Specify a valid reverse strand density file with -r.")
    if not(kw['output']):
        raise Usage("Specify the output file name.")

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
    chrmeta = _get_chrmeta(**kw)
    tfwd = track.track(kw['forward'],format=kw['formatf'],chrmeta=chrmeta)
    trev = track.track(kw['reverse'],format=kw['formatr'],chrmeta=chrmeta)
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

    tout = track.track(kw['output'],fields=fields,chrmeta=chrmeta,info={'datatype':'quantitative'})
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

############## STATS ##############
f = 'stats'
usage[f] = usage['all'] %f +" file [file2 ...]"
description[f] = 'Returns various stats from the scores and feature lengths.'
opts[f] = (("-t", "--format", "File format, if extension not explicit", {'default':None}),
        ("-s", "--selection", "Selection or comma-separated list of chromosome names. Selection \
         can be a dictionary (json) or a string of the form `chr:start-end`.", {'default':None}),
        ("-f", "--fields", "Fields in output", {'default':None}), )

def stats(*args,**kw):
    if len(args) < 1: raise Usage("No input file provided")
    if kw['selection']:
        if kw['selection'].count("{"):
            jsonargs = json.loads(kw['selection'])
            for k,v in jsonargs.iteritems():
                if isinstance(v,basestring): jsonargs[k] = str(v)
            kw['selection'] = dict((str(k),v) for k,v in jsonargs.iteritems())
        elif kw['selection'].count(":"):
            chr,coord = kw['selection'].split(':')
            start,end = coord.split('-')
            kw['selection'] = {'chr':chr,'start':(start,end),'end':(start,end)}
        else:
            kw['selection'] = str(kw['selection']).split(",")
    fields = None
    if kw['fields']:
        fields = str(kw['fields']).split(",")
    if kw['output'] is None:
        output = sys.stdout
    else:
        output = open(kw['output'],'w')
    chrmeta = _get_chrmeta(**kw)
    for infile in args:
        intrack = track.track(infile,format=kw['format'],chrmeta=chrmeta)
        if intrack.info:
            fileinfo = ",".join(["%s=%s" %(k,v) for k,v in intrack.info.iteritems()])
        else: fileinfo = 'None'
        chromlist = ",".join(sorted(intrack.chrmeta.keys())) or "None"
        fields = ",".join(intrack.fields)
        output.write(\
"""*****************************************
File '%s' (%s):
Infos: %s
Chromosomes: %s
Fields: %s
*****************************************\n
""" %(os.path.basename(infile), intrack.format, fileinfo, chromlist, fields))
        track.stats(intrack,out=output,**kw)
        intrack.close()
    output.close()
    return 0

############## CHECK ##############
f = 'check'
usage[f] = usage['all'] %f +" file1 [file2 ...]"
description[f] = 'Checks that a track file is sorted and well-formatted.'
opts[f] = (("-t", "--format", "File format, if extension not explicit", {'default':None}),)

def check(*args,**kw):
    if len(args) < 1: raise Usage("No input file provided")
    if kw['output'] is None: output = sys.stdout
    else: output = open(kw['output'],'w')
    chrmeta = _get_chrmeta(**kw)
    for infile in args:
        intrack = track.track(infile, format=kw['format'], chrmeta=chrmeta)
        cf = track.check_format(intrack, output)
        if cf: output.write("Check format: %s: correct %s format.\n" % (infile,intrack.format))
        co = track.check_ordered(intrack, output)
        if co: output.write("Check order: %s: file is sorted.\n" % (infile,))
        intrack.close()
    return 0

############## SORT ##############
f = 'sort'
usage[f] = usage['all'] %f +" file1 [file2 ...]"
description[f] = 'Sorts a track file. Warning: can use a lot of memory space, according to the file size.'
opts[f] = (("-t", "--format", "File format, if extension not explicit.", {'default':None}),
           ("-x", "--chromosomes", "List of chromosomes in JSON format \
                                    (to order them in a specific way).", {'default':'[]'}),)

def sort(*args,**kw):
    if len(args) < 1: raise Usage("No input file provided")
    chrmeta = _get_chrmeta(**kw)
    for infile in args:
        intrack = track.track(infile,format=kw['format'],chrmeta=chrmeta)
        outname = kw['output'] or intrack.name+'_sorted.'+intrack.format
        outtrack = track.track(outname, chrmeta=intrack.chrmeta)
        instream = intrack.read()
        s = sorted_stream(instream, chrnames=json.loads(kw['chromosomes']))
        outtrack.write(s)
        intrack.close()
    return 0


##################################
############## MAIN ##############
##################################
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
        print >>sys.stderr, '\n',err.msg,'\n'
        print >>sys.stderr, usage[fct]
        if parser: parser.print_help()
        return 1

if __name__ == '__main__':
    sys.exit(main())

