#!/usr/bin/env python

import optparse, json, sys, os
from bbcflib import btrack as track
from bbcflib import genrep

opts = (("-t", "--format", "File format if extension not explicit",{'default':None}),
        ("-s", "--selection", "Selection dictionary or coma-separated list of chromosome names",
         {'default':None}),
        ("-a", "--assembly", "Assembly name or id",{'default':None}),
        ("-f", "--fields", "Fields in output",{'default':None}),
        ("-o", "--output", "Output file (default stdout)",{'default':None}),
        ("-d", "--description", "Only print a description of the file",{'action':"store_true"})
        )

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv = None):
    try:
        usage = "readTrack.py [OPTIONS] file1 [file2 ...]"
        desc = """A generic genomic track data reader."""
        parser = optparse.OptionParser(usage=usage, description=desc)
        for opt in opts:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
        (opt, args) = parser.parse_args()
        if len(args) < 1: 
            raise Usage("No file provided")
        selection = None
        if opt.selection: 
            if opt.selection.count("{"):
                jsonargs = json.loads(opt.selection)
                for k,v in jsonargs.iteritems():
                    if isinstance(v,basestring): jsonargs[k] = str(v)
                selection = dict((str(k),v) for k,v in jsonargs.iteritems())
            else:
                selection = str(opt.selection).split(",")
        fields = None
        if opt.fields:
            fields = str(opt.fields).split(",")
        if opt.output is None: 
            output = sys.stdout
        else: 
            output = open(opt.output,'w')
        for infile in args:
            intrack = track.track(infile,format=opt.format,chrmeta=opt.assembly)
            if opt.description:
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
#                output.write(str(x)+"\n")
        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        sys.exit(1)

if __name__ == '__main__':
    sys.exit(main())

