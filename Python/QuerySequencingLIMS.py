#!/usr/bin/env python
"""
Stand-alone script for DAFLIMS queries
"""

from bbcflib import daflims
from bbcflib.workflows import Usage
from bein import MiniLIMS
from bein.util import use_pickle
from urllib2 import urlopen, URLError
import optparse, sys, json

_usage = __file__+" facility machine run lane [libname] [-o filename] [-m minilims_path]"
_desc = "Retrieve data from sequencing facility lims. Libname is mandatory if multiplexing was used."
_lims = "/data/epfl/bbcf/htsstation/data/mapseq_minilims"

def query(*args,**kwargs):
    facility = args[0]
    params = use_pickle(MiniLIMS(kwargs['minilims']), "global variables")['lims']
    dafl = dict((loc,daflims.DAFLIMS( username=params['user'], password=pwd ))
                for loc,pwd in params['passwd'].iteritems())

    try: 
        info = dafl[facility]._lanedesc(*args)
        print json.dumps(info)
    except ValueError as e:
        if len(args) > 4: 
            msg = "Access not granted to data [%s: %s-%i/%i/%s]" %args
        else:
            msg = "Access not granted to data [%s: %s-%i/%i]" %args
        print msg
        return 1
    try:
        if len(args) > 4:
            libname = args[4]
        else:
            libname = None
        url = dafl[facility]._symlinkname(*args[:4], libname=libname)
        urlok = []
        for k in sorted(url.keys()):
            _ = urlopen(url[k])
            urlok.append(url[k])
        print "\n".join(urlok)
    except (ValueError, URLError) as e:
        print "File not accessible", e
        return 2
    if kwargs.get('output'):
        try:
            filename = dafl[facility]._fetch_symlink(url, kwargs['output'])
            print filename
        except Exception as e:
            print "Problem with downloading %s: %s" %(url,e)
            return 3
    return 0


if __name__ == '__main__':
    try: 
        parser = optparse.OptionParser(usage=_usage, description=_desc)
        parser.add_option("-o","--output",help="download data to file",
                          default=None)
        parser.add_option("-m","--minilims",help="path to an HTSstation minilims",
                          default=_lims)
        (opt, args) = parser.parse_args()
        if len(args) < 4:
            raise Usage("Need at least 'facility', 'machine', 'run', 'lane' informations")
        args[2] = int(args[2])
        args[3] = int(args[3])
        sys.exit(query(*args[:5],minilims=opt.minilims,output=opt.output))
    except Usage as err:
        print >>sys.stderr, '\n',err.msg,'\n'
        parser.print_help()
        sys.exit(10)
