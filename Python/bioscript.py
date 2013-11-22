#!/usr/bin/env python

from bsPlugins import PLUGINS_FILES
import sys, optparse

_usage = "bioscript.py %s [PLUGIN OPTIONS]"
_description = "Command line interface to the bsPlugins library"

def parse_plugin(plugin,name,sub=None):
    usage = _usage%name
    fname = name+"Plugin"
    if not hasattr(plugin,fname):
        if sub is not None:
            fname = sub+"Plugin"
            usage = _usage%(name+" "+sub)
        else:
            allsubs = []
            for k,v in plugin.__dict__.iteritems():
                if getattr(v, 'bs_plugin', '') == 'bs-operation':
                    allsubs.append(k[:-6])
            if len(allsubs) == 1:
                fname = allsubs[0]
            else:
                print "Choose one operation in plugin: ",allsubs
                sys.exit(1)
    desc = getattr(plugin,fname).info.get("description", _description)
    opts = []
    multis = {}
    for param in getattr(plugin,fname).info["in"]:
        addmore = {}
        if param['type'] == 'boolean':
            addmore['action'] = "store_true"
        if param.get('multiple'):
            multis[param['id']] = param['multiple']
        opts.append(("--%s"%param['id'], param['type'], addmore))
    return usage, opts, desc, fname, multis

def main():
    sub = None
    if len(sys.argv) > 1:
        name = sys.argv.pop(1)
        if len(sys.argv) > 1 and not sys.argv[1][0] == '-':
            sub = sys.argv.pop(1)
    else:
        print 'Available plugins are: %s' %", ".join(PLUGINS_FILES)
        return 0
    __import__("bsPlugins."+name)
    plugin = getattr(sys.modules['bsPlugins'],name)
    usage, opts, desc, funcname, multis = parse_plugin(plugin,name,sub)
    parser = optparse.OptionParser(usage=usage, description=desc)
    for opt in opts:
        if len(opt) == 4:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
        elif len(opt) == 3:
            parser.add_option(opt[0],help=opt[1],**opt[2])
    if len(sys.argv) == 1:
        parser.print_help()
        return 0
    (opt, args) = parser.parse_args()
    pl_call = getattr(sys.modules["bsPlugins."+name],funcname)()
    for k, v in multis.iteritems():
        if not v in opt.__dict__: opt.__dict__[v] = {}
        vv = None
        if k in opt.__dict__: vv = opt.__dict__.pop(k)        
        opt.__dict__[v][k] = vv.split(",") if vv else []
    try:
        pl_call(**opt.__dict__)
    except Exception, err:
        print >>sys.stderr, '\n',err,'\n\n'
        parser.print_help()
        return 1

if __name__ == '__main__':
    sys.exit(main())
