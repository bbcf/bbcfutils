#!/usr/bin/env python

import optparse, sys, os, shutil
from bein import MiniLIMS
from bbcflib.common import get_files

#hts_minilims.py -m mapseq -e vwEfIK6vG8iD64kZm5m6 -l -t type=pdf


module_list = ["demulitplexing","mapseq","chipseq","rnaseq","snp","4cseq"]
usage = "hts_minilims.py [OPTIONS]"
description = "Command-line interface to HTSstation minilims: list, copy or link selected files from given execution."

opts = (("-m", "--minilims", "path to personal minilims, or name of an HTSstation module",{'default':''}),
        ("-e", "--execution", "name or number of the execution in the minilims",{'default':0}),
        ("-t", "--tag", "select files having this tag",{'default':None}),
        ("-o", "--output", "local path",{'default':None}),
        ("-l", "--list", "only list files",{'action':"store_true",'default':False}),
        ("-c", "--copy", "copy files to local path", {'action':"store_true",'default':False}),
        ("-s", "--symlink", "create symlinks to files in local path", {'action':"store_true",'default':False}),
        ("", "--basepath","HTS basepath",{'default':"/archive/epfl/bbcf/data/htsstation"}))

class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

def main(argv=None):
    parser = None
    try:
        parser = optparse.OptionParser(usage=usage, description=description)
        for opt in opts:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
        (options, args) = parser.parse_args()
        if options.minilims in module_list:
            M = MiniLIMS(os.path.join(options.basepath,options.minilims+"_minilims"))
        elif os.path.exists(options.minilims):
            M = MiniLIMS(options.minilims)
        else:
            raise Usage("Minilims not found, please specify a path or a module with -m.")
        tags = options.tag
        if tags:
            tags = tags.split(",")
            if tags[0].count("="):
                tags = dict(x.split("=") for x in tags)
            elif tags[0].count(":"):
                tags = dict(x.split(":") for x in tags)
        files = get_files(options.execution,M,select_param=tags)
        fprefix = ''
        if options.list:
            if options.output and os.path.isdir(options.output):
                options.output = os.path.join(options.output,options.execution+".txt")
            outfile = options.output and open(options.output,"w") or sys.stdout
            outfile.write("\t".join(["type","group","name","path","comment"])+"\n")
        else:
            if not(options.output): options.output = "./"
            if not(os.path.isdir(options.output)):
                options.output, fprefix = os.path.split(options.output)
        for t in sorted(files.keys()):
            for k,v in files[t].iteritems():
                fpath = os.path.join(M.file_path,k)
                vv = v.split("[")
                fname = fprefix+vv.pop(0)
                comment = ''
                par_dict = {}
                if vv:
                    vv = vv[0].split("]")
                    par_dict = dict(x.split(":") for x in vv.pop(0).split(","))
                    if vv: comment = vv[0].strip().strip("()")
                if par_dict.get('view') == 'admin': continue
                if options.list: outfile.write("\t".join([t,par_dict.get('groupId',''),fname,fpath,comment])+"\n")
                if options.copy: shutil.copy(fpath, os.path.join(options.output,fname))
                if options.symlink: os.symlink(fpath, os.path.join(options.output,fname))
        if options.list and options.output: outfile.close()
        return 0

    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        if parser: parser.print_help()
        return 1

if __name__ == '__main__':
    sys.exit(main())


