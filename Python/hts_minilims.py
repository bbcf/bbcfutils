#!/usr/bin/env python

import optparse, sys, os, shutil, re
from bein import MiniLIMS
from bbcflib.common import get_files

#hts_minilims.py -m mapseq -e vwEfIK6vG8iD64kZm5m6 -l -t type=pdf


import textwrap
class TextWrapper:
    @staticmethod
    def wrap(text, width=70, **kw):
        if text.count("\n"):
            return [textwrap.wrap(l, width, **kw) for l in text.split("\n")]
        else:
            return textwrap.wrap(text, width, **kw)

    @staticmethod
    def fill(text, width=70, **kw):
        return "\n".join([textwrap.fill(l, width, **kw) for l in text.split("\n")])

optparse.textwrap = TextWrapper()

module_list = ["demultiplexing","mapseq","chipseq","rnaseq","snp","4cseq"]
usage = "hts_minilims.py [OPTIONS]"
description = """Command-line interface to HTSstation minilims: list, copy or link selected files from given execution.
Typical example:
    'hts_minilims.py -m mapseq -e afP4AJAe5EcvEop1EH48 -l -t type=bam'
"""

opts = (("-m", "--minilims", "path to personal minilims, or name of an HTSstation module",{'default':''}),
        ("-e", "--execution", "name or number of the execution in the minilims",{'default':0}),
        ("-t", "--tag", "select files having this tag",{'default':None}),
        ("-o", "--output", "local path",{'default':None}),
        ("-l", "--list", "only list files",{'action':"store_true",'default':False}),
        ("-c", "--copy", "copy files to local path", {'action':"store_true",'default':False}),
        ("-s", "--symlink", "create symlinks to files in local path", {'action':"store_true",'default':False}),
        ("-p", "--programs", "list execution's program arguments and outputs", {'action':"store_true",'default':False}),
        ("-g", "--gdv", "key of a gdv project to send sql files to", {'default':None}),
        ("", "--basepath","HTS basepath",{'default':"/archive/epfl/bbcf/data/htsstation"}),
        ("", "--gdvurl","GDV base url",{'default':None}),
        ("", "--email","GDV user email",{'default':""}),
        ("", "--key","GDV user key",{'default':""}))

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
        if options.gdv:
            if tags: tags['type'] = 'sql'
            else: tags = {'type': 'sql'}
        if options.programs:
            if isinstance(options.execution, basestring):
                exlist = M.search_executions(with_description=options.execution)
                if len(exlist) == 0:
                    exlist = M.search_executions(with_text=options.execution)
                if len(exlist) == 0:
                    raise Usage("Execution with key %s not found in %s." %(options.execution,options.minilims))
                exid = max(exlist)
            else:
                exid = int(options.execution or 0)
            exec_data = M.fetch_execution(exid)['programs']
            outfile = options.output and open(options.output,"w") or sys.stdout
            for prog in exec_data:
                pargs = prog['arguments']
                if tags and all([t not in x for x in pargs for t in tags]):
                    continue
                stout = prog['stdout']
                sterr = prog['stderr']
                if pargs[0] == 'bsub': pargs = str(pargs[-1])
                else: pargs = str(" ".join(pargs))
                outfile.write("\n".join([pargs,stout,'',sterr]))
                outfile.write("\n"+''.join(['-']*40)+"\n")
            outfile.close()
            return 0
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
        if options.gdv:
            gdvpaths = []
            gdvnames = []
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
                if options.gdv:
                    gdvpaths.append(fpath)
                    gdvnames.append(re.sub('\.sql.*','',str(fname)))
        if options.list and options.output: outfile.close()
        if options.gdv:
            from bbcflib import gdv
            gdvurl = options.gdvurl or gdv.default_url
            gdvproject = gdv.get_project(mail=options.email, key=options.key,
                                         project_key=options.gdv)
            if gdvproject.get('project',{}).get('id',0)>0:
                try:
                    tr = gdv.multiple_tracks( mail=options.email, key=options.key,
                                              project_id=gdvproject['project']['id'],
                                              urls=gdvpaths, names=gdvnames, extensions=['sql']*len(gdvpaths),
                                              force=True, serv_url=gdvurl )
                except Exception, err:
                    raise Usage("GDV Tracks Failed: %s\n" %err)
                print """
*********** GDV project at: ***********
%s/public/project?k=%s&id=%s
***************************************
""" %(gdvurl,gdvproject['project']['download_key'],gdvproject['project']['id'])
        return 0

    except Usage, err:
        print >>sys.stderr, "\n"+err.msg+"\n"
        print >>sys.stderr, usage
        if parser: parser.print_help()
        return 1

if __name__ == '__main__':
    sys.exit(main())


