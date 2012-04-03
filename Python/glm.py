#!/usr/bin/env python
"""
...
"""
import os, sys, json
import optparse
from bbcflib import frontend
from bbcflib.common import set_file_descr, unique_filename_in
from bein import execution, MiniLIMS, program
from bein.util import use_pickle


class Usage(Exception):
    def __init__(self,  msg):
        self.msg = msg

@program
def run_glm(data_file, options=[]):
    output_file = unique_filename_in()
    options += ["-o",output_file]
    return {'arguments': ["R","--slave","-f","negbin.test.R","--args",data_file]+options,
            'return_value': output_file}

def main():
    opts = (("-v", "--via", "Run executions using method 'via' (can be 'local' or 'lsf')", {'default': "lsf"}),
            ("-k", "--key", "Alphanumeric key of the new RNA-seq job", {'default': None}),
            ("-d", "--rnaseq_minilims", "MiniLIMS where RNAseq executions and files will be stored.",
                                     {'default': "/srv/rnaseq/public/data/rnaseq_minilims"}),
            ("-w", "--wdir", "Create execution working directories in wdir", {'default': os.getcwd()}),
            ("-c", "--config", "Config file", {'default': None}),
            ("-s", "--sep", "Character separating fields in the input files. Use 'tab' for tab delimiter (not '\t').",
                                     {'default': "tab"}),
            ("--design", "name of the file containing the design matrix (see below).", {'default': None}),
            ("--contrast", "name of the file containing the contrast matrix (see below).", {'default': None}),
           )
    try:
        usage = "run_glm.py data_file [OPTIONS]"
        desc = """  """
        parser = optparse.OptionParser(usage=usage, description=desc)
        for opt in opts:
            if len(opt) == 4:
                parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
            elif len(opt) == 3:
                parser.add_option(opt[0],help=opt[1],**opt[2])
        (opt, args) = parser.parse_args()

        # Treat arguments and options #
        if len(args) == 1: data_file = args[0]
        else: parser.error("glm.py takes only one argument `data_file`.")
        if os.path.exists(opt.wdir): os.chdir(opt.wdir)
        else: parser.error("Working directory '%s' does not exist." % opt.wdir)
        if not opt.rnaseq_minilims: parser.error("Must specify a MiniLIMS to attach to")
        options = []
        if opt.design: options += ["-d",os.path.abspath(opt.design)]
        if opt.contrast: options += ["-c",os.path.abspath(opt.contrast)]
        options += ["-s",opt.sep]

        # Job, MiniLIMS and execution #
        M = MiniLIMS(opt.rnaseq_minilims)
        if opt.key:
            gl = use_pickle( M, "global variables" )
            htss = frontend.Frontend( url=gl['hts_rnaseq']['url'] )
            job = htss.job(opt.key)
            [M.delete_execution(x) for x in M.search_executions(with_description=opt.key,fails=True)]
            description = opt.key
        elif opt.config and os.path.exists(opt.config):
            (job,gl) = frontend.parseConfig(opt.config)
            description = opt.config
        else:
            description = "No description"

        # Program body #
        logfile = open(description+".log",'w')
        with execution(M, description=description, remote_working_directory=opt.wdir) as ex:
            logfile.write("Enter execution. Current working directory: %s \n" % ex.working_directory);logfile.flush()
            print "Current working directory:", ex.working_directory
            output_file = run_glm(ex, data_file, options)
            desc = set_file_descr(output_file,{'step':'2', 'type':'txt',
                                  'comment':'Differential analysis for job %s' % description})
            ex.add(output_file, description=desc)

        logfile.close()
        print json.dumps({'differential_analysis': output_file})

        sys.exit(0)
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2


if __name__ == '__main__':
    sys.exit(main())







