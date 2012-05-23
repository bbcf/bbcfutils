#!/usr/bin/env python

"""

SNP detection workflow.

"""
from bein import execution, MiniLIMS
from bein.util import use_pickle, add_pickle, pause
from bbcflib import daflims, genrep, frontend, gdv, mapseq, common, snp
from bbcflib import email
import sys, os, json, optparse


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv = None):
    opts = (("-v", "--via", "Run executions using method 'via' (can be 'local' or 'lsf')", {'default': "lsf"}),
            ("-k", "--key", "Alphanumeric key of the new RNA-seq job", {'default': None}),
            ("-m", "--mapseq_limspath", "MiniLIMS where a previous Mapseq execution and files has been stored.", 
                                     {'default': "/data/htsstation/mapseq/mapseq_minilims"}),
            ("-w", "--working-directory", "Create execution working directories in wdir",
                                     {'default': os.getcwd(), 'dest':"wdir"}),
            ("-c", "--config", "Config file", {'default': None}),
            ("-d", "--snp_limspath", "MiniLIMS where snp executions and files will be stored.", \
                                     {'default': "/srv/snp/public/data/snp_minilims"}))
    try:
        usage = "run_snp.py [OPTIONS]"
        desc = """........."""
        parser = optparse.OptionParser(usage=usage, description=desc)
        for opt in opts:
            parser.add_option(opt[0],opt[1],help=opt[2],**opt[3])
        (opt, args) = parser.parse_args()

        if os.path.exists(opt.wdir): os.chdir(opt.wdir)
        else: parser.error("Working directory '%s' does not exist." % opt.wdir)
        if not(os.path.exists(opt.snp_limspath)
               and (opt.key != None or (opt.config and os.path.exists(opt.config)))):
            raise Usage("Need a minilims and a job key or a configuration file")
        M = MiniLIMS( opt.snp_limspath )
        if opt.key:
            hts_key = opt.key
            gl = use_pickle( M, "global variables" )
            htss = frontend.Frontend( url=gl['hts_snp']['url'] )
            job = htss.job( hts_key )
            [M.delete_execution(x) for x in M.search_executions(with_description=hts_key,fails=True)]
        elif os.path.exists(opt.config):
            (job,gl) = frontend.parseConfig(opt.config)
            hts_key = job.description
        else:
            raise ValueError("Need either a job key (-k) or a configuration file (-c).")
        mapseq_url = gl.get('hts_mapseq',{}).get('url')
        job.options['ucsc_bigwig'] = True
        if not('create_gdv_project' in job.options):
            job.options['create_gdv_project'] = False
        elif isinstance(job.options['create_gdv_project'],str):
            job.options['create_gdv_project'] = job.options['create_gdv_project'].lower() in ['1','true','t']
        g_rep = genrep.GenRep( gl.get("genrep_url"), gl.get("bwt_root") )
        assembly = genrep.Assembly( assembly=job.assembly_id, genrep=g_rep )
        logfile = open(hts_key+".log",'w')
        debugfile = open(hts_key+".debug",'w')
        debugfile.write(json.dumps(job.options)+"\n\n"+json.dumps(gl)+"\n");debugfile.flush()

        # Program body
        with execution( M, description=hts_key, remote_working_directory=opt.wdir ) as ex:
            (bam_files, job) = mapseq.get_bam_wig_files(ex, job, minilims=opt.mapseq_limspath, 
                                                        hts_url=mapseq_url,
                                                        script_path=gl.get('script_path',''), 
                                                        via=opt.via)
            assert bam_files, "Bam files not found."
            logfile.write("cat genome fasta files\n");logfile.flush()
            genomeRef = snp.untar_genome_fasta(assembly,convert=True)
            logfile.write("done\n");logfile.flush()
            dictPileupFile=dict((chrom,{}) for chrom in genomeRef.keys())
            for idGroup,dictRuns in job.groups.iteritems():
                nbRuns=len(dictRuns["runs"].keys())
                listFormatedFile=[]
                for idRun,dictRun in bam_files[idGroup].iteritems():
                    sampleName=job.groups[idGroup]['name']
                    if(nbRuns>1):
                        sampleName += "_"+dictRun.get('libname',str(idRun))
                        debugfile.write("many runs, need statistics\n");debugfile.flush() 
                    for chrom,ref in genomeRef.iteritems():
                        pileupFilename=common.unique_filename_in()
                        future = snp.sam_pileup.nonblocking( ex, assembly, dictRun["bam"], ref,
                                                             via=opt.via, stdout=pileupFilename )
                        dictPileupFile[chrom][pileupFilename] = (sampleName,future)
            formatedPileupFilename = {}
            for chrom, dictPileup in dictPileupFile.iteritems():
                posAll,parameters = snp.posAllUniqSNP(dictPileup)
                if len(posAll) == 0: continue
                formatedPileupFilename[chrom] = snp.parse_pileupFile(
                    dictPileup, posAll, chrom,
                    minCoverage=parameters[0],
                    minSNP=parameters[1])
            sample_names = dictPileupFile.values()[0].values()
#            output = common.cat(formatedPileupFilename[0],formatedPileupFilename[1:],skip=1)
            output = snp.annotate_snps(formatedPileupFilename,sample_names,assembly)
            description = common.set_file_descr("allSNP.txt",step="SNPs",type="txt")
            ex.add(output[0],description=description)
            description = common.set_file_descr("exonsSNP.txt",step="SNPs",type="txt")
            ex.add(output[1],description=description)
#            codon = snp.synonymous(job,output)
#            ex.add(codon,description=set_file_descr("functionalVariants.txt",
#                                                    step="codon_modification",
#                                                    type="txt"))
        allfiles = common.get_files(ex.id, M)
        logfile.close()
        debugfile.close()
        print json.dumps(allfiles)
        with open(hts_key+".done",'w') as done:
            json.dump(allfiles,done)

        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())


