#!/usr/bin/env python

"""

SNP detection workflow.

run_snp.py -v local -f /Users/julien/Workspace/genomes/sacCer2.tar.gz -c config/snp.config -d snp_minilims

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
    opts = (("-v", "--via", "Run executions locally or using bsub (can be either 'local' or 'lsf')", {'default': "lsf"}),
            ("-k", "--key", "Alphanumeric key of the new SNP job", {'default': None}),
            ("-m", "--mapseq_limspath", "MiniLIMS where a previous Mapseq execution and files has been stored.",
                                     {'default': "/srv/mapseq/public/data/mapseq_minilims"}),
            ("-w", "--working-directory", "Create execution working directories in wdir",
                                     {'default': os.getcwd(), 'dest':"wdir"}),
            ("-c", "--config", "Config file", {'default': None}),
            ("-d", "--snp_limspath", "MiniLIMS where SNP executions and files will be stored.", \
                                     {'default': "/srv/snp/public/data/snp_minilims"}),
            ("-f", "--fasta_path", "Path to a directory containing a fasta file for each chromosome",
                                     {'default':''}),)
    try:
        usage = "run_snp.py [OPTIONS]"
        desc = """Compares sequencing data to a reference assembly to detect SNP."""
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

        if os.path.exists(opt.fasta_path):
            path_to_ref = opt.fasta_path
        elif os.path.exists(assembly.fasta_path()):
            path_to_ref = assembly.fasta_path()

        logfile = open(hts_key+".log",'w')
        debugfile = open(hts_key+".debug",'w')
        debugfile.write(json.dumps(job.options)+"\n\n"+json.dumps(gl)+"\n");debugfile.flush()

        # Program body
        with execution( M, description=hts_key, remote_working_directory=opt.wdir ) as ex:
            print "Current working directory:", ex.working_directory
            (bam_files, job) = mapseq.get_bam_wig_files(ex, job, minilims=opt.mapseq_limspath,
                                                        hts_url=mapseq_url,
                                                        script_path=gl.get('script_path',''),
                                                        via=opt.via)
            assert bam_files, "Bam files not found."
            logfile.write("cat genome fasta files\n");logfile.flush()
            genomeRef = snp.untar_genome_fasta(assembly, path_to_ref, convert=True)
            logfile.write("done\n");logfile.flush()

            # Samtools pileup
            samples = dict((chrom,{}) for chrom in genomeRef.keys()) 
            pileup_dict = dict((chrom,{}) for chrom in genomeRef.keys()) # {'chr.': {}}
            for idGroup,dictRuns in job.groups.iteritems():
                nbRuns = len(dictRuns["runs"].keys())
                listFormattedFile=[]
                for idRun,dictRun in bam_files[idGroup].iteritems():
                    sampleName = job.groups[idGroup]['name']
                    if(nbRuns>1):
                        sampleName += "_"+dictRun.get('libname',str(idRun))
                        debugfile.write("(many runs, need statistics)\n");debugfile.flush()
                    for chrom,ref in genomeRef.iteritems():
                        pileupFilename = common.unique_filename_in()
                        future = snp.sam_pileup.nonblocking( ex, assembly, dictRun["bam"], ref,
                                                             via=opt.via, stdout=pileupFilename )
                        pileup_dict[chrom][pileupFilename] = future # {chr: {filename: future}}
                        samples[chrom][pileupFilename] = sampleName
            sample_names = [samples[pileup_dict.iterkeys().next()][p] \
                            for p in pileup_dict.itervalues().next().iterkeys()]

            # Get & write results
            formattedPileupFilename = {}
            for chrom, dictPileup in pileup_dict.iteritems():
                allSNPpos,parameters = snp.posAllUniqSNP(dictPileup) # {3021:'A'}, (minCoverage, minSNP)
                if len(allSNPpos) == 0: continue
                formattedPileupFilename[chrom] = snp.parse_pileupFile(
                    samples[chrom], allSNPpos, chrom,
                    minCoverage = parameters[0],
                    minSNP = parameters[1])
                #with open(formattedPileupFilename[chrom]) as f:
                #    print f.read()

            # Add exon & codon information
#            output = common.cat(formattedPileupFilename[0],formattedPileupFilename[1:],skip=1)
            outall,outexons = snp.annotate_snps(formattedPileupFilename,sample_names,assembly)
            description = common.set_file_descr("allSNP.txt",step="SNPs",type="txt")
            ex.add(outall,description=description)
            description = common.set_file_descr("exonsSNP.txt",step="SNPs",type="txt")
            ex.add(outexons,description=description)

            # Codon analysis
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


