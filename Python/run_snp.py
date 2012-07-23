#!/usr/bin/env python

"""
SNP detection workflow.
"""
from bein import execution, MiniLIMS
from bein.util import use_pickle, add_pickle, pause
from bbcflib import daflims, genrep, frontend, gdv, mapseq, common, snp
from bbcflib import email
import sys, os, json, optparse, shutil


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
                                     {'default':None}),)
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
        groups = job.groups

        if opt.fasta_path and os.path.exists(opt.fasta_path):
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

            pileup_dict = dict((chrom,{}) for chrom in genomeRef.keys()) # {chr: {}}
            sample_names = []
            for gid, files in bam_files.iteritems():
                sample_name = groups[gid]['name']
                sample_names.append(sample_name)
                runs = [r['bam'] for r in files.itervalues()]
                bam = mapseq.merge_bam(ex,runs)
                pileupFilename = common.unique_filename_in()
                # Samtools pileup
                for chrom,ref in genomeRef.iteritems():
                    future = snp.sam_pileup.nonblocking( ex, assembly, bam, ref,
                                                         via=opt.via, stdout=pileupFilename )
                    pileup_dict[chrom][pileupFilename] = (future,sample_name) # {chr: {filename: (future,name)}}
            chr_filename = {}
            for chrom, dictPileup in pileup_dict.iteritems():
                # Get the results from sam_pileup
                # Write the list of all snps of THIS chromosome, from ALL samples
                allSNPpos,parameters = snp.posAllUniqSNP(dictPileup) # {3021:'A'}, (minCoverage, minSNP)
                if len(allSNPpos) == 0: continue
                # Write results in a temporary file, for this chromosome
                chr_filename[chrom] = snp.write_pileupFile(
                    dictPileup, sample_names, allSNPpos, chrom,
                    minCoverage = parameters[0],
                    minSNP = parameters[1])

            #shutil.copy(chr_filename['chr5'], '../../'+'chr5')
            #shutil.copy(chr_filename['chrV'], '../../'+'yeast_chrV')

            # Add exon & codon information & write the real file
            outall,outexons = snp.annotate_snps(chr_filename,sample_names,assembly,genomeRef)
            description = common.set_file_descr("allSNP.txt",step="SNPs",type="txt")
            ex.add(outall,description=description)
            description = common.set_file_descr("exonsSNP.txt",step="SNPs",type="txt")
            ex.add(outexons,description=description)

            # Create tracks for UCSC and GDV:
            snp.create_tracks(ex,outall,sample_names,assembly)

        allfiles = common.get_files(ex.id, M)
        logfile.close()
        debugfile.close()
        print json.dumps(allfiles)
        with open(hts_key+".done",'w') as done:
            json.dump(allfiles,done)

        if 'email' in gl:
            r = email.EmailReport( sender=gl['email']['sender'],
                                   to=str(job.email).split(','),
                                   subject="snp job "+str(job.description),
                                   smtp_server=gl['email']['smtp'] )
            r.appendBody('''
Your SNP job is finished.

The description was:
'''+str(job.description)+'''
and its unique key is '''+hts_key+'''.

You can retrieve the results at this url:
'''+gl['hts_snp']['url']+"jobs/"+hts_key+"/get_results")
            r.send()
        return 0
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2

if __name__ == '__main__':
    sys.exit(main())


