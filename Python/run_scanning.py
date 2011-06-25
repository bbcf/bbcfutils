#!/usr/bin/env python
from bbcflib.motif              import sqlite_to_false_discovery_rate
from bein                       import MiniLIMS, unique_filename_in, ProgramOutput, program, execution
from bbcflib.genrep             import GenRep
from bbcflib.gdv                import create_gdv_project, add_gdv_track, get_project_id
from bbcflib.frontend           import parseConfig
from bbcflib.common             import ssh_add, scp, normalize_url
from bbcflib.track.format_sql   import Track, new
from os.path                    import basename, expanduser, abspath, normcase, splitext
from sys                        import argv

usage = """%s [-h] [-u via] [-m machine_host][-r remote_path] [-w website] [-p private_key] -f matrix_file -c config_file -d minilims
-h Print this message and exit
-u via Run executions using method 'via' (can be "local" or "lsf")
-m submit track file to selected host machine (e.g sugar)
-r path where put track in host machine
-w website were result will be taken for GDV (e.g http://sugar.epfl.ch)
-p path to private_key for connect to remote host (e.g ~/.ssh/user.pub)
-f dictionary where key is matrix name and value matrix file path (e.g {"Tbf1_snoRNAs":"./Tbf1_snoRNAs.mat"})
-d minilims MiniLIMS where Scanning executions and files will be stored.
-c file Config file
""" %(argv[0])

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv = None):
    genrep              = None
    assembly            = None
    M                   = None
    background          = ""
    matrix              = ""
    original_bed_data   = ""
    random_bed_data     = ""
    original_sql_data   = ""
    random_sql_data     = ""
    track_filtered      = ""
    track_scanned       = ""
    host                = ""
    website             = ""
    track_remote        = ""
    track_web_path      = ""
    via                 = ""
    limspath            = ""
    private_key         = ""
    fdr                 = 0

    try:
        try:
            opts,args = getopt.getopt   (
                                            sys.argv[1:],"hu:m:r:w:p:d:c:"  ,
                                            [
                                                "help","via","machine_host" ,
                                                "remote_path"               ,
                                                "website", "private_key"     ,
                                                "minilims","config"
                                            ]
                                        )
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.exit(0)
            elif o in ("-u", "--via"):
                if a=="local":
                    via = "local"
                elif a=="lsf":
                    via = "lsf"
                else:
                    raise Usage("Via (-u) can only be \"local\" or \"lsf\", got %s." % (a,))
            elif o in ("-w", "--website"):
                website = normalize_url(a)
            elif o in ("-d", "--minilims"):
                limspath = normcase(expanduser(a))
            elif o in ("-m", "--machine_host"):
                host = a
            elif o in ("-r", "--remote_path"):
                track_remote = normcase(expanduser(a))
            elif o in ("-p", "--private_key"):
                private_key = normcase(expanduser(a))
            elif o in ("-f", "--matrix_file"):
                matrix = a
            elif o in ("-c", "--config"):
                config_file = normcase(expanduser(a))
            else:
                raise Usage("Unhandled option: " + o)

        # compute false discovery rate
        with execution(M, description=job.description) as ex:
            # Add user identity
            ssh_add(ex, private_key)

            job, config         = parseConfig(normcase(expanduser(config_file)))
            genrep              = GenRep(config=config)
            assembly            = genrep.assembly(job["assembly_id"])
            background          = genrep.statistics(assembly,output=unique_filename_in)
            M                   = MiniLIMS(limspath)
            if len(job.groups) >2:
                raise ValueError("They are more than 2 group in config file")
            for group in job.groups:
                url = job.groups[group]["url"]
                uri = ""
                if url.startswith("http") or url.startswith("www."):
                    url = normalize_url(url)
                    # download data
                    data    = urllib2.urlopen(url)
                    uri     = unique_filename_in()
                    with open(original_bed_data, "w") as f:
                        f.write(data.read())
                else:
                   uri = normcase(expanduser(uri))
                if job.groups[group]["control"] is True:
                    original_bed_data = uri
                else:
                    random_bed_data = uri

            original_sql_data   = splitext(original_bed_data)[0]+".db"
            random_sql_data     = splitext(random_bed_data)[0]+".db"
            track_filtered      = splitext(original_bed_data)[0]+"_filtered.db"

            # convert to sql
            with Track(original_bed_data, chrmeta=assembly.chromosomes) as track:
                track.convert(original_sql_data, format='sql')
            with Track(random_bed_data, chrmeta=assembly.chromosomes) as track:
                track.convert(random_sql_data, format='sql')

            track_scanned,fdr = sqlite_to_false_discovery_rate(
                                                                ex,
                                                                matrix,
                                                                background,
                                                                genrep,
                                                                assembly.chromosomes,
                                                                original_sql_data,
                                                                random_sql_data,
                                                                threshold=0,
                                                                via=via,
                                                                keep_max_only=True,
                                                                alpha=0.05,
                                                                nb_false_positive_hypotesis=5.0
                                                              )

            # filter track with fdr as treshold
            with new(track_filtered, format="sql", datatype="qualitative") as track_out:
                track_out.meta_track= {"source": basename(original_bed_data)}
                track_out.meta_track.update({"k":"v"})
                with Track(track_scanned, format="sql", chrmeta=assembly.chromosomes) as track_in:
                    for chromosome in track_in.all_chrs:
                        data_list = []
                        for data in track_in.read( {"chr": chromosome, "score": ">= "+str(fdr)}, fields=Track.qualitative_fields ):
                            data_list.append(data)
                        if len(data_list) > 0:
                            track_out.write(chromosome, data_list)
                    #track_out.meta_chr = list(chomosomes_used)

            # TODO if [-m machine_host][-r remote_path] [-w website] [-p private_key]
            # send new track to remote
            scp(ex, track_filtered, track_remote, "jmercier", host)

            # create gdv project
            json        = create_gdv_project(
                                                config["gdv"]["key"], config["gdv"]["email"],
                                                "jonathan", None,
                                                assembly.nr_assembly_id,
                                                config["gdv"]["url"],
                                                public = True
                                            )
            project_id  = get_project_id( json )
            add_gdv_track  (
                                config["gdv"]["key"], config["gdv"]["email"],
                                project_id, track_web_path,
                                gdv_url=config["gdv"]["url"]
                            )
    except Usage, err:
        print >>sys.stderr, err.msg
        print >>sys.stderr, usage
        return 2
if __name__ == '__main__':
    sys.exit(main())
