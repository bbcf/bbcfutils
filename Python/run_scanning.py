#!/usr/bin/env python
from bbcflib.motif              import sqlite_to_false_discovery_rate
from bein                       import MiniLIMS, unique_filename_in, ProgramOutput, program, execution
from bbcflib.genrep             import GenRep
from bbcflib.gdv                import create_gdv_project, add_gdv_track, get_project_id
from bbcflib.frontend           import parseConfig
from bbcflib.common             import scp, normalize_url
from bbcflib.track              import Track, new
from os.path                    import basename, expanduser, abspath, normcase, splitext, isfile, exists
import getopt, sys

usage = """%s [-h] [options] --matrix=/path/to/matrix.mat -c config_file --minilims=minilims_name
-h Print this message and exit
-c --config file Config file
--host submit track file to selected host machine (e.g sugar)
--identity_file public ssh identity file (e.g ~/.ssh/id.rsa)
--matrix file path to matrix
--minilims MiniLIMS where Scanning executions and files will be stored.
--project GDV project name
--remote_path path where put track in host machine
--via via Run executions using method 'via' (can be "local" or "lsf")
--username ssh username for login to host
--website website were result will be taken for GDV (e.g http://sugar.epfl.ch)
""" %(sys.argv[0])

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main(argv = None):
    genrep              = None
    assembly            = None
    M                   = None
    job                 = None
    config              = None
    config_file         = None
    background          = ""
    matrix              = ""
    original_bed_data   = ""
    random_bed_data     = ""
    original_sql_data   = ""
    random_sql_data     = ""
    track_filtered      = ""
    track_scanned       = ""
    project             = ""
    username            = ""
    identity_file       = ""
    host                = ""
    website             = ""
    remote_path         = ""
    track_web_path      = ""
    via                 = ""
    limspath            = ""
    fdr                 = 0
    if argv is None: argv = sys.argv
    try:
        try:
            opts,args = getopt.getopt   (
                                            argv[1:],"hu:c:"  ,
                                            [
                                                "help", "via=", "host="     ,
                                                "remote_path=" , "website=" ,
                                                "minilims=","config="       ,
                                                "matrix=", "username="      ,
                                                "identity_file=", "project="
                                            ]
                                        )
        except getopt.error, msg:
            raise Usage(msg)
        for o, a in opts:
            if o in ("-h", "--help"):
                print __doc__
                print usage
                sys.sys.exit(0)
            elif o == "--via":
                if a=="local":
                    via = "local"
                elif a=="lsf":
                    via = "lsf"
                else:
                    raise Usage("Via (-u) can only be \"local\" or \"lsf\", got %s." % (a,))
            elif o == "--website":
                website = normalize_url(a)
            elif o == "--minilims":
                limspath = normcase(expanduser(a))
            elif o == "--host":
                host = a
            elif o == "--identity_file":
                identity_file = a
            elif o == "--remote_path":
                remote_path = normcase(expanduser(a))
            elif o == "--matrix":
                matrix = {basename(a):normcase(expanduser(a))}
            elif o == "--username":
                username = a
            elif o == "--project":
                project = a
            elif o in ("-c", "--config"):
                config_file = normcase(expanduser(a))
            else:
                raise Usage("Unhandled option: " + o)

        # read config file
        if config_file is None or not exists(config_file) or not isfile(config_file):
            raise Usage("Config file missing")
        else:
            job, config = parseConfig(normcase(expanduser(config_file)))

        genrep              = GenRep(config=config)
        assembly            = genrep.assembly(job.assembly_id)
        M                   = MiniLIMS(limspath)

        # compute false discovery rate
        with execution(M, description=job.description) as ex:
            background = genrep.statistics(assembly,output=unique_filename_in())
            ex.add(background,  description="background:"+background)
            if len(job.groups) >2:
                raise ValueError("They are more than 2 group in config file")
            for group in job.groups:
                if "url" in job.groups[group]["runs"][group]:
                    url = job.groups[group]["runs"][group]["url"]
                    uri = ""
                    if url.startswith("http") or url.startswith("www."):
                        url = normalize_url(url)
                        # download data
                        data    = urllib2.urlopen(url)
                        uri     = unique_filename_in()
                        with open(original_bed_data, "w") as f:
                            f.write(data.read())
                    else:
                       uri = normcase(expanduser(url))
                    if job.groups[group]["control"] is True:
                        original_bed_data = uri
                    else:
                        random_bed_data = uri
                else:
                    random_bed_data = ""

            original_sql_data   = unique_filename_in()
            random_sql_data     = unique_filename_in()
            track_filtered      = unique_filename_in()

            ex.add(original_bed_data,   "bed:"+original_bed_data)
            # convert to sql
            with Track(original_bed_data, chrmeta=assembly.chromosomes) as track:
                track.convert(original_sql_data, format='sql')
                # create random track
                track.shuffle_track(random_sql_data, repeat_number=5)
            ex.add(original_sql_data,   "sql:"+original_sql_data)
            ex.add(random_sql_data,     "sql:"+random_sql_data)
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
            ex.add(track_filtered,      "sql:"+track_filtered)

            # send new track to remote
            if host != "" and remote_path != "" and username != "":
                args = []
                if identity_file != "":
                    args = ["-i " + abspath(expanduser(identity_file)) ]
                source      = abspath(expanduser(track_filtered))
                destination = "%s@%s:%s" %(username, host, remote_path)
                scp(ex, source, destination, args=args)

        # create gdv project
        json        = create_gdv_project(
                                            config["gdv"]["key"], config["gdv"]["email"],
                                            project, None,
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
