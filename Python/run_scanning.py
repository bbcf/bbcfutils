#!/usr/bin/env python

"""
run_scanning is a script for get sequence over represented on  whole genome. The sequences are selected by relation
with a Position Weigth Matrix ( pwm ) and result are filtered by using False discovery method. The results are send
to Genome Data Viewer ( GDV ).
"""

from bbcflib.motif              import sqlite_to_false_discovery_rate
from bein                       import MiniLIMS, unique_filename_in, execution
from bbcflib.genrep             import GenRep
from bbcflib.gdv                import create_gdv_project, add_gdv_track, get_project_id
from bbcflib.frontend           import parseConfig
from bbcflib.common             import scp, normalize_url
from bbcflib.track              import Track, new, load
from os.path                    import basename, expanduser, normcase, isfile, exists
from os                         import sep
from os.path                    import basename, splitext
import getopt, sys, urllib2, logging, sqlite3

USAGE = """%s [-h] [options] --matrix = /path/to/matrix.mat -c config_file --minilims = minilims_name
-h Print this message and exit
-c --config file Config file
--host submit track file to selected host machine (e.g sugar)
--identity_file public ssh identity file (e.g ~/.ssh/id.rsa)
--matrix file path to matrix
--minilims MiniLIMS where Scanning executions and files will be stored.
--project GDV project name
--remote_path path where put track in host machine
--via via Run executions using method "via" (can be "local" or "lsf")
--username ssh username for login to host
--website website were result will be taken for GDV (e.g http://sugar.epfl.ch)
""" % (sys.argv[0])

class Usage(Exception):
    """Display to user all parameter availabe with associated description"""
    def __init__(self, msg):
        super(Usage, self).__init__()
        self.msg = msg

def fix_sqlite_db(database):
    db      = sqlite3.connect(database)
    cursor  = db.cursor()
    chr_name=  cursor.execute("SELECT name FROM chrNames")
    for chromosome in chr_name.fetchall():
        exist = cursor.execute("SELECT COUNT(name) FROM sqlite_master WHERE type='table' AND name='"+chromosome[0]+"';").fetchone()[0]
        if exist == 0:
            cursor.execute("DELETE FROM chrNames WHERE name=?", (chromosome[0],))
    #~ result  = cursor.execute("SELECT name FROM chrNames")
    #~ for i in result.fetchall():
        #~ cursor.execute("CREATE TABLE IF NOT EXISTS '"+i[0]+"'(start INTEGER, end INTEGER, score REAL);")
    db.commit()
    cursor.close()
    db.close()
    track_scanned_signal_tmp = unique_filename_in()
    track_scanned_signal = unique_filename_in()
    # temp result with duplicate element
    with new(track_scanned_signal_tmp, format = "sql", datatype= "quantitative") as t:
        with load(database) as b:
            t.chrmeta = b.chrmeta
            for chrom in b:
                for value in b.read(chrom, fields=['start','end','score', 'strand']):
                    if value[3] == 1:
                        t.write(chrom, (t.chrmeta[chrom]-(value[1]+1), t.chrmeta[chrom]-(value[0]), value[2],))
                    else:
                        t.write(chrom, (value[0:3],))
    # result without duplicate element id duplicate element exist take higher score
    with new(track_scanned_signal, format = "sql", datatype= "quantitative") as t:
        with load(database) as b:
            t.chrmeta = b.chrmeta
    db1     = sqlite3.connect(track_scanned_signal)
    db2     = sqlite3.connect(track_scanned_signal_tmp)
    cursor1 = db1.cursor()
    cursor2 = db2.cursor()
    chr_name=  cursor2.execute("SELECT name FROM chrNames")
    for chromosome in chr_name.fetchall():
        cursor1.execute("CREATE TABLE '"+chromosome[0]+"' (start INTEGER, end INTEGER, score REAL);")
        #~ values = cursor2.execute("SELECT DISTINCT start, end, score FROM '"+chromosome[0]+"'")*
        values = cursor2.execute    ("""
SELECT t.start, t.end, t.score FROM '"""+chromosome[0]+"""' t INNER JOIN (
    SELECT start, end, MAX(score) AS MAXSCORE FROM '"""+chromosome[0]+"""' GROUP BY start
) groupedt ON t.start=groupedt.start AND t.score=MAXSCORE;
                                    """)
        for v in values.fetchall():
            cursor1.execute("INSERT INTO '"+chromosome[0]+"' VALUES (?,?,?) ", (v[0], v[1], v[2],) )
    db1.commit()
    db2.commit()
    cursor1.close()
    cursor2.close()
    db2.close()
    db2.close()
    return track_scanned_signal

def main(argv = None):
    """
    Entry point when program start
    """
    genrep              = None
    assembly            = None
    lims                = None
    job                 = None
    config              = None
    config_file         = None
    background          = ""
    matrix              = ""
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
    result_path         = ""
    track_regions_path  = ""
    via                 = ""
    limspath            = ""
    fdr                 = 0
    runs                = {}
    logging.basicConfig(filename='run_scanning.log',level=logging.INFO)
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt   (
                                            argv[1:],"hu:c:"  ,
                                            [
                                                "help", "via = ", "host = "     ,
                                                "remote_path = " , "website = " ,
                                                "minilims = ","config = "       ,
                                                "matrix = ", "username = "      ,
                                                "identity_file = ", "project = "
                                            ]
                                        )
        except getopt.error, msg:
            raise Usage(msg)
        for option, value in opts:
            if option in ("-h", "--help"):
                print __doc__
                print USAGE
                sys.exit(0)
            elif option == "--via":
                if value == "local":
                    via = "local"
                elif value == "lsf":
                    via = "lsf"
                else:
                    raise Usage("Via (-u) can only be \"local\" or \"lsf\", got %s." % (value,))
            elif option == "--website":
                website = normalize_url(value)
            elif option == "--minilims":
                limspath = normcase(expanduser(value))
            elif option == "--host":
                host = value
            elif option == "--identity_file":
                identity_file = value
            elif option == "--remote_path":
                remote_path = normcase(expanduser(value))
                if not remote_path.endswith(sep):
                    remote_path += sep
            elif option == "--matrix":
                matrix = {basename(value):normcase(expanduser(value))}
            elif option == "--username":
                username = value
            elif option == "--project":
                project = value
            elif option in ("-c", "--config"):
                config_file = normcase(expanduser(value))
            else:
                raise Usage("Unhandled option: " + option)

        # read config file
        if config_file is None or not exists(config_file) or not isfile(config_file):
            raise Usage("Config file missing")
        else:
            job, config = parseConfig(normcase(expanduser(config_file)))

        if project == "":
            project = job.description
        if matrix == "":
            if "matrix" in job.options:
                path = normcase(expanduser(job.options["matrix"]))
                matrix = {basename(path): path}
            else:
                raise Usage("You need give value matrix file ")
        if limspath == "":
            if "minilims" in job.options:
                limspath = job.options["minilims"]
            else:
                raise Usage("You need give value minilims path/name")
        if via == "":
            if "via" in job.options:
                via = job.options["via"]
            else:
                via = "lsf"
        if host == "" and "host" in job.options:
            host = job.options["host"]
        if identity_file == "" and "identity_file" in job.options:
            identity_file = job.options["identity_file"]
        if remote_path == "" and "remote_path" in job.options:
            remote_path = job.options["remote_path"]
        if username == "" and "username" in job.options:
            username = job.options["username"]
        if website == "" and "website" in job.options:
            website = job.options["website"]

        genrep      = GenRep(config = config)
        assembly    = genrep.assembly(job.assembly_id)
        lims        = MiniLIMS(limspath)
        json        = create_gdv_project(
                                            config["gdv"]["key"], config["gdv"]["email"],
                                            project,
                                            assembly.nr_assembly_id,
                                            config["gdv"]["url"],
                                            public = True
                                        )
        project_id  = get_project_id( json )
        # compute false discovery rate
        with execution(lims, description = job.description) as ex:
            background = genrep.statistics  (
                                                assembly,
                                                output = unique_filename_in(),
                                                frequency = True,
                                                matrix_format = True
                                            )
            if len(job.groups) >2:
                raise ValueError("They are more than 2 group in config file")

            for group_number in job.groups:
                group = job.groups[group_number]
                for run_number in group["runs"]:
                    run_iter = job.groups[group_number]["runs"][run_number]
                    if "url" in run_iter:
                        url = run_iter["url"]
                        uri = ""
                        if run_iter["run"] not in runs:
                            runs[run_iter["run"]] = {"name":None, "control":None, "experimental":None}
                        if url.startswith("http") or url.startswith("www."):
                            url = normalize_url(url)
                            # download data
                            data    = urllib2.urlopen(url)
                            uri     = unique_filename_in()
                            with open(uri, "w") as opening_file:
                                opening_file.write(data.read())
                        else:
                            uri = normcase(expanduser(url))
                        if group["control"]:
                            runs[run_iter["run"]]["control"]   = uri
                            runs[run_iter["run"]]["name"]      = basename(uri)
                        else:
                            runs[run_iter["run"]]["experimental"] = uri

            for run in runs:
                current_run         = runs[run]
                original_sql_data   = unique_filename_in()
                random_sql_data     = unique_filename_in()
                track_filtered      = unique_filename_in()
                logging.info( "[%s]" % job.description )
                logging.info( "alias %s => %s" % (current_run["experimental"], track_filtered) )

                # convert data to sql
                with Track(current_run["experimental"], chrmeta = assembly.chromosomes) as track:
                    # Get sqlite file if is not arleady in this format
                    if track.format != "sql" or track.format != "db" or track.format != "sqlite":
                        track.convert(original_sql_data, format = "sql")
                    else:
                        original_sql_data = current_run["experimental"]
                    # Generate a random population from orginal if it is not give from config file
                    if current_run["control"] is None:
                        # create random track
                        track.shuffle_track(random_sql_data, repeat_number = 5)
                    else:
                        with Track(current_run["control"], chrmeta = assembly.chromosomes) as track_random:
                            # Get sqlite file if is not arleady in this format
                            if track_random.format != "sql" or \
                                track_random.format != "db" or \
                                track_random.format != "sqlite":
                                track_random.convert(random_sql_data, format = "sql")
                            else:
                                random_sql_data = current_run["control"]
                track_scanned, fdr, p_value = sqlite_to_false_discovery_rate(
                                                                                ex,
                                                                                matrix,
                                                                                background,
                                                                                genrep,
                                                                                assembly.chromosomes,
                                                                                original_sql_data,
                                                                                random_sql_data,
                                                                                threshold = -100,
                                                                                via = via,
                                                                                keep_max_only = False,
                                                                                alpha = 0.05,
                                                                                nb_sample = 5.0
                                                                            )

                # filter track with fdr as treshold
                with new(track_filtered, format = "sql", datatype = "qualitative") as track_out:
                    chromosome_used     = {}
                    track_out.meta_track = {"source": basename(current_run["experimental"])}
                    track_out.meta_track.update({"k":"v"})
                    with Track(track_scanned, format = "sql", chrmeta = assembly.chromosomes) as track_in:
                        meta = dict([(v["name"], dict([("length", v["length"])])) for v in track_in.chrmeta.values()])
                        for chromosome in track_in.all_chrs:
                            data_list = []
                            for data in track_in.read   (
                                                            {"chr": chromosome, "score": (fdr, sys.maxsize)},
                                                            fields = Track.qualitative_fields
                                                        ):
                                data_list.append(data)
                                chromosome_used[chromosome] = meta[chromosome]
                            if len(data_list) > 0:
                                track_out.write(chromosome, data_list)
                        track_out.chrmeta = chromosome_used
                ex.add(track_filtered,      "sql: filtred %s" % track_filtered)
                logging.info( "scanned: %s" % track_scanned )
                logging.info( "score selected: %f with p: %.3f" % (fdr, p_value) )
                logging.info( "filtred: %s" % track_filtered )


                # fix track
                track_scanned_signal = fix_sqlite_db(track_scanned)
                logging.info( "scanned signal: %s" % track_scanned_signal )
                ex.add(track_scanned_signal, description="%s: sql track signal %s" % (job.description, track_scanned_signal))

                # send filtred track and scanned track to remote
                if host != "" and remote_path != "" and username != "":
                    args = []
                    if identity_file != "":
                        args = ["-i", normcase(expanduser(identity_file)), "-C" ]
                    source_filtred      = normcase(expanduser(track_filtered))
                    source_scanned      = normcase(expanduser(track_scanned_signal))
                    result_destination          = "%s@%s:%s%s%s.db" % (username, host, remote_path, sep, track_filtered)
                    result_path                 = "%s%s%s.db" % (website, sep, track_filtered)
                    track_regions_destination   = "%s@%s:%s%s%s.db" % (username, host, remote_path, sep, track_scanned_signal)
                    track_regions_path          = "%s%s%s.db" % (website, sep, track_scanned_signal)
                    scp(ex, source_filtred, result_destination, args = args)
                    scp(ex, source_scanned, track_regions_destination, args = args)
                else:
                    result_path = track_filtered
                # Send to GDV filtred track
                add_gdv_track  (
                                    config["gdv"]["key"], config["gdv"]["email"],
                                    project_id, result_path,
                                    name    = "filtred_%s" % (splitext( basename( current_run["experimental"] ) )[0]),
                                    gdv_url = config["gdv"]["url"]
                                )
                # Send to GDV scanned track
                add_gdv_track  (
                                    config["gdv"]["key"], config["gdv"]["email"],
                                    project_id, track_regions_path,
                                    name    = "regions_%s" % (splitext( basename( current_run["experimental"] ) )[0]),
                                    gdv_url = config["gdv"]["url"]
                                )
                logging.info( "++++++++++++")
            logging.info( "-------------------END--------------------")
    except Usage, err:
        print >> sys.stderr, err.msg
        print >> sys.stderr, USAGE
        return 2
if __name__ == "__main__":
    sys.exit(main())
