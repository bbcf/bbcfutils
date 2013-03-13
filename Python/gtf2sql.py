#!/usr/bin/env python

from bbcflib import genrep
import os, getopt, sys

opts = dict(getopt.getopt(sys.argv[1:],"d:",[])[0])
basepath = opts.get('-d') or "/data/epfl/bbcf/genrep/nr_assemblies"
basepath += "/%s"
for _a,info in genrep.GenRep().assemblies_available():
    for n in range(100):
        assembly = genrep.Assembly(_a)
        gtf_path = os.path.join(basepath%"gtf","%s_%i.gtf.gz"%(assembly.md5,n))
        if not(assembly.bbcf_valid and os.path.exists(gtf_path)): break
        sql_path = os.path.join(basepath%"annot_tracks","%s_%i.sql"%(assembly.md5,n))
        if os.path.exists(sql_path): continue
        print info, gtf_path, sql_path
        assembly.gtf_to_sql(gtf_path=gtf_path, sql_path=sql_path)
