#!/usr/bin/env python
"""
External script to run microbiome functions
"""

import sys
from bbcflib.microbiome import bam_to_annot_counts, getCountsPerLevel,combine_counts

funcnames = ["bam_to_annot_counts","getCountsPerLevel","combine_counts"]

if len(sys.argv)<2 or sys.argv[1] not in funcnames:
    raise ValueError("Choose one function from %s" %(", ".join(funcnames)))
elif sys.argv[1] == "bam_to_annot_counts":
    if len(sys.argv) < 4:
        raise ValueError("Need at least one bam file and an annotation file")
    bamfiles = sys.argv[2].split(",")
    annotations = sys.argv[3]
    if len(sys.argv) > 4:
        pref_name = sys.argv[4]
    else:
        pref_name = ''
    if len(sys.argv) > 5:
        output = sys.argv[5]
    else:
        output = None
    output = bam_to_annot_counts( bamfiles, annotations, pref_name, output )
elif sys.argv[1] == "getCountsPerLevel":
    if len(sys.argv) < 4:
        raise ValueError("Need at least one file and a level name")
    infile = sys.argv[2]
    level = sys.argv[3]
    if len(sys.argv) > 4:
        output = sys.argv[4]
    else:
        output = None
    output = getCountsPerLevel( infile, level, output )
elif sys.argv[1] == "combine_counts":
    if len(sys.argv) < 5:
        raise ValueError("Need at least one counts file, keys and counts column numbers")
    files = sys.argv[2].split(",")
    idsColsKey = [int(n) for n in sys.argv[3].split(",")]
    idsColsCounts = [int(n) for n in sys.argv[4].split(",")]
    if len(sys.argv) > 5:
        output = sys.argv[5]
    else:
        output = None
    output = combine_counts(counts, idsColsKey, idsColsCounts, output)
