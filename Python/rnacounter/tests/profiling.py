#!/usr/bin/env python
# encoding: utf-8

""" Run

python profiling.py -o zzz -n 1 -c chr18 /archive/epfl/bbcf/jdelafon/test_rnaseq/mefKO.bam /archive/epfl/bbcf/jdelafon/test_rnaseq/mm9_renamed.gtf

python profiling.py -o zzz -n 1 testfiles/gapdhKO.bam testfiles/mm9_3genes_renamed.gtf ;

(need "-o sth" because closing stdout produces a "ValueError: I/O operation or closed file".)
"""

import pstats, cProfile
import pyximport
pyximport.install()
import rnacounter
from rnacounter import rnacounter_main, parse_args
from docopt import docopt

args = docopt(rnacounter.__doc__, version='0.1')
bamname, annotname, options = parse_args(args)

cProfile.runctx("rnacounter_main('%s','%s',%s)" % (bamname,annotname,options),
                globals(), locals(), "stats_profiling.prof")
s = pstats.Stats("stats_profiling.prof")
s.strip_dirs().sort_stats("time").print_stats()




