
Welcome!
========
`Rnacounter` estimates abundances of genes and their different transcripts
from read alignments. Exons and introns can also be quantified.

It provides fast read counting in annotated genomic features as well as a simple,
yet efficient solution to the quantification of isoforms from RNA-seq data.
The method used is described in [<ref>].
A typical run is expected to take less than 2 minutes for a 1Gb BAM file from mouse
RNA sequencing, increasing linearly with the BAM size.

For all these tasks it only requires a BAM file from a read mapping on the genome,
and a single GTF/GFF file describing the exon structure
such as those provided by Ensembl or GenRep.

It is not meant to be used as a library, but through its command-line tool "rnacounter".

The code project is hosted in Github (https://github.com/delafont/rnacounter), GPL-2 licensed.

Usage:
======
See "rnacounter --help" and the tutorial at
http://bbcf.epfl.ch/bbcflib/tutorial_rnacounter.html,
also available in the doc/ folder.

Minimal example::

    rnacounter test.bam test.gtf

Installation:
=============
First ensure that you have numpy installed, or install it with::

    sudo easy_install numpy

Then install rnacounter with easy_install::

    sudo easy_install rnacounter

Or better yet, with pip::

    sudo pip install numpy
    sudo pip install rnacounter

It installs as a standard Python library but includes the executable
and puts it somewhere in your $PATH. Dependencies will be added
automatically.

To uninstall with pip::

    sudo pip uninstall rnacounter

The code is fully compatible with Python 2.7 and Python 3.

Building from source:
=====================
This allows to modify the Cython source code (rnacounter.pyx) before rebuilding.

Clone or download the repository from https://github.com/delafont/rnacounter .

You need cython installed (`pip install cython`).

From where rnacounter.pyx lies (rnacounter/rnacounter/), run::

    sudo python setup.py build_ext

It will recompile to create rnacounter.c, and build it.
Then add the executable (rnacounter/bin/rnacounter) to your $PATH,
or install from the package root (rnacounter/) with::

    sudo python setup.py install

Dependencies:
=============
Tests run with the library versions below, but may work with earlier versions.

* setuptools 7.0+  (installation)
* pysam 0.7.5+     (samtools wrapper)
* numpy 1.6.2+     (efficient numeric arrays)
* scipy 0.9.0+     (NNLS algorithm)
* docopt 0.6.1+    (command-line args parsing)
* cython 0.20+     (translate Python code to C)

Testing:
========
Testing files in the testfiles/ folder:
- gapdhKO.bam: alignment on mm9 with only Gapdh covered.
- mm9_3genes_renamed.gtf: extract of the Ensembl GTF with Gapdh, the gene before and the gene after it.
- mm9_Gapdh_renamed.gtf: extract of the Ensembl GTF with Gapdh only.

Example::

    rnacounter testfiles/gapdhKO.bam testfiles/mm9_3genes_renamed.gtf

The BAM contains 4041 reads all aligning perfectly on Gapdh (ENSMUSG00000057666) exons,
mostly on ENSMUSE00000487077 but also ENSMUSE00000751942 and ENSMUSE00000886744.
Nothing on other exons, which makes it a good example of badly conditioned input data...

The least squares method returns counts on the following transcripts:
ENSMUST00000117757, ENSMUST00000118875, ENSMUST00000147954
and nothing on ENSMUST00000073605, ENSMUST00000144205, ENSMUST00000144588 .

Returns a count of 2456.87 for the gene.

Troubleshooting:
================
Any bug report, usage issue or feature request not listed below can be addressed to
julien.delafontaine@epfl.ch or webmaster.bbcf@epfl.ch .

