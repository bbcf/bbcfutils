
Usage:
======
See "rnacounter --help" and the tutorial at [...<bbcflib tutorials>].


Files description:
==================
* draft_nocython.py:
Python only version (same functionality), for testing purposes. Run as "python draft.py ...".

* rnacounter.pyx:
Cython version that we compile with setup.py.

* setup.py:
The "Makefile" for rnacounter.py.
Produces rnacounter.c, rnacounter.so.

* rnacounter
The (python) executable, imports and runs the C version after compilation.

* benchmark.txt: some execution timings reported, in different conditions.

* tests/:
- test_rnacounter.py: unit tests, run with "nosetests test_rnacounter.py".
- profiling.py: copy of rnacounter with profiling (timing) of internal functions.

* testfiles/:
Folder with testing files, including
- gapdhKO.bam: alignment on mm9 with only Gapdh covered.
- mm9_3genes_renamed.gtf: extract of the Ensembl GTF with Gapdh, the gene before and the gene after.
- mm9_Gapdh_renamed.gtf: extract of the Ensembl GTF with Gapdh only.

* build/: platform-specific compilation stuff.

* backup/:
Contains in particular "rnacounter.pxd":
Redefinitions of rnacounter.py functions headers for faster execution
(specifying C types for the functions' variables before compilation).


To compile and run:
===================

1. Compile rnacounter.pyx:
    python setup.py build_ext --inplace ;

Note: On OSX Mavericks, XCode 5 is bugged and clang requires to add the following flag
to run without raising an error:

    ARCHFLAGS=-Wno-error=unused-command-line-argument-hard-error-in-future \
    python setup.py build_ext --inplace ;

2. Run:
    rnacounter --help
Example:
    rnacounter testfiles/gapdhKO.bam testfiles/mm9_3genes_renamed.gtf


Testing:
=========
Unit tests in folder tests/

The BAM contains 4041 reads all aligning perfectly on Gapdh (ENSMUSG00000057666) exons,
mostly on ENSMUSE00000487077 but also ENSMUSE00000751942 and ENSMUSE00000886744.
Nothing on other exons, which makes it a good example of badly conditioned input data...

The least squares method returns counts on the following transcripts:
ENSMUST00000117757, ENSMUST00000118875, ENSMUST00000147954
and nothing on ENSMUST00000073605, ENSMUST00000144205, ENSMUST00000144588 .

Returns a count of 2459.62 (1091.71 RPK) for the gene.

Reports to "benchmark.txt".
