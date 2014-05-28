
Initial files:
===============

* draft.py:
Python only version, for testing purposes. Run as "python draft.py".

* rnacounter.py:
The one we compile sith setup.py.

* rnacounter.pyx:
The Cython version, still in development.

* rnacounter.pxd:
Redefinitions of rnacounter.py functions headers for faster execution
(specifying C types for the functions' variables before compilation).

* setup.py:
The "Makefile" for Cython.

* rnacounter
The (python) executable, imports and runs the C version after compilation.

* testfiles/:
Folder with testing files, including
- gapdhKO.bam: alignment on mm9 with only Gapdh covered.
- mm9_3genes_renamed.gtf: extract of the Ensembl GTF with Gapdh, the gene before and the gene after.
- mm9_Gapdh_renamed.gtf: extract of the Ensembl GTF with Gapdh only.

* benchmark.txt: some execution timings reported, in different conditions.


After compilation:
===================
* rnacounter.c: C source.
* rnacounter.so: executabe.
* build/: platform-specific compilation stuff.


To compile:
============
1. If there is an update in draft.py to include:
    cp draft.py rnacounter.py ;

2. Compilation:
    python setup.py build_ext --inplace ;

Note: On OSX Mavericks, XCode 5 is bugged and clang requires to add the following
to run without raising an error:

    ARCHFLAGS=-Wno-error=unused-command-line-argument-hard-error-in-future \
    python setup.py build_ext --inplace ;


To run:
========
./rnacounter --help

Example:
./rnacounter testfiles/gapdhKO.bam testfiles/mm9_3genes_renamed.gtf


Testing:
=========
Uni tests to come.

The BAM contains 4041 reads all aligning perfectly on Gapdh (ENSMUSG00000057666) exons,
mostly on ENSMUSE00000487077 but also ENSMUSE00000751942 and ENSMUSE00000886744.
Nothing on other exons, which makes it a good example of badly conditioned input data...

The least squares method returns counts on the following transcripts:
ENSMUST00000117757, ENSMUST00000118875, ENSMUST00000147954
and nothing on ENSMUST00000073605, ENSMUST00000144205, ENSMUST00000144588 .

Returns a count of 2459.62 (1091.71 RPK) for the gene.
When counting a rough 1.0 per read aligned, we get around 7000.

Reports to "benchmark.txt".
