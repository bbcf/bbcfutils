
Initial files:
===============
* draft.py:
Python only version, for testing purposes. Run as "python draft.py".

* rnacounter.py:
Copy of draft.py that we compile with setup.py. Imports functions from rnacounter0.

* rnacounter0.pyx:
Cython functions, to be compiled with setup0.py and then imported by rnacounter.py.

* setup.py:
The "Makefile" for rnacounter.py.
Produces rnacounter.c, rnacounter.so.

* setup0.py:
The "Makefile" for rnacounter0.pyx.
Produces rnacounter0.c, rnacounter0.so.

* rnacounter
The (python) executable, imports and runs the C version after compilation.

* benchmark.txt: some execution timings reported, in different conditions.

* testfiles/:
Folder with testing files, including
- gapdhKO.bam: alignment on mm9 with only Gapdh covered.
- mm9_3genes_renamed.gtf: extract of the Ensembl GTF with Gapdh, the gene before and the gene after.
- mm9_Gapdh_renamed.gtf: extract of the Ensembl GTF with Gapdh only.

* build/: platform-specific compilation stuff.


Backup:
========
* rnacounter.pxd:
Redefinitions of rnacounter.py functions headers for faster execution
(specifying C types for the functions' variables before compilation).
Not used atm.
* diverse working versions of draft.py


To compile and run:
====================
1. Compile rnacounter0.pyx:
    python setup0.py build_ext --inplace ;

2. If there is an update in draft.py to include:
    cp draft.py rnacounter.py ;

3. Compile rnacounter.py:
    python setup.py build_ext --inplace ;

Note: On OSX Mavericks, XCode 5 is bugged and clang requires to add the following
to run without raising an error:

    ARCHFLAGS=-Wno-error=unused-command-line-argument-hard-error-in-future \
    python setup.py build_ext --inplace ;

4. Run:
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
