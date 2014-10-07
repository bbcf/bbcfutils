Using rnacounter to count reads in genomic intervals
====================================================

`rnacounter` estimates abundances of genes and their different transcripts
from read densities. Exons and introns can also be quantified.
It requires a genome-level BAM file and a
GTF/GFF file describing the exon structure, such as those provided by Ensembl or GenRep.
The GTF is assumed to be sorted at least w.r.t. chromosome name,
and the chromosome identifiers in the GTF must be the same as the BAM references.
The method used is described in [<ref>].

Most basic usage::

   rnacounter test.bam test.gtf > counts_table.txt

The command above will create a tab-delimited text file of gene counts, their RPKM
and annotation information sur as the genomic location.
Many options are then available, listed below.

Use `rnacounter join` to merge several output files produced using **the same GTF**,
to create a single table with counts from all samples::

   rnacounter join tab1.txt tab2.txt ... > tab_all.txt


Options:
--------

* :option:`-h`, :option:`--help` and :option:`-v`, :option:`--version`::

    Display information about the program usage / the version currently installed.

* :option:`-s`, :option:`--stranded`::

    If the protocol was strand-specific and this option is provided,
    sense and antisense counts are both reported in two consecutive lines
    with a different tag in the last column.
    They can be split afterwards by piping the result as for instance with
    `... | grep 'antisense'`.
    Using the `--threshold` option together with `--stranded`
    will exclude only elements with both sense and antisense counts under the threshold.

* :option:`-n`, :option:`--normalize`::

    RPKM are automatically calculated together with raw read counts. RPKM are counts
    divided by the length of the transcript as well as by a sample-specific
    normalization constant, usually the total number of aligned reads in the sample (default).
    This value can be changed to a user-defined integer.
    Typically, if you want to compare the same gene in several samples,
    the normalization will cancel out anyway
    and giving `-n 1` will speed up the process since it will skip counting the alignments.
    Some stats programs also require raw counts anyway and do their own normalization.
    To get FPKM instead, see `--fraglength`.

* :option:`-f`, :option:`--fraglength`::

    Since in a transcript of length L there are only L-F+1 different positions where
    a fragment of length F can be cut, one may want to correct for this bias before RPKM
    calculation (then usually called FPKM). Typical fragment lengths are around 350nt;
    default value is 1 (no correction). This is not to be confused with the read length.
    This option can be applied only at the gene- or transcript level.

* :option:`--nh`::

    A flag "NH" can be added to BAM files to indicate the number of times the read
    could be mapped to different locations in the genome. Adding this option
    will take this number into account by adding 1/NH instead of 1 to an exon read count.

* :option:`--noheader`::

    By default the program adds one line with column descriptors on top of the output file.
    For easier piping the result to some other program, one can choose
    not to add the header by adding this option.

* :option:`--exon_cutoff`::

    Often the annotation contains (sometimes artificial) transcript structures that are
    very close to each other and are thus hard to dinstinguish for any model due to
    the read length constraint and lack of coverage on small regions, reducing
    the model's power.
    To address this, one can merge transcripts differing by exonic
    regions of less than that many nucleotides. In the output, only one
    record will be reported, but synonyms will be added in a supplementary column.
    Defaults to read length. Set to 0 to remove transcripts filtering, especially
    with "local" alignments, or to a bigger number to reduce the transcripts variety.

* :option:`--threshold`::

    Features with counts inferior or equal to the given threshold (positive number)
    will not be reported in the ouput. By default everything is reported
    - even with zero counts.

* :option:`--gtf_type`::

    Usually one uses standard (Ensembl etc.) GTF files to count reads in
    exons/genes/transcripts. The only lines of interest are then the ones with
    value "exon" (default) in the 3rd column. If you are counting something else
    or provided your own, differently formatted GTF, with this option you can specify
    the 3rd column value of the lines to consider.

* :option:`--format`::

    One can also give an annotation file in BED format with 4 fields
    (chromosone, start, end, name), in which case each line
    is considered as an independant, disjoint interval with no splicing structure.
    Default is "gtf", can be changed to "bed".
    The 4th column of the BED format (name) must contain *unique* IDs.
    If the input format is "bed", the program cannot know which type of intervals
    is represented, thus will always report them as 'genes' in the output.
    Consistently, it cannot be used in conjunction with the :option:`--type` option.

* :option:`-t`, :option:`--type`::

    The type of feature you want to count reads in. Can be "genes" (default),
    "transcripts", "exons" or "introns".
    One can give multiple comma-separated values, in which case all
    the different features will be mixed in the output but can easily be split
    using the last column tag, as for instance with `... | grep 'exon'`.
    Then if `--method` is specified it must have the same number of values as `type`,
    also as a comma-separated list, or a single one that is applied to all types.

* :option:`-c`, :option:`--chromosomes`::

    Consider only a subset of all chromosomes by providing a comma-separated list
    of chromosome names (that must match those of the GTF and BAM).

* :option:`-o`, :option:`--output`::

    The output is `stdout` by default (output directly to screen), which permits
    redirection to a file. Alternatively one can redirect the standard output to
    a file using this option. If the file name already exists, it will be overwritten.

* :option:`-m`, :option:`--method`::

    Feature counts are inferred from the counts on (slices of) exons
    with the chosen `--method`: "raw" (htseq-count-like) or
    "nnls" (non-negative least squares, see [<ref>]).
    The default is "raw" to not disturb habits, but "nnls" is advised,
    and should be mandatory at the transcripts level (see Example below).


Miscellaneous notes:
--------------------

* Overlapping regions:
  In "raw" counting mode, regions spanned by exon from two or more genes,
  together with the alignements inside these regions, are ignored (ambiguous).
  The "nnls" mode tries to resolve the ambiguity in the same way
  it does for multiple isoforms.

* Multiple alignments:
  Rather than an option/default to remove multiply mapping reads, this filtering
  - if desired - should be done at the mapping step choosing the right parameters,
  or the BAM file can be filtered afterwards. On the contrary if you want to keep
  multiple mapping, you can use the `--nh` option.

* Exons and introns:
  Because annotated exons often overlap a lot, in "raw" mode, "exon" counts are actually
  that of their disjoint slices, and their name in the output table is formatted as
  "exon1|exon2" if a slice is spanned by exon1 and exon2. In "nnls" mode, exon counts
  are inferred from disjoint slices as for genes.

  Intronic regions also annotated as exons in some alternative transcripts are
  ignored whatever the chosen method is. Because they don't have official IDs,
  introns slices are given names following this pattern:
  "<n>I-<gene_id>", if it is the n-th intron of that gene.

* Non-integer counts:
  The fact that some reads cross exon boundaries as well as considering the NH flag
  make the reported number not be integers. They still represent count data and can
  be rounded afterwards if necessary.

* Custom input:
  If your GTF does not represent exons but custom genomic intervals to simply count
  reads in, provide at least a unique `exon_id` in the attributes as a feature name,
  and the type field (column 3) must be set to 'exon' or specified with the
  `--gtf_ftype` option. If not specified, `gene_id`, `transcript_id` and `exon_id`
  will all get the value of `exon_id`.

* Paired-end support:
  At the moment alignments of paired-end reads are not treated specially, i.e.
  all reads are considered as single-end.


Examples:
---------

* Probably the best way to get isoforms counts::

    rnacounter -t transcripts -m nnls --nh -f 350 sample.bam mouse.gtf > transcript_counts.txt

* Compare gene counts between two conditions, HTSeq-like::

    rnacounter -n 1 group1.bam mouse.gtf > gene_counts1.txt
    rnacounter -n 1 group2.bam mouse.gtf > gene_counts2.txt
    rnacounter join gene_counts1.txt gene_counts2.txt > gene_counts.txt

  Then send it to DESeq/EdgeR/whatever other stats program that asks for such a table.


Reference:
----------

<?>

