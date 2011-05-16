import pickle_exon_mapping as p
from Bio import SeqIO
import StringIO
from unittest2 import TestCase, TestSuite, main, TestLoader, skipIf



class TestGenerateMapping(TestCase):
    fasta = StringIO.StringIO(""">YNL228W.1|YNL228W|220646|221422|1
ATGGTTTGGTGTCACTATATTCTTTTGGTCTTGACTTTCTTTCTTTTCACTACGTTTTTC
>YBL095W.1|YBL095W|43274|44086|1
ATGTCCAGAACTATTCCATTTCTATTTAAATTAGTCAACAGGGCAGTAATTTTGCCTACG
>YGR266W.1|YGR266W|1022662|1024767|1
ATGCATGCTACAAACTGGTTCGACGATTGGAACCCAGAAGCTCTTTATAGAGACGATGTC
>YNL228W.2|YNL228W|53342|1312|1
ATGAGTAGTGGA
>YNL228W.3|YNL228W|552343|134|1
TGATG
>YGR266W.2|YGR266W|1232|5252|1
CCC""")

    def test_mapping(self):
        seqs = SeqIO.parse(self.fasta, 'fasta')
        d = p.generate_mapping(seqs)
        gene_labels = ['YNL228W', 'YBL095W', 'YGR266W']
        mapping = [0, 1, 2, 0, 0, 2]
        self.assertEqual(d[0], gene_labels)
        self.assertEqual(d[1], mapping)

main()


    
