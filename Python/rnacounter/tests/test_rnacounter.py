
# Unitesting modules #
try:
    import unittest2 as unittest
    assert unittest
except ImportError:
    import unittest
from numpy.testing import assert_almost_equal

# Nosetest flag #
__test__ = True


from rnacounter import *
import pysam


class Test_Parse_GTF(unittest.TestCase):
    def setUp(self):
        pass
    def test_parse_gtf(self):
        pass


class Test_Cobble(unittest.TestCase):
    # T1  = ===   =======  =  === ===
    # T2     === === ===  === ==   ==
    # T3      =====   =       =     =
    def setUp(self):
        self.exons = []
        self.sts = [0,2,3,4,7, 8, 11,12,16,17,20,21,22,26,27,28]
        self.ens = [1,5,6,9,10,15,14,13,19,18,25,24,23,29,30,31]
        trans = {'T1':(0,2,8,17,20,26), 'T2':(3,7,11,6,21,27), 'T3':(4,8,12,17,22,28)}
        for i in range(len(self.sts)):
            for x in trans:
                if self.sts[i] in trans[x]: t = x; break
            self.exons.append(Exon(chrom='c',start=self.sts[i]*10,end=self.ens[i]*10,gene_id='G',gene_name='g',
                                   name='E%d'%i,strand=-1,transcripts=[t]))
    def test_add_exons(self):
        newexon = self.exons[1] & self.exons[2]
        self.assertListEqual(sorted(list(newexon.transcripts)), ['T1','T2'])
        self.assertEqual(sorted(newexon.name.split('|')), ['E1','E2'])
        self.assertEqual(newexon.gene_name, 'g')
        self.assertEqual(newexon.gene_id, 'G')
        self.assertEqual(newexon.chrom, 'c')
        self.assertEqual(newexon.strand, -1)
    def test_intersect_exons_list(self):
        newexon = intersect_exons_list([self.exons[1]]+self.exons[1:4])
        self.assertListEqual(sorted(list(newexon.transcripts)), ['T1','T2','T3'])
        self.assertListEqual(sorted(newexon.name.split('|')), ['E1','E2','E3'])
        # multiple=True
        newexon = intersect_exons_list([self.exons[1]]+self.exons[1:4], multiple=True)
        self.assertEqual(newexon.name, 'E1|E1|E2|E3')
    def test_cobble(self):
        cobbled = cobble(self.exons)
        csts = [0, 2,3,4,5,6,7,8,9, 10,11,12,13,14, 16,17,18,20,21,22,23,24, 26,27,28,29,30]
        cens = [1, 3,4,5,6,7,8,9,10,11,12,13,14,15, 17,18,19,21,22,23,24,25, 27,28,29,30,31]
        self.assertListEqual([x.start/10 for x in cobbled], csts)
        self.assertListEqual([x.end/10 for x in cobbled], cens)
        self.assertListEqual(sorted(cobbled[3].name.split('|')), ['E1','E2','E3'])
        self.assertListEqual(sorted(cobbled[7].name.split('|')), ['E3','E4','E5'])
        self.assertListEqual(sorted(cobbled[11].name.split('|')), ['E5','E6','E7'])


class Test_Counting(unittest.TestCase):
    def setUp(self):
        pass
    def test_RPK(self):
        cnt = 1000
        length = 20
        norm_cst = 10
        rpk = toRPK(cnt,length,norm_cst)
        self.assertEqual(rpk,5000.)
        cnt2 = fromRPK(rpk,length,norm_cst)
        self.assertEqual(cnt2,cnt)

    def test_count_reads(self):
        # count_reads(exons,ckreads,multiple,stranded)
        pass
    def test_estimate_expression(self):
        #estimate_expression(feat_class, pieces, ids)
        pass


class Test_Executable(unittest.TestCase):
    def setUp(self):
        pass


#----------------------------------------------#
# This code was written by Julien Delafontaine #
# EPFL, BBCF: http://bbcf.epfl.ch/             #
# webmaster.bbcf@epfl.ch                       #
#----------------------------------------------#
