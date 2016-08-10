'''
Created on Mar 18, 2016
@author: eifu
'''
import unittest
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
import histone


class TestHistoneWithDNAModel(unittest.TestCase):
    def testInitiate(self):
        tes1 = histone.HistoneWithDNAModel(position=1)

        self.assertEqual(tes1.position, 1, "position failure")

    def testCpGsite(self):
        test1 = histone.HistoneWithDNAModel(position=-2)
        self.assertEqual(len(test1.CpGislandlist), 4, "something wrong with position -2")

        test2 = histone.HistoneWithDNAModel(position=-1)
        self.assertEqual(len(test2.CpGislandlist), 4, "something wrong with position -1")

        test3 = histone.HistoneWithDNAModel(position=0)
        self.assertEqual(len(test3.CpGislandlist), 2, "something wrong with position 0")

    def testRandomizedFunctions(self):
        test = histone.HistoneWithDNAModel()
        histone.HistoneWithDNAModel.k(test)

    def testVectorize(self):
        test = histone.init_genome_with_dna_model(percentage=100)
        v = histone.vectorize_with_dna_model(test)
        self.assertEqual(v.shape, (4, 81))
        self.assertEqual(sum(v[0]), 81)
        self.assertEqual(sum(v[1]), 0)
        self.assertEqual(sum(v[2]), 0)

        test2 = histone.init_genome_with_dna_model(percentage=0)
        v = histone.vectorize_with_dna_model(test2)
        self.assertEqual(sum(v[0]), 0)
        self.assertEqual(sum(v[1]), 81)
        self.assertEqual(sum(v[2]), 0)


if __name__ == 'main':
    unittest.main()
