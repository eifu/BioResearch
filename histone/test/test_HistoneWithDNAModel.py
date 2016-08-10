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

    def testRandamizedFunctions(self):
        # for _ in range(10):
        test = histone.HistoneWithDNAModel()
        histone.HistoneWithDNAModel.k(test)


if __name__ == 'main':
    unittest.main()
