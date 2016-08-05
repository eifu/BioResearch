'''
Created on May 23, 2016
@author: eifu
'''
import unittest
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
import histone


class TestHistone(unittest.TestCase):
    def testInitiate(self):
        m1 = histone.MHistone(position=1)

        self.assertEqual(m1.status, 'm', 'status failure')
        self.assertEqual(m1.position, 1, "position failure")

    def testIndividualConnection(self):
        m1 = histone.MHistone(position=1)
        m2 = histone.MHistone(position=2, prenode=m1)

        self.assertEqual(m2.prenode.position, 1, 'initiate connection failure')
        self.assertEqual(m1.nextnode, None, 'initiate connection failure 2')

    def testset_adjHistone(self):
        m1 = histone.MHistone(position=1)
        m2 = histone.MHistone(position=2, prenode=m1)
        m1.set_adjhistone(m2)

        self.assertEqual(m2.prenode.position, 1, 'initiate connection failure')
        self.assertNotEqual(m1.nextnode, None, 'initiate connection failure 2')
        self.assertEqual(m1.nextnode, m2, 'initiate connection failure 3')

    def testRandomHistListMethod(self):
        hist_list = histone.init_genome(50, 0, 3, 1)

        self.assertEqual(len(hist_list), 3, "number of elements in list failure")
        self.assertEqual(hist_list[0].nextnode, hist_list[1], 'connection failure 1-1')
        self.assertEqual(hist_list[1].nextnode, hist_list[2], 'connection failure 1-2')
        self.assertEqual(hist_list[2].prenode, hist_list[1], 'connection failure 2-1')
        self.assertEqual(hist_list[1].prenode, hist_list[0], 'connection failure 2-2')

    def testRandomHistListMethod_disconnection(self):
        hist_list = histone.init_genome(50, 0, 3, 1)

        self.assertEqual(hist_list[0].prenode, None, 'disconnection Head and Tail')
        self.assertEqual(hist_list[2].nextnode, None, 'disconnection Head and Tail 2')

        histList2 = histone.init_genome(50, 0, 40, 1)
        self.assertEqual(histList2[0].prenode, None, 'disconnection Head and Tail 3')
        self.assertEqual(histList2[-1].nextnode, None, 'disconnection Head and Tail 4')

    def testbitvec(self):
        hist_list = histone.init_genome()
        bitvec = histone.vectorize(hist_list)
        self.assertTrue(len(bitvec[0]) == 81, 'bitvec num failure')
        self.assertTrue(sum(bitvec[0]) + sum(bitvec[1]) + sum(bitvec[2]) == 81, 'bitvec distri failure')

        hs2 = histone.init_genome(percentage=100)
        bitvec2 = histone.vectorize(hs2)
        self.assertTrue(sum(bitvec2[0]) == 81, 'methylated histone list failure')

        hs3 = histone.init_genome(percentage=0)
        bitvec3 = histone.vectorize(hs3)
        self.assertTrue(sum(bitvec3[1]) == 81, 'unmethylated histone list failure')

    def test_vect(self):
        for _ in range(10):
            hist_list = histone.init_genome(percentage=100)
            bitvec = histone.vectorize(hist_list)
            self.assertEqual(sum(bitvec[0]), 81)


    def testNextGen(self):
        for _ in range(10):
            hist_list = histone.init_genome(percentage=100)
            hist_list2, _, _ = histone.next_genome(hst_list=hist_list,
                                                   window=10,
                                                   k_nuc=1)
            self.assertEqual(hist_list2[40].status, "m")

    def testNextGen2(self):
        m1 = histone.MHistone(position=-2)
        m2 = histone.AHistone(position=-1)
        m3 = histone.AHistone(position=0)
        m4 = histone.AHistone(position=1)
        m5 = histone.AHistone(position=2)
        hist_list = [m1, m2, m3, m4, m5]

        hist_list2 = histone.next_genome(hst_list=hist_list,
                                         window=3,
                                         k_nuc=0.12)

    def testRandom(self):
        for _ in range(10):
            hist_list = histone.init_genome(percentage=0, a_bool=0)
            hist_list2, _, _ = histone.next_genome(hst_list=hist_list,
                                                   window=10,
                                                   k_nuc=0)
            self.assertEqual(hist_list2[40].status, "u")


if __name__ == 'main':
    unittest.main()
