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
        hist_list = histone.init_genome(hst_n=50)

        self.assertEqual(len(hist_list), 50, "number of elements in list failure")
        for i in range(1, 50):
            self.assertEqual(hist_list[i - 1].nextnode, hist_list[i])
            self.assertEqual(hist_list[i].prenode, hist_list[i - 1])

    def testRandomHistListMethod_disconnection(self):
        hist_list = histone.init_genome(50, 0, 3, 1)

        self.assertEqual(hist_list[0].prenode, None, 'disconnection Head and Tail')
        self.assertEqual(hist_list[2].nextnode, None, 'disconnection Head and Tail 2')

        hist_list2 = histone.init_genome(50, 0, 40, 1)
        self.assertEqual(hist_list2[0].prenode, None, 'disconnection Head and Tail 3')
        self.assertEqual(hist_list2[-1].nextnode, None, 'disconnection Head and Tail 4')

    def testbitvec(self):
        hist_list = histone.init_genome()
        bitvec = histone.vectorize(hist_list)
        self.assertTrue(len(bitvec[0]) == 81, 'bitvec num failure')
        self.assertTrue(sum(bitvec[0]) + sum(bitvec[1]) + sum(bitvec[2]) == 81, 'bitvec distri failure')

        hs2 = histone.init_genome(percentage=100)
        bitvec2 = histone.vectorize(hs2)
        self.assertTrue(sum(bitvec2[0]) == 81, 'methylated histone list failure')
        self.assertTrue(sum(bitvec2[0]) + sum(bitvec2[1]) + sum(bitvec2[2]) == 81, 'bitvec distri failure')

        hs3 = histone.init_genome(percentage=0)
        bitvec3 = histone.vectorize(hs3)
        self.assertTrue(sum(bitvec3[1]) == 81, 'unmethylated histone list failure')
        self.assertTrue(sum(bitvec3[0]) + sum(bitvec3[1]) + sum(bitvec3[2]) == 81, 'bitvec distri failure')

    def test_vect(self):
        for _ in range(10):
            hist_list = histone.init_genome(percentage=50)
            bitvec = histone.vectorize(hist_list)
            self.assertEqual(sum(bitvec[2]), 0)

    def testNextGen(self):
        for _ in range(10):
            hist_list = histone.init_genome(percentage=100)
            hist_list2, _, _ = histone.next_genome(hst_list=hist_list,
                                                   window=10,
                                                   k_nuc=1)
            self.assertEqual(hist_list2[40].status, "m")

    def testNextGenome_AllUnmethylated(self):
        for _ in range(10):
            m0 = histone.UHistone(position=-3,
                                  kp2=0,
                                  km=1)
            m1 = histone.AHistone(position=-2,
                                  kp2=0,
                                  km=1,
                                  prenode=m0)
            m0.set_adjhistone(m1)
            m2 = histone.AHistone(position=-1,
                                  kp2=0,
                                  km=1,
                                  prenode=m1)
            m1.set_adjhistone(m2)
            m3 = histone.AHistone(position=0,
                                  kp2=0,
                                  km=1,
                                  prenode=m2)
            m2.set_adjhistone(m3)
            m4 = histone.AHistone(position=1,
                                  kp2=0,
                                  km=1,
                                  prenode=m3)
            m3.set_adjhistone(m4)
            m5 = histone.AHistone(position=2,
                                  kp2=0,
                                  km=1,
                                  prenode=m4)
            m4.set_adjhistone(m5)
            m6 = histone.AHistone(position=3,
                                  kp2=0,
                                  km=1,
                                  prenode=m5)
            m5.set_adjhistone(m6)
            hist_list = [m0, m1, m2, m3, m4, m5, m6]
            for h in hist_list:
                h.set_ka(a_bool=False, k_ace=0)

            hist_list2, _, _ = histone.next_genome(hst_list=hist_list,
                                                   window=2,
                                                   k_nuc=0)
            for h in hist_list:
                self.assertTrue(h.status == "u")

    def testNextGenome_AllAcetylated(self):
        for _ in range(10):
            m0 = histone.AHistone(position=-3,
                                  kp2=0,
                                  km=0)
            m1 = histone.AHistone(position=-2,
                                  kp2=0,
                                  km=0,
                                  prenode=m0)
            m0.set_adjhistone(m1)
            m2 = histone.AHistone(position=-1,
                                  kp2=0,
                                  km=0,
                                  prenode=m1)
            m1.set_adjhistone(m2)
            m3 = histone.AHistone(position=0,
                                  kp2=0,
                                  km=0,
                                  prenode=m2)
            m2.set_adjhistone(m3)
            m4 = histone.AHistone(position=1,
                                  kp2=0,
                                  km=0,
                                  prenode=m3)
            m3.set_adjhistone(m4)
            m5 = histone.AHistone(position=2,
                                  kp2=0,
                                  km=0,
                                  prenode=m4)
            m4.set_adjhistone(m5)
            m6 = histone.AHistone(position=3,
                                  kp2=0,
                                  km=0,
                                  prenode=m5)
            m5.set_adjhistone(m6)
            hist_list = [m0, m1, m2, m3, m4, m5, m6]
            for h in hist_list:
                h.set_ka(a_bool=True, k_ace=1)

            hist_list2, _, _ = histone.next_genome(hst_list=hist_list,
                                                   window=2,
                                                   k_nuc=0)
            for h in hist_list:
                self.assertTrue(h.status == "a")

    def testRandom(self):
        for _ in range(10):
            hist_list = histone.init_genome(percentage=0, a_bool=0)
            hist_list2, _, _ = histone.next_genome(hst_list=hist_list,
                                                   window=10,
                                                   k_nuc=0)
            self.assertEqual(hist_list2[40].status, "u")


if __name__ == 'main':
    unittest.main()
