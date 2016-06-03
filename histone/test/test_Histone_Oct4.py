'''
Created on Mar 18, 2016
@author: eifu
'''
import unittest
import sys
import os

sys.path.append(os.path.join(os.path.dirname(__file__), "..", ".."))
import histone


class TestHistone_Oct4(unittest.TestCase):
    def testInitiate(self):
        tes1 = histone.HistoneOct4(position=1)

        self.assertEqual(tes1.position, 1, "position failure")

        tes2 = histone.HistoneOct4(position=2, prenode=tes1)
        self.assertEqual(tes2.position, 2, "position failure")
        self.assertEqual(tes1.nextnode, None, "failure connection1")
        self.assertEqual(tes2.prenode, tes1, "failure connection2")

    def testInheritance_M_U_A(self):
        tes1 = histone.MHistoneOct4(position=0)
        tes2 = histone.UHistoneOct4(position=-1, nextnode=tes1)
        tes3 = histone.AHistoneOct4(position=-2, nextnode=tes2)

        self.assertEqual(tes1.CpGislandlist, [0, 0], "CpG island list inherited failure ")
        self.assertEqual(tes2.CpGislandlist, [0, 0, 0, 0], "CpG island list inherited failure 2")
        self.assertEqual(tes3.CpGislandlist, [0, 0, 0, 0], "CpG island list inherited failure 3")

        self.assertEqual(tes1.status, 'm', 'status failure1')
        self.assertEqual(tes2.status, 'u', 'status failure2')
        self.assertEqual(tes3.status, 'a', 'status failure3')

    def testRandomGenerator(self):
        hstList = histone.init_genome_oct4(hst_n=10)

        self.assertEqual(len(hstList), 10, "failure create appropriate length of list")
        self.assertEqual(hstList[3].CpGislandlist, [0, 0, 0, 0], "failure set up CpG island1")
        self.assertEqual(hstList[4].CpGislandlist, [0, 0, 0, 0], "failure set up CpG island2")
        self.assertEqual(hstList[5].CpGislandlist, [0, 0], "failure set up CpG island3")
        self.assertEqual(hstList[6].CpGislandlist, [], "failure set up CpG island4")

    def testNextGen(self):
        hstList = histone.init_genome_oct4(percentage=100,
                                                       hst_n=10)
        dict = histone.next_genome_oct4(hstList, 1, 1, 2, 0.01)

        hstList2 = dict['hstL']
        for hst1, hst2 in zip(hstList, hstList2):
            print(str(hst1) + " " + str(hst1.CpGislandlist) + "  --> " + str(hst2) + " " + str(hst2.CpGislandlist))
            # todo it should not do multiple transaction like methylated histone to acetylated histone
        #             self.assertNotEqual(hst2.status, 'a', "something wrong")

    def testvectorize(self):
        histList = histone.init_genome_oct4()
        vect = histone.vectorize_oct4(histList)
        self.assertTrue(len(vect[0]) == 81, 'vect num failure')
        self.assertTrue(sum(vect[0]) + sum(vect[1]) + sum(vect[2]) == 81, 'vect distri failure')
        self.assertTrue(sum(vect[3]) == 0, 'something wrong in CpG island site')

    def testTrackingHist(self):
        hstL = histone.init_genome_oct4(hst_n=10)
        dictH = histone.track_epigenetic_process_oct4(hstL,
                                          20,  # time
                                          1,  # activator
                                          1,  # repressor
                                          0,  # transcription
                                          2,  # window
                                          0.001  # p_off
                                          )
        bit = dictH['vectorize']
        for b in bit:
            print(str(b[0]) + "  " + str(b[1]) + "  " + str(b[2]))

    def testRandom(self):
        histList = histone.init_genome_oct4()
        dictH = histone.next_genome_oct4(histList, 0, 1, 10, p_off=0.01)
        histList2 = dictH["hstL"]
        for hist1, hist2 in zip(histList, histList2):
            print(str(hist1) + "   --> " + str(hist2))
        self.assertEqual(dictH["Eext"], 1, 'R1 does not work correctly')

if __name__ == 'main':
    unittest.main()
