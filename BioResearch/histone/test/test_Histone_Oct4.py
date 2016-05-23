'''
Created on Mar 18, 2016
@author: eifu
'''
import unittest
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__),  "..",".."))
import histone


class TestHistone_Oct4(unittest.TestCase):
    
    def testInitiate(self):
        tes1 = histone.Histone_Oct4(position=1)
        
        self.assertEqual(tes1.position, 1, "position failure")    

        tes2 = histone.Histone_Oct4(position=2,preNode=tes1)
        self.assertEqual(tes2.position, 2,"position failure")
        self.assertEqual(tes1.nextNode, None, "failure connection1")
        self.assertEqual(tes2.preNode,tes1,"failure connection2")
        
    def testInheritance_M_U_A(self):
        tes1 = histone.MHistone_Oct4(position=0)
        tes2 = histone.UHistone_Oct4(position=-1,nextNode=tes1)
        tes3 = histone.AHistone_Oct4(position=-2,nextNode=tes2)
        
        
        
        self.assertEqual(tes1.CpGislandlist, [0,0], "CpG island list inherited failure ")
        self.assertEqual(tes2.CpGislandlist, [0,0,0,0], "CpG island list inherited failure 2")
        self.assertEqual(tes3.CpGislandlist, [0,0,0,0], "CpG island list inherited failure 3")
        
        self.assertEqual(tes1.status, 'm', 'status failue1')
        self.assertEqual(tes2.status, 'u', 'status failue2')
        self.assertEqual(tes3.status, 'a', 'status failue3')

    def testRandomGenerator(self):
        hstList = histone.createRandomHistoneList_Oct4(NUM_OF_HISTONE=10)
        
        self.assertEqual(len(hstList),10,"failure create appropriate length of list")
        self.assertEqual(hstList[3].CpGislandlist, [0,0,0,0], "failure set up CpG island1")
        self.assertEqual(hstList[4].CpGislandlist, [0,0,0,0], "failure set up CpG island2")
        self.assertEqual(hstList[5].CpGislandlist, [0,0], "failure set up CpG island3")
        self.assertEqual(hstList[6].CpGislandlist, [], "failure set up CpG island4")
        
    
    def testNextGen(self):
        hstList = histone.createRandomHistoneList_Oct4(percentage = 100,
                                                       NUM_OF_HISTONE=10)
        ## hstList contains only methylated histones
        
        dict = histone.nextGen_Oct4(hstList,1,1,2)
        
        hstList2 = dict['list']
        for hst1, hst2 in zip(hstList,hstList2):
            print(str(hst1) +" "+ str(hst1.CpGislandlist)+ "  --> " + str(hst2)+ " "+str(hst2.CpGislandlist))
            self.assertNotEqual(hst2.status, 'a', "something wrong")
        
    
    def testbitvec(self):
        histList = histone.createRandomHistoneList_Oct4()
        bitvec = histone.bitvec_Oct4(histList)
        
        
        self.assertTrue(len(bitvec[0])==81, 'bitvec num failure')
        self.assertTrue(sum(bitvec[0])+sum(bitvec[1])+sum(bitvec[2])==81,'bitvec distri failure')
        self.assertTrue(sum(bitvec[3])==0,'something wrong in CpG island site')
        
if __name__ == 'main':
    unittest.main()