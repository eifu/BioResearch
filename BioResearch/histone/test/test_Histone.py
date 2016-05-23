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
        m2 = histone.MHistone(position=2,preNode=m1)
        
        self.assertEqual(m2.preNode.position, 1, 'initiate connection failure')
        self.assertEqual(m1.nextNode, None, 'initiate connection failure 2')            
        
    def testset_adjHistone(self):
        m1 = histone.MHistone(position=1)
        m2 = histone.MHistone(position=2,preNode=m1)
        m1.set_adjHistone(m2)
        
        self.assertEqual(m2.preNode.position, 1, 'initiate connection failure')
        self.assertNotEqual(m1.nextNode, None, 'initiate connection failure 2')   
        self.assertEqual(m1.nextNode, m2, 'initiate connection failure 3')
    
    def testRandomHistListMethod(self):
        histList = histone.createRandomHistoneList(50, 0, 3, 1)
        
        
        self.assertEqual(len(histList), 3, "number of elements in list failure")
        self.assertEqual(histList[0].nextNode, histList[1], 'connection failure 1-1')
        self.assertEqual(histList[1].nextNode, histList[2], 'connection failure 1-2')
        self.assertEqual(histList[2].preNode, histList[1], 'connection failure 2-1')
        self.assertEqual(histList[1].preNode, histList[0], 'connection failure 2-2')
    
    def testRandomHistListMethod_disconnection(self):
        histList = histone.createRandomHistoneList(50, 0, 3, 1)
        
        self.assertEqual(histList[0].preNode,None,'disconnection Head and Tail')
        self.assertEqual(histList[2].nextNode, None, 'disconnection Head and Tail 2')
        
        
        histList2 = histone.createRandomHistoneList(50,0,40,1)
        self.assertEqual(histList2[0].preNode,None,'disconnection Head and Tail 3')
        self.assertEqual(histList2[-1].nextNode, None, 'disconnection Head and Tail 4')

    def testbitvec(self):
        histList = histone.createRandomHistoneList()
        bitvec = histone.bitvec(histList)
        self.assertTrue(len(bitvec[0])==81, 'bitvec num failure')
        self.assertTrue(sum(bitvec[0])+sum(bitvec[1])+sum(bitvec[2])==81,'bitvec distri failure')
        hs2 = histone.createRandomHistoneList(percentage=100)
        bitvec2 = histone.bitvec(hs2)
        self.assertTrue(sum(bitvec2[0])==81,'methylated histone list failure')
        hs3 = histone.createRandomHistoneList(percentage=0)
        bitvec3 = histone.bitvec(hs3)
        self.assertTrue(sum(bitvec3[2])==81,'acetylated histone list failure')

    def testNextGen(self):
        histList = histone.createRandomHistoneList()
        dictH = histone.nextGen(histList, 0, 1, 10)
        histList2 = dictH["list"]
        for hist1, hist2 in zip(histList,histList2):
            print(str(hist1) +"   --> " + str(hist2))
        self.assertEqual(dictH["Eext"], 1, 'R1 does not work correctly')
        
    def testNextGen2(self):
        m1 = histone.MHistone(position=-2)
        m2 = histone.AHistone(position=-1)
        m3 = histone.AHistone(position=0)
        m4 = histone.AHistone(position=1)
        m5 = histone.AHistone(position=2)
        histL = [m1,m2,m3,m4,m5]
        
        dictH = histone.nextGen(histL, 1, 1 , 3)
        histL2 = dictH["list"]
#         self.assertNotEqual(histL2[2].status, 'm', 'K+ will not go beyond one neighbor at one time')
#         self.assertNotEqual(histL2[3].status, 'm', 'K+ will not go beyond one neighbor at one time')
#         self.assertNotEqual(histL2[4].status, 'm', 'K+ will not go beyond one neighbor at one time')
        
        
if __name__ == 'main':
    unittest.main()