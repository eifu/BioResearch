'''
Created on Mar 18, 2016

@author: eifu
'''
import unittest
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
        self.assertEqual(histList[2].preNode, histList[1],'connection failure 2-1')
        self.assertEqual(histList[1].preNode, histList[0],'connection failure 2-2')
    
    def testRandomHistListMethod_disconnection(self):
        histList = histone.createRandomHistoneList(50, 0, 3, 1)
        
        self.assertEqual(histList[0].preNode,None,'disconnection Head and Tail')
        self.assertEqual(histList[2].nextNode, None, 'disconnection Head and Tail 2')
        
        
        histList2 = histone.createRandomHistoneList(50,0,40,1)
        self.assertEqual(histList2[0].preNode,None,'disconnection Head and Tail 3')
        self.assertEqual(histList2[-1].nextNode, None, 'disconnection Head and Tail 4')

if __name__ == 'main':
    unittest.main()