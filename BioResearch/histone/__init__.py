__author__ = 'eifu'


from numpy.random import sample
import numpy as np

class Histone(object):
    """
    Histone is a Python module designed for researching
    methylation modeling. This class allows users to create
    a histone object.
    """


    def __init__(self,position=0,
                 K_PLUS     =   0.176,
                 K_PLUS2    =   0.17,
                 K_MINUS    =   0.117,
                 A_bool     =   False,
                 K_ACE      =   0.12,
                 nextNode=None,
                 preNode=None,
                 percentage=50,
                 copy=False,
                 copy_histone=None):
        """
        Initialize the histone. The position must be
        specified, but K_PLUS, K_MINUS, nextNode, preNode,
        and percentage can be specified. If they are not,
        default values variables will be assigned to those
        instance variables.

        position  ... a position of histone, to differentiate
                  from other histone instances.

        A_bool ... a boolean value, activator of transcription factors
                   if A is true, then the K_ACE is the given paramater.
                   else, K_ACE is 0, which means no histones can be turned to be acetilated.
        preNode ... a neighbor of the instances, judged by position - 1
        nextNode ... a neighbor of the instances, judged by position + 1
        """

        if(copy):   ## copy the info of histone object
            self.position = copy_histone.position

            self.preNode = copy_histone.preNode  # doubly-linked
            if(copy_histone.preNode != None): copy_histone.preNode.nextNode = self 

            self.nextNode = copy_histone.nextNode # doubly-linked
            if(copy_histone.nextNode !=None): copy_histone.nextNode.preNode = self 
                       
                

        else:
            self.position = position
            self.preNode = preNode
            self.nextNode = nextNode
            if(A_bool):
                Histone.K_ACE = K_ACE
            else:
                Histone.K_ACE = 0
            Histone.K_PLUS = K_PLUS
            Histone.K_PLUS2 = K_PLUS2
            Histone.K_MINUS = K_MINUS
            

            
            
    def set_adjHistone(self,nextNode):
        """
        Set the input to self.nextNode, one of its instance variables.
        """
        self.nextNode = nextNode
        nextNode.preNode = self
    
    def set_K_ACE(self,A):
        if(A==True):
            Histone.K_ACE = 0.12 # K_ace default value
        else:
            Histone.K_ACE = 0
    
    def DNAmethylation(self):
        return self
    
    def k_plus(self):
        """
        the unmethylated histone will get methylated by K_PLUS probability
        return methylated object if the histone will get methylated
        """
        return self

    def k_minus(self):
        """
        the methylated histone will be get unmethylated by K_MINUS probability
        or
        the acetilated histone will be get unacetylated by K_MINUS probability
        return unmethylated histone object if the histone will get unmethylated
        """
        
        if(sample()<Histone.K_MINUS):return UHistone(copy=True,copy_histone=self)
        return self


    def k_ace(self):
        """
        the unmethylated hitone will be get acetilated by K_ACE probability
        return acetilated histone if the histone will get acetilated
        """
        return self

    def __lt__(self, other):
        """
        return true when the position of self is smaller than that of other
        """
        return self.position < other.position

    def __str__(self):
        """
        return all info about the instance of this object.
        """
        st=''
        if(self.status=='m'):
            st = "methylated"
        elif(self.status=='a'):
            st = "acetilated"
        else:
            st = "unmethylated"
        sentence = "pos: {}  status: {}\t".format(self.position,st)
        return sentence

    def display(self):
        sentence = "address: {}\nposition: {}\nstatus: {}\n".format(self,self.position,self.status)

        try:
            sentence += "Prev Histone:{}\n".format(self.preNode)
        except(AttributeError):
            sentence += "Prev Histone: None\n"
        try:
            sentence += "Next Histone: {}".format(self.nextNode)
        except(AttributeError):
            sentence += "Next Histone: None"
        finally:
            print(sentence)
            print("###\n")


class MHistone(Histone):
    def __init__(self,**kwarg):
        super().__init__(**kwarg)
        self.status="m"
    
        

class UHistone(Histone):
    def __init__(self,**kwarg):
        super().__init__(**kwarg)
        self.status="u"

    
    def k_plus(self):
        if(self.preNode != None and self.preNode.status == "m" and sample() < Histone.K_PLUS): return MHistone(copy=True,copy_histone=self)
        if(self.nextNode != None and self.nextNode.status == "m" and sample() < Histone.K_PLUS): return MHistone(copy=True,copy_histone=self)
        return self

    def k_minus(self): return self

    def k_ace(self):
        if(sample()<Histone.K_ACE):return AHistone(copy=True,copy_histone=self)
        return self



class AHistone(Histone):
    def __init__(self,**kwarg):
        super().__init__(**kwarg)
        self.status="a"
      
            
    def k_plus(self):
        if(self.preNode != None and self.preNode.status == "m" and sample() < Histone.K_PLUS2):return UHistone(copy=True,copy_histone=self)
        if(self.nextNode != None and self.nextNode.status == "m" and sample() < Histone.K_PLUS2):return UHistone(copy=True,copy_histone=self)
        return self


class Histone_Oct4(Histone):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        if(kwargs.pop('copy',None)==True):
            self.CpGislandlist = kwargs.pop('copy_histone',None).CpGislandlist
        else:
            position = kwargs.pop('position',0)
            if position ==-2:
                self.CpGislandlist = [0,0,0,0]
            elif position == -1:
                self.CpGislandlist = [0,0,0,0]
            elif position == 0:
                self.CpGislandlist = [0,0]
            else:
                self.CpGislandlist = []

class MHistone_Oct4(Histone_Oct4):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.status = "m"
        
    def k_minus(self):
        if(sample()<Histone.K_MINUS):return UHistone_Oct4(copy=True,copy_histone=self)
        return self

class UHistone_Oct4(Histone_Oct4):
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.status = "u"
    def k_plus(self):
        if(self.preNode != None and self.preNode.status == "m" and sample() < Histone.K_PLUS): return MHistone_Oct4(copy=True,copy_histone=self)
        if(self.nextNode != None and self.nextNode.status == "m" and sample() < Histone.K_PLUS): return MHistone_Oct4(copy=True,copy_histone=self)
        return self

    def k_minus(self): return self

    def k_ace(self):
        if(1 in self.CpGislandlist):
            return self
        """
        if histone has CpG island on it, then it wont be accetilated.
        """
        
        if(sample()<Histone.K_ACE):
            return AHistone_Oct4(copy=True,copy_histone=self)
        else:
            return self

    def DNAmethylation(self):
        if( sample() < sum(self.CpGislandlist)*Histone.K_PLUS): 
            return MHistone_Oct4(copy=True,copy_histone=self)
        else:   
            return self
        
class AHistone_Oct4(Histone_Oct4):
    
    def __init__(self,**kwargs):
        super().__init__(**kwargs)
        self.status = "a"
        
    def k_plus(self):
        
        if(self.preNode != None and 
           self.preNode.status == "m" and sample() < Histone.K_PLUS2):
            return UHistone_Oct4(copy=True,copy_histone=self)
        
        if(self.nextNode != None and
            self.nextNode.status == "m" and sample() < Histone.K_PLUS2):
            return UHistone_Oct4(copy=True,copy_histone=self)
        
        return self

def createRandomHistoneList(percentage=50,
                            A=1,
                            NUM_OF_HISTONE=81,
                            K_PLUS = 0.176,
                            K_PLUS2 = 0.17,
                            K_MINUS = 0.117,
                            K_ACE = 0.12):
    """
    percentage ... the probability of having methylated hitone.
    this method returns a list of histone randomly generated with respect to
    the inputs.
    """
    before_promoter = NUM_OF_HISTONE//2
    dstList = [] ## dstList stores histones
    ratio = percentage/100  ## ratio should be float number between 0 and 1

    for i in range(NUM_OF_HISTONE):

        if(sample() < ratio):
            dstList.append(MHistone(position=i-before_promoter,
                                    K_PLUS=K_PLUS,
                                    K_PLUS2=K_PLUS2,
                                    K_MINUS=K_MINUS,
                                    K_ACE=K_ACE,
                                    A_bool=A,
                                    percentage=percentage))
        else:
            dstList.append(AHistone(position=i-before_promoter,
                                     K_PLUS=K_PLUS,
                                     K_PLUS2=K_PLUS2,
                                     K_MINUS=K_MINUS,
                                     K_ACE=K_ACE,
                                     A_bool=A,
                                     percentage=percentage))

        dstList[i-1].set_adjHistone(dstList[i])
    
    dstList[0].preNode = None # disjoint the edge histone to itself.
    return dstList

        
def createRandomHistoneList_Oct4(percentage=50,
                            A=1,
                            NUM_OF_HISTONE=81,
                            BEFORE_PROMOTER=40,
                            K_PLUS = 0.176,
                            K_PLUS2 = 0.17,
                            K_MINUS = 0.117,
                            K_ACE = 0.12):
    """
    this method is a modified version of createRandomHistoneList
    used for Oct4 histones.
    """
    before_promoter = NUM_OF_HISTONE//2
    hstList = [] ## dstList stores histones
    ratio = percentage/100  ## ratio should be float number between 0 and 1

    for i in range(NUM_OF_HISTONE):

        if(sample() < ratio):
            hstList.append(MHistone_Oct4(position=i-before_promoter,
                                         K_PLUS=K_PLUS,
                                         K_PLUS2=K_PLUS2,
                                         K_MINUS=K_MINUS,
                                         K_ACE=K_ACE,
                                         A_bool=A,
                                         percentage=percentage)
                           )
        else:
            hstList.append(AHistone_Oct4(position=i-before_promoter,
                                         K_PLUS=K_PLUS,
                                         K_PLUS2=K_PLUS2,
                                         K_MINUS=K_MINUS,
                                         K_ACE=K_ACE,
                                         A_bool=A,
                                         percentage=percentage)
                           )

        hstList[i-1].set_adjHistone(hstList[i])
    
    hstList[0].preNode = None
    return hstList

def nextGen(histoneList, A, R, window):
    """
    this method takes histone list and returns the next generation out of them.
    """
    result = []
    num_acetylated_in_window = 0
    num_methylated_in_window = 0

    for hist in histoneList:
        
        if(-window//2 <= hist.position and hist.position<= window//2):
            if(hist.status == "a"):
                num_acetylated_in_window += 1
            elif(hist.status == "m"):
                num_methylated_in_window += 1
                
                            
        
        hist = hist.k_minus()
        hist = hist.k_ace()
        hist = hist.k_plus()
        
        result.append(hist)
            
    T = A and (num_acetylated_in_window > 5) 
    """
    WINDOW is size 10(11 histones note that there is E0 between E(-1) and E(1)), 
    so acetylated histones will be dominant if non-acetylated histones are less than 5.
    """
    if(R == 1):    
        Eext = (not T)  
    else:
        Eext = num_methylated_in_window > 2;
        
        
    if(Eext == True):
        result[len(histoneList)//2] = MHistone(copy=True,copy_histone=result[len(histoneList)//2])
        
    return {"list":result,"T":T,"Eext":Eext}

def nextGen_Oct4(histoneList,A,R,window,p_off=0.001179):
    """
    this method takes histone list and returns the next generation of them.
    """
    result = []
    num_acetylated_in_window = 0
    num_methylated_in_window = 0

    for hist in histoneList:
        
        if(-5 <= hist.position and hist.position<= 5):
            if(hist.status == "a"):
                num_acetylated_in_window += 1
            elif(hist.status == "m"):
                num_methylated_in_window += 1
                
            for index in range(len(hist.CpGislandlist)):
                
                if(hist.status=='m' and sample()<0.001133): # p on probability
                    hist.CpGislandlist[index] = 1
                    
                if(sample()<p_off): # p off probability
                    hist.CpGislandlist[index] = 0   
                            
        
        hist = hist.k_minus()
        hist = hist.k_ace()
        hist = hist.k_plus()
        
        hist = hist.DNAmethylation()

        result.append(hist)
            
    T = A and (num_acetylated_in_window > 5) 
    """
    WINDOW is size 10(11 histones note that there is E0 between E(-1) and E(1)), 
    so acetylated histones will be dominant if non-acetylated histones are less than 5.
    """
    if(R == 1):    
        Eext = (not T)  
    else:
        Eext = num_methylated_in_window > 2;
        
        
    if(Eext == True):
        result[len(histoneList)//2] = MHistone_Oct4(copy=True,copy_histone=result[len(histoneList)//2])
        
    return {"list":result,"T":T,"Eext":Eext}

def bitvec(histoneList):
    """
    this method takes a list of histone objects, and returns three dimention numpy array
    in which three bit vectors are stored. 
    """
    v_mlist = np.array([1 if h.status == "m" else 0 for h in histoneList])
    v_ulist = np.array([1 if h.status == "u" else 0 for h in histoneList])
    v_alist = np.array([1 if h.status == "a" else 0 for h in histoneList])
    
    return np.array([v_mlist,
                     v_ulist,
                     v_alist],
                    np.int32)

def bitvec_Oct4(histoneList):
    v_mlist = np.array([1 if h.status == "m" else 0 for h in histoneList])
    v_ulist = np.array([1 if h.status == "u" else 0 for h in histoneList])
    v_alist = np.array([1 if h.status == "a" else 0 for h in histoneList])
    v_cpg = [sum(h.CpGislandlist) for h in histoneList]

    return np.array([v_mlist,
                     v_ulist,
                     v_alist,
                     v_cpg],
                    np.int32)    

def trackingHist(histoneList,time,A,R,T,window):
    for i in range(len(histoneList)):
                histoneList[i].set_K_ACE(A)
    toBeListOfBitVec = []
    toBeListOfT = []
    for _ in range(time):        
        toBeListOfBitVec.append(bitvec(histoneList))
        toBeListOfT.append(T)
        dictH = nextGen(histoneList, A, R,window)
        histoneList = dictH["list"]
        T = dictH["T"]

    return {"bitvec":np.array(toBeListOfBitVec),"histList":histoneList,"TList":toBeListOfT}

def getTimeDiedOut(ls_of_bitvec):
    
    for time in range(len(ls_of_bitvec)):
        if(ls_of_bitvec[time][0][40]==0):
            for else_in_window in range(35,45):
                if(ls_of_bitvec[time][0][else_in_window]==1):
                    break
                return time
    return -1