"""
model4.py

this is an experimental file that uses pandas library and 
plot the histone data according to the environment change.

"""

import histone
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
NUM_OF_HISTONE = 81
BEFORE_PROMOTER= 40
WINDOW = 10
TIME1 = 100
TIME2 = 100

def main():
    #test case (A,R) 1,0 => 0,1
    plt.style.use('ggplot') 
    font = {'family' : 'meiryo'}
    matplotlib.rc('font', **font)
    
    histoneList1 = histone.createRandomHistoneList(50, 1, NUM_OF_HISTONE, BEFORE_PROMOTER)

    tracker = trackingHist(histoneList1, 500, 1,1)
    visualize2(tracker)

def trackingHist(histoneList,time,A,R):
    toBeListOfBitVec = []
    for _ in range(time):        
        toBeListOfBitVec.append(bitvec(histoneList))
        histoneList = nextGen(histoneList, A, R)
    return np.array(toBeListOfBitVec)
    
def nextGen(histoneList,A,R):
    result = [None for _ in range(NUM_OF_HISTONE)]
    start = len(histoneList)//2 - WINDOW//2
    end = len(histoneList)//2 + WINDOW//2
    num_acetylated_in_window = 0
    prev_was_methylated = False
    no_M_in_sequence_in_window = True
    for i in range(len(histoneList)):
        temp_histone = histoneList[i]
        if(start <= i and i<= end):
            if(temp_histone.status == "a"):
                num_acetylated_in_window += 1
                prev_was_methylated = False
            else:
                if(prev_was_methylated == True):
                    no_M_in_sequence_in_window = False
                else:
                    prev_was_methylated = True
        
        temp_histone = temp_histone.k_plus()
        temp_histone = temp_histone.k_ace()
        result[i] = temp_histone.k_minus()
            
    T = A and (num_acetylated_in_window > 5) and no_M_in_sequence_in_window
    """
    WINDOW is size 10(11 histones note that there is E0 between E(-1) and E(1)), 
    so acetylated histones will be dominant if non-acetylated histones are less than 5.
    """
    Eext = ((not T) and (not A)) or R
        
    num_acetylated_in_window = 0
    no_M_in_sequence_in_window = True
    prev_was_methylated = False
    if(Eext == True):
        result[len(histoneList)//2] = histone.MHistone(copy=True,copy_histone=histoneList[len(histoneList)//2])
    return result

def visualize2(list_of_bitvec_list_of_histone):

    for index in range(len(list_of_bitvec_list_of_histone)):
        y_position_m = [i-40 for i in range(NUM_OF_HISTONE) if list_of_bitvec_list_of_histone[index][0][i] == 1]
        x_position_m = np.array([1]*np.sum(list_of_bitvec_list_of_histone[index][0]))
                  
        plt.plot(x_position_m*index,   ## index implies the x coordinates.
                 y_position_m,         
                 ".",color="blue")
        
        y_position_a = [i-40 for i in range(NUM_OF_HISTONE) if list_of_bitvec_list_of_histone[index][1][i] == 1]
        x_position_a = np.array([1]*np.sum(list_of_bitvec_list_of_histone[index][1]))
                  
        plt.plot(x_position_a*index,
                 y_position_a,
                 ".",color="red")
    
    plt.xlim(-0.5,len(list_of_bitvec_list_of_histone)-.5)
    plt.ylim(-40.5,40.5)
    plt.show()

def bitvec(histoneList):
    v_mlist = [1 if h.status == "m" else 0 for h in histoneList]
    v_alist = [1 if h.status == "a" else 0 for h in histoneList]
    v_ulist = [1 if h.status == "u" else 0 for h in histoneList]
    
    return np.array([v_mlist,
                     v_alist,
                     v_ulist],np.int32)
    
    
if __name__ == "__main__":
    main()