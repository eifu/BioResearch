"""
model5.py
"""


import histone
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib
from model4 import getFirstT
NUM_OF_HISTONE = 81
BEFORE_PROMOTER= 40
WINDOW = 10
TIME1 = 1000
TIME2 = 1000

TESTCASE = 200


def main():
    count = 0
    for R in [0,1]:
        for A in [0,1]:
            for secR in [0,1]:
                for secA in [0,1]:
                    submain(R,A,secR,secA)
                    print(count)
                    count += 1

def submain(R,A,secR,secA):
    plt.style.use('ggplot') 
    font = {'family' : 'serif'}
    matplotlib.rc('font', **font)
    
    TList = np.array([0]*(TIME1+TIME2))
    initialT1 = 0
    initialT2 = 0
    for i in range(TESTCASE):
        dictH = TListMaker(R,A,secR,secA)
        initialT1 += dictH["initialT1"]
        initialT2 += dictH["initialT2"]
        TList = TList + dictH["list"]
        print("   " + str(i))
#     print(TList)
    
    

    
    fig = plt.figure()
    plotT(fig, TList)
    
    TList = TList/TESTCASE

    
    plotTsum5days(fig, TList)
    
    plt.suptitle(r"R:{0}, A:{1} $\rightarrow$ R:{2}, A:{3}".format(R,A,secR,secA) + 
                 "\n"+
                 "first half:{0:.1f} second half:{1:.1f}".format(initialT1/TESTCASE,initialT2/TESTCASE))
    
#     plt.show()

    title = "exp3/exp3_200hist_aveT__R{}A{}__R{}A{}.pdf".format(R,A,secR,secA)
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()

def plotT(fig,TList):
    ax = fig.add_subplot(2,1,2)
    ax.plot(np.arange(TIME1+TIME2),TList,"-",color = "red")
    ax.set_ylim(0,TESTCASE)
    ax.set_xlim(-0.5,TIME1+TIME2+0.5)

def TListMaker(R,A,secR,secA):
    A = 0
    R = 1
    
    secA = 1
    secR = 0
    
    T = 0
    
    histoneList1 = histone.createRandomHistoneList(50, A, NUM_OF_HISTONE, BEFORE_PROMOTER)
    dictH = histone.trackingHist(histoneList1, TIME1, A,R,T,WINDOW)
    histL = dictH["histList"]
    TList = dictH["TList"]
    
    dictH2= histone.trackingHist(histL,TIME2,secA,secR,TList[-1],WINDOW)
    histL = dictH2["histList"]
    TList2 = dictH2["TList"]
    
    initialT1 = getFirstT(TList)
    initialT2 = getFirstT(TList2)
    
    return {"list":np.array(TList+TList2),"initialT1":initialT1,"initialT2":initialT2}
    
def plotTsum5days(fig,TList):
    ax = fig.add_subplot(2,1,1)
    days = []
    acc = 0
    for time in range(TIME1-24,TIME1): ## day1
        acc = acc + TList[time]
    days.append(acc)
    
    acc = 0
    for time in range(TIME1,TIME1+24): ## day2
        acc = acc + TList[time]
    days.append(acc)    

    acc = 0
    for time in range(TIME1+24,TIME1+48): ## day3
        acc = acc + TList[time]
    days.append(acc)    

    acc = 0
    for time in range(TIME1+48,TIME1+72): ## day4
        acc = acc + TList[time]
    days.append(acc)    
    
    acc = 0
    for time in range(TIME1+72,TIME1+96): ## day5
        acc = acc + TList[time]
    days.append(acc)    

    ax.bar(np.arange(5),days,color="red")
    ax.set_xticklabels(("day1","day2","day3","day4","day5"))
    ax.set_ylim(0,24)
    ax.set_yticks((0,4,8,12,16,20,24))
    
    
def getFirstT(list_of_T):
    for i in range(len(list_of_T)):
        if(list_of_T[i]==1):
            return i
    return 0
    


if __name__ == "__main__":
    main()