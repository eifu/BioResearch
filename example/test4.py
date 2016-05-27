"""
model4.py

this is an experimental file that uses pandas library and 
plot the histone data according to the environment change.

"""

import histone
from histone.figure import plotHist
from histone.figure import plotWindow
from histone.figure import plotT
from histone.figure import plotStatistics
from histone.figure import plotStatistics2
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib
NUM_OF_HISTONE = 81
WINDOW = 10
TIME1 = 1000
TIME2 = 1000
DELTA = 5


def main():
#     count = 0
#     for R in [0,1]:
#         for A in [0,1]:
#             for secR in [0,1]:
#                 for secA in [0,1]:
#                     submain(R, A, secR, secA)
#                     print(count)
#                     count += 1
#
# def submain(R,A,secR,secA):
    #test case (A,R) 1,0 => 0,1
    R = 0
    A = 1
    secR = 1
    secA =  0
    
    T = 0
    plt.style.use('ggplot') 
    font = {'family' : 'sans-serif'}
    matplotlib.rc('font', **font)
    
    histoneList1 = histone.init_genome(percentage=50, a_bool=A, hst_n=NUM_OF_HISTONE)

    fig = plt.figure()
    
    dictH = histone.track_epigenetic_process(hst_list=histoneList1,
                                             time=TIME1,
                                             a_bool=A,
                                             r_bool=R,
                                             t_bool=T)
    tracker = dictH["vectorize"]
    hstL = dictH["hstL"]
    TList = dictH["TList"]
    
    dictH2 = histone.track_epigenetic_process(hst_list=hstL,
                                              time=TIME2,
                                              a_bool=secA,
                                              r_bool=secR,
                                              t_bool=TList[-1]
                                              )
    tracker2 = dictH2["vectorize"]
    histL = dictH2["hstL"]
    TList2 = dictH2["TList"]
    
    finalTracker = np.concatenate((tracker,tracker2))
    plotHist(fig,finalTracker)
    plotT(fig, TList+TList2)
    plotWindow(fig, finalTracker)
    dictStat1 = plotStatistics(fig, tracker)
    dictStat2 = plotStatistics2(fig,tracker2)
    plotTsum5days(fig, TList+TList2)
    
    plt.suptitle(r"R:{0}, A:{1} $\rightarrow$ R:{2}, A:{3}".format(R,A,secR,secA)
                  + "\n" + 
                  r"AM:{0:.3f} SD:{1:.3f} $\rightarrow$ AM:{2:.3f} SD:{3:.3f}".format(dictStat1["AM"], dictStat1["SD"],dictStat2["AM"], dictStat2["SD"] ))
#     plt.show()
    title = "fig_test4/test4_R{}A{}__R{}A{}.pdf".format(R,A,secR,secA)
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()


def plotTsum5days(fig,TList):
    ax = fig.add_subplot(6,1,5)
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
    
    firstT1 = getFirstT(TList[:TIME1])
    if(firstT1 == -1):
        firstT1 = "none"
    
    firstT2 = getFirstT(TList[TIME1:])
    if(firstT2 == -1):
        firstT2 = "none"
        
    ax.set_xlabel("initial T in the first half is " + str(firstT1)+"\n initial T in the second half is " + str(firstT2))


def getFirstT(TList):
    for i,t in enumerate(TList):
        if t == 1:
            return i
    return -1

if __name__ == "__main__":
    main()
