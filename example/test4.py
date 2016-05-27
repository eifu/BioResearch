"""
model4.py

this is an experimental file that uses pandas library and 
plot the histone data according to the environment change.

"""

import histone
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import matplotlib
import math
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
    
    dictH = histone.track_epigenetic_process(hst_list=histoneList1,time=TIME1,a_bool=A,r_bool=R,t_bool=T)
    tracker = dictH["vectorize"]
    hstL = dictH["hstL"]
    TList = dictH["TList"]
    
    dictH2 = histone.track_epigenetic_process(hst_list=hstL,time=TIME2,a_bool=secA,r_bool=secR,t_bool=TList[-1])
    tracker2 = dictH2["vectorize"]
    histL = dictH2["hstL"]
    TList2 = dictH2["TList"]
    
    finalTracker = np.concatenate((tracker,tracker2))
    plotHist(fig,finalTracker)
    plotT(fig, TList+TList2)
    plotWindow(fig, finalTracker)
    dictStat1 = plotStatistics(fig, finalTracker)
    dictStat2 = plotStatistics2(fig,finalTracker)
    plotTsum5days(fig, TList+TList2)
    
    plt.suptitle(r"R:{0}, A:{1} $\rightarrow$ R:{2}, A:{3}".format(R,A,secR,secA)
                  + "\n" + 
                  r"AM:{0:.3f} SD:{1:.3f} $\rightarrow$ AM:{2:.3f} SD:{3:.3f}".format(dictStat1["AM"], dictStat1["SD"],dictStat2["AM"], dictStat2["SD"] ))
#     plt.show()
    title = "fig_test4/test4_R{}A{}__R{}A{}.pdf".format(R,A,secR,secA)
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()

def plotStatistics(fig,l_of_bitvec):
    bx = fig.add_subplot(9,4,22)

    count_m = [0 for _ in range(WINDOW+1)]
    count = 0
    for h in range(NUM_OF_HISTONE//2- WINDOW//2,NUM_OF_HISTONE//2+ WINDOW//2+1):
        for time in range(TIME1//2,TIME1,DELTA):
            if(l_of_bitvec[time][0][h] == 1):
                count_m[count] = count_m[count] + 1

        count = count + 1
    xaxis = [i for i in range(-WINDOW//2,WINDOW//2+1)]
    bx.barh(xaxis,count_m,align="center")
    bx.set_xticks([])
    bx.set_yticks([])


    print(count_m)
    acc = 0
    for i in range(WINDOW+1):
        acc = acc + count_m[i]
    AM = acc / len(count_m)

    acc =0
    for i in range(WINDOW+1):
        acc = acc + (count_m[i] - AM)**2
    SD = math.sqrt(acc/len(count_m))
    return {"AM":AM,"SD":SD}



def plotStatistics2(fig,l_of_bitvec):
    bx = fig.add_subplot(9,4,24)

    count_m = [0 for _ in range(WINDOW+1)]
    count = 0
    for h in range(NUM_OF_HISTONE//2- WINDOW//2,NUM_OF_HISTONE//2+ WINDOW//2+1):
        for time in range(TIME1+(TIME2//2),TIME1+TIME2,DELTA):
            if(l_of_bitvec[time][0][h] == 1):
                count_m[count] = count_m[count] + 1

        count = count + 1
    xaxis = [i for i in range(-WINDOW//2,WINDOW//2+1)]
    bx.barh(xaxis,count_m,align="center")
    bx.set_xticks([])
    bx.set_yticks([])

    acc = 0
    for i in range(WINDOW+1):
        acc = acc + count_m[i]
    AM = acc / len(count_m)

    acc =0
    for i in range(WINDOW+1):
        acc = acc + (count_m[i] - AM)**2
    SD = math.sqrt(acc/len(count_m))


    return {"AM":AM,"SD":SD}

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
    

def plotT(fig,list_of_T):
    ax = fig.add_subplot(9,1,4)
    time = [i for i in range(len(list_of_T))]
    ax.plot(time, list_of_T,"-", color="red")
    ax.set_yticks([])
    ax.set_ylim(-0.5,1.5)
    ax.set_xlim(-0.5,TIME1+TIME2+0.5)
    ax.set_xticks([])


def getFirstT(list_of_T):
    for i in range(len(list_of_T)):
        if(list_of_T[i]==1):
            return i
    return -1


def plotWindow(fig, l_of_bitvec):
    ax = fig.add_subplot(9,1,5)
    for time in range(len(l_of_bitvec)):
        y_position_m = [i-40 for i in range(NUM_OF_HISTONE//2-WINDOW//2,NUM_OF_HISTONE//2+WINDOW//2) if l_of_bitvec[time][0][i]==1]
        x_position_m = np.array([1]*len(y_position_m))
        
        ax.plot(x_position_m*time,
                y_position_m,
                ",",color="blue")
        
        y_position_a = [i-40 for i in range(NUM_OF_HISTONE//2-WINDOW//2,NUM_OF_HISTONE//2+WINDOW//2) if l_of_bitvec[time][2][i]==1]
        x_position_a = np.array([2]*len(y_position_a))
        
        ax.plot(x_position_a*time,
                y_position_a,
                ",",color="red")
    
    ax.set_xlim(-0.5,len(l_of_bitvec))
    ax.set_ylim(-5,5)
    ax.set_yticks([])
    ax.set_xticks([])


def plotHist(fig,list_of_bitvec_of_histoneList):
    
    ax = fig.add_subplot(3,1,1)
    for time in range(len(list_of_bitvec_of_histoneList)):
        y_position_m = [i-40 for i in range(NUM_OF_HISTONE) if list_of_bitvec_of_histoneList[time][0][i] == 1]
        x_position_m = np.array([1]*np.sum(list_of_bitvec_of_histoneList[time][0]))
                  
        ax.plot(x_position_m*time,   ## time implies the x coordinates.
                 y_position_m,         
                 ",",color="blue")
        
        y_position_a = [i-40 for i in range(NUM_OF_HISTONE) if list_of_bitvec_of_histoneList[time][2][i] == 1]
        x_position_a = np.array([1]*np.sum(list_of_bitvec_of_histoneList[time][2]))
                  
        ax.plot(x_position_a*time,
                 y_position_a,
                 ",",color="red")
    
    ax.set_xlim(-0.5,len(list_of_bitvec_of_histoneList)-0.5)
    ax.set_ylim(-40.5,40.5)

if __name__ == "__main__":
    main()
