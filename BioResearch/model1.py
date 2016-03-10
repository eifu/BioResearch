# model1.py

import histone, math
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
"""
this file is from proj5_4.py, with new histone file, 
which checks the 50% domincance
of acetilated histones in the window.

"""

##############################

NUM_OF_HISTONE = 81
BEFORE_PROMOTER = 40
WINDOW = 10
TIME1 = 50
TIME2 = 150
delta = 5
##############################

def submain1(fig,trackerList,R,A):
    """
    this function is to create a graph of all histones status
    """
    bx = fig.add_subplot(311)
    bx.tick_params(left ="off",labelleft="off")

    for h in range(NUM_OF_HISTONE):
        trackerM = [i for i in range(TIME1+TIME2) if trackerList[h][i] == "m"]
        y = [h for i in range(len(trackerM))]
        bx.plot(trackerM,y,",",color = "blue")
        trackerA = [i for i in range(TIME1+TIME2) if trackerList[h][i] == "a"]
        y = [h for i in range(len(trackerA))]
        bx.plot(trackerA,y,",",color = "red")

    bx.set_xlim(-0.5,TIME1+TIME2)
    bx.set_ylabel("histones' status")

    bx.set_ylim(-1,NUM_OF_HISTONE+0.5)
    bx.set_title(r"Percentage of histones: $R={}, A={}$".format(R,A))

def submain2(fig,trackerList):
    """
    this function shows the histones status within the window(= -5 ~ +5 )
    
    """
    bx = fig.add_subplot(613)
    bx.tick_params(left ="off",
                   labelleft="off",
                   bottom="off",
                   labelbottom="off")

    # time = [i for i in range(TIME1,TIME1+TIME2,delta)]
    for h in range(NUM_OF_HISTONE//2- WINDOW//2,NUM_OF_HISTONE//2+ WINDOW//2 + 1):
        trackerM = [i for i in range(TIME1+TIME2) if (trackerList[h][i] == "m" and i >= TIME1) ]
        y = [h - NUM_OF_HISTONE//2 for i in range(len(trackerM))]
        bx.plot(trackerM,y,".",color = "blue")
        trackerA = [i for i in range(TIME1+TIME2) if (trackerList[h][i] == "a" and i>= TIME1) ]
        y = [h - NUM_OF_HISTONE//2 for i in range(len(trackerA))]
        bx.plot(trackerA,y,".",color = "red")
        # print(h- NUM_OF_HISTONE//2 )

    bx.set_xlim(0,TIME1+TIME2)
    bx.set_ylim (-WINDOW//2,WINDOW//2)

def submain3(fig,trackerList):
    """
    this function creates a bar graph of histones within the window
    """
    bx = fig.add_subplot(614)
    time = [i for i in range(TIME1,TIME1+TIME2,delta)]

    list_a1,list_a2 = turnTrackerlistToList_a1a2(trackerList)

    bx.plot(time, list_a1,"-",color="blue")
    bx.plot(time, list_a2,"-",color="red")

    # bx.set_xlim(TIME1,TIME1+TIME2)
    bx.set_ylim(0,10)
    bx.set_ylabel("percentage\n in\n  windown")
    bx.set_xlim(0,TIME1+TIME2)

    bx.grid(True,axis='y')

def submain4(fig,time,TEextTrackerList):
    """
    this function is to create a graph of T status
    """
    cx = fig.add_subplot(615)
    # T_List = T_list_maker(trackerList)
    T_List = [TEextTrackerList[i][0] for i in range(TIME1+TIME2)]
    cx.plot(time,T_List,"-",color="red",drawstyle="steps")
    #cx.scatter(time,T_List)
    cx.set_xlabel("time")
    cx.set_xlim(-0.5,TIME1+TIME2)
    cx.set_ylabel("T")
    cx.set_yticks((0,1))
    cx.set_yticklabels(("Off","On"))
    cx.set_ylim(-0.5,1.5)


def submain5(fig,trackerList):
    bx = fig.add_subplot(616)

    count_m = [0 for i in range(WINDOW+1)]
    count = 0
    for h in range(NUM_OF_HISTONE//2- WINDOW//2,NUM_OF_HISTONE//2+ WINDOW//2+1):
        for t in range(TIME1,TIME1+TIME2,delta):
            if(trackerList[h][t] == "m"):
                count_m[count] = count_m[count] + 1

        count = count + 1
    xaxis = [i for i in range(-WINDOW//2,WINDOW//2+1)]
    bx.barh(xaxis,count_m,align="center")


    #print(count_m)
    acc = 0
    for i in range(WINDOW+1):
        acc = acc + count_m[i]
    AM = acc / len(count_m)

    acc =0
    for i in range(WINDOW+1):
        acc = acc + (count_m[i] - AM)**2
    SD = math.sqrt(acc/len(count_m))

    textstr = '$\mu=%.2f$\n$\sigma=%.2f$'%(AM, SD)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    bx.text(0.75, 0.9, textstr,transform=bx.transAxes,fontsize=10,
        verticalalignment='top', bbox=props)
    # bx.set_xlabel(r"$AM: {}$, $SD: {}$".format(round(AM, 3),round(SD, 3)))


def submain(A,R):
    histoneList = histone.createRandomHistoneList(A=A)
    T = 1
    Eext = 0
    percentage = 50
    trackerList,TEextTrackerList  = histone.trackingHistone(histoneList=histoneList,R=R,A=A,T=T,Eext=Eext,TIME=TIME1+TIME2)

    return trackerList, TEextTrackerList

def main():
    A = 1
    R = 1
    trackerList, TEextTrackerList = submain(A, R)
    fig= pyplot.figure()
    time = np.linspace(0,TIME1+TIME2-1,TIME1+TIME2)


    submain1(fig, trackerList,R,A)
    submain2(fig, trackerList)
    submain3(fig, trackerList)
    submain4(fig, time, TEextTrackerList)
    submain5(fig, trackerList)

    # title = "exp1/exp1_A{}_R{}_num{}.pdf".format(A,R,num)
    # pp = PdfPages(title)
    # pp.savefig(fig)
    # pp.close()

    pyplot.show()


# def main():
#     for i in range(20):
#         for A in range(2):
#             for R in range(2):
#                 main1(i,A,R)
#         print(i)

def turnTrackerlistToList_a1a2(trackerList):
    list_a1 = [0 for i in range(TIME1,TIME1+TIME2,delta)] # list of # of methylated histones
    list_a2 = [0 for i in range(TIME1,TIME1+TIME2,delta)] # list of # of acetilated histones
    for i in range(NUM_OF_HISTONE//2- WINDOW//2,NUM_OF_HISTONE//2 + WINDOW//2 +1):
        counter = 0
        for t in range(TIME1,TIME1+TIME2,delta):
            if(trackerList[i][t] =="m"):
                list_a1[counter] += 1
            elif(trackerList[i][t] =="a"):
                list_a2[counter] += 1
            counter += 1
    return list_a1, list_a2

if __name__ == "__main__":
    main()
