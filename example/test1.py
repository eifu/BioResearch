# model2.py

import histone
import math
from numpy.random import sample
import matplotlib.pyplot as pyplot
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
"""
this file is from nidek3.py, with new histone file, 
which creates a 4 by 4 plots graphs of every case 
of A and R.

"""

##############################

NUM_OF_HISTONE = 81
BEFORE_PROMOTER = 40
WINDOW = 10
TIME1 = 1000
TIME2 = 1000
delta = 5

##############################


def main():
    testnum = 0
    fig= pyplot.figure()
    # for A in range(2):
    #     for R in range(2):
    #         submain(fig,A,R,testnum)
    #         testnum += 1
    #         print(A*50+R*25)

    submain(fig,1,1,1)

    # fig.text(0.03,0.19,"$R+,A+$")
    # fig.text(0.03,0.39,"$R-,A+$")
    # fig.text(0.03,0.59,"$R+,A-$")
    # fig.text(0.03,0.79,"$R-,A-$")
    #
    # fig.text(0.19,0.95,"$R-,A-$")
    # fig.text(0.39,0.95,"$R+,A-$")
    # fig.text(0.59,0.95,"$R-,A+$")
    # fig.text(0.79,0.95,"$R+,A+$")
    
    
    # title = "figure1/fig1_histones.pdf"
    # pp = PdfPages(title)
    # pp.savefig(fig)
    # pp.close()

    pyplot.show()


def submain(fig,A,R,num):
    testnum = num * 4
    for secA in range(2):
        for secR in range(2):
            testnum += 1
            subsubmain(fig,A,R,secA,secR,testnum)
            print("   done")

def subsubmain(fig,A,R,secondA,secondR,testnum):
    dicth = setHistones(A,R)
    vectorize, TTrackerList = dicth['vectorize'], dicth['TList']
    time = np.linspace(0,TIME1+TIME2-1,TIME1+TIME2)
    
#     submain1(fig, trackerList,R,A,secondR,secondA)
#     submain4(fig,time,finalTEextTrackerList,avg_first_T_on,testnum)
    submain1(fig,vectorize,R,A,testnum)

#     pyplot.show()
    
def setHistones(A,R):
    histoneList = histone.init_genome(a_bool=A)
    T = 0
    Eext = 0
    return histone.track_epigenetic_process(hst_list=histoneList,
                                            a_bool=A,
                                            r_bool=R,
                                            t_bool=T,
                                            time=TIME1
                                            )
    
def submain1(fig,vectorize,R,A,num):
    """
    this function is to create a graph of all histones status
    """
    bx = fig.add_subplot(4,4,num)
    bx.tick_params(left ="off",labelleft="off")

    for h in range(NUM_OF_HISTONE):
        trackerM = [i for i in range(TIME1) if vectorize[i][0][h] == 1]
        y = [h for i in range(len(trackerM))]
        bx.plot(trackerM,y,",",color = "blue")
        trackerA = [i for i in range(TIME1) if vectorize[i][2][h] == 1]
        y = [h for i in range(len(trackerA))]
        bx.plot(trackerA,y,",",color = "red")

    bx.set_xlim(-0.5,TIME1+TIME2)
    bx.set_xticks([0,TIME1,TIME2])
#     bx.set_ylabel("histones' status")

    bx.set_ylim(-1,NUM_OF_HISTONE+0.5)
#     bx.set_title(r"Percentage of histones: $R={}, A={} \to R={}, A = {}$".format(R,A,secondR,secondA))

def submain4(fig,time,TEextTrackerList, ave,num):
    """
    this function is to create a graph of T status
    """
    cx = fig.add_subplot(4,4,num)
    # T_List = T_list_maker(trackerList)
    T_List = [TEextTrackerList[i][0] for i in range(TIME1+TIME2)]
    cx.plot(time,T_List,"-",color="red",drawstyle="steps")
    #cx.scatter(time,T_List)
    cx.set_xlabel("time")
    cx.set_xticks((0,TIME1,TIME2))
    cx.set_xlim(-0.5,TIME1+TIME2)
    cx.set_yticks([])
#     cx.set_yticklabels(("O%","100%"))
    cx.set_ylim(-0.5,1+.5)
    textstr = 'average:\n %.2f'%(ave)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    cx.text(0.5, 0.9, textstr,transform=cx.transAxes,fontsize=10,
        verticalalignment='top', bbox=props)

def index_First_T_on(list):
    for i in range(len(list)):
        if(list[i][0] != 0):
            return i
    return 0


        

if __name__ == "__main__":
    main()

