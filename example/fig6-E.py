"""
figure 6 C, E

"""

import histone
import histone.figure as figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib
import os


NUM_OF_HISTONE = 81
WINDOW = 10
TIME1 = 504  # 3 week in hour
TIME2 = 504  # 3 week in hour
DELTA = 1

NUMEXAMPLE = 5


def main():
    variation_list_vecgenetimeseries = np.zeros((24, NUMEXAMPLE, TIME2, 3, NUM_OF_HISTONE))
    k_minus = 0.11

    variation_of_kp = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25, 0.3]
    for i, kp in enumerate(variation_of_kp):
        variation_list_vecgenetimeseries[i] = submain(kp, k_minus)
        print("done", kp)

    # TODO write this method in histone.figure
    dump3dim = np.zeros((24 * TIME2, NUM_OF_HISTONE))
    for var, oneversion_list_vecgenetimeseries in enumerate(variation_list_vecgenetimeseries):
        for vecgenetimeseries in oneversion_list_vecgenetimeseries:
            for t in range(TIME2):
                dump3dim[TIME2*var + t] += vecgenetimeseries[t][0]

    filename = "dumpdata__k-{}__{}examples.csv".format(k_minus, NUMEXAMPLE)
    with open(filename, 'wb') as f:
        # f.write(bytes(bytes(i) for i in range(-40, 41)))
        # np.savetxt(f,label,delimiter=',', newline='\n')
        np.savetxt(f,
                   dump3dim,
                   fmt='%d',
                   delimiter=',',
                   newline='\n')

    # day8 = 24*8
    # dumpChart = np.zeros((NUM_OF_HISTONE,24))
    # for pos, _ in enumerate(dumpChart):
    #     # one_pos_after8day = np.zeros(25)
    #     # dumpChart[pos][0] = pos - 40
    #     for var, onevariation_list_vecgenetimeseries in enumerate(variation_list_vecgenetimeseries):
    #         for ex, gene_seq_timeseries in enumerate(onevariation_list_vecgenetimeseries):
    #             dumpChart[pos][var] += gene_seq_timeseries[day8][0][pos]
    #
    #
    # dumpChart2 = np.zeros((NUM_OF_HISTONE,25))
    # for pos, d in enumerate(dumpChart):
    #     dumpChart2[pos][0] = pos-40
    #     dumpChart[pos] =  dumpChart[pos] * 100 / NUMEXAMPLE
    #     dumpChart2[pos][1:] = dumpChart[pos]
    #
    ## TODO write this method in histone.figure
    # # dumpChart = dumpChart * 100 / NUMEXAMPLE
    # # print(dumpChart)
    # filename = "dump_k+~pos_chart_k-{}__{}examples.csv".format(k_minus,NUMEXAMPLE)
    # with open(filename, 'wb') as f:
    #     # f.write(bytes(bytes(i) for i in range(-40, 41)))
    #     # np.savetxt(f,label,delimiter=',', newline='\n')
    #     np.savetxt(f,
    #                dumpChart2,
    #                fmt='%d',
    #                delimiter=',',
    #                newline='\n')

    # plt.style.use('ggplot')
    # font = {'family': 'sans-serif'}
    # matplotlib.rc('font', **font)
    # fig = plt.figure()



    # figure.figure6c_and_6e(fig, variation_list_vecgenetimeseries)
    #
    # plt.show()
    #
    # title = "fig6-CE/fig__{}examples__k-{}.pdf".format(NUMEXAMPLE, k_minus)
    # pp = PdfPages(title)
    # pp.savefig(fig)
    # pp.close()

# TODO write this method in histone.figure
def dumpout(path_num,k_plus,k_minus):
    result = np.zeros((TIME2,81))
    for i in range(NUMEXAMPLE):
        fname = "__"+str(path_num+2)+"__k+" + str(k_plus)+"__k-"+str(k_minus)+"/example"+str(i)+".csv"
        hst_timeseries = histone.read_hstcsv(filename=fname,time=TIME2)


        for t, hst_seq in enumerate(hst_timeseries):
            result[t] += hst_seq[0]

    result_percentage = result/NUMEXAMPLE*100

    savename = "__"+str(path_num+2)+"__k+" + str(k_plus)+"__k-"+str(k_minus)+"/figEC_percentage__{}examples__k+{}_k-{}.csv".format(NUMEXAMPLE,k_plus,k_minus)
    with open(savename, 'wb') as f:
        np.savetxt(f,
                   result_percentage,
                   fmt='%d',
                   delimiter=',',
                   newline='\n')

# TODO write this method in histone.figure
def dumpout2(path_num,k_plus,k_minus):
    result = np.zeros((TIME2,81))
    for i in range(NUMEXAMPLE):
        fname = "__"+str(path_num+2)+"__k+" + str(k_plus)+"__k-"+str(k_minus)+"/example"+str(i)+".csv"
        hst_timeseries = histone.read_hstcsv(filename=fname,time=TIME2)


        for t, hst_seq in enumerate(hst_timeseries):
            result[t] += hst_seq[0]

    result_percentage = result/NUMEXAMPLE*100

    savename = "__"+str(path_num+2)+"__k+" + str(k_plus)+"__k-"+str(k_minus)+"/figEC_percentage__{}examples__k+{}_k-{}.csv".format(NUMEXAMPLE,k_plus,k_minus)
    with open(savename, 'wb') as f:
        np.savetxt(f,
                   result_percentage,
                   fmt='%d',
                   delimiter=',',
                   newline='\n')


def submain( k_plus, k_minus):
    one_variation = np.zeros((NUMEXAMPLE, TIME2, 3, NUM_OF_HISTONE))
    # os.mkdir("__"+str(path_num+2)+"__k+" + str(k_plus)+"__k-"+str(k_minus))
    for ex in range(NUMEXAMPLE):
        one_variation[ex] = subsubmain(k_plus, k_minus)
        print(ex)
    # dumpout(path_num, k_plus, k_minus)
    return one_variation


def subsubmain(k_plus, k_minus):
    R = 0
    A = 1
    secR = 1
    secA = 1

    T = 0

    histoneList1 = histone.init_genome(percentage=50,
                                       a_bool=A,
                                       hst_n=NUM_OF_HISTONE,
                                       kp=k_plus,
                                       kp2=k_plus,
                                       km=k_minus,
                                       ka=k_minus
                                       )

    dictH = histone.track_epigenetic_process(hst_list=histoneList1,
                                             time=TIME1,
                                             a_bool=A,
                                             r_bool=R,
                                             t_bool=T
                                             )
    # tracker = dictH["vectorize"]
    hstL = dictH["hstL"]
    TList = dictH["TList"]

    dictH2 = histone.track_epigenetic_process(hst_list=hstL,
                                              time=TIME2,
                                              a_bool=secA,
                                              r_bool=secR,
                                              t_bool=TList[-1]
                                              )
    tracker2 = dictH2["vectorize"]

    # finalTracker = np.concatenate((tracker, tracker2))
    # histone.save_hst_timeseries(tracker2,
    #                             "__" + str(path_num + 2) + "__k+" + str(k_plus) + "__k-" + str(k_minus) + "/example" + str(ex) + ".csv")

    return tracker2


if __name__ == "__main__":
    main()
