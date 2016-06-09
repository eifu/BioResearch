"""
figure 6 C, E

"""
import histone
import histone.figure as figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os



NUM_OF_HISTONE = 81
WINDOW = 10
TIME1 = 504  # 3 week in hour
TIME2 = 504  # 3 week in hour
DELTA = 1

NUMEXAMPLE = 5


def main():
    k_minus = 0.11
    kp_list = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25, 0.3]
    kp_n = len(kp_list)
    kplist_samplelist_genets= np.zeros((kp_n, NUMEXAMPLE, TIME2, 3, NUM_OF_HISTONE))

    for i, kp in enumerate(kp_list):
        kplist_samplelist_genets[i] = submain(kp, k_minus)
        print("done", kp)

    compressed = histone.compress_kplist_samplelist_hstseqts(kplist_samplelist_genets)
    os.mkdir("__k-"+str(k_minus)+"/")
    filename3d = "__k-"+str(k_minus)+"/"+"dumpdata3d__k-{}__{}examples.csv".format(k_minus, NUMEXAMPLE)
    filename2d = "__k-"+str(k_minus)+"/"+"dumpdata2d__pos__k-{}__{}examples.csv".format(k_minus, NUMEXAMPLE)

    histone.write_dump3d_kp_time_hst(compressed, filename3d, TIME2)
    day8 = 8 * 24
    histone.write_dump2d_pos_kp(kplist_samplelist_genets, filename2d, day8)  # for plot2


def submain( k_plus, k_minus):
    one_variation = np.zeros((NUMEXAMPLE, TIME2, 3, NUM_OF_HISTONE))
    # os.mkdir("__"+str(path_num+2)+"__k+" + str(k_plus)+"__k-"+str(k_minus))
    for ex in range(NUMEXAMPLE):
        one_variation[ex] = subsubmain(k_plus, k_minus)
        print(ex, k_minus)
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
