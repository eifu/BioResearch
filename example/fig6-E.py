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

NUM_OF_HISTONE = 81
WINDOW = 10
TIME1 = 504  # 3 week in hour
TIME2 = 504  # 3 week in hour
DELTA = 1

NUMEXAMPLE = 10


def main():
    fig = plt.figure()
    variation_list_vecgenetimeseries =[]
    k_minus = 0.1

    variation_list_vecgenetimeseries.append(submain(0.0001, k_minus))
    variation_list_vecgenetimeseries.append(submain(0.001, k_minus))

    for k_plus in np.arange(0.01, 0.21, 0.01):
        variation_list_vecgenetimeseries.append(submain(k_plus, k_minus))
        print("done  ", k_plus )
    variation_list_vecgenetimeseries.append(submain(0.25, k_minus))
    variation_list_vecgenetimeseries.append(submain(0.3,k_minus))
    print(len(variation_list_vecgenetimeseries))
    # list_after8days = list([after8days])
    # print(list_after8days,len(list_after8days))

    figure.figure6c_and_6e(fig, variation_list_vecgenetimeseries)

    plt.show()

    title = "fig6-CE/fig_test1__{}examples__k-{}.pdf".format(NUMEXAMPLE,k_minus)
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()


def submain(k_plus,k_minus):
    count = 0
    one_variation = np.zeros((NUMEXAMPLE,TIME1+TIME2,3,NUM_OF_HISTONE))
    for i in range(NUMEXAMPLE):
        one_variation[i] = subsubmain(count, k_plus, k_minus)
        count += 1

    return one_variation


def subsubmain(count, k_plus, k_minus):
    R = 0
    A = 1
    secR = 1
    secA = 1

    T = 0
    plt.style.use('ggplot')
    font = {'family': 'sans-serif'}
    matplotlib.rc('font', **font)

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

    finalTracker = np.concatenate((tracker, tracker2))

    # list_tracker.append(finalTracker)

    print(count)

    return finalTracker



if __name__ == "__main__":
    main()
