import histone
import histone.figure as figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib

NUM_OF_HISTONE = 81
WINDOW = 10
TIME1 = 504  # 3 week
TIME2 = 504  # 3 week
DELTA = 1


def main():
    k_plus = 0.2
    k_minus = 0.117
    k_nuc = 0.6
    trackerlist = []

    # for c in range(50):
    #     trackerlist = submain(trackerlist, k_plus, k_minus)
    #     print(c)

    # fig = plt.figure()
    # histone.figure.kinetic_model(fig, trackerlist)


    # plt.show()

# def submain(trackerlist, k_plus, k_minus):
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
                                             t_bool=T,
                                             K_NUC=k_nuc
    )
    tracker = dictH["vectorize"]
    hstL = dictH["hstL"]
    TList = dictH["TList"]

    dictH2 = histone.track_epigenetic_process(hst_list=hstL,
                                              time=TIME2,
                                              a_bool=secA,
                                              r_bool=secR,
                                              t_bool=TList[-1],
                                              K_NUC=0.2
    )
    tracker2 = dictH2["vectorize"]

    finalTracker = np.concatenate((tracker, tracker2))
    finalTList = np.concatenate((TList, dictH2["TList"]))

    fig = plt.figure()
    histone.figure.sequence(fig, finalTracker, 3, 1, 1)
    histone.figure.transcription(fig, finalTList, 9, 1, 4)
    histone.figure.window(fig, finalTracker, 9, 1, 5)
    histone.figure.m_stat(fig, tracker, 9, 4, 22)
    histone.figure.m_stat(fig, tracker2, 9, 4, 24)
    plt.show()

    title = "test.pdf"
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()
    



if __name__ == "__main__":
    main()
