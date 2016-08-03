import histone
import histone.figure
from time import strftime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import numpy as np

NUM_OF_HISTONE = 81
WINDOW = 10
TIME1 = 3 * 7 * 24  # 3 week in hours
TIME2 = 3 * 7 * 24  # 3 week in hours

k_plus = 0.145
k_minus = 0.145
k_ace = 0

R = 0
A = 1
secR = 1
secA = 1


def main():
    for k_nuc in np.arange(0.05, 0.96, 0.05):
        submain(k_nuc)
        print(k_nuc)

def submain(k_nuc):
    # k_nuc = 0.95
    percent = 100

    k_nuc_default = 0
    histoneList1 = histone.init_genome(percentage=percent,
                                       a_bool=A,
                                       hst_n=NUM_OF_HISTONE,
                                       kp=k_plus,
                                       kp2=k_plus,
                                       km=k_minus,
                                       ka=k_ace
                                       )

    T = 0
    dictH = histone.track_epigenetic_process(hst_list=histoneList1,
                                             time=TIME1,
                                             a_bool=A,
                                             r_bool=R,
                                             t_bool=T,
                                             ace_prob=k_ace,
                                             nuc_prob=k_nuc_default
                                             )
    tracker = dictH["vectorize"]
    hstL = dictH["hstL"]
    TList = dictH["TList"]
    dictH2 = histone.track_epigenetic_process(hst_list=hstL,
                                              time=TIME2,
                                              a_bool=secA,
                                              r_bool=secR,
                                              t_bool=TList[-1],
                                              ace_prob=k_ace,
                                              nuc_prob=k_nuc
                                              )
    tracker2 = dictH2["vectorize"]

    finalTracker = np.concatenate((tracker, tracker2))
    finalTList = np.concatenate((TList, dictH2["TList"]))

    # basic setup for matplotlib
    plt.style.use('ggplot')
    font = {'family': 'sans-serif'}
    matplotlib.rc('font', **font)

    fig = plt.figure()

    histone.figure.sequence(fig, finalTracker, 3, 1, 1, kace=k_ace, knuc=k_nuc_default, kace2=k_ace, knuc2=k_nuc)
    histone.figure.window(fig, finalTracker, 6, 1, 3)

    k_nuc_default = True
    histoneList1 = histone.init_genome(percentage=percent,
                                       a_bool=A,
                                       hst_n=NUM_OF_HISTONE,
                                       kp=k_plus,
                                       kp2=k_plus,
                                       km=k_minus,
                                       ka=k_ace
                                       )

    T = 0
    dictH = histone.track_epigenetic_process(hst_list=histoneList1,
                                             time=TIME1,
                                             a_bool=A,
                                             r_bool=R,
                                             t_bool=T,
                                             ace_prob=k_ace,
                                             nuc_prob=k_nuc_default
                                             )
    tracker = dictH["vectorize"]
    hstL = dictH["hstL"]
    TList = dictH["TList"]
    dictH2 = histone.track_epigenetic_process(hst_list=hstL,
                                              time=TIME2,
                                              a_bool=secA,
                                              r_bool=secR,
                                              t_bool=TList[-1],
                                              ace_prob=k_ace,
                                              nuc_prob=k_nuc
                                              )
    tracker2 = dictH2["vectorize"]

    finalTracker = np.concatenate((tracker, tracker2))
    finalTList = np.concatenate((TList, dictH2["TList"]))

    histone.figure.sequence(fig, finalTracker, 3, 1, 3, kace=k_ace, knuc=k_nuc_default, kace2=k_ace, knuc2=k_nuc)
    histone.figure.window(fig, finalTracker, 6, 1, 4)

    # plt.show()

    # to save a figure to pdf file
    title = "example/figure/knuc/k+{}_k-{}/kace{}_knuc{}_percent{}.pdf".format(round(k_plus,4), round(k_minus,4), round(k_ace,4), round(k_nuc, 4),
                                                                               percent)
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()


if __name__ == "__main__":
    main()
