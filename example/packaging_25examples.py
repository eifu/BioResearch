import histone
import histone.figure
from time import strftime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
from matplotlib import cm
import numpy as np

NUM_OF_HISTONE = 81
WINDOW = 10
TIME1 = 3 * 7 * 24  # 3 week in hours
TIME2 = 3 * 7 * 24  # 3 week in hours

k_plus = 0.145
k_minus = 0.145

percent = 100

a_bool = 1
a_bool2 = 1

k_ace = 0
k_ace2 = 0

N = 25


def main():
    p_off_list = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
    # p_off_list = [0.5]
    for p_off in p_off_list:
        submain(p_off)


def submain(p_off):
    # basic setup for matplotlib
    plt.style.use('ggplot')
    font = {'family': 'sans-serif'}
    matplotlib.rc('font', **font)

    kn2_list = np.arange(0.025, 0.125, 0.025)
    variation = len(kn2_list)
    cm_subsection = np.linspace(0.0, 1.0, 2 * variation)
    colors = [cm.jet(x) for x in cm_subsection]
    hours = np.arange(TIME1 + TIME2)
    fig1 = plt.figure()
    container = np.zeros((variation, N, TIME1 + TIME2))

    for i_kn2, kn2 in enumerate(kn2_list):

        ax = fig1.add_subplot(variation, 1, i_kn2 + 1)
        if i_kn2 == 0:
            ax.set_title("kn1: 1, p off: {}".format(p_off))
        ax.set_xlim(0, TIME1 + TIME2)
        ax.set_ylim(-0.3, 1.3)
        ax.set_xticks(np.arange(0, TIME1 + TIME2, 24 * 7))
        if i_kn2 == 3:
            ax.set_xticklabels([str(i // (24 * 7) + 1) + " week" for i in range(0, TIME1 + TIME2, 24 * 7)])
        else:
            ax.set_xticklabels([])
        ax.set_yticks([])
        ax.set_ylabel("k2: {}".format(round(kn2, 4)))

        for i in range(N):
            container[i_kn2][i] = subsubmain(p_off, kn2)
            print(p_off, i_kn2, i)
            ax.fill_between(hours, 0,  container[i_kn2][i], alpha=1/N, color=colors[i_kn2])

    # to save a figure to pdf file
    title = "packaging_{}examples_p_off{}.pdf".format(N, p_off)
    pp = PdfPages(title)
    pp.savefig(fig1)
    pp.close()

    for i_kn2, kn2 in enumerate(kn2_list):
        fig2 = plt.figure()

        for i in range(N):
            ax = fig2.add_subplot(N, 1, i + 1)
            if i == 0:
                ax.set_title("kn1: 1, kn2: {}, p off: {}".format(round(kn2, 4),
                                                                 round(p_off, 4)))
            ax.set_xlim(0, TIME1 + TIME2)
            ax.set_ylim(-0.3, 1.3)
            ax.set_xticks(np.arange(0, TIME1 + TIME2, 24 * 7))
            if i == N - 1:
                ax.set_xticklabels([str(i // (24 * 7) + 1) + " week" for i in range(0, TIME1 + TIME2, 24 * 7)])
            else:
                ax.set_xticklabels([])
            ax.set_yticks([])
            ax.fill_between(hours, 0, container[i_kn2][i], color=colors[i_kn2])
            # ax.fill_between(hours, )
        # to save a figure to pdf file
        title = "packaging_25examples_p_off{}_kn2_{}.pdf".format(round(p_off, 4),
                                                   round(kn2, 4))
        pp = PdfPages(title)
        pp.savefig(fig2)
        pp.close()


def subsubmain(p_off, kn2):
    k_nuc = 1

    histone_list1 = histone.init_genome_with_dna_model(percentage=percent,
                                                       a_bool=a_bool,
                                                       hst_n=NUM_OF_HISTONE,
                                                       kp=k_plus,
                                                       kp2=k_plus,
                                                       km=k_minus,
                                                       ka=k_ace
                                                       )

    histone_list1[38].set_cpg_island_list(4)  # turn -2 pos on
    histone_list1[39].set_cpg_island_list(4)  # turn -1 pos on
    histone_list1[40].set_cpg_island_list(2)  # turn 0 pos on

    t_bool = 0
    p_bool = True

    dict1 = histone.track_epigenetic_process_with_dna_model(hst_list=histone_list1,
                                                            time=TIME1,
                                                            a_bool=a_bool,
                                                            t_bool=t_bool,
                                                            p_bool=p_bool,
                                                            ace_prob=k_ace,
                                                            nuc_prob=k_nuc,
                                                            p_off=p_off
                                                            )
    tracker = dict1["vectorize"]
    hst1 = dict1["hstL"]
    t_list1 = dict1["TList"]
    p_list1 = dict1["PList"]

    dict2 = histone.track_epigenetic_process_with_dna_model(hst_list=hst1,
                                                            time=TIME2,
                                                            a_bool=a_bool2,
                                                            t_bool=t_list1[-1],
                                                            p_bool=p_list1[-1],
                                                            ace_prob=k_ace2,
                                                            nuc_prob=kn2,
                                                            p_off=p_off
                                                            )
    tracker2 = dict2["vectorize"]
    t_list2 = dict2["TList"]
    p_list2 = dict2["PList"]

    # final_tracker = np.concatenate((tracker, tracker2))
    # final_t_list = np.concatenate((t_list1, t_list2))
    final_p_list = np.concatenate((p_list1, p_list2))

    return final_p_list


if __name__ == "__main__":
    main()
