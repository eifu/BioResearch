import histone
import histone.figure
from time import strftime
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import numpy as np

"""
demonstlation of histone library.
In this file, we demonstlate the histone code without DNA methylation model.
As are imported above, you have to have two libraries, numpy, matplotlib on
your local computer.

This code outputs a figure that has a list of histone in y-axis, and
time series in x-axis.

NUM_OF_HISTONE ... global variable, a length of histones that you keep track.
WINDOW ... global variable, a length of locus in the list of histones.
TIME1 ... global variable, a time length of how long histones are held before
the environment changes.
TIME2 ... global variable, a time length of how long histones are held after
the environment changes.
"""

NUM_OF_HISTONE = 81
WINDOW = 10
TIME1 = 3 * 7 * 24  # 3 week in hours
TIME2 = 3 * 7 * 24  # 3 week in hours
K_PLUS = 0.145
K_MINUS = 0.145

PERCENT = 100

NUC_PROB = 1
ACE_PROB = 0.1
NUC_PROB2 = 0
ACE_PROB2 = 0.1

P_OFF = 0.01


def main():


    """
    init_genome

    initializes a list of histone objects.
    if you want to have a list of histones with
    specific characteristics, you can customize it with
    keyword arguments.
    input:
    (list of keyword arguments)
    - percentage .. probability of methylated histones.
    - hst_n .. number of histones in a list
    - kp .. probability of k plus
    - km .. probability of k minus
    output:
    a list of histone objects
    """
    histone_list1 = histone.init_genome_with_dna_model(percentage=PERCENT,
                                                       hst_n=NUM_OF_HISTONE,
                                                       kp=K_PLUS,
                                                       km=K_MINUS,
                                                       )

    histone_list1[38].set_cpg_island_list(4)  # turn -2 pos on
    histone_list1[39].set_cpg_island_list(4)  # turn -1 pos on
    histone_list1[40].set_cpg_island_list(2)  # turn 0 pos on

    """
    track_epigenetic_process

    returns a dictionary of the result of the records under
    a specific environment. you can setup the environment
    with keyword arguments.
    input:
    (list of keyword arguments)
    - hst_list .. the initial condition of list of histones
    - time .. the time length of tracking
    - ace_prob .. K-ace probability.
    - nuc_prob .. K-nuc probability.
    - p_off .. probability of CpG island getting unmethylated.
    output:
    dictionary object that contains three data,
    key: vectorize, value: a record of timeseries of list of histones.
    key: hstL, value: a final status of list of histones
    key: TList, value: a record of timeseries of whether transcription happens.
    """
    dict1 = histone.track_epigenetic_process_with_dna_model(hst_list=histone_list1,
                                                            time=TIME1,
                                                            ace_prob=ACE_PROB,
                                                            nuc_prob=NUC_PROB,
                                                            p_off=P_OFF
                                                            )
    tracker = dict1["vectorize"]
    hst1 = dict1["hstL"]
    t_list1 = dict1["TList"]
    p_list1 = dict1["PList"]

    dict2 = histone.track_epigenetic_process_with_dna_model(hst_list=hst1,
                                                            time=TIME2,
                                                            ace_prob=ACE_PROB2,
                                                            nuc_prob=NUC_PROB2,
                                                            p_off=P_OFF
                                                            )
    tracker2 = dict2["vectorize"]
    t_list2 = dict2["TList"]
    p_list2 = dict2["PList"]

    final_tracker = np.concatenate((tracker, tracker2))
    final_t_list = np.concatenate((t_list1, t_list2))
    final_p_list = np.concatenate((p_list1, p_list2))

    # basic setup for matplotlib
    plt.style.use('ggplot')
    font = {'family': 'sans-serif'}
    matplotlib.rc('font', **font)

    fig = plt.figure()

    title = "test_{}\nk+{}_k-{}\nkace1:{} knuc1:{} -> kace2:{} knuc2:{}".format(
        strftime("%Y_%m%d_%H%M"), K_PLUS, K_MINUS, ACE_PROB, NUC_PROB, ACE_PROB2, NUC_PROB2)
    fig.suptitle(title, size=8)

    histone.figure.sequence(fig, final_tracker, 3, 1, 1)
    histone.figure.window_with_dna_model(fig, final_tracker, 6, 1, 3)
    histone.figure.package_with_dna_model(fig, final_p_list, final_tracker, 6, 1, 4)
    histone.figure.transcription(fig, final_t_list, 6, 1, 5)
    histone.figure.m_stat(fig, tracker, 6, 4, 22, start_time_ratio=0.5, end_time_ratio=1)
    histone.figure.m_stat(fig, tracker2, 6, 4, 23, start_time_ratio=0, end_time_ratio=0.5)
    histone.figure.m_stat(fig, tracker2, 6, 4, 24, start_time_ratio=0.5, end_time_ratio=1)
    plt.show()

    # to save a figure to pdf file
    fig.tight_layout()
    fig.subplots_adjust(top=0.88)

    path = "demo/result/data_{}__k+{}k-{}__knuc1_{}__kace1_{}____knuc2_{}__kace2_{}.pdf".format(
        strftime("%Y_%m%d_%H%M"), K_PLUS, K_MINUS, NUC_PROB, ACE_PROB, NUC_PROB2, ACE_PROB2)

    pp = PdfPages(path)
    pp.savefig(fig)
    pp.close()


if __name__ == "__main__":
    main()
