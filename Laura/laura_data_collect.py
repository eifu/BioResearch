import histone
import histone.io as io
import numpy as np
import os

TIME1 = 504  # 3 week in hour
TIME2 = 504  # 3 week in hour

HST_N = 81
EXAMPLE_N = 5000

K_PLUS = 0.145
K_MINUS = 0.145

if not os.path.exists("Laura/data"):
    os.mkdir('Laura/data')

if not os.path.exists('Laura/data/withNUC'):
    os.mkdir('Laura/data/withNUC')

dir1 = 'Laura/data/withNUC/data{}/'.format(EXAMPLE_N)
if not os.path.exists(dir1):
    os.mkdir(dir1)

"""
dumpdata2d__k+0.145k-0.145_5000examples.csv
this is 5000 examples of the main histone status data.
the main histone status data is characterized by vectorized function. Read the function in __init__.py in histone.

final_hst_list_k+0.145k-0.145_5000examples.csv
this is 5000 examples of the all locus info. 1 is M. 0 is U. -1 is A.

packaging__k+0.145k-0.145_5000examples.csv
the packaging data from the beginning to the end to experiment. by default, it is 6 weeks, with 3 week - 3 week.

transcription__k+0.145k-0.145_5000examples.csv
the transcription data from the beginning to the end to experiment. by default, it is 6 weeks, with 3 week - 3 week.


Play this with changing the number of EXAMPLE_N so that you can see how the data is stored.
"""

def main():
    ka1 = 0
    ka2 = 0
    kn1 = 0
    kn2 = 1

    dir2 = dir1 + "kn{}ka{}_kn{}ka{}/".format(round(kn1, 4),
                                              round(ka1, 4),
                                              round(kn2, 4),
                                              round(ka2, 4))
    if not os.path.exists(dir2):
        os.mkdir(dir2)

        one_var_tracker, one_var_all_week_hst, one_var_pack, one_var_t = submain(kn1, ka1, kn2, ka2)
        print("kn:{}, ka:{} ,done km:{}, kp:{}".format(round(kn2, 4),
                                                       round(ka2, 4),
                                                       round(K_MINUS, 4),
                                                       round(K_PLUS, 4)))

        # for tracker info
        filename2d = dir2 + "dumpdata2d__k+{}k-{}_{}examples.csv".format(round(K_MINUS, 4),
                                                                         round(K_PLUS, 4),
                                                                         EXAMPLE_N)
        compressed = io.compress_onekp_samplelist_hstseqts(one_var_tracker)
        io.write_dump2d_onekp_time_hst(compressed, filename2d, TIME1+TIME2)

        # for final histone list info
        filename2d = dir2 + "final_hst_list_k+{}k-{}_{}examples.csv".format(round(K_MINUS, 4),
                                                                            round(K_PLUS, 4),
                                                                            EXAMPLE_N)
        compressed = io.compress_all_week_hst_locus_vec(one_var_all_week_hst)
        io.write_dump2d_final_hst_list(compressed, filename2d, HST_N)

        # for packaging info
        filename2d = dir2 + "packaging__k+{}k-{}_{}examples.csv".format(round(K_MINUS, 4),
                                                                        round(K_PLUS, 4),
                                                                        EXAMPLE_N)
        compressed = io.compress_packaging_samplelist(one_var_pack)
        io.write_dump2d_onekp_time_hst(compressed, filename2d, TIME1+TIME2)

        # for transcription info
        filename2d = dir2 + "transcription__k+{}k-{}_{}examples.csv".format(round(K_MINUS, 4),
                                                                            round(K_PLUS, 4),
                                                                            EXAMPLE_N)
        compressed = io.compress_packaging_samplelist(one_var_t)
        io.write_dump2d_onekp_time_hst(compressed, filename2d, TIME1+TIME2)


def submain(kn1, ka1, kn2, ka2):
    one_var_m = np.zeros((EXAMPLE_N, TIME1+TIME2, 3, HST_N))
    one_var_hst_list = np.zeros((EXAMPLE_N, TIME1+TIME2, 11))  # week 3 histone list at locus
    one_var_pack = np.zeros((EXAMPLE_N, TIME1+TIME2))
    one_var_t = np.zeros((EXAMPLE_N, TIME1+TIME2))
    for ex in range(EXAMPLE_N):
        one_var_m[ex], one_var_hst_list[ex], one_var_pack[ex], one_var_t[ex] = subsubmain( k_nuc1=kn1, k_ace1=ka1,
                                                                           k_nuc2=kn2, k_ace2=ka2)

        print("kn1:{}, ka2:{}, -> kn2:{}, ka2:{}, km:{}, kp:{}, complete {}%".format(round(kn1, 4),
                                                                                     round(ka1, 4),
                                                                                     round(kn2, 4),
                                                                                     round(ka2, 4),
                                                                                     round(K_MINUS, 4),
                                                                                     round(K_PLUS, 4),
                                                                                     round(ex * 100 / EXAMPLE_N, 4)
                                                                                     )
              )
    return one_var_m, one_var_hst_list, one_var_pack, one_var_t


def subsubmain( k_nuc1, k_ace1, k_nuc2, k_ace2):

    initial_hst_list = histone.init_genome(percentage=50,
                                           hst_n=HST_N,
                                           kp=K_PLUS,
                                           km=K_MINUS,
                                           )

    dict1 = histone.track_epigenetic_process(hst_list=initial_hst_list,
                                             time=TIME1,
                                             ace_prob=k_ace1,
                                             nuc_prob=k_nuc1
                                             )
    tracker1 = dict1["vectorize"]
    p_list1 = dict1["PList"]
    t_list1 = dict1["TList"]
    hst_list = dict1["hstL"]

    dict2 = histone.track_epigenetic_process(hst_list=hst_list,
                                             time=TIME2,
                                             ace_prob=k_ace2,
                                             nuc_prob=k_nuc2
                                             )
    tracker2 = dict2["vectorize"]
    p_list2 = dict2["PList"]
    t_list2 = dict2["TList"]

    tracker = np.concatenate((tracker1, tracker2))
    p_list = np.concatenate((p_list1, p_list2))
    t_list = np.concatenate((t_list1, t_list2))

    locus_hst_list = np.zeros((TIME1+TIME2, 11))

    for i, hst in enumerate(tracker):
        locus_hst_list[i] = hst[0][35:46] - hst[2][35:46]

    return tracker, locus_hst_list, p_list, t_list


if __name__ == "__main__":
    main()
