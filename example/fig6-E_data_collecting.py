"""
figure 6 C, E

"""
import histone
import histone.io as io
import numpy as np
import os

TIME1 = 504  # 3 week in hour
TIME2 = 504  # 3 week in hour

HST_N = 81

example_n = 5

if not os.path.exists("example/data"):
    os.mkdir('example/data')

if not os.path.exists('example/data/withNUC'):
    os.mkdir('example/data/withNUC')

dir1 = 'example/data/withNUC/data{}/'.format(example_n)
if not os.path.exists(dir1):
    os.mkdir(dir1)


def main():
    kn2_list = np.arange(0, 0.5, 0.025)
    for kn2 in kn2_list:
        submain1(kn2)


def submain1(kn2):
    # kn2 = 1
    ka1 = 0.1
    ka2 = 0.1
    kn1 = 1

    kmkppair = [(0.145, 0.145)]

    dir2 = dir1 + "kn{}ka{}_kn{}ka{}/".format(round(kn1, 4),
                                              round(ka1, 4),
                                              round(kn2, 4),
                                              round(ka2, 4))
    if not os.path.exists(dir2):
        os.mkdir(dir2)

    for km, kp in kmkppair:
        one_var_tracker, one_var_last_week_hst, one_var_pack = submain(kp, km, kn1, ka1, kn2, ka2)
        print("kn:{}, ka:{} ,done km:{}, kp:{}".format(round(kn2, 4),
                                                       round(ka2, 4),
                                                       round(km, 4),
                                                       round(kp, 4)))

        # for tracker info
        filename2d = dir2 + "dumpdata2d__k+{}k-{}_{}examples.csv".format(round(km, 4),
                                                                         round(kp, 4),
                                                                         example_n)
        compressed = io.compress_onekp_samplelist_hstseqts(one_var_tracker)
        io.write_dump2d_onekp_time_hst(compressed, filename2d, TIME2)

        # for final histone list info
        filename2d = dir2 + "final_hst_list_k+{}k-{}_examples.csv".format(round(km, 4),
                                                                          round(kp, 4),
                                                                          example_n)
        compressed = io.compress_last_week_hst_vec(one_var_last_week_hst)
        io.write_dump2d_final_hst_list(compressed, filename2d, HST_N)

        # for packaging info
        filename2d = dir2 + "packaging__k+{}k-{}_{}examples.csv".format(round(km, 4),
                                                                        round(kp, 4),
                                                                        example_n)
        compressed = io.compress_packaging_samplelist(one_var_pack)
        io.write_dump2d_onekp_time_hst(compressed, filename2d, TIME2)


def submain(k_plus, k_minus, kn1, ka1, kn2, ka2):
    one_var_m = np.zeros((example_n, TIME2, 3, HST_N))
    one_var_hst_list = np.zeros((example_n, 24*7,  HST_N)) # week 3 histone list
    one_var_pack = np.zeros((example_n, TIME2))
    for ex in range(example_n):
        one_var_m[ex], one_var_hst_list[ex], one_var_pack[ex] = subsubmain(k_plus, k_minus, k_nuc1=kn1, k_ace1=ka1,
                                                                           k_nuc2=kn2, k_ace2=ka2)
        print("kn:{}, ka:{}, km:{}, kp:{}, example number:{}".format(round(kn2, 4),
                                                                     round(ka2, 4),
                                                                     round(k_minus, 4),
                                                                     round(k_plus, 4),
                                                                     ex))
    return one_var_m, one_var_hst_list, one_var_pack


def subsubmain(k_plus, k_minus, k_nuc1, k_ace1, k_nuc2, k_ace2):
    R = 0
    A = 1
    secR = 1
    secA = 1

    T = 0
    P = True

    initial_hst_list = histone.init_genome(percentage=50,
                                           a_bool=A,
                                           hst_n=HST_N,
                                           kp=k_plus,
                                           kp2=k_plus,
                                           ka=k_ace1,
                                           km=k_minus,
                                           )

    dict1 = histone.track_epigenetic_process(hst_list=initial_hst_list,
                                             time=TIME1,
                                             a_bool=A,
                                             t_bool=T,
                                             p_bool=P,
                                             ace_prob=k_ace1,
                                             nuc_prob=k_nuc1
                                             )
    # tracker = dictH["vectorize"]
    hst_list = dict1["hstL"]
    t_list = dict1["TList"]
    p_list = dict1["PList"]

    dict2 = histone.track_epigenetic_process(hst_list=hst_list,
                                             time=TIME2,
                                             a_bool=secA,
                                             t_bool=t_list[-1],
                                             p_bool=p_list[-1],
                                             ace_prob=k_ace2,
                                             nuc_prob=k_nuc2
                                             )
    tracker2 = dict2["vectorize"]
    p_list2 = dict2["PList"]

    week3_hst_list = np.zeros((24*7, HST_N))
    for i, hst in enumerate(tracker2[2*24*7:]):
        week3_hst_list[i] = hst[0] - hst[2]

    return tracker2, week3_hst_list, p_list2


if __name__ == "__main__":
    main()
