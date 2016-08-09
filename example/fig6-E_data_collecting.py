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

example_n = 1000

if not os.path.exists("example/data"):
    os.mkdir('example/data')

if not os.path.exists('example/data/withNUC'):
    os.mkdir('example/data/withNUC')

dir = 'example/data/withNUC/data{}/'.format(example_n)
if not os.path.exists(dir):
    os.mkdir(dir)


def main():
    kn2_list = np.arange(0, 0.5, 0.025)
    for kn2 in kn2_list:
        submain1(kn2)


def submain1(kn2):
    ka1 = 0
    ka2 = 0
    kn1 = 1


    kmkppair = [(0.145, 0.145)]

    dir2 = dir + "kn{}ka{}_kn{}ka{}/".format(round(kn1, 4),
                                             round(ka1, 4),
                                             round(kn2, 4),
                                             round(ka2, 4))
    if not os.path.exists(dir2):
        os.mkdir(dir2)

    for km, kp in kmkppair:
        onekp_samplelist_genets, packaginglist = submain(kp, km, kn1, ka1, kn2, ka2)
        print("kn:{}, ka:{} ,done km:{}, kp:{}".format(round(kn2, 4),
                                                       round(ka2, 4),
                                                       round(km, 4),
                                                       round(kp, 4)))
        print(packaginglist,packaginglist.shape)
        compressed = io.compress_onekp_samplelist_hstseqts(onekp_samplelist_genets)

        filename2d = dir2 + "dumpdata2d__k+{}k-{}_{}examples.csv".format(round(km, 4),
                                                                         round(kp, 4),
                                                                         example_n)
        io.write_dump2d_onekp_time_hst(compressed, filename2d, TIME2)

        filename2d = dir2 + "packaging__k+{}k-{}_{}examples.csv".format(round(km, 4),
                                                                        round(kp, 4),
                                                                        example_n)
        compre = io.compress_packaging_samplelist(packaginglist)
        io.write_dump2d_onekp_time_hst(compre, filename2d, TIME2)


def submain(k_plus, k_minus, kn1, ka1, kn2, ka2):
    one_var_m = np.zeros((example_n, TIME2, 3, HST_N))
    one_var_pack = np.zeros((example_n, TIME2))
    for ex in range(example_n):
        one_var_m[ex], one_var_pack[ex] = subsubmain(k_plus, k_minus, k_nuc1=kn1, k_ace1=ka1, k_nuc2=kn2, k_ace2=ka2)
        print("kn:{}, ka:{}, km:{}, kp:{}, example number:{}".format(round(kn2, 4),
                                                                     round(ka2, 4),
                                                                     round(k_minus, 4),
                                                                     round(k_plus, 4),
                                                                     ex))
    return one_var_m, one_var_pack


def subsubmain(k_plus, k_minus, k_nuc1, k_ace1, k_nuc2, k_ace2):
    R = 0
    A = 1
    secR = 1
    secA = 1

    T = 0
    P =True

    histoneList1 = histone.init_genome(percentage=50,
                                       a_bool=A,
                                       hst_n=HST_N,
                                       kp=k_plus,
                                       kp2=k_plus,
                                       ka=k_ace1,
                                       km=k_minus,
                                       )

    dictH = histone.track_epigenetic_process(hst_list=histoneList1,
                                             time=TIME1,
                                             a_bool=A,
                                             t_bool=T,
                                             p_bool=P,
                                             ace_prob=k_ace1,
                                             nuc_prob=k_nuc1
                                             )
    # tracker = dictH["vectorize"]
    hstL = dictH["hstL"]
    TList = dictH["TList"]
    PList = dictH["PList"]

    dictH2 = histone.track_epigenetic_process(hst_list=hstL,
                                              time=TIME2,
                                              a_bool=secA,
                                              t_bool=TList[-1],
                                              p_bool=PList[-1],
                                              ace_prob=k_ace2,
                                              nuc_prob=k_nuc2
                                              )
    tracker2 = dictH2["vectorize"]
    p_list2 = dictH2["PList"]

    return tracker2, p_list2


if __name__ == "__main__":
    main()
