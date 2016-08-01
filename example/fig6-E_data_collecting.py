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
dir = 'example/data50_withNUC/'


def main():
#     kn_list = [i for i in np.arange(0.1,0.6,0.1)]
#     ka_list = [i for i in np.arange(0.1,0.6,0.1)]
#     for kn in kn_list:
#         for ka in ka_list:
#             os.mkdir(dir + "kn{}__ka{}/".format(kn, ka))
#             smain(kn,ka)
#
# def smain(KN,KA):
    kn = 0.1
    ka = 0.1

    km_list = [i for i in np.arange(0.05, 0.21, 0.05)]
    kp_list = [i for i in np.arange(0.05, 0.21, 0.05)]

    for km in km_list:
        os.mkdir(dir+"kn{}__ka{}/".format(KN,KA) + "__k-{}/".format(round(km, 4)))
        for i, kp in enumerate(kp_list):
            onekp_samplelist_genets = submain(kp, km, kn)
            print("done", kp)

            compressed = io.compress_onekp_samplelist_hstseqts(onekp_samplelist_genets)

            filename2d = dir + "kn{}__ka{}/".format(KN,KA) + "__k-{}/".format(round(km, 4)) \
                         + "dumpdata2d__k+{}__{}examples.csv".format(round(kp, 4), example_n)
            io.write_dump2d_onekp_time_hst(compressed, filename2d, TIME2)


def submain(k_plus, k_minus, k_nuc):
    one_variation = np.zeros((example_n, TIME2, 3, HST_N))
    for ex in range(example_n):
        one_variation[ex] = subsubmain(k_plus, k_minus, k_nuc)
        print(ex, k_minus)
    return one_variation


def subsubmain(k_plus, k_minus, k_nuc, k_ace):
    R = 0
    A = 1
    secR = 1
    secA = 1

    T = 0

    histoneList1 = histone.init_genome(percentage=50,
                                       a_bool=A,
                                       hst_n=HST_N,
                                       kp=k_plus,
                                       kp2=k_plus,
                                       ka=k_ace,
                                       km=k_minus,
                                       )

    dictH = histone.track_epigenetic_process(hst_list=histoneList1,
                                             time=TIME1,
                                             a_bool=A,
                                             r_bool=R,
                                             t_bool=T,
                                             K_NUC=k_nuc
                                             )
    # tracker = dictH["vectorize"]
    hstL = dictH["hstL"]
    TList = dictH["TList"]

    dictH2 = histone.track_epigenetic_process(hst_list=hstL,
                                              time=TIME2,
                                              a_bool=secA,
                                              r_bool=secR,
                                              t_bool=TList[-1],
                                              K_NUC=k_nuc
                                              )
    tracker2 = dictH2["vectorize"]

    # finalTracker = np.concatenate((tracker, tracker2))
    # histone.save_hst_timeseries(tracker2,
    #                             "__" + str(path_num + 2) + "__k+" + str(k_plus) + "__k-" + str(k_minus) + "/example" + str(ex) + ".csv")

    return tracker2


if __name__ == "__main__":
    main()
