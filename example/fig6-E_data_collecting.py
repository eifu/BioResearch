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
dir = 'example/data{}_withNUC/'.format(example_n)
if not os.path.exists(dir):
    os.mkdir(dir)


def main():
    kn = 0.1
    ka = 0.1
    if not os.path.exists(dir+"kn{}__ka{}/".format(round(kn,4), round(ka,4))):
        os.mkdir(dir + "kn{}__ka{}/".format(round(kn,4), round(ka,4)))

    km_list = [i for i in np.arange(0.05, 0.21, 0.05)]
    kp_list = [i for i in np.arange(0.05, 0.21, 0.05)]

    for km in km_list:
        if not os.path.exists(dir+"kn{}__ka{}/".format(kn,ka) + "__k-{}/".format(round(km, 4))):
            os.mkdir(dir+"kn{}__ka{}/".format(kn,ka) + "__k-{}/".format(round(km, 4)))
        for i, kp in enumerate(kp_list):
            onekp_samplelist_genets = submain(kp, km, kn, ka)
            print("done", kp)

            compressed = io.compress_onekp_samplelist_hstseqts(onekp_samplelist_genets)

            filename2d = dir + "kn{}__ka{}/".format(kn,ka) + "__k-{}/".format(round(km, 4)) \
                         + "dumpdata2d__k+{}__{}examples.csv".format(round(kp, 4), example_n)
            io.write_dump2d_onekp_time_hst(compressed, filename2d, TIME2)


def submain(k_plus, k_minus, k_nuc, k_ace):
    one_variation = np.zeros((example_n, TIME2, 3, HST_N))
    for ex in range(example_n):
        one_variation[ex] = subsubmain(k_plus, k_minus, k_nuc, k_ace)
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
                                             K_ACE=k_ace,
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
                                              K_ACE=k_ace,
                                              K_NUC=k_nuc
                                              )
    tracker2 = dictH2["vectorize"]

    return tracker2


if __name__ == "__main__":
    main()
