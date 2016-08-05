import histone
import histone.io as io
import numpy as np
import os

TIME1 = 504  # 3 week in hour
TIME2 = 504  # 3 week in hour

HST_N = 81

example_n = 100

if not os.path.exists("example/data"):
    os.mkdir("example/data")

if not os.path.exists("example/data/withNUC"):
    os.mkdir("example/data/withNUC")

dir = 'example/data/withNUC/data{}_packagedata/'.format(example_n)
if not os.path.exists(dir):
    os.mkdir(dir)


def main():
    kn2_list = np.arange(0, 0.4, 0.025)
    ka2_list = np.arange(0, 0.4, 0.025)

    kn1 = 1
    ka1 = 0

    km = 0.145
    kp = 0.145

    for kn2 in kn2_list:
        for ka2 in ka2_list:
            samplelist = submain(kp, km, kn1, ka1, kn2, ka2)
            print("kn:{}, ka:{} ,done km:{}, kp:{}".format(round(kn2, 5),
                                                           round(ka2, 5),
                                                           round(km, 5),
                                                           round(kp, 5)))

            filename = dir + "package_data__kn{}ka{}_kn{}ka{}__k+{}k-{}_{}examples.csv".format(round(kn1, 5),
                                                                                               round(ka1, 5),
                                                                                               round(kn2, 5),
                                                                                               round(ka2, 5),
                                                                                               round(km, 5),
                                                                                               round(kp, 5), example_n)
            with open(filename, 'wb') as f:
                np.savetxt(f,
                           samplelist,
                           fmt='%d',
                           delimiter=',',
                           newline='\n')


def submain(k_plus, k_minus, kn1, ka1, kn2, ka2):
    one_variation = np.zeros((example_n, 4))
    for ex in range(example_n):
        one_variation[ex] = subsubmain(k_plus, k_minus, k_nuc1=kn1, k_ace1=ka1, k_nuc2=kn2, k_ace2=ka2)
        print(
            "kn1: {}, ka1: {},  kn2: {}, ka2: {}  example number:{}".format(round(kn1, 4),
                                                                            round(ka1, 4),
                                                                            round(kn2, 4),
                                                                            round(ka2, 4),
                                                                            ex))
    return one_variation


def subsubmain(k_plus, k_minus, k_nuc1, k_ace1, k_nuc2, k_ace2):
    R = 0
    A = 1
    secR = 1
    secA = 1

    T = 0
    P = True

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
                                             r_bool=R,
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
                                              r_bool=secR,
                                              t_bool=TList[-1],
                                              p_bool=PList[-1],
                                              ace_prob=k_ace2,
                                              nuc_prob=k_nuc2
                                              )

    PList2 = dictH2["PList"]
    finalPList = np.concatenate((PList, PList2))

    pp = 0
    pu = 0
    up = 0
    uu = 0
    for i, b in enumerate(finalPList[TIME1:]):
        if finalPList[TIME1 + i - 1] == True and b == True:
            pp += 1
        elif finalPList[TIME1 + i - 1] == True and b == False:
            pu += 1
        elif finalPList[TIME1 + i - 1] == False and b == True:
            up += 1
        else:
            uu += 1

    return [pp, pu, up, uu]


if __name__ == "__main__":
    main()
