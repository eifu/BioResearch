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

# p_off = 0.01

if not os.path.exists("example/data"):
    os.mkdir('example/data')

if not os.path.exists('example/data/dna_model'):
    os.mkdir('example/data/dna_model')

dir1 = 'example/data/dna_model/data{}/'.format(example_n)
if not os.path.exists(dir1):
    os.mkdir(dir1)


def main():
    kn1_list = [0, 1]
    kn2_list = np.arange(0, 0.5, 0.025)
    p_off_list = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
    for p_off in p_off_list:
        for kn1 in kn1_list:
            for kn2 in kn2_list:
                submain1(kn1, kn2, p_off)


def submain1(kn1, kn2, p_off):
    ka1 = 0
    ka2 = 0

    km_kp_pair = [(0.145, 0.145)]

    dir2 = dir1 + "kn{}ka{}_kn{}ka{}/".format(round(kn1, 4),
                                              round(ka1, 4),
                                              round(kn2, 4),
                                              round(ka2, 4))
    if not os.path.exists(dir2):
        os.mkdir(dir2)

    for km, kp in km_kp_pair:
        one_var_tracker, one_var_week3_hst, one_var_pack, one_var_cpg_sum = submain(kp, km, kn1, ka1, kn2, ka2, p_off)
        print("kn:{}, ka:{} ,done km:{}, kp:{}  poff".format(round(kn2, 4),
                                                             round(ka2, 4),
                                                             round(km, 4),
                                                             round(kp, 4),
                                                             round(p_off, 4)))

        # for tracker info
        filename2d = dir2 + "dumpdata2d_p_off{}_k+{}k-{}_{}examples.csv".format(round(p_off, 4),
                                                                                round(kp, 4),
                                                                                round(km, 4),
                                                                                example_n)
        compressed = io.compress_onekp_samplelist_hstseqts(one_var_tracker)
        io.write_dump2d_onekp_time_hst(compressed, filename2d, TIME2)

        # for final histone list info
        filename2d = dir2 + "final_hst_list_p_off{}_k+{}k-{}_{}examples.csv".format(round(p_off, 4),
                                                                                    round(kp, 4),
                                                                                    round(km, 4),
                                                                                    example_n)
        compressed = io.compress_last_week_hst_vec(one_var_week3_hst)
        io.write_dump2d_final_hst_list(compressed, filename2d, HST_N)

        # for packaging info
        filename2d = dir2 + "packaging_p_off{}_k+{}k-{}_{}examples.csv".format(round(p_off, 4),
                                                                               round(kp, 4),
                                                                               round(km, 4),
                                                                               example_n)
        compressed = io.compress_packaging_samplelist(one_var_pack)
        io.write_dump2d_onekp_time_hst(compressed, filename2d, TIME2)

        # for cpg sum info
        filename2d = dir2 + "cpg_sum_p_off{}_k+{}k-{}_{}examples.csv".format(round(p_off, 4),
                                                                             round(kp, 4),
                                                                             round(km, 4),
                                                                             example_n)
        print(one_var_cpg_sum)
        compressed = io.compress_cpg_samplelist(one_var_cpg_sum)
        io.write_dump2d_cpg_sum(compressed, filename2d)


def submain(k_plus, k_minus, kn1, ka1, kn2, ka2, p_off):
    one_var_m = np.zeros((example_n, TIME2, 4, HST_N))
    one_var_hst_list = np.zeros((example_n, 24 * 7, 11))  # week 3 histone list
    one_var_pack = np.zeros((example_n, TIME2))
    one_var_cpg_sum = np.zeros((example_n, TIME1 + TIME2))
    for ex in range(example_n):
        one_var_m[ex], one_var_hst_list[ex], one_var_pack[ex], one_var_cpg_sum[ex] = subsubmain(k_plus=k_plus,
                                                                                                k_minus=k_minus,
                                                                                                k_nuc1=kn1,
                                                                                                k_ace1=ka1,
                                                                                                k_nuc2=kn2,
                                                                                                k_ace2=ka2,
                                                                                                p_off=p_off)

        print("kn1:{}, ka2:{}, -> kn2:{}, ka2:{}, km:{}, kp:{}, complete {}%".format(round(kn1, 4),
                                                                                     round(ka1, 4),
                                                                                     round(kn2, 4),
                                                                                     round(ka2, 4),
                                                                                     round(k_minus, 4),
                                                                                     round(k_plus, 4),
                                                                                     round(ex * 100 / example_n, 4)
                                                                                     )
              )
    return one_var_m, one_var_hst_list, one_var_pack, one_var_cpg_sum


def subsubmain(k_plus, k_minus, k_nuc1, k_ace1, k_nuc2, k_ace2, p_off):
    A = 0
    secA = 0

    T = 0
    P = True

    initial_hst_list = histone.init_genome_with_dna_model(percentage=50,
                                                          a_bool=A,
                                                          hst_n=HST_N,
                                                          kp=k_plus,
                                                          kp2=k_plus,
                                                          ka=k_ace1,
                                                          km=k_minus,
                                                          )

    if k_nuc1 == 1:
        initial_hst_list[38].set_cpg_island_list(4)
        initial_hst_list[39].set_cpg_island_list(4)
        initial_hst_list[40].set_cpg_island_list(2)

    dict1 = histone.track_epigenetic_process_with_dna_model(hst_list=initial_hst_list,
                                                            time=TIME1,
                                                            a_bool=A,
                                                            t_bool=T,
                                                            p_bool=P,
                                                            ace_prob=k_ace1,
                                                            nuc_prob=k_nuc1,
                                                            p_off=p_off
                                                            )
    tracker = dict1["vectorize"]
    hst_list = dict1["hstL"]
    t_list = dict1["TList"]
    p_list = dict1["PList"]
    cpg_sum_list = np.zeros(TIME1 + TIME2)
    for i, hst in enumerate(tracker):
        cpg_sum_list[i] = sum(hst[3][35:46])

    dict2 = histone.track_epigenetic_process_with_dna_model(hst_list=hst_list,
                                                            time=TIME2,
                                                            a_bool=secA,
                                                            t_bool=t_list[-1],
                                                            p_bool=p_list[-1],
                                                            ace_prob=k_ace2,
                                                            nuc_prob=k_nuc2,
                                                            p_off=p_off
                                                            )
    tracker2 = dict2["vectorize"]
    p_list2 = dict2["PList"]

    week3_hst_list = np.zeros((24 * 7, 11))
    for i, hst in enumerate(tracker2[2 * 24 * 7:]):
        week3_hst_list[i] = hst[0][35:46] - hst[2][35:46]

    for i, hst in enumerate(tracker2):
        cpg_sum_list[i + TIME1] = sum(hst[3][35:46])

    return tracker2, week3_hst_list, p_list2, cpg_sum_list


if __name__ == "__main__":
    main()
