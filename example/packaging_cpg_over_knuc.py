import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm


example_n = 100
dir1 = "example/data/dna_model/data{}/".format(example_n)
TIME1 = 504
TIME2 = 504

def main():

    p_off_list = [1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001]
    for p_off in p_off_list:
        submain(p_off)
        print("done {}".format(p_off))

def submain(p_off):

    kn1_list = [0, 1]
    kn2_list = np.arange(0, 0.5, 0.025)

    ka1 = 0
    ka2 = 0

    kp = 0.145
    km = 0.145

    plt.style.use('ggplot')
    font = {'family': 'sans-serif'}
    matplotlib.rc('font', **font)
    
    fig = plt.figure()


    variation = len(kn2_list)
    cm_subsection = np.linspace(0.0, 1.0, 2 * variation)
    colors = [cm.jet(x) for x in cm_subsection]


    for i_kn1, kn1 in enumerate(kn1_list):

        ax = fig.add_subplot(4, 1,  i_kn1 * 2 + 1)
        if i_kn1 == 0:
            ax.set_ylabel("p off: {}".format(round(p_off, 4)), fontsize=8)
        ax.set_xlim(0, TIME1 +TIME2)
        ax.set_xticks(np.arange(0, TIME1+TIME2, 24 * 7))
        ax.set_xticklabels([str(i//(24*7) + 1) + "week" for i in range(0, TIME1+TIME2, 24*7)])
        ax.set_ylim(0, 10*example_n)
        ax.set_yticks(np.arange(0, 1001, 250))
        ax.set_yticklabels([])
        
        plot_list = []
        bx = fig.add_subplot(4, 1,  i_kn1 * 2 + 2)
        bx.set_xlim(0, TIME1 + TIME2)
        bx.set_xticks([i for i in range(0, TIME1+TIME2, 24*7)])
        bx.set_xticklabels([str(i//(24*7) + 1) + "week" for i in range(0, TIME1+TIME2, 24 * 7)])
        bx.set_ylim(-5, 100.5)
        bx.set_yticks(np.arange(0, 101, 25))
        bx.set_yticklabels(str(i) + "%" for i in range(0, 101, 25))
        
            
        for i_kn2, kn2 in enumerate(kn2_list):

            if i_kn2 == 0:
                if kn1 == 1:
                    ax.set_title("kn1: On")
                else:
                    ax.set_title("kn1: Off")
                    
            dir2 = dir1 + "kn{}ka{}_kn{}ka{}/".format(round(kn1, 4),
                                                      round(ka1, 4),
                                                      round(kn2, 4),
                                                      round(ka2, 4))
            data_cpg_filepath = dir2 + "cpg_sum_p_off{}_k+{}k-{}_{}examples.csv".format(round(p_off, 4),
                                                                                        round(kp, 4),
                                                                                        round(km, 4),
                                                                                        example_n)
            data_cpg = np.genfromtxt(data_cpg_filepath, delimiter=',')
            t = np.arange(TIME1 + TIME2)
            ax.plot(t, data_cpg,
                    color=colors[i_kn2])
                
            data_pack_filepath = dir2 + "packaging_p_off{}_k+{}k-{}_{}examples.csv".format(round(p_off,4),
                                                                                           round(kp, 4),
                                                                                           round(km, 4),
                                                                                           example_n)
            
                
            data_pack = np.genfromtxt(data_pack_filepath, delimiter=',', dtype=np.int8)
            t = np.arange(TIME1, TIME1+TIME2)
            p = bx.plot(t, data_pack,
                            '-',
                            color = colors[i_kn2],
                            label='knuc2: {}'.format(round(kn2, 4)),
                )
            plot_list.append(p[0])


        l = plt.legend(handles=plot_list[::4], loc=1, fontsize=6)
        bx = plt.gca().add_artist(l)

    title = "packaging_over_poff{}_over_knuc_withcolors.pdf".format(p_off)
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()

if __name__ == "__main__":
    main()
