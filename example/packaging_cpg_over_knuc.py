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

    for i_p_off, p_off in enumerate(p_off_list):
        for i_kn1, kn1 in enumerate(kn1_list):
            plot_list = []
            ax = fig.add_subplot(7,2,i_p_off * 2 + i_kn1+1)
            ax.set_xlim(TIME1, TIME1 + TIME2)
            ax.set_xticks([i for i in range(TIME1, TIME1+TIME2, 24*7)])
            ax.set_xticklabels([str(i/(24*7)) + "week" for i in range(TIME1, TIME1+TIME2, 24 * 7)])
            ax.set_ylim(-0.5, 100.5)
            ax.set_yticks(np.arange(0, 101, 20))
            ax.set_yticklabels(str(i) + "%" for i in range(0, 101, 20))
            if i_kn1 == 0:
                ax.set_ylabel("p off: {}".format(round(p_off, 4)), fontsize=4)
            
            for i_kn2, kn2 in enumerate(kn2_list): 
                if i_p_off == 0 and i_kn2 == 0:
                    ax.set_title("kn1: {}".format(kn1))
                
                dir2 = dir1 + "kn{}ka{}_kn{}ka{}/".format(round(kn1, 4),
                                                        round(ka1, 4),
                                                        round(kn2, 4),
                                                        round(ka2, 4))
                
                data_pack_filepath = dir2 + "packaging_p_off{}_k+{}k-{}_{}examples.csv".format(round(p_off,4),
                                                                                               round(kp, 4),
                                                                                               round(km, 4),
                                                                                               example_n)
            
                
                data_pack = np.genfromtxt(data_pack_filepath, delimiter=',', dtype=np.int8)
                t = np.arange(TIME1, TIME1+TIME2)
                p = ax.plot(t, data_pack, color = colors[i_kn2] , label='knuc2: {}'.format(round(kn2, 4)))
                plot_list.append(p[0])
            l = plt.legend(handles=plot_list[::4], loc=1, fontsize=6)
            ax = plt.gca().add_artist(l)
                
    plt.show()


    title = "test_dna_model_over_poff_packaging_over_knuc_withcolors.pdf"
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()

if __name__ == "__main__":
    main()
    
        
        
