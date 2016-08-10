import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

dir1 = "example/data/withNUC/"
example_n = 1000


def main():
    ka_list = [0, 0.025, 0.05, 0.075, 0.1]
    for ka in ka_list:
        submain(ka)


def submain(ka):
    kp = 0.145
    km = 0.145

    kn1 = 1
    ka1 = ka

    kn2_list = np.arange(0, 0.5, 0.025)
    variation = len(kn2_list)
    ka2 = ka
    cm_subsection = np.linspace(0.0, 1.0, variation)
    colors = [cm.jet(x) for x in cm_subsection]

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    plot_list1 = []

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    plot_list2 = []

    for i, kn2 in enumerate(kn2_list):
        # read data
        print("read data...")
        filename = dir1 + "data{0:d}/kn{1}ka{2}_kn{3}ka{4}/packaging__k+{5}k-{6}_{7:d}examples.csv".format(example_n,
                                                                                                           kn1,
                                                                                                           ka1,
                                                                                                           round(kn2,
                                                                                                                 4),
                                                                                                           ka2,
                                                                                                           kp,
                                                                                                           km,
                                                                                                           example_n)
        data = np.genfromtxt(filename)
        x = np.arange(504)

        p1 = ax1.plot(x, data, color=colors[i])
        plot_list1.append(p1[0])

        p1 = ax2.plot(x, data, color=colors[i])
        plot_list2.append(p1[0])

    ax1.set_yscale("log")
    ax1.set_xlim(0, 504)
    ax1.set_xticks(np.arange(0, 504, 168))
    ax1.set_xticklabels([str(i + 1) + " week" for i in range(3)])
    plt.legend(plot_list1[::4], ["k nuc {}".format(round(i, 4)) for i in np.arange(0, 0.5, 0.1)], fontsize=6)
    title = "packaging_over_knuc_withcolors__log__ka{}.pdf".format(ka)
    pp = PdfPages(title)
    pp.savefig(fig1)
    pp.close()

    ax2.set_ylim(0, example_n + 1)
    ax2.set_yticks(np.arange(0, example_n + 1, example_n // 10))
    ax2.set_yticklabels([str(i) + "%" for i in range(0, 101, 10)])
    ax2.set_xlim(0, 504)
    ax2.set_xticks(np.arange(0, 504, 168))
    ax2.set_xticklabels([str(i + 1) + " week" for i in range(3)])
    plt.legend(plot_list2[::4], ["k nuc {}".format(round(i, 4)) for i in np.arange(0, 0.5, 0.1)], fontsize=6)
    title = "packaging_over_knuc_withcolors__ka{}.pdf".format(ka)
    pp = PdfPages(title)
    pp.savefig(fig2)
    pp.close()


if __name__ == "__main__":
    main()
