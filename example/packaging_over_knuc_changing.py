import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy.optimize import curve_fit

dir1 = "example/data/dna_model/"
example_n = 100
kp = 0.145
km = 0.145

kn1 = 1


def func(t, a, b, c):
    return a * np.exp(-b * t) + c


def main():
    poff_list = [0.0, 0.0025, 0.005, 0.0075]
    ka_list = [0]
    for poff in poff_list:
        for ka in ka_list:
            submain(poff, ka)


def submain(poff, ka):
    # ka = 0
    ka1 = ka
    ka2 = ka

    kn2_list = np.arange(0, 0.5, 0.025)
    # kn2_list=[0.0]
    variation = len(kn2_list)
    cm_subsection = np.linspace(0.0, 1.0, 2 * variation)
    colors = [cm.jet(x) for x in cm_subsection]

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    plot_list1 = []
    plot_list2 = []
    # plot_list_curve = []

    for i, kn2 in enumerate(kn2_list):
        # read data
        print("read data...")
        filename = dir1 + "data{0:d}/kn{1}ka{2}_kn{3}ka{4}/packaging_p_off{5}_k+{6}k-{7}_{8:d}examples.csv".format(
            example_n,
            kn1,
            ka1,
            round(kn2,
                  4),
            ka2,
            poff,
            kp,
            km,
            example_n)
        data = np.genfromtxt(filename)
        x = np.arange(504)

        p1 = ax1.plot(x, data,
                      color=colors[i],
                      label="k nuc {}".format(round(kn2,4)))
        plot_list1.append(p1[0])
        try:
            fit_param, fit_covariance = curve_fit(func, x, data)
            p1 = ax1.plot(x, func(x, fit_param[0], fit_param[1], fit_param[2]),
                          color=colors[variation + i],
                          label=fit_param)
            plot_list2.append(p1[0])
            # plot_list_curve.append(fit_paramj
        except RuntimeError:
            print("error curve fit failure")
    ax1.set_yscale("log")

    ax1.set_xlim(0, 504)
    ax1.set_xticks(np.arange(0, 504, 168))
    ax1.set_xticklabels([str(i + 1) + " week" for i in range(3)])

    plt.legend(plot_list1[::4], ["k nuc {}".format(round(i, 4)) for i in np.arange(0, 0.5, 0.1)], fontsize=6)
    title = "test_dna_model_poff{}_packaging_over_knuc_withcolors__log__ka{}.pdf".format(poff, ka)
    pp = PdfPages(title)
    pp.savefig(fig1)
    pp.close()

    ax1.set_yscale("linear")
    ax1.set_ylim(0, example_n + 1)
    ax1.set_yticks(np.arange(0, example_n + 1, example_n // 10))
    ax1.set_yticklabels([str(i) + "%" for i in range(0, 101, 10)])
    l = plt.legend(handles=plot_list1[::4], loc=1, fontsize=6)
    ax = plt.gca().add_artist(l)
    plt.legend(handles=plot_list2, loc=2, fontsize=4)
    title = "test_dna_model_poff{}_packaging_over_knuc_withcolors__ka{}.pdf".format(poff, ka)
    pp = PdfPages(title)
    pp.savefig(fig1)
    pp.close()


if __name__ == "__main__":
    main()
