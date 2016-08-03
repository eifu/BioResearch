import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

example_n = 300


def main():


    km = 0.145
    kp = 0.145

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    ka1 = 0
    ka2 = 0
    kn_list = np.arange(0.05, 0.96, 0.05)
    for i_kn, kn2 in enumerate(kn_list):
        kn1 = 0
        dir2 = "example/data{}_withNUC/kn{}ka{}_kn{}ka{}/".format(example_n,round(kn1, 4), round(ka1, 4), round(kn2, 4),
                                                                  round(ka2, 4))

        filename = dir2 + "dumpdata2d__k+{}k-{}_{}examples.csv".format(km, kp, example_n)

        data2d = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')

        acc = 0
        for row in data2d:
            for el in row[35:46]:
                acc += el

        ave = acc * 100 / (504 * 11 * example_n)

        p1 = ax.plot(i_kn,ave,".", color="blue")

        kn1 = 1
        dir2 = "example/data{}_withNUC/kn{}ka{}_kn{}ka{}/".format(example_n,round(kn1, 4), round(ka1, 4), round(kn2, 4),
                                                                  round(ka2, 4))

        filename = dir2 + "dumpdata2d__k+{}k-{}_{}examples.csv".format(km, kp, example_n)

        data2d = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')

        acc = 0
        for row in data2d:
            for el in row[35:46]:
                acc += el

        ave = acc * 100 / (504 * 11 * example_n)

        p2 = ax.plot(i_kn,ave,".", color="red")

    ax.set_xticks(np.arange(19))
    ax.set_xticklabels(kn_list)

    ax.set_ylim(-5, 105)
    ax.set_yticks(np.arange(0,101,20))
    ax.set_yticklabels((str(i)+"%" for i in range(0,101,20)))
    plt.legend((p1[0], p2[0]), ('initial knuc 0', 'initial knuc 1'))
    plt.show()

    title = "example/figure/test2.pdf"


    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()

if __name__ == "__main__":
    main()