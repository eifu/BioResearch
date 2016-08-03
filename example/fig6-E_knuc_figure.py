import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

example_n = 100
dir1 = "example/data{}_withNUC_strongmemory/".format(example_n)


def main():
    km = 0.145
    kp = 0.145

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    ka1 = 0
    ka2 = 0
    kn_list = np.arange(0.001, 0.051, 0.001)
    for i_kn, kn2 in enumerate(kn_list):
        kn1 = 0
        dir2 = dir1 +"kn{}ka{}_kn{}ka{}/".format(round(kn1, 4), round(ka1, 4), round(kn2, 4),
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
        dir2 = dir1 + "kn{}ka{}_kn{}ka{}/".format(round(kn1, 4), round(ka1, 4), round(kn2, 4),
                                                                  round(ka2, 4))

        filename = dir2 + "dumpdata2d__k+{}k-{}_{}examples.csv".format(km, kp, example_n)

        data2d = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')

        acc = 0
        for row in data2d:
            for el in row[35:46]:
                acc += el

        ave = acc * 100 / (504 * 11 * example_n)

        p2 = ax.plot(i_kn,ave,".", color="red")

    ax.set_xticks(np.arange(len(kn_list)))
    ax.set_xticklabels(kn_list)

    ax.set_ylim(-5, 105)
    ax.set_yticks(np.arange(0,101,20))
    ax.set_yticklabels((str(i)+"%" for i in range(0,101,20)))
    plt.legend((p1[0], p2[0]), ('initial knuc 0', 'initial knuc 1'))
    plt.show()

    title = "example/figure/test_strongmemory_locus.pdf"


    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()

if __name__ == "__main__":
    main()