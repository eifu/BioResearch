import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

example_n = 300
dir1 = "example/data/withNUC/data{}/".format(example_n)
hst_n = 81
kp = 0.145
km = 0.145
kn2_list = np.arange(0, 0.2, 0.01)  # length 20


def main():
    for week in range(0, 7):
        submain(week)
        print(week, " week is done")


def submain(week):

    ka = 0

    plt.style.use('ggplot')
    font = {'family': 'sans-serif'}
    matplotlib.rc('font', **font)
    fig1 = plt.figure()

    for i_kn2, kn2 in enumerate(kn2_list):

        kn1 = 0
        ax = fig1.add_subplot(20, 2, 2 * i_kn2 + 1)
        if i_kn2 == 0:
            ax.set_title("kn1: 0")
        ax.set_yticks([])
        if i_kn2 % 4 == 0:
            ax.set_ylabel(round(kn2, 4))
        ax.set_xticks([])
        filename = dir1 + "kn{}ka{}_kn{}ka{}/".format(round(kn1, 4),
                                                      round(ka, 4),
                                                      round(kn2, 4),
                                                      round(ka, 4)) + "final_hst_list_k+{}k-{}_{}examples.csv".format(
            kp,
            km, example_n)

        data = np.genfromtxt(filename, delimiter=',')
        data = data.reshape((example_n, 24 * 7 * 6, 11))

        container = np.zeros(12)  # from 0 to 11. total 12 paterns
        for one_sample in data:
            for _, row in enumerate(one_sample[24 * 7 * week: 24 * 7 * (week+1)]):
                container[sum(row)] += 1
        ax.set_ylim(0, max(container))
        ax.bar([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], container)

        kn1 = 1
        ax = fig1.add_subplot(20, 2, 2 * i_kn2 + 2)
        if i_kn2 == 0:
            ax.set_title("kn1: 1")
        ax.set_yticks([])
        if i_kn2 % 4 == 0:
            ax.set_ylabel(round(kn2, 4))
        ax.set_xticks([])
        filename = dir1 + "kn{}ka{}_kn{}ka{}/".format(round(kn1, 4),
                                                      round(ka, 4),
                                                      round(kn2, 4),
                                                      round(ka, 4)) + "final_hst_list_k+{}k-{}_{}examples.csv".format(
            kp,
            km,
            example_n)

        data = np.genfromtxt(filename, delimiter=',')

        data = data.reshape((example_n, 24 * 7 * 6, 11))

        container = np.zeros(12)  # from 0 to 11. total 12 paterns
        for one_sample in data:

            for _, row in enumerate(one_sample[24 * 7 * week: 24 * 7 * (week+1)]):
                container[sum((row + row * row) / 2)] += 1

        ax.set_ylim(0, max(container))
        ax.bar([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11], container)

        print(i_kn2 / 20)

    title = "week{}_histogram__ka{}__k+{}k-{}_N{}.pdf".format(week,round(ka, 4),
                                                             round(kp, 4),
                                                             round(km, 4),
                                                             example_n)
    pp = PdfPages(title)
    pp.savefig(fig1)
    pp.close()


if __name__ == "__main__":
    main()
