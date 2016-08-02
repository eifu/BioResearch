import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

example_n = 200

dir = "example/data500_withNUC/"


def main():
    ka_list = [0.1, 0.2]
    kn_list = [0.1, 0.2, 0.3]

    # ka_list = [0.1]
    # kn_list = [0.1]

    # ka = 0.1
    # kn = 0.1

    kmkppair = [(0.05, 0.11), (0.05, 0.12), (0.06, 0.14), (0.06, 0.15), (0.07, 0.16), (0.07, 0.17), (0.08, 0.19),
                (0.08, 0.2), (0.11, 0.25), (0.14, 0.3)]

    # kmkppair = [(0.05,0.11)]

    # km = 0.05
    # kp = 0.11

    plt.style.use('ggplot')
    fig = plt.figure()

    for i_kmkp, (km, kp) in enumerate(kmkppair):
        ax = fig.add_subplot(3,4, i_kmkp + 1, projection="3d")

        x = np.arange(3)
        y = np.arange(2)
        # hist, _, _ = np.histogram2d(x, y, bins=[3,3])
        # print(hist,hist.shape, len(hist))
        hist_acc = np.zeros((2,3))

        for i_ka, ka in enumerate(ka_list):
            for i_kn, kn in enumerate(kn_list):
                filename = dir + "kn{}__ka{}/".format(kn, ka) + "__k-{}/dumpdata2d__k+{}__{}examples.csv".format(km, kp,
                                                                                                             example_n)

                data2d = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')

                acc = 0
                for row in data2d:
                    for el in row:
                        acc += el

                ave = acc / (504 * 81 * example_n)

                hist_acc[i_ka][i_kn] = ave

        elements = (len(x)) * (len(y))
        xpos, ypos = np.meshgrid(x, y)

        xpos = xpos.flatten()
        ypos = ypos.flatten()
        zpos = np.zeros(elements)
        dx = 0.5 * np.ones_like(zpos)
        dy = dx.copy()
        dz_acc = hist_acc.flatten()

        ax.bar3d(xpos, ypos, zpos, dx, dy, dz_acc, color='b')
        ax.set_xlabel('kn {}'.format(round(kn,4)))
        ax.set_ylabel('ka {}'.format(round(ka,4)))
        ax.text2D(0.05, 0.95, "km:{}, kp:{}".format(km,kp), transform=ax.transAxes)

    plt.show()


if __name__ == "__main__":
    main()