import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

example_n = 200

dir = "example/data200_withNUC/"


def main():
    ka_list = np.arange(0.1,0.91,0.05)
    kn_list = np.arange(0.1,0.91,0.05)

    km = 0.08
    kp = 0.19

    plt.style.use('ggplot')
    fig = plt.figure()

    ax = fig.add_subplot(1,1, 1, projection="3d")

    x = np.arange(len(kn_list))
    y = np.arange(len(ka_list))
    hist_acc = np.zeros((17,17))  # (y, x) order

    for i_ka, ka in enumerate(ka_list):
        for i_kn, kn in enumerate(kn_list):
            filename = dir + "kn{}__ka{}/".format(round(kn,4), round(ka,4)) + "__k-{}/dumpdata2d__k+{}__{}examples.csv".format(round(km,4), round(kp,4),
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
    ax.set_xlabel('kn ')
    ax.set_xticks(np.arange(20))
    ax.set_xticklabels(round(i,4) for i in np.arange(0.1, 0.91, 0.05))
    ax.set_ylabel('ka ')
    ax.set_yticks(np.arange(20))
    ax.set_yticklabels(round(i,4) for i in np.arange(0.1,0.91,0.05))
    ax.text2D(0.05, 0.95, "km:{}, kp:{}".format(km,kp), transform=ax.transAxes)

    plt.show()


if __name__ == "__main__":
    main()