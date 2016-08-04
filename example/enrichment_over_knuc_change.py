import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


dir1 = "example/data/withNUC/"
example_n = 500


def main():

    kp = 0.145
    km = 0.145

    kn1 = 1
    ka1 = 0

    kn2 = 0
    ka2 = 0

    # read data
    print("read data...")
    filename = dir1 + "data{}/kn{}ka{}_kn{}ka{}/dumpdata2d__k+{}k-{}_{}examples.csv".format(example_n,
                                                                                        round(kn1, 4),
                                                                                        round(ka1, 4),
                                                                                        round(kn2, 4),
                                                                                        round(ka2, 4),
                                                                                        round(kp, 4),
                                                                                        round(km, 4),
                                                                                        example_n)
    hist = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')

    # create a graph
    print("create graph...")
    plt.style.use('ggplot')
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection="3d")

    x = np.arange(81)
    y = np.arange(24 * 7 * 3)

    elements = (len(x)) * (len(y))
    xpos, ypos = np.meshgrid(x, y)

    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)
    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz_acc = hist.flatten()


    ax.set_xlabel("histone pos")
    ax.set_xticks((0, 35, 40, 45, 80))
    ax.set_xticklabels((-40, -5, 0, 5, 40))
    ax.set_ylabel("time")
#    ax.set_yscale("log")
#    ax.set_zscale("log")
    ax.text2D(0.05, 0.95, "km:{}, kp:{}".format(km, kp), transform=ax.transAxes)
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz_acc, color='b')
    plt.show()
    # save a graph
    # print("save graph...")


if __name__ == "__main__":
    main()
