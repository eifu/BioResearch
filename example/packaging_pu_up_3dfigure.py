import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D

example_n = 500

dir = "example/data/withNUC/data{}/packagedata/".format(example_n)


def main():
    km = 0.145
    kp = 0.145

    kn1 = 0
    ka1 = 0

    kn2_list = np.arange(0, 0.4, 0.025)
    ka2_list = np.arange(0, 0.4, 0.025)

    pu_matrix = np.zeros((16, 16))
    up_matrix = np.zeros((16, 16))

    for i_kn2, kn2 in enumerate(kn2_list):
        for i_ka2, ka2 in enumerate(ka2_list):
            filename = dir + "package_data__kn{}ka{}_kn{}ka{}__k+{}k-{}_{}examples.csv".format(kn1,
                                                                                               ka1,
                                                                                               round(kn2, 4),
                                                                                               round(ka2, 4),
                                                                                               kp,
                                                                                               km,
                                                                                               example_n)

            data = np.genfromtxt(filename,
                                 skip_header=0,
                                 skip_footer=0,
                                 delimiter=',',
                                 dtype=np.int64
                                 )
            sum_pp_and_pu = np.array([0, 0])
            sum_up_and_uu = np.array([0, 0])

            for row in data:
                sum_pp_and_pu += row[0:2]
                sum_up_and_uu += row[2:4]

            tmp = sum_pp_and_pu[1] / (sum_pp_and_pu[0] + sum_pp_and_pu[1])
            if tmp != 0:

                pu_matrix[i_kn2][i_ka2] = -np.log(tmp)
            else:
                pu_matrix[i_kn2][i_ka2] = 0

            tmp = sum_up_and_uu[0] / (sum_up_and_uu[0] + sum_up_and_uu[1])
            if tmp != 0:

                up_matrix[i_kn2][i_ka2] = -np.log(sum_up_and_uu[0] / (sum_up_and_uu[0] + sum_up_and_uu[1]))
            else:
                up_matrix[i_kn2][i_ka2] = 0
            print("kn{0:.3f}ka{1:.3f}, pu{2:.3f},up{3:.3f}".format(kn2, ka2, sum_pp_and_pu[1] / (
            sum_pp_and_pu[0] + sum_pp_and_pu[1]), sum_up_and_uu[0] / (sum_up_and_uu[0] + sum_up_and_uu[1])))
    x = np.arange(16)
    y = np.arange(16)
    elements = (len(x)) * (len(y))
    xpos, ypos = np.meshgrid(x, y)
    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)
    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111, projection="3d")
    ax1.set_xlabel("k ace")
    ax1.set_ylabel("k nuc")
    ax1.set_xticks(np.arange(16))
    ax1.set_xticklabels(kn2_list, fontsize=4)
    ax1.set_yticks(np.arange(16))
    ax1.set_yticklabels(ka2_list, fontsize=4)
    dz_pu = pu_matrix.flatten()
    ax1.bar3d(xpos, ypos, zpos, dx, dy, dz_pu, color='b')
    ax1.set_label("p->u")
    title = "fig6_pu.pdf"
    pp = PdfPages(title)
    pp.savefig(fig1)
    pp.close()

    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111, projection="3d")
    ax2.set_xlabel("k ace")
    ax2.set_ylabel("k nuc")
    ax2.set_xticks(np.arange(16))
    ax2.set_xticklabels(ka2_list, fontsize=4)
    ax2.set_yticks(np.arange(16))
    ax2.set_yticklabels(kn2_list, fontsize=4)
    ax2.set_label("u->p")
    dz_up = up_matrix.flatten()
    ax2.bar3d(xpos, ypos, zpos, dx, dy, dz_up, color='b')

    plt.show()

    title = "fig6_up.pdf"
    pp = PdfPages(title)
    pp.savefig(fig2)
    pp.close()


if __name__ == "__main__":
    main()
