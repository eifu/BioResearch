import histone
import histone.io as io
import histone.figure as figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm
from scipy import stats
from matplotlib import lines


example_n = 5000
dir = "data"+str(example_n)




def main():
    kp_list = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25, 0.3]
    km_list = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25, 0.3]
    plt.style.use('ggplot')
    fig = plt.figure()
    ax = fig.add_subplot(221, projection="3d")

    bx = fig.add_subplot(222, projection='3d')
    # bx.set_zlim3d(0.1,0.3)
    bx.set_xlabel("k+")
    bx.set_ylabel("k-")
    cx = fig.add_subplot(223, projection='3d')
    ccx = fig.add_subplot(224, projection='3d')

    # x, y = np.random.rand(2, 100) * 4
    x = np.arange(24)
    y = np.arange(24)
    hist, xedges, yedges = np.histogram2d(x, y,bins=[24,24])

    hist_a0 = np.zeros_like(hist)
    hist_b = np.zeros_like(hist)
    hist_g = np.zeros_like(hist)
    hist_cc = np.zeros_like(hist)

    # print(hist, hist.shape)

    km_list = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25,0.3]
    for i_km,km in enumerate(km_list):
        filename = dir+"/linregress/fig6_k-{}__regress_amax__gamma__beta.csv".format(round(km,4))
        # print("reading dump file .. ")
        data2d = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')
        # print('done reading dump file..')

        # print(i)

        gamma_01to03 = []

        for ii, h in enumerate(hist_a0[i_km]):
            # hist_b[i][ii] =data2d[ii][2]
            if 0.1 < data2d[ii][1]  < 0.3:
                hist_g[i_km][ii] = data2d[ii][1]
                gamma_01to03.append(ii)

        for ii in gamma_01to03:
            hist_a0[i_km][ii] = data2d[ii][0] # store amax to the graph
            hist_b[i_km][ii] = - data2d[ii][2]
            hist_cc[i_km][ii] = data2d[ii][3]

        print("i_km", i_km, "km:", km_list[i_km], [(ii, round(kp_list[ii],4)) for ii in gamma_01to03])

    # data2d = np.genfromtxt("fig6_test6_k-0.11__regress_amax__gamma__beta.csv", skip_header=0, skip_footer=0, delimiter=',')
    # for ii,h in enumerate(hist_g[11]):
    #     hist_a[11][ii] = data2d[ii][0]
    #     hist_g[11][ii] = data2d[ii][1]
    #     hist_b[11][ii] = data2d[ii][2]

    elements = (len(x)) * (len(y))
    xpos, ypos = np.meshgrid(x, y)

    # xpos = xpos.flatten()
    # ypos = ypos.flatten()
    zpos = np.zeros(elements)

    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()

    # ax_s = ax.plot_surface(xpos,ypos,hist_a,rstride=1, cstride=1)
    # bx.plot_surface(xpos,ypos,hist_g,rstride=1, cstride=1)
    # cx.plot_surface(xpos,ypos,hist_b,rstride=1, cstride=1)
    # ccx.plot_surface(xpos,ypos,hist_cc,rstride=1, cstride=1)


    ll = []

    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)
    dz_a = hist_a0.flatten()

    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz_a, color='b')

    dz_g = hist_g.flatten()
    bx.bar3d(xpos, ypos, zpos, dx, dy, dz_g, color='b')

    dz_b = hist_b.flatten()
    cx.bar3d(xpos, ypos, zpos, dx, dy, dz_b, color='b')

    dz_cc = hist_cc.flatten()
    ccx.bar3d(xpos,ypos, zpos, dx,dy,dz_cc, color='b')

    plt.show()

    title = "fig6_histogram_test.pdf"
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()


if __name__ == "__main__":
    main()
