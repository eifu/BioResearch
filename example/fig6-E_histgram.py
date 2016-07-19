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


example_n = 500
dir = "data"+str(example_n)

def main():

    plt.style.use('ggplot')
    fig = plt.figure()
    ax = fig.add_subplot(221, projection="3d")
    bx = fig.add_subplot(222, projection='3d')
    cx = fig.add_subplot(223, projection='3d')

    # x, y = np.random.rand(2, 100) * 4
    x = np.arange(24)
    y = np.arange(24)
    hist, xedges, yedges = np.histogram2d(x, y,bins=[24,24])

    hist_a = np.zeros_like(hist)
    hist_b = np.zeros_like(hist)
    hist_g = np.zeros_like(hist)


    km_list = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25]
    # for i,km in enumerate(km_list):
    #     filename = dir+"/__k-"+str(km)+"/fig6_regress_amax__gamma__beta.csv"
    #     print("reading dump file .. ")
    #     data2d = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')
    #     print('done reading dump file..')
    #
    #     print(i)
    #     for ii, h in enumerate(hist_a[i]):
    #         # hist_a[i][ii] = data2d[ii][0] # store amax to the graph
    #        hist_g[i][ii] = (-1)*data2d[ii][1]
    #        # hist_b[i][ii] = data2d[ii][2]

    data2d = np.genfromtxt("fig6_test6_k-0.11__regress_amax__gamma__beta.csv", skip_header=0, skip_footer=0, delimiter=',')
    for ii,h in enumerate(hist_g[11]):
        hist_a[11][ii] = data2d[ii][0]
        hist_g[11][ii] = data2d[ii][1]
        hist_b[11][ii] = data2d[ii][2]

    elements = (len(x)) * (len(y))
    xpos, ypos = np.meshgrid(x, y)

    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)

    print(xpos.shape,ypos.shape,zpos.shape)
    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz_a = hist_a.flatten()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz_a, color='b', zsort='average')

    dz_g = hist_g.flatten()
    bx.bar3d(xpos, ypos, zpos, dx, dy, dz_g, color='b', zsort='average')

    dz_b = hist_g.flatten()
    cx.bar3d(xpos, ypos, zpos, dx, dy, dz_b, color='b', zsort='average')

    plt.show()

    title = "fig6_histogram_test.pdf"
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()


if __name__ == "__main__":
    main()
