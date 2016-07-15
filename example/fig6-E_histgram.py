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


def main():

    plt.style.use('ggplot')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")

    km_list = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25]
    # for km in km_list:
    #
    #     filename = "fig6_test6_k-0.11__regress_amax__gamma__beta.csv"
    #     print("reading dump file .. ")
    #     data2d = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')
    #     print('done reading dump file..')

    filename = "fig6_test6_k-0.11__regress_amax__gamma__beta.csv"
    print("reading dump file .. ")
    data2d = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')
    print('done reading dump file..')

    # x, y = np.random.rand(2, 100) * 4
    x = np.arange(24)
    y = np.arange(24)
    hist, xedges, yedges = np.histogram2d(x, y,bins=[24,24])

    hist = np.zeros_like(hist)
    print(hist)

    print(len(data2d))

    for i, h in enumerate(hist[11]):
        hist[11][i] = data2d[i][0]

    print(hist)
    print(xedges)
    print(yedges)
    elements = (len(x)) * (len(y) )
    xpos, ypos = np.meshgrid(x , y )

    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)

    print(xpos.shape,ypos.shape,zpos.shape)
    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = hist.flatten()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')

    plt.show()


if __name__ == "__main__":
    main()
