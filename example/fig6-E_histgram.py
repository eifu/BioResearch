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
    for km in km_list:

        filename = str(km)+"/fig6_regress_amax__gamma__beta.csv"
        print("reading dump file .. ")
        data2d = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')
        print('done reading dump file..')

    x, y = np.random.rand(2, 100) * 4
    hist, xedges, yedges = np.histogram2d(x, y, bins=4)

    elements = (len(xedges) - 1) * (len(yedges) - 1)
    xpos, ypos = np.meshgrid(xedges[:-1] + 0.25, yedges[:-1] + 0.25)

    xpos = xpos.flatten()
    ypos = ypos.flatten()
    zpos = np.zeros(elements)
    dx = 0.5 * np.ones_like(zpos)
    dy = dx.copy()
    dz = hist.flatten()

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')

    plt.show()


if __name__ == "__main__":
    main()
