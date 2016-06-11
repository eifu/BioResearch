"""
figure 6 C, E

"""
import histone
import histone.io as io
import histone.figure as figure
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D
import numpy as np



NUM_OF_HISTONE = 81
WINDOW = 10
TIME1 = 504  # 3 week in hour
TIME2 = 504  # 3 week in hour
DELTA = 1

NUMEXAMPLE = 5


def main():
    k_minus = 0.3

    filename = "__k-{}/dumpdata3d__k-{}__{}examples.csv".format(k_minus, k_minus, NUMEXAMPLE)
    read_data = io.read_dump(filename,TIME2)
    plt.style.use('ggplot')
    fig = plt.figure()
    figure.figure6c_and_6e(fig, read_data, NUMEXAMPLE)

    plt.show()

    # title = "fig6-CE/fig__{}examples__k-{}.pdf".format(NUMEXAMPLE, k_minus)
    # pp = PdfPages(title)
    # pp.savefig(fig)
    # pp.close()

if __name__ == "__main__":
    main()
