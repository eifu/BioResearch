import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from mpl_toolkits.mplot3d import Axes3D


dir1 = "example/data/withNUC/"
example_n = 1000


def main():

    kp = 0.145
    km = 0.145

    kn1 = 1
    ka1 = 0

    kn2 = 0
    ka2 = 0

    # read data
    print("read data...")
    filename = dir1 + "data{0:d}/kn{1}ka{2}_kn{3}ka{4}/packaging__k+{5}k-{6}_{7:d}examples.csv".format(example_n,
                                                                                           kn1,
                                                                                           ka1,
                                                                                           kn2,
                                                                                           ka2,
                                                                                           kp,
                                                                                           km,
                                                                                           example_n)
    data = np.genfromtxt(filename)
    x = np.arange(504)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x,data)
    # ax.set_ylim(0, example_n+1)
    # ax.set_yticks(np.arange(0,example_n+1,example_n//10))
    # ax.set_yticklabels([str(i)+"%" for i in range(0, 101, 10)])
    ax.set_yscale("log")
    ax.set_xlim(0, 504)
    ax.set_xticks(np.arange(0, 504, 168))
    ax.set_xticklabels([str(i+1) + " week" for i in range(3)])

    plt.show()

    title = "packaging_over_knuc.pdf"
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()

if __name__ == "__main__":
    main()
