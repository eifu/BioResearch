import numpy as np
import matplotlib.pyplot as plt
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
    plt.plot(x,data)
    plt.show()

if __name__ == "__main__":
    main()
