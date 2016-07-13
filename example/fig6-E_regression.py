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
from matplotlib import cm
from scipy import stats


TIME1 = 504  # 3 week in hour
TIME2 = 504  # 3 week in hour

# def main():
#     km_list = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25, 0.3]
#     for i,km in enumerate(km_list):
#         submain(round(km,4))
#         print(i)

def main():
    
    example_n = 500
    k_minus = 0.11

    filename = "__k-{}/dumpdata3d__k-{}__{}examples.csv".format(k_minus, k_minus, example_n)
    read_data = io.read_dump3d_kp_time_hst(filename, TIME2)

    plt.style.use('ggplot')
    fig = plt.figure()
    # ax = fig.add_subplot(1,1,1)

    hours = np.arange(8 * 24)
    print(read_data.shape)
    container = np.zeros((24,8*24))
    for i, oneversion in enumerate(read_data):

        enrichment = np.zeros(8*24) # enrichment is the timeseries of enrichment of locus of histone methylation
        for t, mseq in enumerate(oneversion[:24*8]):
            enrichment[t] = sum(mseq[35:46])
        enrichment /= example_n # divided by example_n because mseq is the sum of sample lists
        # enrichment /= 11
        # print(enrichment)
        container[i] = enrichment

    cm_subsection = np.linspace(0.0, 1.0, 24)
    colors = [cm.jet(x) for x in cm_subsection]
    # reverse
    for i, enrichment in enumerate(container):
        ax = fig.add_subplot(4,6,i+1)
        # print(max(enrichment))
        container[i] = max(enrichment) - enrichment
        # print(container[i])
        ax.plot(hours, container[i], color=colors[i])
        a0 = max(enrichment)
        cmax = 0
        amax = 0
        gamma = 0
        for a in np.arange(a0+0.0001, a0*5, 0.01):
            z = np.array([np.log((a-enrichment[i])/a) for i in range(8*24)])

            slope, intercept, r_value, p_value, std_err = stats.linregress(hours, z)
            if r_value > cmax:
                cmax = r_value
                amax = a
                gamma = slope
        ax.set_ylim(0, 11)

        print("cmax {}   amax {}   gamma {}".format(cmax,amax,gamma))





    # figure.figure6c_and_6e(fig, read_data, example_n)

    # plt.show()

    title = "fig6-CE/fig__regress__{}examples__k-{}.pdf".format(example_n, k_minus)
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()

if __name__ == "__main__":
    main()
