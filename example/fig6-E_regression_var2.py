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
from matplotlib import lines

TIME1 = 504  # 3 week in hour
TIME2 = 504  # 3 week in hour


# def main():
#     km_list = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25, 0.3]
#     for i,km in enumerate(km_list):
#         submain(round(km,4))
#         print(i)
#
# def submain(k_minus):
#     example_n = 500
def main():
    k_minus = 0.11

    example_n = 500
    filename = "__k-{}/dumpdata3d__k-{}__{}examples.csv".format(k_minus, k_minus, example_n)
    read_data = io.read_dump3d_kp_time_hst(filename, TIME2)

    plt.style.use('ggplot')
    fig = plt.figure()

    hours = np.arange(8 * 24)
    print(read_data.shape)
    container = np.zeros((24, 8 * 24))

    result = np.zeros((24,3),dtype=float)  # 24 versions of k+ with amax, gamma, betta

    for i, oneversion in enumerate(read_data):

        enrichment = np.zeros(8 * 24)  # enrichment is the timeseries of enrichment of locus of histone methylation
        for t, mseq in enumerate(oneversion[:24 * 8]):
            enrichment[t] = sum(mseq[35:46])
        enrichment /= example_n  # divided by example_n because mseq is the sum of sample lists
        # enrichment /= 11
        # print(enrichment)
        container[i] = enrichment

    # colors is a list of color data used for figures.
    colors = [cm.jet(x) for x in np.linspace(0.0, 1.0, 24)]

    for i, enrichment in enumerate(container):
        ax = fig.add_subplot(4, 6, i + 1)

        ax.plot(hours, container[i], "-", color=colors[i])
        a0 = max(enrichment) + 0.0001
        amax = max(enrichment) + 0.0001

        rev_enrichment = [np.log((amax - e) / amax) for e in enrichment]
        ran = list(range(50))
        gamma, beta, cmax, _, _ = stats.linregress(ran, rev_enrichment[:50])
        ax.plot(hours, rev_enrichment, ".-", color=colors[i])
        func = lambda x: gamma * x + beta
        estimated_y = [func(x) for x in hours]
        ax.plot(hours, estimated_y, '--', color='gray')


        en50 = enrichment[:50]  # experimentally set up 50 first examples

        for a in np.arange(amax, amax * 4, 0.001):

            z = np.array([np.log((a - en) / a) for i, en in enumerate(en50)])

            slope, intercept, r_value, _, _ = stats.linregress(ran, z)
            if r_value > cmax:
                cmax = r_value
                amax = a
                gamma = slope
                beta = intercept
        ax.set_xticks([0, 50, 100, 150, 200])
        ax.set_xticklabels((0, 50, 100, 150, 200), fontsize=4)
        ax.set_yticks([0])
        ax.set_yticklabels([0], fontsize=4)

        final_z = np.array([np.log((amax - en) / amax) for i, en in enumerate(enrichment)])
        ax.plot(hours, final_z, "-", color=colors[-i])

        func = lambda x: gamma * x + beta
        estimated_y = [func(x) for x in hours]
        ax.plot(hours, estimated_y, '--', color='black')
        print("cmax {}  a0 {} amax {}   gamma {}  beta {}".format(cmax,a0, amax, gamma, beta))

        ax.set_title(
            "cmax {} amax {}\n gamma {} beta{}".format(round(cmax, 4), round(amax, 4), round(gamma, 4), round(beta, 4)),
            fontsize=4)

        result[i][0] = amax
        result[i][1] = gamma
        result[i][2] = beta

    # figure.figure6c_and_6e(fig, read_data, example_n)

    # plt.show()


    title = "fig6-CE/fig___regress__{}examples__k-{}.pdf".format(example_n, k_minus)
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()


    with open("fig6___{}example__k-{}.csv".format(example_n,k_minus), 'wb') as f:
        np.savetxt(f,
                   result,
                   fmt='%f',
                   delimiter=',',
                   newline='\n')


if __name__ == "__main__":
    main()
