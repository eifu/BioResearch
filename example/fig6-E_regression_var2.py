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
#     km_list = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25]
#     for i,km in enumerate(km_list):
#         print(km)
#         submain(round(km,4))
#         print(i)
#
# def submain(k_minus):
def main():
    k_minus = 0.11
    kp_list = [0.0001, 0.001] + [i for i in np.arange(0.01, 0.21, 0.01)] + [0.25,0.3]

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
        ax = fig.add_subplot(4, 6, 12 (i//6) + (i % 6) + 1)
        ax.set_xticks([0, 50, 100, 150, 200])
        ax.set_xticklabels((0, 50, 100, 150, 200), fontsize=4)
        ax.set_yticks([0])
        ax.set_yticklabels([0], fontsize=4)

        a0 = max(enrichment) + 0.0001

        bx = fig.add_subplot(8, 6, 12 * (i//6) + (i%6) + 7 )
        bx.set_xticks([])
        bx.set_yticks([])

        bx.plot(hours, container[i], "-", color=colors[i])

        rev_enrichment = [np.log((a0 - e) / a0) for e in enrichment]
        ax.plot(hours, rev_enrichment, ".-", color=colors[i])

        # # range 50 -> gray
        # range50 = list(range(50))
        # gamma50, beta50, _, _, _ = stats.linregress(range50, rev_enrichment[:50])
        #
        # func = lambda x: gamma50 * x + beta50
        # estimated_y = [func(x) for x in hours]
        # ax.plot(hours, estimated_y, '--', color='gray')
        #
        # # range 100 -> black
        # range100 = list(range(100))
        # gamma100, beta100, _, _, _ = stats.linregress(range100, rev_enrichment[:100])
        #
        # func = lambda x: gamma100 * x + beta100
        # estimated_y = [func(x) for x in hours]
        # ax.plot(hours, estimated_y, '--', color='black')

        r_list =[]

        range_list = [25, 50, 75, 100]
        co = ["#d3d3d3","#939393","#545454","#151515"]
        for index, r in enumerate(range_list):
            ran = list(range(r))
            g, b, r, _, _ = stats.linregress(ran, rev_enrichment[:r])
            func = lambda x: g * x + b
            estimated_y = [func(x) for x in hours]
            ax.plot(hours, estimated_y, '--', color=co[index])

            en_re = [a0 * (1 - np.exp(g*t + b)) for t in hours]
            bx.plot(hours, en_re, "--", color=co[index])

            r_list.append(r)

        # for a in np.arange(amax, amax * 10, 0.01):
        #
        #     z = np.array([np.log((a - en) / a) for i, en in enumerate(en50)])
        #
        #     slope, intercept, r_value, _, _ = stats.linregress(ran, z)
        #     if r_value > cmax:
        #         cmax = r_value
        #         amax = a
        #         gamma = slope
        #         beta = intercept


        # bx = fig.add_subplot(8, 6, 12 * (i//6) + (i%6) + 7 )
        # bx.set_xticks([])
        # bx.set_yticks([])
        #
        #
        # bx.plot(hours, container[i], "-", color=colors[i])

        # en_50 = [a0 * (1 - np.exp(gamma50*t + beta50)) for t in hours ]
        # bx.plot(hours, en_50, "--", color="gray")
        #
        #
        # bn_100 = [a0 * (1 - np.exp(gamma100*t + beta100)) for t in hours ]
        # ax.plot(hours, en_100, "--", color='black')



        # final_z = np.array([np.log((amax - en) / amax) for i, en in enumerate(enrichment)])
        # ax.plot(hours, final_z, "-", color=colors[-i])

        # func = lambda x: gamma * x + beta
        # estimated_y = [func(x) for x in hours]
        # ax.plot(hours, estimated_y, '--', color='black')
        # print("k+ {}  a0 {} amax {}   gamma {}  beta {}".format(kp_list[i],a0, 0, gamma50, beta50))

        t = "k+ {} \n 25:{} 50:{}\n 75:{} 100{}".format(kp_list[i],
                                                        round(r_list[0],3),
                                                        round(r_list[1],3),
                                                        round(r_list[2],3),
                                                        round(r_list[3],3))
        ax.set_title(t, fontsize=4)

        # result[i][0] = amax
        # result[i][1] = gamma
        # result[i][2] = beta

    title = "fig6_test5_regress___k-{}.pdf".format(k_minus)
    pp = PdfPages(title)
    pp.savefig(fig)
    pp.close()

    # with open("fig6_test_regress_amax__gamma__beta.csv".format(k_minus), 'wb') as f:
    #     np.savetxt(f,
    #                result,
    #                fmt='%f',
    #                delimiter=',',
    #                newline='\n')


if __name__ == "__main__":
    main()
