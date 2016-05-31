import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
import numpy as np


def sequence(fig, vectgene_timeseries):
    hst_n = len(vectgene_timeseries[0][0])

    ax = fig.add_subplot(3, 1, 1)
    for time in range(len(vectgene_timeseries)):
        y_position_m = [i - 40 for i in range(hst_n) if vectgene_timeseries[time][0][i] == 1]
        x_position_m = np.array([1] * np.sum(vectgene_timeseries[time][0]))

        ax.plot(x_position_m * time,  # time implies the x coordinates.
                y_position_m,
                ",", color="blue")

        y_position_a = [i - 40 for i in range(hst_n) if vectgene_timeseries[time][2][i] == 1]
        x_position_a = np.array([1] * np.sum(vectgene_timeseries[time][2]))

        ax.plot(x_position_a * time,
                y_position_a,
                ",", color="red")

    ax.set_xlim(-0.5, len(vectgene_timeseries) - 0.5)
    ax.set_ylim(-40.5, 40.5)
    ax.set_xticks([w for w in range(0, 2017, 168)])
    ax.set_xticklabels(("week" + str(w) for w in range(1, len(vectgene_timeseries) // 168 + 1)))


def window(fig, vectorizedgenome_timeseries):
    hst_n = len(vectorizedgenome_timeseries[0][0])
    w = 10

    ax = fig.add_subplot(9, 1, 5)
    for time in range(len(vectorizedgenome_timeseries)):
        y_position_m = [i - 40 for i in range(hst_n // 2 - w // 2, hst_n // 2 + w // 2 + 1)
                        if vectorizedgenome_timeseries[time][0][i] == 1]
        x_position_m = np.array([1] * len(y_position_m))

        ax.plot(x_position_m * time,
                y_position_m,
                ",", color="blue")

        y_position_a = [i - 40 for i in range(hst_n // 2 - w // 2, hst_n // 2 + w // 2 + 1)
                        if vectorizedgenome_timeseries[time][2][i] == 1]
        x_position_a = np.array([1] * len(y_position_a))

        ax.plot(x_position_a * time,
                y_position_a,
                ",", color="red")

        ax.set_xlim(-0.5, len(vectorizedgenome_timeseries))
        ax.set_ylim(-5, 5)
        ax.set_yticks([])
        ax.set_xticks([])


def transcription(fig, TList):
    total_time = len(TList)
    ax = fig.add_subplot(9, 1, 4)
    time = [i for i in range(len(TList))]
    ax.plot(time, TList, "-", color="red")
    ax.set_yticks([])
    ax.set_ylim(-0.5, 1.5)
    ax.set_xlim(-0.5, total_time + 0.5)
    ax.set_xticks([])


def m_stat(fig, vectorizedgenome_timeseries):
    bx = fig.add_subplot(9, 4, 22)

    w = 10
    hst_n = len(vectorizedgenome_timeseries[0][0])
    time = len(vectorizedgenome_timeseries)
    delta = 5

    count_m = [0 for _ in range(w + 1)]

    count = 0
    for h in range(hst_n // 2 - w // 2, hst_n // 2 + w // 2 + 1):
        for time in range(time // 2, time, delta):
            if vectorizedgenome_timeseries[time][0][h] == 1:
                count_m[count] += 1

        count += 1
    xaxis = [i for i in range(-w // 2, w // 2 + 1)]
    bx.barh(xaxis, count_m, align="center")
    bx.set_xticks([])
    bx.set_yticks([])

    acc = 0

    for each in count_m:
        acc += each
    AM = acc / len(count_m)

    acc = 0

    for each in count_m:
        acc += (each - AM) ** 2
    SD = math.sqrt(acc / len(count_m))
    return {"AM": AM, "SD": SD}


def m_stat2(fig, vectorizedgenome_timeseries):
    bx = fig.add_subplot(9, 4, 24)

    w = 10
    count_m = [0 for _ in range(w + 1)]
    hst_n = len(vectorizedgenome_timeseries[0][0])
    time = len(vectorizedgenome_timeseries)
    delta = 5

    count = 0
    for h in range(hst_n // 2 - w // 2, hst_n // 2 + w // 2 + 1):
        for time in range(time // 2, time, delta):
            if vectorizedgenome_timeseries[time][0][h] == 1:
                count_m[count] += 1

        count += 1
    xaxis = [i for i in range(-w // 2, w // 2 + 1)]
    bx.barh(xaxis, count_m, align="center")
    bx.set_xticks([])
    bx.set_yticks([])

    acc = 0
    for each in count_m:
        acc += each
    AM = acc / len(count_m)

    acc = 0
    for each in count_m:
        acc += (each - AM) ** 2

    SD = math.sqrt(acc / len(count_m))

    return {"AM": AM, "SD": SD}


def cc(arg):
    return colorConverter.to_rgba(arg, alpha=0.6)


def dynamic_change(fig, list_vectorized_gene_timeseries, delta=2):
    hst_n = len(list_vectorized_gene_timeseries[0][0][0])  # default 81
    w = 10  # window default 10
    time = len(list_vectorized_gene_timeseries[0])
    example_n = len(list_vectorized_gene_timeseries)

    table_m = []
    table_m_am = []
    table_m_sd = []
    table_a = []
    table_a_am = []
    table_a_sd = []

    days = [0, 3, 5]
    for day in days:
        container_m = np.array([0 for _ in range(0, hst_n, delta)])
        container_a = np.array([0 for _ in range(0, hst_n, delta)])

        for vectorized_gene_timeseries in list_vectorized_gene_timeseries:
            container_m += vectorized_gene_timeseries[time // 2 + day * 24][0][0::delta]
            container_a += vectorized_gene_timeseries[time // 2 + day * 24][2][0::delta]

        m_am = container_m / example_n
        a_am = container_a / example_n
        container_m_err = np.array([0.0 for _ in range(0, hst_n, delta)])
        container_a_err = np.array([0.0 for _ in range(0, hst_n, delta)])
        for vectorized_gene_timeseries in list_vectorized_gene_timeseries:
            container_m_err += pow(abs(m_am - vectorized_gene_timeseries[time // 2 + day * 24][0][0::delta]), 2)
            container_a_err += pow(abs(a_am - vectorized_gene_timeseries[time // 2 + day * 24][2][0::delta]), 2)
        m_sd = np.sqrt(container_m_err / example_n)
        a_sd = np.sqrt(container_a_err / example_n)
        table_m.append(container_m)
        table_m_am.append(m_am)
        table_m_sd.append(m_sd)
        table_a.append(container_a)
        table_a_am.append(a_am)
        table_a_sd.append(a_sd)

    print("am")
    for am in table_a_am:
        print(am)
    print("sd")
    for sd in table_a_sd:
        print(sd)

    x = np.arange(0, hst_n, delta)
    pos_m = [3, 5, 7]
    pos_a = [4, 6, 8]
    for i in [0, 1, 2]:
        bx = fig.add_subplot(4, 2, pos_m[i])
        bx.fill_between(x, table_m[i], 0, color="blue")
        bx.errorbar(x, table_m[i], fmt=',-', elinewidth=0.2, capsize=1, label="data", mfc="lightblue",linestyle="none", color="blue", yerr=table_m_sd[i], ecolor="black")

        bx.set_yticks([i for i in range(0, example_n + 1, example_n // 5)])
        bx.set_yticklabels(("0%", "20%", "40%", "60%", "80%", "100%"))
        bx.set_ylim(0, example_n + 1)
        bx.set_xticks([0, 35, 40, 45, 81])
        bx.set_xticklabels((-40, -5, 0, 5, 40))
        bx.set_ylabel("day " + str(days[i]))
        bx.legend(loc='upper right', fontsize="x-small")

        bx = fig.add_subplot(4, 2, pos_a[i])
        bx.fill_between(x, table_a[i], 0, color="red")
        bx.errorbar(x, table_a[i], fmt=',-', elinewidth=0.2, capsize=1, label="data", mfc="tomato", linestyle="none", color="red", yerr=table_a_sd[i], ecolor="black")

        bx.set_yticks([i for i in range(0, example_n + 1, example_n // 5)])
        bx.set_yticklabels(("0%", "20%", "40%", "60%", "80%", "100%"))
        bx.set_ylim(0, example_n + 1)
        bx.set_xticks([0, 35, 40, 45, 81])
        bx.set_xticklabels((-40, -5, 0, 5, 40))
        bx.legend(loc='upper right', fontsize='x-small')

    ax = fig.add_subplot(4, 2, 1, projection='3d')
    zs = [0, 3, 5]
    verts = []
    for i in [0, 1, 2]:
        ys = table_m[i]
        ys[0], ys[-1] = 0, 0
        verts.append(list(zip(x, ys)))

    poly = PolyCollection(verts, facecolors=[cc('b'), cc('b'), cc('b')])

    poly.set_alpha(0.7)
    ax.add_collection3d(poly, zs=zs, zdir='y')

    ax.set_xlabel('X genome')
    ax.set_xlim3d(0, 81)
    ax.set_xticks([0, 35, 40, 45, 81])
    ax.set_xticklabels((-40, -5, 0, 5, 40))
    ax.set_ylabel('Y days')
    ax.set_ylim3d(-1, 6)
    ax.set_yticks([0, 3, 5])
    ax.set_yticklabels((0, 3, 5))
    ax.set_zlabel('Z freq')
    ax.set_zlim3d(0, example_n)
    ax.set_zticks([i for i in range(0, example_n + 1, example_n // 5)])
    ax.set_zticklabels(("0%", "20%", "40%", "60%", "80%", "100%"))

    ax = fig.add_subplot(4, 2, 2, projection='3d')
    verts = []
    for i in [0, 1, 2]:
        ys = table_a[i]
        ys[0], ys[-1] = 0, 0
        verts.append(list(zip(x, ys)))

    poly = PolyCollection(verts, facecolors=[cc('r'), cc('r'), cc('r')])

    poly.set_alpha(0.4)
    ax.add_collection3d(poly, zs=zs, zdir='y')

    ax.set_xlabel('X genome')
    ax.set_xlim3d(0, 81)
    ax.set_xticks([0, 35, 40, 45, 81])
    ax.set_xticklabels((-40, -5, 0, 5, 40))
    ax.set_ylabel('Y days')
    ax.set_ylim3d(-1, 6)
    ax.set_yticks([0, 3, 5])
    ax.set_yticklabels((0, 3, 5))
    ax.set_zlabel('Z freq')
    ax.set_zlim3d(0, example_n)
    ax.set_zticks([i for i in range(0, example_n + 1, example_n // 5)])
    ax.set_zticklabels(("0%", "20%", "40%", "60%", "80%", "100%"))
