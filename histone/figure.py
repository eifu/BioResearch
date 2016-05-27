import numpy as np
import math


def plotHist(fig, vectorizedgenome_timeseries):

    hst_n = len(vectorizedgenome_timeseries[0][0])

    ax = fig.add_subplot(3, 1, 1)
    for time in range(len(vectorizedgenome_timeseries)):
        y_position_m = [i - 40 for i in range(hst_n) if vectorizedgenome_timeseries[time][0][i] == 1]
        x_position_m = np.array([1] * np.sum(vectorizedgenome_timeseries[time][0]))

        ax.plot(x_position_m * time,  # time implies the x coordinates.
                y_position_m,
                ",", color="blue")

        y_position_a = [i - 40 for i in range(hst_n) if vectorizedgenome_timeseries[time][2][i] == 1]
        x_position_a = np.array([1] * np.sum(vectorizedgenome_timeseries[time][2]))

        ax.plot(x_position_a * time,
                y_position_a,
                ",", color="red")

    ax.set_xlim(-0.5, len(vectorizedgenome_timeseries) - 0.5)
    ax.set_ylim(-40.5, 40.5)


def plotWindow(fig, vectorizedgenome_timeseries):

    hst_n = len(vectorizedgenome_timeseries[0][0])
    w = 10

    ax = fig.add_subplot(9, 1, 5)
    for time in range(len(vectorizedgenome_timeseries)):
        print(hst_n // 2 - w // 2, hst_n // 2 + w // 2)
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

def plotT(fig, TList):

    total_time = len(TList)
    ax = fig.add_subplot(9,1,4)
    time = [i for i in range(len(TList))]
    ax.plot(time, TList, "-", color="red")
    ax.set_yticks([])
    ax.set_ylim(-0.5,1.5)
    ax.set_xlim(-0.5,total_time+0.5)
    ax.set_xticks([])


def plotStatistics(fig, vectorizedgenome_timeseries):
    bx = fig.add_subplot(9,4,22)

    w = 10
    hst_n = len(vectorizedgenome_timeseries[0][0])
    time = len(vectorizedgenome_timeseries)
    delta = 5

    count_m = [0 for _ in range(w+1)]

    count = 0
    for h in range(hst_n//2 - w//2,hst_n//2 + w//2+1):
        for time in range(time//2, time, delta):
            if vectorizedgenome_timeseries[time][0][h] == 1:
                count_m[count] += 1

        count += 1
    xaxis = [i for i in range(-w//2,w//2+1)]
    bx.barh(xaxis, count_m, align="center")
    bx.set_xticks([])
    bx.set_yticks([])

    acc = 0

    for each in count_m:
        acc += each
    AM = acc / len(count_m)

    acc = 0

    for each in count_m:
        acc += (each - AM)**2
    SD = math.sqrt(acc/len(count_m))
    return {"AM": AM, "SD": SD}


def plotStatistics2(fig, vectorizedgenome_timeseries):
    bx = fig.add_subplot(9,4,24)

    w = 10
    count_m = [0 for _ in range(w+1)]
    hst_n = len(vectorizedgenome_timeseries[0][0])
    time = len(vectorizedgenome_timeseries)
    delta = 5

    count = 0
    for h in range(hst_n//2- w//2, hst_n//2+ w//2+1):
        for time in range(time//2, time, delta):
            if vectorizedgenome_timeseries[time][0][h] == 1:
                count_m[count] += 1

        count += 1
    xaxis = [i for i in range(-w//2, w//2+1)]
    bx.barh(xaxis, count_m, align="center")
    bx.set_xticks([])
    bx.set_yticks([])

    acc = 0
    for each in count_m:
        acc += each
    AM = acc / len(count_m)

    acc =0
    for each in count_m:
        acc += (each - AM)**2

    SD = math.sqrt(acc/len(count_m))

    return {"AM": AM, "SD": SD}


def plotStatistics2(fig, vectorizedgenome_timeseries):
    bx = fig.add_subplot(9,4,24)

    hst_n = len(vectorizedgenome_timeseries[0][0])
    time = len(vectorizedgenome_timeseries)
    w = 10
    delta = 5

    count_m = [0 for _ in range(w+1)]

    count = 0
    for h in range(hst_n//2- w//2, hst_n//2+ w//2+1):
        for time in range(time//2, time, delta):
            if vectorizedgenome_timeseries[time][0][h] == 1:
                count_m[count] += 1

        count += 1
    xaxis = [i for i in range(-w//2,w//2+1)]
    bx.barh(xaxis, count_m, align="center")
    bx.set_xticks([])
    bx.set_yticks([])

    acc = 0
    for each in count_m:
        acc += each
    AM = acc / len(count_m)

    acc =0
    for each in count_m:
        acc += (each - AM)**2

    SD = math.sqrt(acc/len(count_m))

    return {"AM":AM,"SD":SD}