import numpy as np

def save_hst_timeseries(hst_timeseries, filename):
    time = len(hst_timeseries)
    kind = len(hst_timeseries[0])
    hst_n = len(hst_timeseries[0][0])

    new_array = np.zeros((time, hst_n))
    for i, t in enumerate(hst_timeseries):
        new_array[i] = t[0] - t[2]

    print('creating new array...')
    with open(filename, 'wb') as f:
        np.savetxt(f,
                   new_array,
                   fmt='%d',
                   delimiter=',',
                   newline='\n')
    print('created and saved')


def read_hstcsv(filename, time, kind=3, hst_n=81):
    print("reading arrays from file")
    data = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')
    hst_timeseries = np.zeros((time, kind, hst_n))
    print('finish reading file')
    for t, compressed_array in enumerate(data):
        hst_timeseries[t][0] = (compressed_array + compressed_array * compressed_array) / 2
        hst_timeseries[t][1] = compressed_array * compressed_array
        hst_timeseries[t][2] = pow((compressed_array - compressed_array * compressed_array) / 2, 2)
    return hst_timeseries

#
def compress_kplist_samplelist_hstseqts(kplist_samplelist_hstseqts):
    kp_n = len(kplist_samplelist_hstseqts)
    # example_n = len(kplist_samplelist_hstseqts[0])
    time = len(kplist_samplelist_hstseqts[0][0])
    hst_n = len(kplist_samplelist_hstseqts[0][0][0][0])
    compressed = np.zeros((kp_n, time, hst_n))
    for var, oneversion_list_vecgenetimeseries in enumerate(kplist_samplelist_hstseqts):
        for _, vecgenetimeseries in enumerate(oneversion_list_vecgenetimeseries):
            for t in range(time):
                compressed[var][t] += vecgenetimeseries[t][0]

    return compressed

def compress_onekp_samplelist_hstseqts(onekp_samplelist_hstseqts):
    time = len(onekp_samplelist_hstseqts[0])
    hst_n = len(onekp_samplelist_hstseqts[0][0][0])
    compressed = np.zeros((time, hst_n))
    for _, vecgenetimeseries in enumerate(onekp_samplelist_hstseqts):
        for t in range(time):
            compressed[t] += vecgenetimeseries[t][0]

    return compressed

def compress_last_week_hst_vec(last_week_hst_vec):
    locus_n = 11
    example_n = len(last_week_hst_vec)
    return last_week_hst_vec.reshape(example_n * 24 * 7, locus_n)

def compress_packaging_samplelist(packaging_list):
    time = len(packaging_list[0])
    compressed = np.zeros(time)
    for _, onesample in enumerate(packaging_list):
        for t in range(time):
            compressed[t] += onesample[t]
    return compressed

def write_dump2d_pos_kp(kplist_samplelist_hstseqts, filename, time_inhour):
    kp_n = len(kplist_samplelist_hstseqts)
    example_n = len(kplist_samplelist_hstseqts[0])
    hst_n = len(kplist_samplelist_hstseqts[0][0][0][0])
    dump_table = np.zeros((kp_n+1, hst_n))
    dump_table[0] = [i for i in range(-40, 41)]

    for kp, samplelist_vecgenetimeseries in enumerate(kplist_samplelist_hstseqts):
        for _, vecgene_timeseries in enumerate(samplelist_vecgenetimeseries):
            dump_table[kp+1] += vecgene_timeseries[time_inhour][0]

    dump_table[1:] = dump_table[1:] * 100 / example_n

    with open(filename, 'wb') as f:
        np.savetxt(f,
                   dump_table.transpose(),
                   fmt='%d',
                   delimiter=',',
                   newline='\n')


def write_dump3d_kp_time_hst(data3d, filename, time, kp_n=24, hst_n=81):
    """
    :param data3d: kp list ~> time series ~> hst seq
    """
    data2d = data3d.reshape(kp_n * time, hst_n) # reshape to 2d for saving data
    with open(filename, 'wb') as f:
        np.savetxt(f,
                   data2d,
                   fmt='%d',
                   delimiter=',',
                   newline='\n')

def write_dump2d_onekp_time_hst(data2d, filename, time, kp_n=24, hst_n=81):
    """
    :param data3d: kp list ~> time series ~> hst seq
    """
    # data2d = data3d.reshape(time, hst_n) # reshape to 2d for saving data
    with open(filename, 'wb') as f:
        np.savetxt(f,
                   data2d,
                   fmt='%d',
                   delimiter=',',
                   newline='\n')

def write_dump2d_final_hst_list(data2d, filename, hst_n):
    with open(filename, "wb") as f:
        np.savetxt(f,
                   data2d,
                   fmt='%d',
                   delimiter=',',
                   newline='\n')

def read_dump3d_kp_time_hst(filename, time, kp_n=24, hst_n=81):
    data2d = np.genfromtxt(filename, skip_header=0, skip_footer=0, delimiter=',')
    print('reading dump file..')
    return data2d.reshape(kp_n, time, hst_n)  # convert to 3d


def read_dump2d_onekp_time_hst(filename):
    return np.genfromtxt(filename, skip_footer=0,skip_header=0,delimiter=',')