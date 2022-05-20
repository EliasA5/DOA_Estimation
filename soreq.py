####################################################################################################
#                                             soreq.py                                             #
####################################################################################################
#                                                                                                  #
# Authors: J. M.                                                                                   #
#                                                                                                  #
# Created: 28/02/22                                                                                #
#                                                                                                  #
# Purpose: Create a dataset from the SOREQ data recordings.                                        #
#                                                                                                  #
####################################################################################################



#*************#
#   imports   #
#*************#
import glob, os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tarfile

from functools import reduce

from obspy.core.inventory import read_inventory
from obspy.core.stream import read
from obspy.geodetics.base import gps2dist_azimuth
from obspy.geodetics import degrees2kilometers
from obspy.core.utcdatetime import UTCDateTime

from scipy.io import savemat

from tqdm import tqdm


#****************#
#   load files   #
#****************#
arrival = pd.read_csv('./geres.arrival.csv')
assoc = pd.read_csv('./geres.assoc.csv')
origin = pd.read_csv('./geres.origin.csv')

# concatenate with other csv files...
arrival = pd.concat([arrival, pd.read_csv('./geres.arrival-306-346.csv')], ignore_index=True)
assoc = pd.concat([assoc, pd.read_csv('./geres.assoc-306-346.csv')], ignore_index=True)
origin = pd.concat([origin, pd.read_csv('./geres.origin-306-346.csv')], ignore_index=True)

# print number of entries in both arrival and assoc for each day
days, nums = np.unique(arrival.JDATE[arrival.ARID.isin(assoc.ARID.unique())].values,
                       return_counts=True)
print([str(day) + ': ' + str(nums[i]) for i, day in enumerate(days)])

inv = read_inventory('GERES.xml')


#****************************#
#   loop through all files   #
#****************************#
files = [file for file in glob.glob("*.tar.gz")]
for file in tqdm(files):

    # load and clean stream
    stream = read(file)
    ref_sensor = inv.get_channel_metadata(stream[0].get_id())
    stream.merge(method=0, fill_value='latest')

    removed = []
    # delete the sensors with no metadata
    for i, tr in enumerate(stream):
        try:
            meta = inv.get_channel_metadata(tr.get_id())
            removed.append(i)
        except Exception:
            stream.remove(stream[i])
    print(removed)
    print('remaining sensors: ', stream.count())

    # # remove (sensitivity and) the response of the sensors from the data
    # stream.remove_sensitivity(inv)
    # stream.remove_response(inventory=inv, output="VEL")
    # stream.filter('lowpass_cheby_2', freq=10)

    # stream.plot()

    # # get angle to north
    # north = {'latitude': 90, 'longitude': 0, 'elevation': 0}  # north
    # print(gps2dist_azimuth(ref_sensor['latitude'], ref_sensor['longitude'],
    #                        north['latitude'], north['longitude'])[1])

    # initialize
    doa = []
    err = []
    sig = []
    slo = []
    stv = []
    inc = []
    snr = []

    # get matches and corresponding info for each file
    pre = 40*60   # number of samples before the event
    snapshots = 40*240   # number of (time) samples in total
    for arid, time in zip(arrival.ARID, arrival.LDDATE):

        if arid in assoc.ARID.values:
            if str(arrival.loc[arrival['ARID'] == arid, 'JDATE'].values[0])[-3:] == file[11:14]:

                # set to true to visualize the slice of the stream
                if False:
                    cut = stream.slice(UTCDateTime(time) - pre / 40,
                                       UTCDateTime(time) + (snapshots - 1 - pre) / 40)
                    if cut.count():
                        print(str(arrival.loc[arrival['ARID'] == arid, 'JDATE'].values[0])[-3:],
                              file[11:14])
                        print(cut.count())
                        cut.plot()

                # get time and slice stream
                distances = np.zeros((len(stream), 3))
                traces = np.zeros((len(stream), snapshots))
                for i, tr in enumerate(stream):
                    cut = tr.slice(UTCDateTime(time) - pre / 40,
                                   UTCDateTime(time) + (snapshots - 1 - pre) / 40)

                    if len(cut.data):

                        # build sensor map <-> distances (for steering vectors)
                        try: meta = inv.get_channel_metadata(tr.get_id())
                        except Exception: print(tr.get_id())
                        distances[i, :] = [gps2dist_azimuth(ref_sensor['latitude'], ref_sensor['longitude'],
                                                            meta['latitude'], ref_sensor['longitude'])[0],
                                           gps2dist_azimuth(ref_sensor['latitude'], ref_sensor['longitude'],
                                                            ref_sensor['latitude'], meta['longitude'])[0],
                                           (ref_sensor['elevation'] - meta['elevation'])]

                        # get recordings and pad if necessary
                        traces[i] = np.pad(cut.data, (0, snapshots - len(cut.data)), constant_values=0)

                # if stream contains data at timeslot -> save
                if np.any(traces):
                    doa.append(assoc.SEAZ[np.where(assoc.ARID == arid)[0]].values)
                    err.append(arrival.DELAZ[np.where(arrival.ARID == arid)[0]].values)
                    sig.append(traces)
                    slo.append(1/(1000*degrees2kilometers(1/arrival.SLOW[np.where(arrival.ARID == arid)[0]].values)))
                    stv.append(distances)
                    inc.append(arrival.EMA[np.where(assoc.ARID == arid)[0]].values)
                    snr.append(arrival.SNR[np.where(assoc.ARID == arid)[0]].values)

    #print(file, len(doa))

    # save data set in matFiles folder
    matRoot = 'matFiles2'
    try:
        os.mkdir(matRoot)
    except Exception:
        pass

    mdict = {}
    mdict['data'] = sig
    mdict['distances'] = stv
    mdict['doa'] = doa
    mdict['err'] = err
    mdict['slow'] = slo
    mdict['removed'] = removed
    mdict['incidence'] = inc
    mdict['snr'] = snr

    savemat(os.path.join(matRoot,'mdict_' + file[:-7] + '.mat'), mdict)
