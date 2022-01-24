from obspy import *
from obspy.geodetics.base import gps2dist_azimuth
import tarfile
import numpy as np
import glob, os
import pandas as pd
from functools import reduce
from scipy.io import savemat

Arrival=pd.read_csv('./GERES.leb.arrival.csv',header=0)
Assoc=pd.read_csv('./GERES.leb.assoc.csv',header=0)
Detection=pd.read_csv('./GERES.detection.csv',header=0)
Origin=pd.read_csv('./GERES.leb.origin.csv',header=0)
inv = read_inventory('GERES.xml')

Arrival_TIME_numpy = Arrival.TIME.to_numpy()
Assoc_TIME_numpy = np.array(list(map(lambda x: UTCDateTime(x).timestamp, Arrival.LDDATE)))
Detection_TIME_numpy = Detection.TIME.to_numpy()
Origin_TIME_numpy = Origin.TIME.to_numpy()

files = [file for file in glob.glob("*.tar.gz")]
delete = False #set to true to remove extracted files automatically
filesToClean = []
ps = []
for file in files:
    #extract tar.gz file into ./tmp/(last 3 numbers of file name)
    with tarfile.open(file, 'r') as tar:
        tarmembers = tar.getmembers()
        p = tarmembers[0].path[-3:]
        ps.append(p)
        for member in tarmembers:
            member.path = member.path[member.path.rfind('/')+1:]
        tar.extractall(path = './tmp/'+p, members=tarmembers)

    
    #get extracted file names, used to read the traces, and cleanup at the end of the script
    filesToClean.extend([file for file in glob.glob("./tmp/"+p+"/IM*")])
    stream = read("./tmp/"+p+"/IM*", format="mseed")
    print(stream.print_gaps())
    #some sensors have alot of gaps, we can fix this by using stream.merge()
    #https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.merge.html#obspy.core.stream.Stream.merge
    
    stream.merge(method=0, fill_value='latest')

    distances = np.empty([len(stream), 3])
    ref_sensor = inv.get_channel_metadata(stream[0].get_id())
    idx = []
    for i,tr in enumerate(stream):
        try:
            meta = inv.get_channel_metadata(tr.get_id())
            distances[i,:] = [gps2dist_azimuth(ref_sensor['latitude'], ref_sensor['longitude'], meta['latitude'], ref_sensor['longitude'])[0], gps2dist_azimuth(ref_sensor['latitude'], ref_sensor['longitude'], ref_sensor['latitude'], meta['longitude'])[0], (ref_sensor['elevation']-meta['elevation'])]
        except Exception:
            #some sensors have no metadata i.e we can't know thier distance
            idx.append(i)

    #we delete the sensors with no metadata
    distances = np.delete(distances, idx, 0)
    for i in idx:
        stream.remove(stream[i])


    indexes = []
    #another threshhold for the lines below which is the length of the signal in seconds
    thresh = len(stream[0].times())/40
    #start time of the first trace time signal
    reference_start_time = stream[0].times("timestamp")[0]
    for i in range(0, len(Arrival_TIME_numpy)):
        # #threshhold for difference at the start time of the trace, if less than 1 hour (3600) seconds then we have a match
        # if(abs(Arrival_TIME_numpy[i] - reference_start_time) < 3600): 
        #     indexes.append(i)
        #if an arrival happened sometime after the signal started, we got a match
        if(Arrival_TIME_numpy[i] > reference_start_time and Arrival_TIME_numpy[i] - reference_start_time < thresh):
            indexes.append(i)
    Arrival_Matches = Arrival.loc[indexes, :]
    Assoc_Matches = Assoc.loc[indexes, :] #assoc and arrival have the same indexes
    Origin_Matches = Origin.loc[Origin.ORID.isin(Assoc_Matches.ORID), :] #getter for origin

    #whatever you want to do with the data from here

    #you can export data to csv file to load into matlab by using:
    matchRoot = 'Matches'
    try:
        os.mkdir(matchRoot)
    except Exception:
        pass
    Arrival_Matches.to_csv(os.path.join(matchRoot, 'Arrival_Matches_' + p + '.csv'), index=False)
    Assoc_Matches.to_csv(os.path.join(matchRoot,'Assoc_Matches_' + p + '.csv'), index=False)
    Origin_Matches.to_csv(os.path.join(matchRoot,'Origin_Matches_' + p + '.csv'), index=False)
    #you can remove the response of the sensors from the data by using:
    stream.remove_response(inventory=inv,output="VEL")
    #note that each .mat file is around 500 mb
    matRoot = 'matFiles'
    try:
        os.mkdir(matRoot)
    except Exception:
        pass
    # #save the distance vectors
    # distDict = {}
    # distDict['distances'] = distances
    # #or
    # #distDict['distances_' + p] = distances
    # savemat(os.path.join(matRoot,'distances_' + p + '.mat'), distDict)
     
    # #save the trace waveforms (data)
    # traceDataDict = {}
    # traceDataDict['data'] = np.array(list(map(lambda x: x.data, stream)), dtype=object)
    # savemat(os.path.join(matRoot,'traceData_' + p + '.mat'), traceDataDict)
    #can also be saved as one dictionary
    mdict = {}
    mdict['distances'] = distances
    mdict['data'] = np.array(list(map(lambda x: x.data, stream)), dtype=object)
    savemat(os.path.join(matRoot,'mdict_' + p + '.mat'), mdict)

    #until here
    
#cleanup files, set to false by default
if(delete):
    for file in filesToClean:
        os.remove(file)
    for p in ps:
        os.rmdir('./tmp/'+p)
    os.rmdir('./tmp')