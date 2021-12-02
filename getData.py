

from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read_inventory, read
from functools import reduce
from obspy.geodetics.base import gps2dist_azimuth
from scipy.io import savemat
import numpy as np
from random import randrange
from datetime import timedelta, datetime

def random_date(start, end):
    delta = end - start
    int_delta = (delta.days * 24 * 60 * 60) + delta.seconds
    random_second = randrange(int_delta)
    return start + timedelta(seconds=random_second)

client = Client("GFZ")

random_time = False #change to True to get a random time, picks a time between start_time and end_time
#write stations in the following format [networkID, Station Names]
stations = [["IS", "MMA*"], ["IS", "MMB*"], ["IS", "MMC*"]] 

start_time = datetime(year = 2012, month = 1, day = 1) 
end_time = datetime(year = 2019, month = 12, day = 30)

if(random_time):
    t0 = UTCDateTime(random_date(start_time, end_time).isoformat())
else:
    t0 = UTCDateTime("2019-01-19T00:00:00.000")

mins_4 = lambda x: x + 4 * 60
day_1 = lambda x: x + 24 * 60 * 60

def query_data(stations, starttime, endtime, level="response", attach_response="True"):
    inv = query_inventories(stations)
    return query_data_with_inv(inv, starttime, endtime, level, attach_response), inv

def query_data_with_inv(inv, starttime, endtime, level="response", attach_response="True"):
    return client.get_waveforms_bulk(make_bulk(inv, starttime, endtime), level= level, attach_response= attach_response)

def make_bulk(inv, starttime, endtime):
    channels = inv.get_contents().get('channels')
    bulk = list(map(lambda x: x.split('.'), channels))
    list(map(lambda x: x.extend([starttime, endtime]), bulk))
    return bulk

def query_inventory(network, station, level = "response"):
    return client.get_stations(network=network, station=station, level=level)

def query_inventories(stations, level = "response"):
    return reduce(lambda x,y: x+y, map(lambda station: query_inventory(*station, level), stations))

choose_time = lambda starttime, endtime: query_data(stations, starttime, endtime)
stream, inv = choose_time(t0, mins_4(t0))
temp = stream.pop(0)
stream.insert(3, temp)
temp = None

stream.remove_response(output="VEL")
matrix = np.empty((len(stream), 9600))

channel_metadata = {'format' : ["latitude", "longitude", "elevation", "local_depth", "azimuth", "dip"]}
distances = np.empty([len(stream), 3])
ref_sensor = inv.get_channel_metadata(stream[0].get_id())

for i,tr in enumerate(stream):
    if(tr.stats.channel != 'BHZ'):
        tr.filter("highpass", freq=1.0)
    meta = inv.get_channel_metadata(tr.get_id())
    distances[i,:] = [gps2dist_azimuth(ref_sensor['latitude'], ref_sensor['longitude'], meta['latitude'], ref_sensor['longitude'])[0], gps2dist_azimuth(ref_sensor['latitude'], ref_sensor['longitude'], ref_sensor['latitude'], meta['longitude'])[0], (ref_sensor['elevation']-meta['elevation'])]
    channel_metadata[tr.stats.station] = list(meta.values())
    matrix[i] = tr.data[(tr.times(reftime=t0) >= 0)][0:9600]
    # remove_response(output = "VEL")

mdict = {}
mdict['data'] = matrix
mdict['r_m'] = distances
mdict['meta'] = channel_metadata
savemat("data.mat", mdict)