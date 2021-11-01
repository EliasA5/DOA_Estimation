

from obspy.clients.fdsn import Client
from obspy import UTCDateTime, read_inventory, read
from functools import reduce
from obspy.core.inventory.inventory import read_inventory
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

start_time = datetime(year = 2012, month = 1, day = 1) #.isoformat()
end_time = datetime(year = 2019, month = 12, day = 30)

t0 = UTCDateTime("2019-01-19T00:00:00.000")
#t0 = UTCDateTime(random_date(start_time, end_time).isoformat())


stations = [["IS", "MMA*"], ["IS", "MMB*"], ["IS", "MMC*"]] #[networkID, Station Names]

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
stream.remove_response(output="VEL")
matrix = np.empty((len(stream), 9600))

for i,tr in enumerate(stream):
    if(tr.stats.channel != 'BHZ'):
        tr.filter("highpass", freq=1.0)
    matrix[i] =  tr.data[(tr.times(reftime=t0) >= 0)][0:9600]
    # remove_response(output = "VEL")

mdict = {'time': t0}
mdict['data'] = matrix
savemat("data.mat", mdict)