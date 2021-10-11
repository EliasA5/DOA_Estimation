

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from functools import reduce
from scipy.io import savemat
import numpy as np

client = Client("GFZ")

t0 = UTCDateTime("2014-09-12T00:00:00.0000Z")

stations = [["IS", "MMA*"], ["IS", "MMB*"], ["IS", "MMC*"]] #[networkID, Station Names]

mins_4 = lambda x: x + 4 * 60
day_1 = lambda x: x + 24 * 60 * 60

def query_data(stations, starttime, endtime, level="response"):
    inv = query_inventories(stations)
    return client.get_waveforms_bulk(make_bulk(inv, starttime, endtime), attach_response= True), inv

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
#stream.remove_response(output="VEL")
matrix = np.empty((len(stream), 9600))

# for i,tr in enumerate(stream):
#     matrix[i] =  tr.data[(tr.times(reftime=t0) >= 0)][0:9600]
#     # remove_response(output = "VEL")

# mdict = {'time': t0}
# mdict['data'] = matrix
# savemat("data.mat", mdict)