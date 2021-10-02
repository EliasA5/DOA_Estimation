

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from functools import reduce
from scipy.io import savemat

client = Client("GFZ")

t0 = UTCDateTime("2014-02-27T06:45:00.000")
t1 = t0 + 4 * 60

stations = [["IS", "MMA*"], ["IS", "MMB*"], ["IS", "MMC*"]]

add_4_mins = lambda x: x + 4 * 60

def query_data(network, station, starttime, endtime, level="response"):
    inv = client.get_stations(network=network, station=station, level=level).get_contents().get('channels')
    bulk = list(map(lambda x: x.split('.'), inv))
    list(map(lambda x: x.extend([starttime, endtime]), bulk))
    return client.get_waveforms_bulk(bulk)

def wrapper(func, arg_list):
    return lambda starttime: func(*arg_list, starttime, add_4_mins(starttime))

choose_time = lambda starttime: reduce(lambda x,y: x+y, map(lambda station: wrapper(query_data, station)(starttime), stations))
stream = choose_time(t0)
#print(stream)

for i, tr in enumerate(stream):
    mdict = {'station': tr.id}
    mdict['data'] = tr.data
    savemat("data-%d.mat" % i, mdict)