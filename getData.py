

from obspy.clients.fdsn import Client
from obspy import UTCDateTime

client = Client("GFZ")

t0 = UTCDateTime("2014-02-27T06:45:00.000")
t1 = t0 + 4 * 60

#st = client.get_waveforms("IU", "ANMO", "00", "LHZ", t, t + 60 * 60)

inventoryA = client.get_stations(network="IS", station="MMA*", level="response")

inventoryB = client.get_stations(network="IS", station="MMB*", level="response")

inventoryC = client.get_stations(network="IS", station="MMC*", level="response")

channel_seedsA = inventoryA.get_contents().get('channels')

channel_seedsB = inventoryB.get_contents().get('channels')

channel_seedsC = inventoryC.get_contents().get('channels')

result = list(map(lambda x: x.split('.'), channel_seedsA))
list(map(lambda x: x.extend([t0, t1]), result))
st = client.get_waveforms_bulk(result)