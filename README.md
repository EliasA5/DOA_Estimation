

Seismic data extracor, includes waveforms, metadata (latitude, longitute, elevation, local_depth, azimuth, dip) and sensor distance.

stations: a list of lists for extracting data from specific stations, each inner list must be of the form [networkID, Station Names].
random_time: set True to extract the data from a random time between start_time, end_time.



To use the script download obspy, geographiclib:
https://github.com/obspy/obspy/wiki#installation
https://geographiclib.sourceforge.io/

our implementation is tested on python 3.8.12



