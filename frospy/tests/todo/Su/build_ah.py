from frospy.util.data_request import process_downloader_data
from frospy.util.read import read_inventory as read_inv
from frospy.util.read import read_std_cat

inv = None #read_inv('/data/talavera/data/300hrs/rawZ/060994A/inv.xml')
#cat = read_std_cat('060994A')
cat = read_std_cat('100494B')
components = ["Z"]
sampling_rate = 0.1
localfolder = '/data/talavera/data/300hrs/raw'
inspected_only = 'both'
rotate_traces = False
rm_response = True
rm_tidal = True
keep_longest_traces = False
cut_to_same_length = False
starttime_check=True
remtidah='/net/home/deuss/bin/remtidah'

process_downloader_data(inv, cat, components, sampling_rate, localfolder,
                        inspected_only,
                        rotate_traces,
                        rm_response,
                        rm_tidal,
                        keep_longest_traces,
                        cut_to_same_length,
			starttime_check,
			remtidah)

