from frospy.util.data_request import downloader
from frospy.util.read import read_std_cat

inv = read_std_inv()
#inv = inv.select(station='ANTO')
#inv = None
#cat_eq = read_std_cat('060994A')
cat_eq = read_std_cat('100494B')
length_in_days = 13
components = ["Z"]
sampling_rate = 0.1
localfolder = '/data/talavera/data/300hrs/raw'

downloader(inv, cat_eq, components, sampling_rate,
           localfolder, length_in_days,
           do_process_data=False, verbose=False)
