# Allows to run spectrum and saves data in ASCII
from frospy.preprocessing.spectrum import spectrum
from frospy.preprocessing.spectrum import taper
from frospy.preprocessing.spectrum import printw
import sys

event = sys.argv[1] # Event name

# self
spectrum(data='/net/home/talavera/radial-inversion/broadband-syn/0-7.5mHz/%s.ahx'%(event),
syn=['/net/home/talavera/radial-inversion/broadband-syn/0-7.5mHz/%s.ahx.syn'%(event),
'/net/home/talavera/radial-inversion/broadband-syn/7.5-8.5mHz/%s.ahx.syn'%(event),
'/net/home/talavera/radial-inversion/broadband-syn/8.5-9mHz/%s.ahx.syn'%(event),
'/net/home/talavera/radial-inversion/broadband-syn/9-10mHz/%s.ahx.syn'%(event),
'/net/home/talavera/radial-inversion/01s00/01s00/%s/%s.ahx.syn.fil' %(event,event),
'/net/home/talavera/radial-inversion/02s00/02s00/%s/%s.ahx.syn.fil' %(event,event),
'/net/home/talavera/radial-inversion/03s00/03s00/%s/%s.ahx.syn.fil' %(event,event),
'/net/home/talavera/radial-inversion/04s00/04s00/%s/%s.ahx.syn.fil' %(event,event),
'/net/home/talavera/radial-inversion/05s00/05s00/%s/%s.ahx.syn.fil' %(event,event),
'/net/home/talavera/radial-inversion/06s00/06s00/%s/%s.ahx.syn.fil' %(event,event),
'/net/home/talavera/radial-inversion/07s00/07s00/%s/%s.ahx.syn.fil' %(event,event),
'/net/home/talavera/radial-inversion/08s00/08s00/%s/%s.ahx.syn.fil' %(event,event),
'/net/home/talavera/radial-inversion/09s00/09s00/%s/%s.ahx.syn.fil' %(event,event),
'/net/home/talavera/radial-inversion/10s00/10s00/%s/%s.ahx.syn.fil' %(event,event)],
syn_label=['0-7.5mhz', '7.5-8.5mHz', '8.5-9mHz','9-10mHz', '01s00', '02s00', '03s00','04s00', '05s00', '06s00', '07s00', '08s00', '09s00', '10s00'],
tw=[10,60], fw=[0.5, 10],
show_modes=['1S0', '4S2', '5S2', 
            '2S0', '7S2', '6S2', 
            '3S0', '9S2', '8S2', '5S7', 
            '4S0', '10S2', '11S2', 
            '5S0', '13S2', '9S7', 
            '6S0', '16S2', '15S2', 
            '7S0', '18S2', 
            '8S0', '20S2', '21S2',
            '9S0', '22S2',
            '10S0','25S2'])
