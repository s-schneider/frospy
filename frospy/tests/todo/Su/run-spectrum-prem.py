# Allows to run spectrum and saves data in ASCII
from frospy.preprocessing.spectrum import spectrum
from frospy.preprocessing.spectrum import taper
from frospy.preprocessing.spectrum import printw
import sys

event = sys.argv[1] # Event name

# self
spectrum(data='/net/home/talavera/radial-inversion/broadband-syn/0-3mHz/%s.ahx'%(event),
syn=['/net/home/talavera/radial-inversion/broadband-syn/0-3mHz/%s.ahx.syn'%(event), #s20rts
'/net/home/talavera/radial-inversion/broadband-syn/3-4.5mHz/%s.ahx.syn'%(event), #s20rts
'/net/home/talavera/radial-inversion/broadband-syn/4.5-5.5mHz/%s.ahx.syn'%(event), #s20rts
'/net/home/talavera/radial-inversion/broadband-syn/5.5-6.5mHz/%s.ahx.syn'%(event), #s20rts
'/net/home/talavera/radial-inversion/broadband-syn/6.5-7.5mHz/%s.ahx.syn'%(event), #s20rts
'/net/home/talavera/radial-inversion/broadband-syn/7.5-8.5mHz/%s.ahx.syn'%(event), #s20rts
'/net/home/talavera/radial-inversion/broadband-syn/8.5-9mHz/%s.ahx.syn'%(event), #s20rts
'/net/home/talavera/radial-inversion/broadband-syn/9-10mHz/%s.ahx.syn'%(event), #s20rts

'/net/home/talavera/radial-inversion/broadband-syn/0-3mHz/prem/%s.ahx.syn'%(event), #prem
'/net/home/talavera/radial-inversion/broadband-syn/3-4.5mHz/prem/%s.ahx.syn'%(event), #prem
'/net/home/talavera/radial-inversion/broadband-syn/4.5-5.5mHz/prem/%s.ahx.syn'%(event), #prem
'/net/home/talavera/radial-inversion/broadband-syn/5.5-6.5mHz/prem/%s.ahx.syn'%(event), #prem
'/net/home/talavera/radial-inversion/broadband-syn/6.5-7.5mHz/prem/%s.ahx.syn'%(event), #prem
'/net/home/talavera/radial-inversion/broadband-syn/7.5-8.5mHz/prem/%s.ahx.syn'%(event), #prem
'/net/home/talavera/radial-inversion/broadband-syn/8.5-9mHz/prem/%s.ahx.syn'%(event), #prem
'/net/home/talavera/radial-inversion/broadband-syn/9-10mHz/prem/%s.ahx.syn'%(event)], #prem

syn_label=['0-3mHz-s20rts', '3-4.5mHz-s20rts', '4.5-5.5mHz-s20rts', '5.5-6.5mHz-s20rts', '6.5-7.5mHz-s20rts', '7.5-8.5mHz-s20rts', '8.5-9mHz-s20rts', '9-10mHz-s20rts',
           '0-3mHz-prem', '3-4.5mHz-prem', '4.5-5.5mHz-prem', '5.5-6.5mHz-prem', '6.5-7.5mHz-prem', '7.5-8.5mHz-prem', '8.5-9mHz-prem', '9-10mHz-prem'],
show_modes=True)
