import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LinearSegmentedColormap
from matplotlib import ticker
from frospy.util.read import read_modes
from frospy.util.read import read_std_cat
from frospy.preprocessing.spectrum import read_data, read_syn, get_times,taper, fourier_transform
from frospy.preprocessing.spectrum import _plot_modes_in_ax, get_mode
from frospy.preprocessing.data_correction import get_baz_synseis
from frospy.preprocessing.spectrum import read_segment_file
from obspy.geodetics import gps2dist_azimuth
#from operator import itemgetter

def misfit(X, Xsyn, s, e):
    d2 = sum(abs(X[s:e+1])**2.)
    rmf = sum((abs(X[s:e+1]) - abs(Xsyn[s:e+1]))**2.) / d2
    cmf = sum(abs(X[s:e+1] - Xsyn[s:e+1])**2.) / d2

    return rmf, cmf # real and complex

event='060994A'
data='PREM/%s.ahx'%event
syn=['PREM/%s.ahx.syn'%event, # self
     'PREM/%s.ahx.prem.syn'%event, # coupled 
    # 'IC/%s.zero.ahx.syn.fil'%event, # group S20RTS
     'S20RTS/%s.ahx.s20rts.syn'%event,
    # 'IC/%s.woodhouse.ahx.syn.fil'%event,
    # 'IC/%s.romanowicz.ahx.syn.fil'%event,
    # 'IC/%s.tromp.ahx.syn.fil'%event,
    # 'IC/%s.begtramp.ahx.syn.fil'%event,
     'Qtest/%s.ahx.syn'%event,
     'HT/%s.ahx.syn'%event]#,
    # 'DE/%s.ahx.syn'%event,
    # 'QM1/%s.ahx.syn'%event] 
syn_label=['PREM-self',
           'PREM-cross', 
         #  'S20RTS-group', 
           'S20RTS',
         #  'IC Woodhouse', 
         #  'IC Romanowicz', 
         #  'IC Tromp', 
         #  'IC BegTramp',
           'Qtest',
           'HT']#, 
         #  'DE', 
         #  'QM1']
segment = read_segment_file('segments/%s.dat'%event)
#stations = [s for s in segment]; stations.sort()
stations = [s for s in segment]#.iterkeys()];
#mode=['${}_1 S_0$']; fw = [1.618, 1.641]; tw_range = range(80,120,5)
#mode=['${}_2 S_0$']; fw = [2.489, 2.534]; tw_range = range(65,110,5)
mode=['${}_3 S_0$']; fw = [3.2550, 3.280]; tw_range = range(95,100,5); tw1=35
#mode=['${}_4 S_0$']; fw = [4.085, 4.125]; tw_range = range(50,95,5); tw1=25
#mode=['${}_5 S_0$']; fw = [4.872, 4.895]; tw_range = range(60,110,5)
#mode=['${}_6 S_0$']; fw = [5.724, 5.759]; tw_range = range(50,110,5)

cat = read_std_cat(cmt_id=event)
st, st_work, syn = read_data(data, syn)
st_syn, st_syn_work = read_syn(syn, st, syn_label) 
color = plt.cm.get_cmap("rainbow", len(st_syn))

tw_shift = []; j = 0   

for t in tw_range:
    tw=[tw1,t]
    taper_shape='hanning'
    Htstart = tw[0] * 3600.
    Htend = tw[1] * 3600.
    Htmid = (tw[1] - tw[0])/2.    
    Htstart_org, Htend_org, Htmid_org = Htstart, Htend, Htmid

    for sta in stations: 
    #for sta in  ['CMO']:
        # startup
        Fxx = []; Fxxsyn = []
        xp = []; xpsyn = []
        freq = [] 
        mf=[[],[]]

        # from segment file
        wstart = segment[sta]['fw1']
        wend = segment[sta]['fw2']
        xscale = segment[sta]['weight']

        # data   
        tr_data = st.select(station=sta)[0]
        print 'before', tr_data
        times = get_times(tr_data, Htstart, Htend, Htstart_org, Htend_org)
        tstart, nstart, nend, t, Htstart, Htend, Htmid, Terr = times
        trfill = taper(nstart, nend, tr_data, taper_shape) # taper
        f, FT, XP, delomeg = fourier_transform(tr_data, Htstart, Htend, nstart, nend, taper_shape, 1) # Fouriertransform and filter
        freq.append(f)
        Fxx.append(FT)
        xp.append(XP)
        
        #setting limits
        startlabel = []; endlabel = [];
        startlabel.append(int(np.round(wstart * 2. * np.pi/(1000. * delomeg))) + 1)
        endlabel.append(int(np.round(wend * 2. * np.pi/(1000. * delomeg))) + 1)
        startlabel=startlabel[0]
        endlabel=endlabel[0]
        
        # syn
        for k in range(len(st_syn)):
            tr_syn = st_syn[k].select(station=sta)[0]
            print 'before', tr_syn
            trsynfill = taper(nstart, nend, tr_syn, taper_shape)
            void, FT, XP, void = fourier_transform(tr_syn, Htstart, Htend, nstart, nend, taper_shape, 1)
            Fxxsyn.append(FT[startlabel:endlabel+1])
            xpsyn.append(XP[startlabel:endlabel+1])
            rmf, cmf = misfit(Fxx[0]/xscale, FT/xscale, startlabel, endlabel)
            mf[0].append((1-rmf)*100); mf[1].append((1-cmf)*100)
        
        print ''
        #station metadata
        s_lat = tr_data.stats.ah.station.latitude
        s_lon = tr_data.stats.ah.station.longitude
        e_lat = cat[0].origins[0].latitude
        e_lon = cat[0].origins[0].longitude
        baz = get_baz_synseis(e_lat, e_lon, s_lat, s_lon)
        d_km = gps2dist_azimuth(s_lat, s_lon, e_lat, e_lon)[0]/1000.
        meta = [s_lat, s_lon, baz, d_km]
        
        tw_shift.append([])
        tw_shift[j].append(sta) # station
        tw_shift[j].append(meta) # lat, lon, back azimuth and epicentral distance
        tw_shift[j].append(tw) # time window
        tw_shift[j].append(mf) # real and complex misfit
        tw_shift[j].append(Fxx[0][startlabel:endlabel+1]) # Data
        tw_shift[j].append(Fxxsyn)
        j = j + 1  
                
#plotting
            
station = 'ADK'
one_sta = [ts for ts in tw_shift if ts[0]==station]
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot(1,1,1)
dlabel = '%s' % one_sta[0][0]
ax.plot(abs(one_sta[0][4]),
        linestyle='solid', color='black', linewidth=0.75, label=dlabel)
for k in range(len(st_syn)):        
    ax.plot(abs(one_sta[0][5][k]),
            linestyle='dashed', color=color(k), linewidth=0.75, label=syn_label[k])
box = ax.get_position() # Shrink current axis's height by 10% on the bottom
ax.set_position([box.x0, box.y0 + box.height * 0.15,
                 box.width, box.height * 1.5])
modedata = read_modes()#; modes = modedata
modes=get_mode(['3s0', '8S2', '9S2'], modedata)
f_lim = [f[startlabel], f[endlabel+1]]
_plot_modes_in_ax(ax, modes, f_lim)
ax.set_ylabel('Amplitude')
ax.set_xlabel('frequency (mHz)')
ax.legend()
mf = ticker.ScalarFormatter(useMathText=True)
mf.set_powerlimits((-2,2)) # for exponents greater than +-2
plt.gca().yaxis.set_major_formatter(mf)