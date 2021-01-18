# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['ps.useafm'] = True
matplotlib.rcParams['pdf.use14corefonts'] = True
matplotlib.rcParams['text.usetex'] = True

fig = plt.figure()
fig.set_size_inches(7.5, 6)
ax = fig.add_subplot(1, 1, 1)


a = 0.4
y = np.linspace(1,50,100)
x=np.exp(-a*y)

# Right plot / left yaxis
ax.plot(x, y, color="r")
ax.fill_between(x, y, 1, color="r")

ax.set_xlabel('Number of worldwide earthquakes per year', fontsize=16)
ax.set_ylabel('Moment\n Magnitude', rotation=0, y=1.02, fontsize=14)
ax.axhspan(ymin=1, ymax=2,  xmin=-1, xmax=1, lw=1, color="lightgray", zorder=0)
ax.axhspan(ymin=3, ymax=4,  xmin=-1, xmax=1, lw=1, color="lightgray", zorder=0)
ax.axhspan(ymin=5, ymax=6,  xmin=-1, xmax=1, lw=1, color="lightgray", zorder=0)
ax.axhspan(ymin=7, ymax=8,  xmin=-1, xmax=1, lw=1, color="lightgray", zorder=0)
ax.axhspan(ymin=9, ymax=10, xmin=-1, xmax=1, lw=1, color="lightgray", zorder=0)
ax.set_xlim([-0.65,0.65])
ax.set_ylim([1.01,10.5])

ax.yaxis.set_major_locator(plt.MaxNLocator(11))
ax.spines['bottom'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.tick_params(axis='x', which='both', bottom=False,labelsize=12)
ax.tick_params(axis='y', which='both',labelsize=12)
ax.set_xticklabels([])
ax.set_title("Notable earthquakes", ha="right", x=0.45, 
             fontsize=16, fontweight="bold")



# Left plot / right yaxis
ax1 = ax.twinx()
ax1.plot(-x, y, color="r")
ax1.fill_between(-x, y, 1, color="r")
ax1.scatter(-np.exp(-a*9.1), 9.1, s=35, color="k", zorder=50) #Tohoku, 2011
ax1.text(-np.exp(-a*9.1)*1.5, 9.1, s="Japan, 2011", ha='right')
ax1.scatter(-np.exp(-a*8.8), 8.8, s=35, color="k", zorder=50) #Chile, 2010
ax1.text(-np.exp(-a*8.8)*1.5, 8.8, s="Chile, 2010", ha='right')
ax1.scatter(-np.exp(-a*8.3), 8.3, s=35, color="k", zorder=50) #Bolivia, 1994
ax1.text(-np.exp(-a*8.3)*1.5, 8.3,  s="Bolivia, 1994", ha='right') 
ax1.scatter(-np.exp(-a*7.6), 7.6, s=35, color="k", zorder=50) #Fiji, 1994
ax1.text(-np.exp(-a*7.6)*1.5, 7.6,  s="Fiji, 1994", ha='right') 
ax1.scatter(-np.exp(-a*6.3), 6.3, s=35, color="k", zorder=50) #Managua, 1972
ax1.text(-np.exp(-a*6.3)*1.2, 6.3,  s="Nicaragua, 1972", ha='right') 
ax1.scatter(np.exp(-a*9), 9, s=35, color="k", zorder=50) #Huracane
ax1.text(np.exp(-a*9)*1.6, 9, s="Category 5 huracane", ha='left')
ax1.scatter(np.exp(-a*8.6), 8.6, s=35, color="k", zorder=50) #Krakatoa, Indonesia
ax1.text(np.exp(-a*8.6)*1.5, 8.6, s="Krakatoa volcanic eruption", ha='left')
ax1.scatter(np.exp(-a*8.05), 8.05, s=35, color="k", zorder=50) #Atomic test
ax1.text(np.exp(-a*8.05)*1.5, 8.05, s="Worlds's largest nuclear test", ha='left')
ax1.scatter(np.exp(-a*7), 7, s=35, color="k", zorder=50) #Tropical cyclo
ax1.text(np.exp(-a*7)*1.3, 7, s="Tropical cyclon", ha='left')
ax1.scatter(np.exp(-a*6.25), 6.25, s=35, color="k", zorder=50) #Hiroshima
ax1.text(np.exp(-a*6.25)*1.3, 6.25, s="Hiroshima atomic bomb", ha='left')
ax1.scatter(np.exp(-a*5), 5, s=35, color="k", zorder=50) #Twin towers
ax1.text(np.exp(-a*5)*1.3, 5, s="9/11 terrorist attack", ha='left')
ax1.scatter(np.exp(-a*4.7), 4.7, s=35, color="k", zorder=50) #Tornado
ax1.text(np.exp(-a*4.7)*1.2, 4.7, s="Average tornado", ha='left')
ax1.scatter(np.exp(-a*3.3), 3.3, s=35, color="k", zorder=50) #Lighting bolt
ax1.text(np.exp(-a*3.3)*1.1, 3.3, s="Large lighting bolt", ha='left')

elabel = ["1,000,000", "100,000", "12,000", "2,000", "200", "20", "3", "$<$1"]
for i,label in zip(np.arange(2,10), elabel):
    if i == 9:
        ax1.text(x=0, y=i, s=label, ha='center', va='center', zorder=20, 
                 fontsize=10)
    else:
        ax1.hlines(y=i, xmin=-np.exp(-a*i), xmax=np.exp(-a*i), 
                   lw=1, color="gray", zorder=10)
        ax1.text(x=0, y=i, s=label, ha='center', va='center', zorder=20, 
                 fontsize=12, bbox=dict(facecolor='r', edgecolor="None"))


yticks = ["", " 56", 
          "1.8$\cdot 10^{3}$",  " 56$\cdot 10^{3}$", 
          "1.8$\cdot 10^{6}$",  " 56$\cdot 10^{6}$",  
          "1.8$\cdot 10^{9}$",  " 56$\cdot 10^{9}$", 
          "1.8$\cdot 10^{12}$", " 56$\cdot 10^{12}$", ]
ax1.set_yticks(np.linspace(1, len(yticks), len(yticks)))
ax1.set_yticklabels(yticks, rotation='horizontal', ha='left', fontsize=12)
ax1.set_ylabel('Energy\n Release', rotation=0, y=1.1, ha="right",fontsize=14)
ax1.text(x=0.68, y=0, s="explosive\nequivalent\nin (kg)", fontsize=12)
ax1.set_ylim([1.01,10.5])

ax1.yaxis.set_major_locator(plt.MaxNLocator(11))
ax1.spines['bottom'].set_visible(False)
ax1.spines['top'].set_visible(False)
ax1.tick_params(axis='x', which='both', bottom=False)
ax1.set_xticklabels([])

ax1.set_title("Energy equivalents", ha="left", x=0.55, 
              fontsize=16, fontweight="bold") #.expandtabs()
plt.savefig('mag_v_freq.png',
              dpi=400,
              bbox_inches='tight', 
              pad_inches=0.0,
              )
