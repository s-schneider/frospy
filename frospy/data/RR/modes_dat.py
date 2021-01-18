import nmpy
from frospy.core.modes import read as read_modes
from frospy.core.modes import Modes
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import glob
import os

folder = [
    '/quanta1/home/simons/splitting/modes/00t05',
    '/quanta1/home/simons/splitting/modes/00t06',
    '/quanta1/home/simons/splitting/modes/00t07',
    '/quanta1/home/simons/splitting/modes/00t08',
    '/quanta1/home/simons/splitting/modes/00t09',
    '/quanta1/home/simons/splitting/modes/00t22',
    '/quanta1/home/simons/splitting/modes/01t01',
    '/quanta1/home/simons/splitting/modes/01t02',
    '/quanta1/home/simons/splitting/modes/01t03',
    '/quanta1/home/simons/splitting/modes/01t04',
    '/quanta1/home/simons/splitting/modes/01t05',
    '/quanta1/home/simons/splitting/modes/01t06',
    '/quanta1/home/simons/splitting/modes/01t07',
    '/quanta1/home/simons/splitting/modes/01t08',
    '/quanta1/home/simons/splitting/modes/01t09',
    '/quanta1/home/simons/splitting/modes/02t08'
]

r = ['0 t  5',
     '0 t  6',
     '0 t  7',
     '0 t  8',
     '0 t  9',
     '0 t 10',
     '0 t 11',
     '0 t 12',
     '0 t 13',
     '0 t 14',
     '0 t 15',
     '0 t 16',
     '0 t 17',
     '0 t 18',
     '0 t 19',
     '0 t 20',
     '0 t 21',
     '0 t 22',
     '1 t  1',
     '1 t  2',
     '1 t  3',
     '1 t  4',
     '1 t  5',
     '1 t  6',
     '1 t  7',
     '1 t  8',
     '1 t  9',
     '2 t  2',
     '2 t  4',
     '2 t  8']

r_self = [
     '0 t  5',
     '0 t  6',
     '0 t  7',
     '0 t  8',
     '0 t  9',
     '0 t 22',
     '1 t  1',
     '1 t  2',
     '1 t  3',
     '1 t  4',
     '1 t  5',
     '1 t  6',
     '1 t  7',
     '1 t  8',
     '1 t  9',
     '2 t  8']

m = {
 '00t05': 100,
 '00t06': 100,
 '00t07': 100,
 '00t08': 100,
 '00t09': 100,
 '00t22': 100,
 '01t01': 100,
 '01t02': 100,
 '01t03': 100,
 '01t04': 100,
 '01t05': 100,
 '01t06': 100,
 '01t07': 100,
 '01t08': 100,
 '01t09': 100,
 '02t08': 100}

m_cc = {'00t11-00t10': 1, '00t12-00t11': 1}

modes_dir = '/tmp/eejit_simons/splitting/modes'
modes_dir = '/quanta1/home/simons/splitting/modes/'

m_simon = []
my_modes = glob.glob("%s/?????/segments_T/" % modes_dir)
for f in my_modes:

    mode = f.split('self_coupling')[0].split('/')[-3]
    if len(glob.glob(f)):
        print(mode)
        m_simon += [mode]


allm = read_modes()
M = Modes()
Msimon = Modes()
r = r_self
r = [''.join(x.split()) for x in r]
for x in r:
    M += allm.select(name=x)

for x in m_simon:
    Msimon += allm.select(name=x)

allmt = allm.select(mtype='T')
allmt.sort(['l'])

n0 = 0
x = []
x_highlight = []
y = []
y_highlight = []
l_print = False

x_simon = []
y_simon = []

fig, ax = plt.subplots()
fig.set_size_inches(7, 12)

for m in allmt:
    if m.l > 40:
        continue
    if m.n == n0:
        if l_print is False:
            ax.text(m.l-2, m.freq, m.n, fontsize=12)
            l_print = True
        y.append(m.freq)
        x.append(m.l)
        if m in M:
            x_highlight.append(m.l)
            y_highlight.append(m.freq)
        if m in Msimon:
            x_simon.append(m.l)
            y_simon.append(m.freq)
    else:
        n0 = m.n
        ax.plot(x, y, 'k-', linewidth=0.9)
        ax.scatter(x, y, s=8, c='k')
        x = []
        y = []
        l_print = False

ax.plot(x, y, 'k-', linewidth=0.9)
ax.scatter(x, y, s=8, c='k')

ax.scatter(x_highlight, y_highlight, s=80, facecolors='none', edgecolors='k',
           marker='s')

ax.scatter(x_simon, y_simon, s=80, facecolors='none', edgecolors='r',
           marker='s')


ax.set_xlabel('angular order l', fontsize=22)
ax.set_ylabel('frequency (mHz)', fontsize=22)
ml = MultipleLocator(1)
ax.yaxis.set_minor_locator(ml)
ax.xaxis.set_minor_locator(ml)

# plt.tick_params(axis='both', which='minor', labelsize=15)
plt.show()
