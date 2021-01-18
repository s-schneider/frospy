# Normal Modes Toolbox
Small toolbox to help with spectra and data-preparation for splitting function measurements.
## Content
 * [Features and commands]( #features-and-commands)
 * [Installation](#installation)
 * [Run](#run)
 * [Update](#update)
 * [Troubleshooting](#troubleshooting)
 * [Usage](#usage)
 * [Commands and kwargs](#commands-and-kwargs)
   * [Picking Frequency Windows](#picking-frequency-windows)
   * [Timewindow and tapershape test](#timewindow-and-tapershape-test)


## Features and commands

| Command | Description |
| :- |:-|
|bts    |    Back to the start    |
|d    |    Delete current station    |
|fsep    |    Alter the amount of shown modes    |
|fsep    |    frequency resolution of mode-lines in plot    |
|fw    |    re-enter the frequency window    |
|goto    |    Jump to Station/Channel by index    |
|help    |    display this window    |
|list or l    |     list all stations and indices    |
|list all or la    |    list all station with meta informaton    |
|load segments or ls    |    Load and plot segment file    |
|mfac    |    Multiplication factor for synthetics    |
|modes    |    Display Mode-frequencies    |
|next, nsp    |    go to next station    |
|p all    |    Pick freq-windows, set for all stations    |
|p    |    Pick freq-windows for segment-file    |
|print fpeaks    |    write max freqs to file    |
|print segments    |    write segments to file ($IDdat)    |
|print station    |    Print all stations to 'stationdat'    |
|printw all    |    write curr win for all stations to file    |
|printw    |    write curr win for curr station to file    |
|psp    |    go to previous station    |
|qcycle    |    sets tw and fw for a given mode    |
|quit    |    exit spectrum    |
|reset    |    Reload original files    |
|save    |    Save current stream    |
|set segments or ss |   Set loaded segments as picks    |
|search    |    search for station attribut    |
|taper    |    set shape of taper window in time-domain    |
|taper test    |    testing all available taper shapes    |
|tw    |    re-enter the time-window    |
|tw test    |    testing different timewindows around qcycle    |
|unload segments    |    Unload segment file, remove plotted lines    |

For most of the [commands](#commands) keyword-arguments (kwargs) are allowed


## Installation
I recommend to run this toolbox using anaconda.

[Follow the instructions on their homepage](https://www.anaconda.com/download/)

[How to manage your environments](https://conda.io/docs/user-guide/tasks/manage-environments.html)

Create a new environment:
```
$ conda create -n nmpy python=3.7
$ conda activate nmpy
(nmpy) $
```

Install required packages
```
(nmpy) $ ./install_requirements
```

Install nmpy
```
(nmpy) $ git clone git@git.science.uu.nl:DeepEarth-UU/modes/nmPy.git
(nmpy) $ cd into nmPy
(nmpy) $ chmod 755 INSTALL.sh
(nmpy) $ chmod 755 update.sh
(nmpy) $ ./INSTALL.sh 1
```
The number following INSTALL.sh defines the following (default is 1):
* - 1 installs nmpy with a startup file for ipython
* - 0 installs nmpy without startup file

start ipython

voilá

## Run
To run spectrum enter the following line
```
(nmpy) $ ipython run_spectrum.py
```

## Update

To update the current branch, run  
`./update.sh`

To update other branches, run  
`./update.sh Branch_Name`

## Troubleshooting
If you have issues with the environment, e.g. you get an `ImportError` for 
`import obspy`, make sure all installations are removed outside the environment.

```
(nmpy) $ source deactivate nmpy
$ conda uninstall PACKAGE
$ source activate nmpy
(nmpy) $  conda install PACKAGE
```


If you see a message like
```
On branch development
Your branch is up-to-date with 'origin/development'.

Changes not staged for commit:
  (use "git add <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	modified:   nmpy/util/geodesy.so

no changes added to commit (use "git add" and/or "git commit -a")
```
although it is in `.gitignore` you need to clean up your tracked files in git.  
Use the following commands:
```
git rm -r --cached .
git add .
git commit -m ".gitignore is now working"
git push
```
## Usage
```
from frospy.spectrum.app import spectrum
data = 'PATH-TO-DATA-FILE'
syn = 'PATH-TO-SYNTHETICS-FILE'
spectrum(data, syn)
```

### Commands and kwargs
After a command is entered, some additional information or keywords may be
needed as a keyword argument (kwarg). E.g. after entering `modes`, the name of
the modes is asked:  
Example 1
```
Waiting for user input (type help for all commands)
 -->  modes
Showing modes? (all / none / names (separated by space), e.g. 0S15) -->  0s15
```

The additional keyword can be added after the command as well,
separated by ` - `, having the same outcome as Example 1:  
Example 2
```
Waiting for user input (type help for all commands)
 -->  modes - 0s15
 ```
 The same goes for `tw`, `fw` etc.

### Picking Frequency Windows
 To pick a frequency window for a segment file enter:
 ```
Waiting for user input (type help for all commands)
 -->  p
 ```
 A window will pop up, in which you can pick 2 frequencies by clicking on the
 data. If you already now your frequencies you can enter
 ```
Waiting for user input (type help for all commands)
 -->  p - f1 f2
 ```
 The same goes for `pa`, which sets/calculates the frequency picks, timewindow
 and weighting for all stations, and `pf`, which does it only for the following
 stations.

### Timewindow and tapershape test
* To test a range of timewindows you may enter `tw test`, which will ask for a
mode to calculate the qcycle 'qT' of that it. 2 new plot windows pop up with
timewindows in a range of [1, qT/2] to [1, 3/2*qT] for data and synthetics.

* To test a range of taper shapes you may enter `taper test`, calculating
the spectrum in the current time- and frequency-window for the following
functions ([more documentation here][1]):
  * hanning - window
  * hamming - window
  * blackman - window
  * kaiser - window
  * bartlett - window

[1]:https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.hanning.html
