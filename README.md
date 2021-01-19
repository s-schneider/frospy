# Free Oscillation Toolbox - FrosPy
Small toolbox to help with visualization of spectra and data-handling for splitting function measurements.
## Content
 * [Toroidal Mode Splittingfunctions](#toroidal-mode-splittingfunctions)
 * [Core Classes]( #core-classes)
 * [Installation](#installation)
 * [Update](#update)
 * [Usage](#usage)
 * [Spectrum](#spectrum)
 * [Run](#run)
 * [Commands and kwargs](#commands-and-kwargs)
 * [Picking Frequency Windows](#picking-frequency-windows)
 * [Timewindow and tapershape test](#timewindow-and-tapershape-test)

## Toroidal Mode Splittingfunctions

The splittingfunctions of our [toroidal mode overtone paper](https://doi.org/10.1093/gji/ggaa567) can be found here:
  [Toroidal mode overtones](https://github.com/s-schneider/frospy/tree/main/frospy/data/SAS/cst-coef-T.dat)



## Core Classes

* Modes / Mode
* Spectrum
* Pick / Segment
* Splittingfunc / Set
* Setup

## Installation
I recommend to run this toolbox using anaconda.

[Follow the instructions on their homepage](https://www.anaconda.com/download/)

[How to manage your environments](https://conda.io/docs/user-guide/tasks/manage-environments.html)

Create a new environment:
```
$ conda create -n frospy_env python=3.7
$ conda activate frospy_env
(frospy_env) $
```

Install required packages
```
(frospy_env) $ ./install_requirements
```

Install nmpy
clone this repo
```
(frospy_env) $ cd into frospy
(frospy_env) $ chmod 755 install.sh
(frospy_env) $ chmod 755 update.sh
(frospy_env) $ ./install.sh
```

start ipython

voilá

## Update

To update the current branch, run  
`./update.sh`

To update other branches, run  
`./update.sh Branch_Name`

## Usage
```
from frospy.spectrum.app import spectrum
data = 'PATH-TO-DATA-FILE'
syn = 'PATH-TO-SYNTHETICS-FILE'
spectrum(data, syn)
```

## Spectrum
Spectrum is an example tool that uses frospy's core classes to plot normal mode frequency spectra.

## Run
To run spectrum enter the following line
```
(nmpy) $ ipython run_spectrum.py
```

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
