# Free Oscillation Toolbox - FrosPy
Small toolbox to help with visualization of spectra and data-handling for splitting function measurements.

## Content
 * [Wiki](https://github.com/s-schneider/frospy/wiki/Home:-Free-Oscillation-Toolbox---FrosPy)
 * [Toroidal Mode Splittingfunctions](#toroidal-mode-splittingfunctions)
 * [Installation](#installation)
 * [Update](#update)
 * [Usage](#usage)

## Splittingfunction Coefficients

The splitting functions of [Schneider & Deuss, 2020](https://doi.org/10.1093/gji/ggaa567) and [Talavera-Soza & Deuss, 2020](https://doi.org/10.1093/gji/ggaa499) can be found here:
http://www.frospy.org



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
