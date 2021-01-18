4# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np

mypwd='/net/home/talavera/radial-inversion/03s00/09s02'
n = int(open('%s/matrix.dat'%mypwd).readline().rstrip())
matrix=np.loadtxt('%s/matrix.dat'%mypwd, skiprows=1).view(dtype=np.complex128).reshape((n,-1))

rr=np.loadtxt('%s/rr.dat'%mypwd).view(dtype=np.complex128).reshape((n,-1))
rr1=np.loadtxt('%s/rr1.dat'%mypwd).view(dtype=np.complex128).reshape((n,-1))
rr1_test=np.conjugate(rr).T

I = np.dot(rr1,rr)
I_test = np.dot(rr1_test,rr)
trace = np.trace(I_test) - np.trace(I)