#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 21 11:55:14 2019

@author: talavera
"""


from frospy.util.converter import convert_AB_to_cst
import numpy as np

ReA = [-1036]
imA = [-56]

A2 = [1076, -94, -587]
A4 = [773, 330, -141, 252, -190]
AA = [A2, A4]

B2 = [106, -757]
B4 = [-210, -134, 234, -19]
BB = [B2, B4]

w_center = 4042.58e-6
cst = []

void, mcst = convert_AB_to_cst(ReA, [], 0)
cst.extend(mcst * w_center)
        
void, mcst = convert_AB_to_cst(imA, [], 0)
cst.extend(mcst * w_center)
for A,B,deg in zip(AA,BB,[2,4]):
    void, mcst = convert_AB_to_cst(A, B, deg)
    cst.extend(mcst * w_center)

#print cst

# errors

ReA_err = [22]
imA_err = [16]

A2_err = [48, 83, 95]
A4_err = [51, 86, 90, 79, 40]
AA_err = [A2_err, A4_err]

B2_err = [73, 84]
B4_err = [70, 84, 77, 41]
BB_err = [B2_err, B4_err]

w_center = 4042.58e-6
cst_err = []

for coeff, err in zip(ReA, ReA_err):
    err = coeff+err
void, mcst = convert_AB_to_cst([err], [], 0)
cst_err.extend(mcst * w_center)
   
for coeff, err in zip(imA, imA_err):
    err = coeff+err     
void, mcst = convert_AB_to_cst([err], [], 0)
cst_err.extend(mcst * w_center)

i = 2
for A,Ae,B,Be,deg in zip(AA,AA_err,BB,BB_err,[2,4]):
    errA = []; errB = []
    for coeff, err in zip(A, Ae):
        errA.append(coeff+err)
    for coeff, err in zip(B, Be):
        errB.append(coeff+err)
#    print(A,Ae,errA,B,Be,errB, deg)
    void, mcst = convert_AB_to_cst(errA, errB, deg)
    cst_err.extend(mcst * w_center)

#print cst_err

err = []
for c, ce in zip(cst, cst_err):
    err.append(abs(c-ce))
  
for c, ce in zip(cst, err):
    print c, ce
    
#s = 4
#ss = 0
#tt = 0
#for s in range(0,s+1,2):
#    for t in range(0, s+1):
#        if t == 0:
#            tt = tt + t
#            c = np.sqrt(4. * np.pi) * w*(A[ss])
#            mcst.append(c)
#            print 'Re c%s%s = %s; A%s%s = %s '%(s,t,c,
#                                           s,t,A[ss])
#            if s == 0:
#                c = np.sqrt(4. * np.pi) * w*(imA[0])
#                mcst.append(c)
#                print 'Im c%s%s = %s; im(A%s%s) = %s'%(ss,t,c,
#                                                   s,t,imA[0])
#        elif t > 0:
#            c = (-1)**t * np.sqrt(2. * np.pi) *w* (A[ss] - 1j*B[tt-1])
#            mcst.append(c.real)
#            mcst.append(c.imag)
#            
#            print 'Re c%s%s = %s; A%s%s = %s; B%s%s = %s, %s, %s'%(s,t,c.real,
#                                                           s,t,A[ss],
#                                                           s,t,B[tt-1], ss, tt-1)
#            print 'Im c%s%s = %s; A%s%s = %s; B%s%s = %s'%(s,t,c.imag,
#                                                           s,t,A[ss],
#                                                           s,t,B[tt-1])
#        ss += 1
#        