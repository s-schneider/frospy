#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 17 17:04:29 2017

@author: karaoglu
"""

import numpy as np
import struct

MATRIX_FILE = "matrix.dat"

###################
# READ ATAB FILES #
###################

with open(MATRIX_FILE,mode='rb') as fp:
    filecontent = fp.read()
    nelem = struct.unpack("i", filecontent[:4])[0]
   
    matrix = np.zeros([nelem,nelem],dtype=np.complex_)
    matrix_re = np.zeros([nelem,nelem],dtype=float)
    matrix_im = np.zeros([nelem,nelem],dtype=float)
#    matrix_list = [[], []]
    icnt = 0
    for i in range(nelem):
        for j in range(nelem):                
            realpart = struct.unpack("d", filecontent[(4+icnt*8):(4+(icnt+1)*8)])[0]
            #print realpart
            
            icnt = icnt + 1
            
            imagpart = struct.unpack("d", filecontent[(4+icnt*8):(4+(icnt+1)*8)])[0]
            icnt = icnt + 1
            
            matrix[j,i] = np.complex(realpart,imagpart)
            matrix_re[j,i] = realpart
            matrix_im[j,i] = imagpart
#            matrix_list[0].append(realpart)
#            matrix_list[1].append(imagpart)

#np.savetxt('matrix.dat.ascii',matrix,fmt='%16.4e',delimiter='\n', 
#           header=str(nelem), comments='')

np.savetxt('matrix.dat.ascii', matrix.view(float).reshape(-1, 2) ,
           fmt='% 20.4e', header=str(nelem), comments='')


