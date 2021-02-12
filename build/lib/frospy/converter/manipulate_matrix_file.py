# -*- coding: utf-8 -*-
"""
Created on Sat May 19 15:50:29 2018

Read the matrix file in one format (binary or ascii format) and print it
in the other format.

@author: karaoglu
"""

from abc import ABCMeta, abstractmethod
import sys
import struct
import numpy as np


class Matrix(object):
    """docstring for Matrix."""
    __metaclass__ = ABCMeta

    def __init__(self, matrixfile):
        self.matrixfile = matrixfile
        self.matrixrank = None
        self.matrix = None
        self.filetype = self.setfiletype()

    @abstractmethod
    def setfiletype(self):
        pass

    @abstractmethod
    def loadmatrix(self):
        pass

    @abstractmethod
    def writematrix(self):
        pass


class BinMatrix(Matrix):

    def setfiletype(self):
        return 'binary'

    def loadmatrix(self):
        fp = open(self.matrixfile, 'rb')

        filecontent = fp.read()
        fp.close()

        self.matrixrank = struct.unpack("i", filecontent[:4])[0]
        self.matrix = np.zeros([self.matrixrank, self.matrixrank],
                               dtype=np.complex_)
        icnt = 0
        for i in range(self.matrixrank):
            for j in range(self.matrixrank):
                cont = filecontent[(4+icnt*8):(4+(icnt+1)*8)]
                realpart = struct.unpack("d", cont)[0]
                icnt = icnt + 1
                imagpart = struct.unpack("d", cont)[0]
                icnt = icnt + 1
                self.matrix[j, i] = np.complex(realpart, imagpart)

    def writematrix(self):
        fp = open('ascii_matrix.dat', 'w')
        fp.write("%d\n" % self.matrixrank)

        for i in range(self.matrixrank):
            for j in range(self.matrixrank):
                Re = np.real(self.matrix[j, i])
                Im = np.imag(self.matrix[j, i])
                fp.write("%16.4e %16.4e\n" % (Re, Im))

        fp.close()


class AsciiMatrix(Matrix):

    def setfiletype(self):
        return 'ascii'

    def loadmatrix(self):
        fp = open(self.matrixfile, 'r')
        self.matrixrank = int(fp.readline())
        filecontent = np.loadtxt(fp)
        fp.close()

        self.matrix = np.zeros([self.matrixrank, self.matrixrank],
                               dtype=np.complex_)
        mat = (filecontent[:, 0] + np.complex(0, 1)*filecontent[:, 1])
        self.matrix = mat.reshape([self.matrixrank, self.matrixrank])

    def writematrix(self):
        fp = open('binary_matrix.dat', 'wb')
        fp.write(struct.pack('i', self.matrixrank))

        for i in range(self.matrixrank):
            for j in range(self.matrixrank):
                fp.write(struct.pack('d', np.real(self.matrix[i, j])))
                fp.write(struct.pack('d', np.imag(self.matrix[i, j])))

        fp.close()


def main():

    if(len(sys.argv) < 2):
        raise RuntimeError("Provide the matrix file as an argument")

    try:
        mymatrix = BinMatrix(sys.argv[1])
        mymatrix.loadmatrix()
        mymatrix.writematrix()
        print("Input: Binary matrix file\nOutput: Ascii file")
    except Exception:
        mymatrix = AsciiMatrix(sys.argv[1])
        mymatrix.loadmatrix()
        mymatrix.writematrix()
        print("Input: Ascii file\nOutput: Binary file")


if __name__ == "__main__":
    main()
