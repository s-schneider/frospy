import numpy as np


def convert_AB_to_cst(A, B, s):
    """
    Converter as introduced in:
    Tromp and Zanzerkia (1995): Toroidal Splitting Observations

    converts coefficients A and B of a degree s to splitting function
    coefficients of degree s

    type A, B: list
    type s: list

    returns:

    cst: array with t value on the x-axis and complex c_st on y-axis
         e.g. for s = 2
         cst[0] = [ -2, -1, 0, 1, 2 ]
         cst[1] = [ c_(2)(-2), c_(2)(-1), c_(2)(0), c_(2)(1), c_(2)(2) ]

    mcst: array containing A,B of degree s converted to c_st ordered as
          c_0t
          Re(c_2t)
          Im(c_2t)
          Re(c_4t)
          Im(c_4t)
          ...

    e.g.:

    from frospy.util.read import read_AB
    from frospy.converter.AB_to_cst import convert_AB_to_cst

    AB = read_AB("AB_s2_tromp_zanzerkia.dat")
    A = AB['0T4'].A
    B = AB['0T4'].B
    s = 2 # only returns values for one degree at a time
    cst, mcst = convert_AB_to_cst(A, B, s)

    mcst results should be multiplied by the center frequency
    """

    cst = []
    mcst = []
    for t in range(-s, s+1):
        if t < 0:
            c = np.sqrt(2. * np.pi) * (A[abs(t)] + 1j*B[abs(t)-1])
        elif t == 0:
            c = np.sqrt(4. * np.pi) * (A[t])
            mcst.append(c)
        elif t > 0:
            c = (-1)**t * np.sqrt(2. * np.pi) * (A[t] - 1j*B[t-1])
            mcst.append(c.real)
            mcst.append(c.imag)

        cst.append([t, c])

    cst = np.array(cst, dtype='object').transpose()
    mcst = np.array(mcst)
    return cst, mcst
