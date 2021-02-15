from frospy.splitting.load import loadmodel
from obspy.core import AttribDict

def test_loadmodel1():
    SF = loadmodel(modes=['0T4'], format='RR')
    test_csts = AttribDict({'0': [ 1.1910889878085518 ],
                            '2': [-0.24,  0.54, -0.55, -0.87, -0.11]})
    for mode, coeff in SF.cst.items():
        assert mode == '0T4'
        for deg, c in coeff.items():
            for i, _c in enumerate(c):
                assert _c == test_csts[deg][i]
