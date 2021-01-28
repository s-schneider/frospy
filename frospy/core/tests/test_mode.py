from frospy.core.modes import Mode
from obspy.core import AttribDict

def test_mode():
    m = Mode(n=0, type='S', l=2, freq=0.309, Q=509.684)

    assert m.n == 0
    assert m.type == 'S'
    assert m.l == 2
    assert m.freq == 0.309
    assert m.Q == 509.684

def test_modes():
    assert 0 == 0
