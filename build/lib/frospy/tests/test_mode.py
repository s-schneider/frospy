from frospy.core.modes import Mode, Modes, format_name

def test_mode():
    m = Mode(n=0, type='S', l=2, freq=0.309, Q=509.684)

    assert m.n == 0
    assert m.type == 'S'
    assert m.l == 2
    assert m.freq == 0.309
    assert m.Q == 509.684

def test_modes():
    M = Modes([Mode(), Mode(), Mode()])
    assert len(M) == 3

    m1 = Modes([Mode(), Mode(), Mode()])
    assert len(m1) == 3

    m2 = Modes([Mode(), Mode()])
    assert len(m2) == 2

    M = m1 + m2
    assert len(M) == 5

def test_mode_name_formatting():
    s = format_name('0s2', 2)
    assert s == '00s02'
