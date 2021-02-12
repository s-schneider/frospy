#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from future.builtins import *  # NOQA @UnusedWildImport

from frospy.core.spectrum.spectrum import Stats


# x = spectrum(data, tw=[2, 10], fw=[1, 2], nowindow=True, station='TAM')


class TestSpectrum:

    def test_stats(self):
        x = Stats({})
        assert len(x.defaults) == 8
        assert len(x.defaults['tw']) == 0
        assert x.defaults['taper'] == 'hanning'
        assert x.defaults['station'] is None
        assert x.defaults['record'] is None
        assert x.defaults['times'] is None
        assert x.defaults['delomeg'] == 0
        assert x.defaults['origin_Tdiff'] == 0

    def test_load_cmt(self):
        return

    def test_load_segments(self):
        return

    def test_search(self):
        return

    def test_fpeaks(self):
        return


class TestUserinput():
    pass
# if __name__ == '__main__':
    # unittest.main(defaultTest='suite')
