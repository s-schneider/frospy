from frospy.postprocessing.plot import plot_all_data
# from frospy.postprocessing.selection import segment_removal
import matplotlib.pyplot as plt
from frospy.tests import data as testdata


def test_plot_and_select():

    sp = None
    comp = 'T'
    rundir = "%s/segments_T" % (testdata.__path__.__dict__['_path'][0])
    data_path = "%s" % (testdata.__path__.__dict__['_path'][0])
    fw = [4.82, 4.84]
    tw = [3, 40]
    threshold = 2.0e-3

    sp = plot_all_data(rundir=rundir, comp=comp, seg_suffix='segment',
                       data_path=data_path, spectra=sp, weight_fw=fw,
                       select_by='peak', threshold=threshold,
                       save_removed_segments=True, tw=tw)
    plt.close()
    # segment_removal(rundir, segments=sp.del_segments, verbose=True,
    #                 seg_suffix='segment', print_seg=False)
    assert len(sp.all_segments) == 4
    assert len(sp.del_segments) == 0
