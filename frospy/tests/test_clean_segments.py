from frospy.scripts.segments import clean_segments
from frospy.tests import data as testdata


def test_clean_segments():
    path = "%s/deg8" % (testdata.__path__.__dict__['_path'][0])
    out = clean_segments.main(rundir=path, damp=0.0001, print_seg=False)

    assert len(out) == 71
