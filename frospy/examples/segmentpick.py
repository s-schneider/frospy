# The following example demonstrates how to create Pick object of
# a frequency window for a given Spectrum object.
# An example file is provided in the repository.
# In this example replace the 'PATH/TO/REPO' with the actual path
# to where you cloned/downloaded the repository.

from obspy import read
st = read('/PATH/TO/REPO/frospy/data/examples/031111B.ahx')
print(st)

from frospy import Spectrum
from frospy.spectrum.controllers import get_pick

spec = Spectrum(st[0], 5, 50)
print(spec)

pick = get_pick(spec, [1.025, 1.05], event='031111B', weighting='sum')
print(pick)

# The metadata of the Pick can be accessed vie the stats keyword:
print(pick.stats)

# Multiple Pick objects can be grouped in a Segment object,
# which is demonstrated in the following example:

from frospy import Segment
Seg = Segment()
Seg += pick
print(Seg)
