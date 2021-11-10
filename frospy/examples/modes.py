# To import the mode catalog and select mode  you can use the followin
# lines of code, which creates a Modes object, containing only the mode :

from frospy.core.modes import read as read_modes
allmodes = read_modes()
modes = allmodes.select(name='0S2')
mode = modes[0]

# printing the overtone number, type and angular degree of the mode
print(mode.n, mode.type, mode.l)

# printing the PREM center frequency and Q value and Q-cycle
print(mode.freq, mode.Q)
print(mode.qcycle)

# It is also possible to create new Mode objects and
# give your own values as an input:

from frospy.core.modes import Mode

mode = Mode(n=0, type='S', l=2, freq=0.309, Q=509)
print(mode.n, mode.type, mode.l)
print(mode.freq, mode.Q)
print(mode.qcycle)
