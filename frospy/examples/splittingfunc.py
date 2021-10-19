import matplotlib.pyplot as plt
from frospy import load
SF = load(modes='0T4', format='SAS')

# To access all cst, dst and corresponding errors
print(SF.cst)
print(SF.cst_errors)
print(SF.dst)
print(SF.dst_errors)

# To access a particular mode and degree type
# (Both the mode and degree identifiers are strings)
print(SF.cst["0T4"]['2'])
print(SF.cst_errors["0T4"]['2'])
print(SF.get_fQ(mode_name="0T4"))


# To plot the splitting function coefficients together with their
# corresponding uncertainties per degree
SF.plot()

# To plot the associated splitting function map
SF.plot_map()

# The minimum and maximum degree plotted can be indicated
# for both plots using smin and smax.
SF.plot(smin=2, smax=2)
SF.plot_map(smax=2)

# To write the coefficients of a mode to
# an external file use (the default format is pickle)
SF.write(filename=modes.out, format="dat")


# To import a Set of modes and plot it as a branch
from frospy import load
from frospy.plot.branch import branch

SF = load(format='SAS', modes=['0T4', '0T5', '0T6', '0T7', '0T8'], return_set=True)
# Center frequencies
branch(SF_in=SF, degree=0, plot_f_center=True, mode='t', n=0)
plt.show()

#Q values
branch(SF_in=SF, degree=0, plot_Q=True)
plt.show()

# Degree 2 coefficients
branch(SF_in=SF, degree=2, mode='t', n=0)
plt.show()

# To add more splittingfunctions to Set
SF += load(modes='0T9', format='SAS')
SF.append(load(modes='0T10', format='SAS'))

# Writing the whole Set to file
from frospy.core.splittingfunc.write import write
write(SF)
