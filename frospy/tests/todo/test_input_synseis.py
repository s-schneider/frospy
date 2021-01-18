"""
Trial for a testing script, to make sure, the input files are correct
for mdcplmrho methods.
If an error occurs the function returns a non zero value to the system
Simon Schneider, 2017
"""
import sys

infile = str(sys.argv[1])
try:
    verbose = str(sys.argv[2])
    if verbose in ['true', 'True', 'verbose', 'Verbose']:
        verbose = True
    else:
        verbose = False
except Exception:
    verbose = False


def test_input_synseis(input, verbose):
    ifile = [line.strip() for line in open(input)]
    retCode = 0
    errlog = []
    # Check for meaning of first integer, whats the link between this number
    # and the modes?
    if int(ifile[0]) > len(ifile):
        retCode = 5
        errlog.append('Wrong number of coupled modes in input file')

    for i, line in enumerate(ifile[1:]):
        if i == len(ifile) - 2 and not line:
            continue
        if len(line) != 8:
            retCode = 5
            m = 'Line %i: Mode-name "%s", must have length of 8' % (i+1, line)
            errlog.append(m)

        elements = line.split()
        if len(elements[0]) != 3:
            retCode = 5
            m = 'Line %i: %s in "%s", must be 3 digits'
            m = m % (i+1, elements[0], line)
            errlog.append()

        if len(elements[2]) != 2:
            retCode = 5
            m = 'Line %i: %s in "%s", must be 2 digits'
            m = m % (i+1, elements[2], line)
            errlog.append(m)

        if elements[1].islower():
            retCode = 5
            m = 'Line %i: %s in "%s", must be upper case'
            m = m % (i+1, elements[1], line)
            errlog.append(m)

    if retCode == 5:
        print('\n Error in file: %s' % input)
        if verbose:
            for err_msg in errlog:
                print(err_msg)
        raise SystemExit(retCode)


if __name__ == '__main__':
    test_input_synseis(infile, verbose)
