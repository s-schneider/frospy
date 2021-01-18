from frospy.util.base import sort_human
from frospy.util.base import split_digit_nondigit


def write(SF, format='txt'):
    """
    format = 'dat' -> Arwens format, that is downloadable
    http://www.geo.uu.nl/~deuss/img/cst-coef.dat
    """
    if format == 'dat':
        msg = ""
        for sfunc in SF:
            for mode in sfunc.cst:
                n, name, l = split_digit_nondigit(mode)
                msg += "{}  {}\n".format(n, l)
                degs = sort_human(list(sfunc.cst[mode].keys()))
                for deg in degs:
                    c = sfunc.cst[mode][deg]
                    c_err = sfunc.cst_errors[mode][deg]['uncertainty']

                    if deg == '0':
                        # msg += "{}\n".format(deg)
                        msg += "{:>6.2f} {:>6.2f}\n".format(c[0], sfunc.dst[mode]['0'][0])
                        msg += "{:>6.2f} {:>6.2f}".format(c_err[0], sfunc.dst_errors[mode]['0']['uncertainty'][0])
                    else:
                        cf = ["{:>6.2f}".format(x) for x in c]
                        cf_errors = ["{:>6.2f}".format(x) for x in c_err]
                        # msg += "{}\n".format(deg)
                        msg += ' '.join(cf)
                        msg += '\n'
                        msg += ' '.join(cf_errors)
                    msg += '\n'


        text_file = open("cst-coef.dat", "wt")
        n = text_file.write(msg)
        text_file.close()
    return
