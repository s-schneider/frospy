from frospy.util.base import sort_human
from frospy.util.base import split_digit_nondigit


def write(SF, filename='cst-coef', format='dat', coupling='self',
          verbose=False):
    """
    format = 'dat' -> Arwens format, that is downloadable
    http://www.geo.uu.nl/~deuss/img/cst-coef.dat
    """
    if format == 'dat':
        msg = create_cst_content(SF, coupling, verbose)
        text_file = open("{}.{}".format(filename, format), "wt")
        _n = text_file.write(msg)
        text_file.close()
    return


def create_cst_content(SF, coupling, verbose=False):
    if coupling.lower() in ['self', 'sc']:
        cmod = 'sc'
        msg = ("Self-coupled splitting function coefficients in micro Hz. "
               "The first line contains the n, name, l of\nthe mode."
               "Then there are two lines for each angular order s of the "
               "splitting function, starting with s=0.\n"
               "The first of these contains the c_st coefficients and the "
               "second line the corresponding error values.\n"
               "The coefficients are ordered "
               "c_s0, Re(c_s1), Im(c_s1), Re(c_s2), Im(c_s2), ... etc.\n")

    if coupling.lower() in ['cross', 'cc', 'gc']:
        cmod = 'cc'
        msg = ("Cross-coupled splitting function coefficients in micro Hz. "
               "The first line contains the n1, type, l1 of "
               "the first mode and the n2, type, l2 of the second mode. "
               "The second line contains the minimum and "
               "maximum angular order s of the splitting function "
               "coefficients.\nThen there are two lines for each angular "
               "order of the splitting function. The first of these"
               "contains the c_st coefficients and the second line the "
               "corresponding error values.\n The coefficients are ordered "
               "c_s0, Re(c_s1), Im(c_s1), Re(c_s2), Im(c_s2), ... etc.\n")
    if coupling.lower() == 'all':
        cmod = 'all'
        msg = ("Splitting function coefficients in micro Hz. "
               "The first line contains the n1, type, l1 of "
               "the first mode and, if applicable, the n2, type, l2 of the "
               "second mode. "
               "The second line contains the minimum and "
               "maximum angular order s of the splitting function "
               "coefficients.\nThen there are two lines for each angular "
               "order of the splitting function. The first of these"
               "contains the c_st coefficients and the second line the "
               "corresponding error values.\n The coefficients are ordered "
               "c_s0, Re(c_s1), Im(c_s1), Re(c_s2), Im(c_s2), ... etc.\n")
    for sfunc in SF:
        for mode in sfunc.cst:
            if verbose:
                print(mode)
            if cmod == 'sc':
                if '-' in mode:
                    continue
                n, name, l = split_digit_nondigit(mode)
                msg += "{} {} {}\n".format(n, name, l)
            elif cmod == 'cc':
                if '-' not in mode:
                    continue
                ccdegs = list(sfunc.cst[mode].keys())
                sdn = split_digit_nondigit(mode)
                n1, name1, l1, _, n2, name2, l2 = sdn[:]
                msg += "{} {} {}  {} {} {}\n".format(n1, name1, l1, n2,
                                                     name2, l2)
                msg += "{} {}\n".format(ccdegs[0], ccdegs[-1])
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
    return msg
