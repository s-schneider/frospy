import sqlite3
import getpass
# from .postgresql import load_connection_info, create_db, create_table


def write_cst(SF, db_path, model, author=None, lcut='all', verbose=False):

    try:  # See if its iterable otherwise convert to list
        iter(SF)
    except TypeError:
        SF = [SF]

    if author is None:
        author = getpass.getuser()

    values = []
    for s in SF:
        damp = s.stats.damp
        for kind in ['cst', 'dst']:
            sf = getattr(s, kind)
            sf_errors = getattr(s, kind + '_errors')
            for mode, degrees in sf.items():
                for deg, cst in degrees.items():
                    # convert csts to strings
                    cst = [str(x) for x in cst]
                    if len(sf_errors) > 0:
                        if deg == '0':
                            ce = sf_errors[mode][deg]['uncertainty'][0]
                            ce_u = sf_errors[mode][deg]['upper_uncertainty'][0]
                            ce_l = sf_errors[mode][deg]['lower_uncertainty'][0]
                            ce_c = sf_errors[mode][deg]['confidence_level'][0]
                            cm = "%s %s %s %s" % (str(ce), str(ce_u),
                                                  str(ce_l), str(ce_c))
                            cst_errors = [cm]
                        else:
                            ce = sf_errors[mode][deg]['uncertainty']
                            ce_u = sf_errors[mode][deg]['upper_uncertainty']
                            ce_l = sf_errors[mode][deg]['lower_uncertainty']
                            ce_c = sf_errors[mode][deg]['confidence_level'][0]

                            cst_errors = []
                            for i, e in enumerate(ce):
                                try:
                                    cm = "%s %s %s %s" % (str(e),
                                                          str(ce_u[i]),
                                                          str(ce_l[i]),
                                                          str(ce_c))
                                except TypeError:
                                    cm = "%s %s %s %s" % (str(e),
                                                          str(ce_u),
                                                          str(ce_l),
                                                          str(ce_c))
                                cst_errors += [cm]

                            # cst_errors = [str(c) for c in ce]
                        cst_str = []

                        if verbose:
                            print(mode, kind, deg)
                        for c, ce in zip(cst, cst_errors):
                            if verbose:
                                msg = "cst, err, upper, lower, conf:", c, ce
                                print(msg)
                            cst_str.append("%s %s" % (c, ce))
                    else:
                        cst_str = ['%s 0 0 0 0' % c for c in cst]

                    if '-' in mode:
                        mtype = 'CC'
                    elif 'T' in mode.upper():
                        mtype = 'T'
                    elif 'S' in mode.upper():
                        mtype = 'S'

                    v = (mode, kind, deg, ', '.join(cst_str),
                         lcut, damp, author, mtype, 0)
                    values.append(v)




    db = sqlite3.connect(db_path)
    c = db.cursor()

    schema = "CREATE TABLE if not exists '{table}'"
    schema += " (mode text, kind text, deg text, coeffs text, "
    schema += "lcut text, damp real, author text, type text, id int)"
    table = schema.format(table=model)
    c.execute(table)

    # If necessary replace csts in database?
    schema = "CREATE UNIQUE INDEX if not exists 'uidx' ON '{table}'"
    schema += " (mode, kind, author, deg);"
    table = schema.format(table=model)
    c.execute(table)

    exe = "REPLACE INTO '{table}' VALUES (?,?,?,?,?,?,?,?,?)"
    exe = exe.format(table=model)
    c.executemany(exe, values)
    db.commit()
    db.close()
    return


# def write_cst2postgresql(SF, db_path, model, author=None, lcut='all',
#                          verbose=False):
#
#     try:  # See if its iterable otherwise convert to list
#         iter(SF)
#     except TypeError:
#         SF = [SF]
#
#     if author is None:
#         author = getpass.getuser()
#
#     values = []
#     for s in SF:
#         damp = s.stats.damp
#         for kind in ['cst', 'dst']:
#             sf = getattr(s, kind)
#             sf_errors = getattr(s, kind + '_errors')
#             for mode, degrees in sf.items():
#                 for deg, cst in degrees.items():
#                     # convert csts to strings
#                     cst = [str(x) for x in cst]
#                     if len(sf_errors) > 0:
#                         if deg == '0':
#                             ce = sf_errors[mode][deg]['uncertainty'][0]
#                             ce_u = sf_errors[mode][deg]['upper_uncertainty'][0]
#                             ce_l = sf_errors[mode][deg]['lower_uncertainty'][0]
#                             ce_c = sf_errors[mode][deg]['confidence_level'][0]
#                             cm = "%s %s %s %s" % (str(ce), str(ce_u),
#                                                   str(ce_l), str(ce_c))
#                             cst_errors = [cm]
#                         else:
#                             ce = sf_errors[mode][deg]['uncertainty']
#                             ce_u = sf_errors[mode][deg]['upper_uncertainty']
#                             ce_l = sf_errors[mode][deg]['lower_uncertainty']
#                             ce_c = sf_errors[mode][deg]['confidence_level'][0]
#
#                             cst_errors = []
#                             for i, e in enumerate(ce):
#                                 try:
#                                     cm = "%s %s %s %s" % (str(e),
#                                                           str(ce_u[i]),
#                                                           str(ce_l[i]),
#                                                           str(ce_c))
#                                 except TypeError:
#                                     cm = "%s %s %s %s" % (str(e),
#                                                           str(ce_u),
#                                                           str(ce_l),
#                                                           str(ce_c))
#                                 cst_errors += [cm]
#
#                             # cst_errors = [str(c) for c in ce]
#                         cst_str = []
#
#                         if verbose:
#                             print(mode, kind, deg)
#                         for c, ce in zip(cst, cst_errors):
#                             if verbose:
#                                 msg = "cst, err, upper, lower, conf:", c, ce
#                                 print(msg)
#                             cst_str.append("%s %s" % (c, ce))
#                     else:
#                         cst_str = ['%s 0 0 0 0' % c for c in cst]
#
#                     if '-' in mode:
#                         mtype = 'CC'
#                     elif 'T' in mode.upper():
#                         mtype = 'T'
#                     elif 'S' in mode.upper():
#                         mtype = 'S'
#
#                     v = (mode, kind, deg, ', '.join(cst_str),
#                          lcut, damp, author, mtype, 0)
#                     values.append(v)
#
#
#
#
#     db = sqlite3.connect(db_path)
#     c = db.cursor()
#
#     schema = "CREATE TABLE if not exists '{table}'"
#     schema += " (mode text, kind text, deg text, coeffs text, "
#     schema += "lcut text, damp real, author text, type text, id int)"
#     table = schema.format(table=model)
#     c.execute(table)
#
#     # If necessary replace csts in database?
#     schema = "CREATE UNIQUE INDEX if not exists 'uidx' ON '{table}'"
#     schema += " (mode, kind, author, deg);"
#     table = schema.format(table=model)
#     c.execute(table)
#
#     exe = "REPLACE INTO '{table}' VALUES (?,?,?,?,?,?,?,?,?)"
#     exe = exe.format(table=model)
#     c.executemany(exe, values)
#     db.commit()
#     db.close()
#     return


def _write_cst_coeffs(cst, dst, db_path, model, author=None, lcut='all',
                      verbose=False):

    if author is None:
        author = getpass.getuser()

    values = []

    damp = 1
    for val, kind in zip([cst, dst], ['cst', 'dst']):
        for mode, degrees in val.items():
            for deg, _cst in degrees.items():
                # convert csts to strings
                _cst = [str(x) for x in _cst]
                cst_str = ['%s 0 0 0 0' % c for c in _cst]
                v = (mode, kind, deg, ', '.join(cst_str),
                     lcut, damp, author)
                values.append(v)
    db = sqlite3.connect(db_path)
    c = db.cursor()

    schema = "CREATE TABLE if not exists '{table}'"
    schema += " (mode text, kind text, deg text,"
    schema += "coeffs text, lcut text, damp real, author text)"
    table = schema.format(table=model)
    c.execute(table)

    # If necessary replace csts in database?
    schema = "CREATE UNIQUE INDEX if not exists 'uidx' ON '{table}'"
    schema += " (mode, kind, author, deg);"
    table = schema.format(table=model)
    c.execute(table)

    exe = "REPLACE INTO '{table}' VALUES (?,?,?,?,?,?,?)"
    exe = exe.format(table=model)
    c.executemany(exe, values)
    db.commit()
    db.close()
    return


def _write_db(segment, db_path, verbose=False):
    values = {}
    tables = []
    for pick in segment:
        chan = pick.stats.channel
        if chan not in values:
            values[chan] = []
            tables.append(chan)

        if pick.stats.snr is None:
            snr = 'None'
        else:
            snr = pick.stats.snr

        if pick.stats.modes == '':
            modes = ''
        else:
            modes = ','.join(pick.stats.modes.names)

        v = (pick.station, pick.fw1, pick.fw2, pick.tw1, pick.tw2, pick.weight,
             pick.stats.event, modes, snr,
             pick.stats.author, pick.stats.data_origin)

        values[chan].append(v)

    db = sqlite3.connect(db_path)
    c = db.cursor()

    for t in tables:
        if verbose is True:
            print(t)
        schema = "CREATE TABLE if not exists '{table}"
        schema += "' (station text, fw1 real, fw2 real, tw1 real, tw2 real, "
        schema += "weight real, event text, modes text, snr real, "
        schema += "author text, data_origin text)"
        table = schema.format(table=t)

        c.execute(table)
        exe = "INSERT INTO '{table}' VALUES (?,?,?,?,?,?,?,?,?,?,?)"
        exe = exe.format(table=t)
        if verbose is True:
            print(values[t])
        c.executemany(exe, values[t])
    db.commit()
    db.close()
    return
