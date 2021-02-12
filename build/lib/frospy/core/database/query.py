from __future__ import print_function
import sqlite3
from frospy.util.base import is_number
import sys


def write_misfit2db(db_path, table, event, station, tw1, tw2, fw1, fw2, mf,
                    damp, iter_no):

    value = (event, station, tw1, tw2, fw1, fw2, mf, damp, iter_no)
    db = sqlite3.connect(db_path)
    c = db.cursor()

    create = "CREATE TABLE if not exists "+table
    create += " (event text, station text, tw1 real, tw2 real,"
    create += "fw1 real, fw2 real, misfit real, damp real, iter int)"

    c.execute(create)
    exe = "INSERT INTO "+table+" VALUES (?,?,?,?,?,?,?,?,?)"
    c.execute(exe, value)
    db.commit()
    db.close()
    return


def db_query(db_path, table, key=None, value=None, select='*',
             condition_dict=None, condition_dict_or=None):
    """
    SQL: SELECT select FROM table WHERE key=val

    db_path : path to .db file
    table : name of the table
    key : key value, e.g. station, event or misfit
    value : value of that key, e.g. key event, value 060994A
    """

    db = sqlite3.connect(db_path)
    c = db.cursor()
    val = (value, )
    if key is not None and value is not None:
        select_msg = "SELECT %s FROM '%s' WHERE %s=?" % (select, table, key)
        c.execute(select_msg, val)

    elif condition_dict is not None:
        query = "SELECT %s FROM '%s' WHERE " % (select, table)
        val = ()
        for key, value in condition_dict.items():
            query += key+"=? AND "
            val += (value, )
        query = query[:-4]
        c.execute(query, val)

    elif condition_dict_or is not None:
        query = "SELECT %s FROM '%s' WHERE " % (select, table)
        val = ()
        for key, value in condition_dict_or.items():
            for v in value:
                query += key+"=? OR "
                val += (v, )
        query = query[:-4]
        print(query, val)
        c.execute(query, val)

    else:
        ex_msg = "SELECT %s FROM '%s';" % (select, table)
        c.execute(ex_msg)
    fetch = c.fetchall()
    db.close()
    return fetch


def sql_query(db_path, query):
    db = sqlite3.connect(db_path)
    c = db.cursor()
    c.execute(query)
    fetch = c.fetchall()
    db.close()
    return fetch


def get_db_tables(db_path, pragmas=False):
    db = sqlite3.connect(db_path)
    c = db.cursor()
    c.execute("SELECT name FROM sqlite_master WHERE type='table';")
    table = c.fetchall()
    if pragmas is True:
        table = dict([[x[0], ''] for x in table])
        for fetch in table.keys():
            c.execute('PRAGMA table_info("%s");' % fetch)
            p = []
            for row in c.fetchall():
                p.append(row)
            table[fetch] = p
    else:
        table = [x[0] for x in table]
    db.close()
    return table


def segments_query(db_path, channel, modes, min_snr=None, author=None,
                   event=None, station=None):

    db = sqlite3.connect(db_path)
    c = db.cursor()
    val = ()

    query = "SELECT * FROM '%s' WHERE " % (channel)
    if type(modes) == list:
        query += '('
        for m in modes:
            query += "modes=? OR "
            val += (m, )
        query = query[:-4] + ") "
    else:
        query += "modes=?"
        val += (modes, )

    if min_snr is not None:
        query += " AND snr >=?"
        val += (min_snr, )

    if author is not None:
        query += " AND author =?"
        val += (author, )

    if event is not None:
        if type(event) is list:
            query += '('
            for e in event:
                query += "events=? OR "
                val += (e, )
            query = query[:-4] + ") "
        else:
            query += " AND event=?"
            val += (event, )

    if station is not None:
        if type(station) is list:
            query += '('
            for s in station:
                query += "station=? OR "
                val += (e, )
            query = query[:-4] + ") "
        else:
            query += " AND station =?"
            val += (station, )

    query += ";"

    c.execute(query, val)
    fetch = c.fetchall()
    db.close()
    return fetch


def cst_query(db_path, model, modes):
    if type(modes) is not list:
        modes = [modes]

    db = sqlite3.connect(db_path)
    c = db.cursor()
    val = ()

    query = "SELECT * FROM '%s' WHERE " % (model)

    for mode in modes:
        query += " mode=? OR "
        val += (mode, )

    query = query[:-4] + ";"

    c.execute(query, val)
    fetch = c.fetchall()
    db.close()

    return fetch


def update_tables(db_path, verbose=False):
    """
    Updates inversion.db files, that don't have unique indices
    to ones with them
    """
    try:
        db = sqlite3.connect(db_path)
        c = db.cursor()
        tables = get_db_tables(db_path, pragmas=True)
        if verbose is True:
            print("Tables in database %s\n" % tables.keys())

        ta = tables.keys()
        ta.sort()
        for t in ta:
            p = tables[t]
            if verbose is True:
                print("%s: Creating transaction" % t, end='\r')
                sys.stdout.flush()
            c.execute('PRAGMA foreign_keys=off')
            c.execute('BEGIN TRANSACTION')
            c.execute("ALTER TABLE '%s' RENAME TO temp_table" % t)

            if verbose is True:
                print("%s: Creating new table" % t, end='\r')
                sys.stdout.flush()
            create = "CREATE TABLE if not exists '%s' (" % t
            for entry in p:
                create += "%s %s," % (entry[1], entry[2])
            create = create[:-1] + ');'
            c.execute(create)

            if verbose is True:
                print("%s: Creating unique indices" % t, end='\r')
                sys.stdout.flush()
            if is_number(t) or t == 'initial':
                create = "CREATE UNIQUE INDEX "
                create += "if not exists 'uidx' ON '%s' (iterno);" % t
                c.execute(create)
            else:
                create = "CREATE UNIQUE INDEX "
                create += "if not exists 'uidx' ON "
                create += "'%s' (event, station, damp, iter);" % t
                c.execute(create)

            if verbose is True:
                print("%s: Merging tables" % t, end='\r')
                sys.stdout.flush()

            create = "INSERT INTO '%s' (" % t
            for entry in p:
                create += "%s," % (entry[1])
            create = create[:-1] + ')'
            # c.execute(create)

            create += " SELECT "
            for entry in p:
                create += "%s," % (entry[1])
            create = create[:-1]
            create += " FROM temp_table"
            c.execute(create)

            if verbose is True:
                print("%s: Dropping old table" % t, end='\r')
                sys.stdout.flush()
            c.execute("DROP TABLE temp_table")
            db.commit()
            c.execute('PRAGMA foreign_keys=off')
            if verbose is True:
                s = "%s:          Merge complete\n" % t
                # s = '{0:0<55}'.format(s)
                print(s, end='\r')
                sys.stdout.flush()

        db.close()
    except Exception as e:
        print("\n%s" % e)
        db.close()
    return
