# -*- coding: utf-8 -*-
"""
Module for handling nmPy Modes objects.

:copyright:
    Simon Schneider (s.a.schneider@uu.nl)
:license:
    GNU Lesser General Public License, Version 3
    (https://www.gnu.org/copyleft/lesser.html)
"""

from future.utils import native_str
from obspy.core.util import AttribDict
from copy import deepcopy
import fnmatch
from frospy import data as frospydata
from frospy.util.base import split_digit_nondigit
import json
import os
import pickle
import sys


def read(ifile=None, format=None, modenames=None):
    """
    modenames: string or list of strings
    """
    if modenames is not None:
        if type(modenames) is not str:
            # if type(modenames) is not unicode:
            if type(modenames) is not list:
                raise IOError('modenames is neither string nor list')
        # Put in correct format
        if type(modenames) is str:  # or type(modenames) is unicode:
            modenames = [modenames]
        for i, m in enumerate(modenames):
            m = split_digit_nondigit(m)
            m = ''.join([str(int(m[0])), m[1].upper(), str(int(m[2]))])
            modenames[i] = m

    if format is None and ifile is not None:
        # try to guess format from file extension
        _, format = os.path.splitext(ifile)
        format = format[1:]

    if format == 'pickle':
        return _read_pickle(ifile)

    elif format == 'json':
        modes = Modes()
        with open(ifile) as json_file:
            data = json.load(json_file)

        for values in data['modes']:
            modes += Mode(header=values)
        return modes

    else:
        path = frospydata.__path__[0] + "/AD/modes-full.json"
        modes = Modes()
        with open(path) as json_file:
            data = json.load(json_file)

        for values in data['modes']:
            if modenames is not None:
                if values['name'] not in modenames:
                    continue
            modes += Mode(header=values)

        # Old .dat file, use .json instead
        # path = frospydata.__path__[0] + "/AD/modes-full.dat"
        # modesl = np.loadtxt(path, dtype={'names': ('n', 'name', 'l',
        #                                  'frequency', 'Q', 'sensitivity'),
        #                     'formats': (np.int, '|S15', np.int, np.float,
        #                                 np.float, '|S15')},
        #                     delimiter=', ').tolist()
        #
        # modes = Modes()
        # for entry in modesl:
        #     name = str(entry[0]) + str(entry[1]).upper() + str(entry[2])

        #     mode = Mode({'n': entry[0],
        #                  'type': str(entry[1]).upper(),
        #                  'l': entry[2],
        #                  'freq': entry[3],
        #                  'Q': entry[4],
        #                  'sens': entry[5],
        #                  'name': name,
        #                  })
        #     modes += mode
        return modes


class Mode(AttribDict):
    """
    :type n: integer
    :param n: overtone number of mode
    :type type: string
    :param type: S or T for Spheroidal or Toroidal
    :type l: integer
    :param l: angular order
    :type name: string
    :param name: name of the mode
    :type sens: string
    :param sens: sensitivity of the mode
    :type freq: float
    :param freq: frequency in PREM
    :type Q: float
    :param Q: Q in PREM
    :type qcycle: float
    :param qcycle: Q-cycle of mode in hours

    .. rubric:: Example

    1. Qcycle test
        >>> mode = Mode({'freq': 1.0, 'Q': 1.0})
        >>> mode.qcycle
        0.6111111111111112
    """

    readonly = ['qcycle']
    defaults = {
        'n': -1,
        'type': '',
        'l': -1,
        'name': '',
        'sens': '',
        'freq': 0.0,
        'Q': 0.0,
        'qcycle': 0.0,
        'avg_tw': None
    }

    _refresh_keys = {'n', 'type', 'l', 'name', 'sens', 'freq', 'Q', 'avg_tw'}

    def __init__(self, header={}, **args):
        if header == {}:
            if len(args) != 0:
                header = args.copy()
                header['name'] = "{}{}{}".format(header['n'], header['type'].upper(), header['l'])
        super(Mode, self).__init__(header)

    def __setitem__(self, key, value):
        """
        """
        if key in self._refresh_keys:
            if key in ['n', 'l']:
                value = int(value)
            elif key in ['type', 'name', 'sens']:
                value = str(value)
            elif key in ['freq', 'Q']:
                value = float(value)
            elif key in ['avg_tw']:
                if value is not None:
                    value = [float(x) for x in value]
            # equivalent to AttribDict.__setitem__(self, key, value)
            super(Mode, self).__setitem__(key, value)
            try:
                Q = self.Q
                f = self.freq / 1000.
                qc = 1.1*(Q/f) / 3600.
                self.__dict__['qcycle'] = qc
            except ZeroDivisionError:
                self.__dict__['qcycle'] = 0.0

    __setattr__ = __setitem__

    def __str__(self, extended=False):
        """
        Return better readable string representation of Stats object
        """
        _pretty_str = '%s |' % (self.name)
        if extended is True:
            _pretty_str += ' n: %2d | type: %s | l: %2d |' % (self.n,
                                                              self.type,
                                                              self.l)
        _pretty_str += ' freq: %.3f | Q: %7.3f' % (self.freq, self.Q)
        if extended is True:
            _pretty_str += ' | sensitivity: %4s |' % (self.sens)
            _pretty_str += ' qcycle: %.3f |' % (self.qcycle)
            if self.avg_tw is not None:
                _str = ' Avg. Timewindow: %.1f - %.1f' % (self.avg_tw[0],
                                                          self.avg_tw[1])
                _pretty_str += _str
        return _pretty_str

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))

    def folder_name(self):
        return "%02d%s%02d" % (self.n, self.type.lower(), self.l)


class Modes(object):
    """
    Container for frospy.core.modes.Mode AttribDicts
    """
    defaults = {
        'names': []
    }

    _refresh_keys = {'names'}

    def __init__(self, modes=None):
        self.modes = []
        self.names = []
        if isinstance(modes, Mode):
            modes = [modes]
        if modes is not None:
            self.modes.extend(modes)
            for m in modes:
                self.names += [m.name]

    def __str__(self, extended=False):
        out = str(len(self.modes)) + ' Mode(s) in Modes:\n'
        if len(self.modes) <= 20 or extended is True:
            out = out + "\n".join([_i.__str__() for _i in self])
        else:
            out = out + "\n" + self.modes[0].__str__() + "\n" + \
                '...\n(%i other modes)\n...\n' % (len(self.modes) - 2) + \
                self.modes[-1].__str__() + '\n\n[Use "print(' + \
                'Modes.__str__(extended=True))" to print all modes]'
        return out

    def _repr_pretty_(self, p, cycle):
        p.text(self.__str__(extended=p.verbose))

    def __nonzero__(self):
        """
        A Modes object is considered zero if contains no Mode.
        """
        return bool(len(self.modes))

    def __len__(self):
        """
        Return number of Mode in the Modes object.

        .. rubric:: Example

        >>> modes = Modes([Mode(), Mode(), Mode()])
        >>> len(modes)
        3
        """
        return len(self.modes)

    def __iter__(self):
        """
        Return a robust iterator for Modes.modes.

        Doing this it is safe to remove modes from Modes inside of
        for-loops using Modes's :meth:`~frospy.core.modes.Modes.remove`
        method. Actually this creates a new iterator every time a mode is
        removed inside the for-loop.

        """
        return list(self.modes).__iter__()

    def __add__(self, other):
        """
        Add two Modes or a Modes with a single Mode.

        :type other: :class:`~frospy.core.modes.Modes` or
            :class:`~frospy.core.modes.Mode`
        :param other: Modes or Mode object to add.
        :rtype: :class:`~frospy.core.modes.Modes`
        :returns: New Modes object containing references to the Mode of the
            original Modes

        .. rubric:: Examples

        1. Adding two Modes

            >>> st1 = Modes([Mode(), Mode(), Mode()])
            >>> len(st1)
            3
            >>> st2 = Modes([Mode(), Mode()])
            >>> len(st2)
            2
            >>> Modes = st1 + st2
            >>> len(Modes)
            5

        2. Adding Modes and Mode

            >>> Modes2 = st1 + Mode()
            >>> len(Modes2)
            4
        """
        if isinstance(other, Mode):
            other = Modes([other])
        if not isinstance(other, Modes):
            raise TypeError
        modes = self.modes + other.modes
        return self.__class__(modes=modes)

    def __iadd__(self, other):
        """
        Add two modes with self += other.

        .. rubric:: Example

        >>> modes = Modes([Mode(), Mode(), Mode()])
        >>> len(modes)
        3

        >>> modes += Modes([Mode(), Mode()])
        >>> len(modes)
        5

        >>> modes += Mode()
        >>> len(modes)
        6
        """
        if isinstance(other, Mode):
            other = Modes([other])
        if not isinstance(other, Modes):
            raise TypeError
        self.extend(other.modes)
        return self

    def __getitem__(self, index):
        """
        __getitem__ method of frospy.Modes objects.

        :return: Mode objects
        """
        return self.modes.__getitem__(index)

    def append(self, mode):
        """
        Append a single Mode object to the current Modes object.
        """
        if isinstance(mode, Mode):
            self.modes.append(mode)
        else:
            msg = 'Append only supports a single mode object as an argument.'
            raise TypeError(msg)
        return self

    def copy(self):
        return deepcopy(self)

    def extend(self, mode_list):
        """
        Extend the current Modes object with a list of Mode objects.
        """
        if isinstance(mode_list, list):
            for _i in mode_list:
                # Make sure each item in the list is a trace.
                if not isinstance(_i, Mode):
                    msg = 'Extend only accepts a list of Mode objects.'
                    raise TypeError(msg)
            self.modes.extend(mode_list)
            for m in mode_list:
                self.names += [m.name]
        elif isinstance(mode_list, Modes):
            self.modes.extend(mode_list.modes)
            for m in mode_list:
                self.names += [m.name]
        else:
            msg = 'Extend only supports a list of Mode objects as argument.'
            raise TypeError(msg)
        return self

    def remove(self, mode):
        """
        Remove the first occurrence of the specified Mode object in the
        Modes object. Passes on the remove() call to self.modes.
        """
        self.modes.remove(mode)
        self.names.remove(mode.name)
        return self

    def select(self, name=None, mtype=None, sens=None, overtone=None,
               angular_order=None):
        """
        Return new Modes object only with these modes that match the given
        criteria (e.g. all modes with ``type="S"``).
        """

        if name:
            # Checking for right format: Digit Letter Digit
            name = format_name(name)
        modes = []
        for mode in self:
            # skip trace if any given criterion is not matched
            if name is not None:
                if not fnmatch.fnmatch(str(mode.name.upper()),
                                       str(name).upper()):
                    continue
            if mtype is not None:
                if not fnmatch.fnmatch(mode.type.upper(),
                                       mtype.upper()):
                    continue
            if sens is not None:
                if not fnmatch.fnmatch(mode.sens.upper(),
                                       sens.upper()):
                    continue
            if overtone is not None:
                if mode.n != overtone:
                    continue
            if angular_order is not None:
                if mode.l != angular_order:
                    continue

            modes.append(mode)
        return self.__class__(modes=modes)

    def sort(self, keys=['freq', 'n', 'l', 'Q', 'qcycle',
                         'type', 'sens', 'name', 'cc'], reverse=False):
        """
        Sort the modes in the Modes object.
        """
        # check if list
        msg = "keys must be a list of strings. Always available items to " + \
            "sort after: \n'freq', 'n', 'l', 'Q', 'qcycle', 'type', 'sens'" + \
            " 'name', 'cc'"
        if not isinstance(keys, list):
            raise TypeError(msg)
        # Loop over all keys in reversed order.
        for _i in keys[::-1]:
            if _i == 'cc':
                self.modes.sort(
                    key=lambda x: tuple(map(str, x.name[2:].split("-")))
                                )
            elif _i == 'sens':
                self.modes.sort(
                    key=lambda x: getattr(x, _i), reverse=(not reverse)
                                )
            else:
                self.modes.sort(key=lambda x: getattr(x, _i), reverse=reverse)
        return self

    def write(self, fname, overwrite=False, format=None):
        try:
            if not overwrite:
                i = 0
                while True:
                    msg = ''
                    f, ext = os.path.splitext(fname)
                    if os.path.exists(fname):
                        msg += '\033[93mMode-file exist,'
                        msg += 'not overwriting\033[0m'
                        if i == 0:
                            f = "%s_%s" % (f, str(i))
                        i += 1
                        a = "_%s" % str(i-1)
                        b = "_%s" % str(i)
                        f = f.replace(a, b)
                        fname = ''.join([f, ext])
                    else:
                        print(msg)
                        break
            filename = fname
            if format is None:
                # try to guess format from file extension
                _, format = os.path.splitext(filename)
                format = format[1:]

            if format == 'pickle':
                _write_pickle(self, filename)
            else:
                with open(filename, 'w') as fh:
                    # for each stationname 'key' fw and tw is written
                    for mode in self:
                        fh.write("%s\t " % str(mode.name))
                        fh.write("%d\t%s\t%d\t" % (mode.n, mode.type, mode.l))
                        fh.write("%f\t%f\t" % (mode.freq, mode.Q))
                        fh.write("%f\t%s\n" % (mode.qcycle, mode.sens))
                msg = "\033[92mMode-file written to %s\033[0m" % filename
                print(msg)
        except IOError:
            msg = "\033[91mCan't save file\n"
            msg += "Error message: %s\033[0m" % sys.exc_info()[1]
            print(msg)
        return


def _read_pickle(filename, **kwargs):
    """
    Read and return Modes from pickled Modes file.

    .. warning::
        This function should NOT be called directly, it registers via the
        nmPy :func:`~frospy.core.modes.read` function, call this instead.

    :type filename: str
    :param filename: Name of the pickled Modes file to be read.
    :rtype: :class:`~frospy.core.modes.Modes`
    :return: A Modes object.
    """
    kwargs = {}

    if isinstance(filename, (str, native_str)):
        with open(filename, 'rb') as fp:
            return pickle.load(fp, **kwargs)
    else:
        return pickle.load(filename, **kwargs)


def _write_pickle(mode, filename, protocol=2, **kwargs):
    """
    Write a Python pickle of current modes.

    .. note::
        Writing into PICKLE format allows to store additional attributes
        appended to the current Modes object or any contained Mode.

    .. warning::
        This function should NOT be called directly, it registers via the
        the :meth:`~frospy.core.modes.Modes.write` method of an
        nmPy :class:`~frospy.core.modes.Modes` object, call this instead.

    :type mode: :class:`~frospy.core.modes.Modes`
    :param mode: The Modes object to write.
    :type filename: str
    :param filename: Name of file to write.
    :type protocol: int, optional
    :param protocol: Pickle protocol, defaults to ``2``.
    """
    if isinstance(filename, (str, native_str)):
        with open(filename, 'wb') as fp:
            pickle.dump(mode, fp, protocol=protocol)
    else:
        pickle.dump(mode, filename, protocol=protocol)


def format_name(mode, number_of_digits=1):
    def name(mode):
        m = split_digit_nondigit(mode)
        for i, x in enumerate(m):
            if x.isdigit() and int(x) != 0:
                N = len(x.lstrip('0'))
                if number_of_digits < N:
                    m[i] = x.lstrip('0')
                else:
                    m[i] = '0' * abs(N - number_of_digits) + x.lstrip('0')
            elif x.isdigit() and int(x) == 0:
                m[i] = '0' * number_of_digits
        return ''.join(m)

    if '-' in mode:
        _mode = mode.split('-')
        _name = []
        for x in _mode:
            _name += [name(x)]
        return '-'.join(_name)
    else:
        return name(mode)
