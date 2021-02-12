from future.utils import native_str
import pickle


def write_pickle(object, filename, protocol=2, **kwargs):
    """
    Write a Python pickle of current segments.

    .. note::
        Writing into PICKLE format allows to store additional attributes
        appended to the current Segment object or any contained Pick.

    .. warning::
        This function should NOT be called directly, it registers via the
        the :meth:`~frospy.core.segment.Segment.write` method of an
        nmPy :class:`~frospy.core.segment.Segment` object, call this instead.

    :type object: any python type
    :param object: The object to write.
    :type filename: str
    :param filename: Name of file to write.
    :type protocol: int, optional
    :param protocol: Pickle protocol, defaults to ``2``.
    """
    if isinstance(filename, (str, native_str)):
        with open(filename, 'wb') as fp:
            pickle.dump(object, fp, protocol=protocol)
    else:
        pickle.dump(object, filename, protocol=protocol)
