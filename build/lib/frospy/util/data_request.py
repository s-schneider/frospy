from __future__ import absolute_import, print_function
import obspy
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Stream
from obspy.geodetics import locations2degrees, gps2dist_azimuth
from obspy.taup import TauPyModel
from obspy.core.event import Catalog, Event, Magnitude, Origin
from obspy.core.event.event import EventDescription
from obspy.clients.fdsn.client import FDSNNoDataException

import sys
from frospy.util.array_util import (center_of_gravity, attach_network_to_traces,
                                  attach_coordinates_to_traces,
                                  geometrical_center)
from frospy.util.read import read_std_cat, read_std_inv, read_st
from frospy.preprocessing.data_correction import (remove_response,
                                                rotate, rotate_cmt,
                                                check_delta_start_origin,
                                                check_stream_stat_doublettes,
                                                select_stations)
from frospy.preprocessing.inspect import remove_delta_pulses
from frospy.util.read import read_inventory as read_inv
from frospy.util.base import (get_chan_from_stream, cut2shortest)
from frospy.util.ah import attach_ah_header
import glob
import shutil
import os
import time
from distutils.dir_util import mkpath

try:
    import instaseis
except ImportError:
    pass


def blockPrint():
    sys.stdout = open(os.devnull, 'w')


# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def downloader(cat, components, sampling_rate, localfolder,
               length_in_days, inv=None, do_process_data=False,
               verbose=False):
    """
    param inv:
    type  inv:

    param cat: catalog of events or cmt_id-string of the event
    type  cat:

    param components: List of Channel names, e.g. ['Z', 'N', 'E']
    type  components: list of strings

    param sampling_rate: sampling_rate in Hz
    type  sampling_rate: float

    param localfolder:
    type  localfolder:

    param length_in_days:
    type  length_in_days:
    """
    def blockPrint():
        sys.stdout = open(os.devnull, 'w')

    # Restore
    def enablePrint():
        sys.stdout = sys.__stdout__

    # Check for right channel, depending on sampling_rate
    if sampling_rate == 0.1:
        channel_key = 'VH'
    elif sampling_rate == 1.:
        channel_key = 'LH'
    elif sampling_rate == 10.:
        channel_key = 'MH'

    if components != '*':
        if type(components) != list:
            msg = "Components must be list of strings e.g.: ['BHZ', 'BHN']"
            print(msg)
            return
        else:
            numcomps = len(components)
            comp_tmp = str()
            component_folder = []
            for _i, comp in enumerate(components):
                if _i == 0:
                    comp_tmp += '%s%s' % (channel_key, comp)
                else:
                    comp_tmp += ',%s%s' % (channel_key, comp)
                # Create subfolder with channel name
                component_folder.append('%s%s' % (channel_key, comp))

            components = comp_tmp

    msg = "Checking inventory and catalog"
    print(msg)

    if type(cat) == str:
        cat = read_std_cat(cat)

    numevents = len(cat)
    numstations = 0
    if inv is not None:
        for net in inv:
            numstations += len(net)

        msg = "Looking for %i events and %i stations\n" % (numevents,
                                                           numstations)
        print(msg)
    else:
        msg = "Looking for Arwens stations"
        print(msg)
        inv = read_std_inv()

    msg = "Data will be saved in %s\n\n" % localfolder
    print(msg)
    print('Downloading data:\n')
    for ev_num, event in enumerate(cat):
        ev_num += 1
        cmt_id = None
        for desc in event.event_descriptions:
            if len(desc.text) == 7:
                if desc.text[0].isdigit() and desc.text[-1].istitle():
                    cmt_id = desc.text
        if cmt_id is None:
            cmt_id = str(event.origins[0].time)
        path = "%s/%s/" % (localfolder, cmt_id)
        cat = obspy.core.event.Catalog()
        cat.append(event)
        start = event.origins[0].time
        end = event.origins[0].time + 3600 * 24 * length_in_days

        pathsinglefiles = "%s/singlefiles/" % path

        folders = [pathsinglefiles]
        for comp in folders:
            mkpath(comp)

        os.chdir(path)
        if not verbose:
            enablePrint()
            print('%i/%i\t%s\n' % (ev_num, numevents, event.short_str()))
            blockPrint()

        out = data_request('IRIS', inv=inv, cat=cat, channels=components,
                           savefile='station', normal_mode_data=True,
                           record_startt=start, record_endt=end,
                           file_format='pickle')
        streamall, inv, cat = out[:]
        print(streamall)
        inv.write('inv.xml', format='STATIONXML')
        cat.write('cat.xml', format='QUAKEML')

        files = glob.glob('*.pickle')
        for f in files:
            print(path + f, pathsinglefiles + f)
            shutil.move(path + f, pathsinglefiles + f)

        if do_process_data:
            if not verbose:
                enablePrint()
                print('Processing data...\n\n')
                blockPrint()
            process_data(path, inv, event, cmt_id, numcomps, sampling_rate,
                         component_folder, localfolder)
        if not verbose:
            enablePrint()

    return


def process_downloader_data(inv, cat, components, sampling_rate, localfolder,
                            inspected_only=False,
                            remove_deltas=False,
                            rotate_traces=True,
                            rm_response=True,
                            rm_tidal=True,
                            keep_longest_traces=True,
                            cut_to_same_length=False,
                            starttime_check=True,
                            mean_snr=None,
                            remtidah='~/bin/remtidah'):
    """
    remtidah = '~/bin/remtidah'

    from frospy.util.data_request import process_downloader_data
    from frospy.util.read import read_std_cat

    inv = None
    cat = read_std_cat('011012A')
    components = ["N", "E", "Z"]
    sampling_rate = 0.1
    localfolder = '/data/simons/3-comp-data/120hrs/VH_raw'
    inspected_only = False
    remove_deltas = True
    rotate_traces = True
    rm_response = True
    rm_tidal = True
    keep_longest_traces = True
    cut_to_same_length = False

    process_downloader_data(inv, cat, components, sampling_rate, localfolder,
                            inspected_only,
                            remove_deltas,
                            rotate_traces,
                            rm_response,
                            rm_tidal,
                            keep_longest_traces,
                            cut_to_same_length)

    """

    # if sampling_rate == 0.1:
    #     channel_key = 'VH'
    # elif sampling_rate == 1.:
    #     channel_key = 'LH'
    # elif sampling_rate == 10.:
    #     channel_key = 'MH'
    #
    # if components != '*':
    #     if type(components) != list:
    #         msg = "Components must be list of strings e.g.: ['BHZ', 'BHN']"
    #         print(msg)
    #         return
    #     else:
    #         components = [channel_key + c for c in components]

    if type(cat) == str:
        cat = read_std_cat(cat)

    for ev_num, event in enumerate(cat):
        ev_num += 1
        cmt_id = None
        for desc in event.event_descriptions:
            if len(desc.text) == 7:
                if desc.text[0].isdigit() and desc.text[-1].istitle():
                    cmt_id = desc.text
        if cmt_id is None:
            cmt_id = str(event.origins[0].time)
        path = "%s/%s/" % (localfolder, cmt_id)
        cat = obspy.core.event.Catalog()
        cat.append(event)

        pathsinglefiles = "%ssinglefiles/" % path

        folders = [pathsinglefiles]
        for comp in folders:
            if not os.path.exists(comp):
                mkpath(comp)

        os.chdir(path)
        with open('error.log', 'a') as logfile:
            logmsg = "\n##### Processing Log at %s  #####\n\n"
            logfile.write(logmsg % time.ctime())

        if not inv:
            msg = "Reading inventory"
            print(msg)
            inv = read_inv('inv.xml')

        # for _i, c in enumerate(components):
        #     if rotate_traces:
        #         if c.endswith('N'):
        #             components[_i] = c.replace('N', 'R')
        #         elif c.endswith('E'):
        #             components[_i] = c.replace('E', 'T')
        #
        #     subfolder = '%s/%s' % (localfolder, components[_i])
        #     if not os.path.exists(subfolder):
        #         mkpath(subfolder)

        pathsinglefiles = path + "/singlefiles/"

        if not len(glob.glob("%s/*.pickle" % pathsinglefiles)) == 0:
            files = glob.glob('%s/*.pickle' % (path))
            for f in files:
                file = f.split('/')[-1]
                shutil.move(path + file, pathsinglefiles+file)

        files = glob.glob('%s/*.pickle' % (pathsinglefiles))

        if len(files) == 0:
            msg = 'No files found'
            print(msg)
            return

        msg = 'Processing:\n'
        msg += '|--- Inspected files ---|'
        print(msg)

        if starttime_check is True:
            starttime_check = 3
        if inspected_only is True or inspected_only == 'both':
            process_data(files, localfolder, inv, event, cmt_id,
                         sampling_rate,
                         inspected_only,
                         remove_delta_pulses,
                         rotate_traces,
                         rm_response,
                         rm_tidal,
                         keep_longest_traces,
                         cut_to_same_length,
                         starttime_check,
                         remtidah=remtidah)

        if inspected_only is False or inspected_only == 'both':
            for j, entry in enumerate(files):
                files[j] = entry.replace('inspected.pickle', 'pickle')
            fname = '%s_raw' % cmt_id
            print('')
            msg = '\n|--- Corresponding raw files ---|'
            print(msg)
            process_data(files, localfolder, inv, event, cmt_id,
                         sampling_rate,
                         inspected_only,
                         remove_deltas,
                         rotate_traces,
                         rm_response,
                         rm_tidal,
                         keep_longest_traces=True,
                         cut_to_same_length=True,
                         starttime_check=starttime_check,
                         remtidah=remtidah,
                         filename=fname,
                         mean_snr=mean_snr)

    return


def process_data(files, savefolder, inv, event, cmt_id,
                 sampling_rate,
                 inspected_only=False,
                 remove_deltas=False,
                 rotate_traces=True,
                 rm_response=True,
                 rm_tidal=True,
                 keep_longest_traces=True,
                 cut_to_same_length=True,
                 starttime_check=3,
                 remtidah='~/bin/remtidah',
                 filename='auto',
                 mean_snr=None):

    """
    Data processing routine for a given structure. More documentation will
    follow

    remtidah = '~/bin/remtidah'
    """
    def update_msg(success_total, failed_total,  filetotal):
        pmsg = 'Processing Data [\033[92m succeed\033[0m/ '
        pmsg += '\033[91mfailed\033[0m/'
        pmsg += ' total]'
        if failed_total != 0:
            pmsg += ' [\033[92m %i\033[0m/ \033[91m%i\033[0m/ %i]:'
        else:
            pmsg += ' [\033[92m %i\033[0m/ %i/ %i]:'

        pmsg = pmsg % (success_total, failed_total,  filetotal)
        return pmsg

    def no_data_left(st, f, log):
        if len(st) == 0:
            msg = '-----------------------------------------------------------'
            msg += "File: %s\n" % f
            msg += "%s\n" % log
            msg += "\nError message: \n%s\nNo Data left in stream\n" % f

            with open('error.log', 'a') as logfile:
                logfile.write(msg)
            return True
        return False

    # Sorting data by channel
    filetotal = len(files)
    failed_total = 0
    success_total = 0
    for _filenum, f in enumerate(files):
        try:
            ofile = '%s%s.ahx' % (cmt_id, success_total)
            st = read_st(f)
            attach_ah_header(st, event, inv, add_response=False)

            pmsg = update_msg(success_total, failed_total,  filetotal)
            readmsg = pmsg + '%40s' % 'Read file'
            log = readmsg
            print('%s' % readmsg, end='\r')
            sys.stdout.flush()

            st = check_stream_stat_doublettes(st)

            if starttime_check:
                cmsg = pmsg + '%40s' % 'Checking starttime of file'
                log += "\n%s" % cmsg
                print('%s' % cmsg, end='\r')
                st = check_delta_start_origin(st, event, starttime_check)
            if cut_to_same_length:
                cmsg = pmsg + '%40s' % 'Checking length of file'
                log += "\n%s" % cmsg
                print('%s' % cmsg, end='\r')
                st = cut2shortest(st)

            if remove_deltas:
                cmsg = pmsg + '%40s' % 'Removing delta pulses'
                log += "\n%s" % cmsg
                print('%s' % cmsg, end='\r')
                for i, tr in enumerate(st):
                    st[i] = remove_delta_pulses(tr, verbose=False)

            if rotate_traces:
                cmt_path = '//nfs/stig/simons/alldata/cmts'
                cmt_file = "%s/%s.cmt" % (cmt_path, cmt_id)
                if os.path.exists(cmt_file):
                    event = read_std_cat(cmt_file)[0]
                    attach_ah_header(st, event, inv, add_response=False)
                    st = rotate_cmt(st, cmt_file=cmt_file,
                                    correction='geocentric', verbose=False)

                    rotmsg = pmsg + '%40s' % 'Rotating Data (local cmt)'
                else:
                    event = event
                    st = rotate(st, inv=inv, event=event, verbose=False,
                                correction='geocentric')

                    rotmsg = pmsg + '%40s' % 'Rotating Data (globalcmt)'
                    with open('error.log', 'a') as logfile:
                        logmsg = "\n%s\n%s" % (f, rotmsg)
                        logfile.write(logmsg)

                print('%s' % rotmsg, end='\r')
                sys.stdout.flush()
                log += "\n%s" % rotmsg

            if no_data_left(st, f, log):
                failed_total += 1
                continue

            if rm_tidal:
                rmtidemsg = pmsg + '%40s' % 'Removing tidal signal'
                log += "\n%s" % rmtidemsg
                print('%s' % rmtidemsg, end='\r')
                sys.stdout.flush()

                ofile_tmp = cmt_id + '%s.wtide.ahx' % _filenum
                st.write(ofile_tmp, format='AH')
                st_old = st.copy()
                os.system('%s %s %s' % (remtidah, ofile_tmp, ofile))
                st = read_st(ofile)
                for tr, tr_old in zip(st, st_old):
                    tr.stats = tr_old.stats
                os.system('rm %s %s' % (ofile_tmp, ofile))

            if rm_response:
                rmsg = pmsg + '%40s' % 'Removing response'
                log += "\n%s" % rmsg
                print('%s' % rmsg, end='\r')
                sys.stdout.flush()
                st = remove_response(st, sampling_rate=sampling_rate,
                                     pre_filt=[0.00015, 0.0002, .5, .51])

                st.sort(['channel'])
                st.sort(['station'])

            if no_data_left(st, f, log):
                failed_total += 1
                continue

            savemsg = pmsg + '%40s' % 'Saving ah-file'
            log += "\n%s" % savemsg
            print('%s' % savemsg, end='\r')
            sys.stdout.flush()

            for c in get_chan_from_stream(st):
                st_select = st.select(channel=c)
                file_path = "%s/%s/%s" % (savefolder, c, cmt_id)
                if not os.path.exists(file_path):
                    os.makedirs(file_path)
                st_path = "%s/%s-%s.ahx" % (file_path, cmt_id, _filenum)
                st_select.write(st_path, format='AH')
            success_total += 1

        except Exception:
            if glob.glob('*.ahx'):
                os.system('rm *.ahx')
            failed_total += 1
            msg = '-----------------------------------------------------------'
            msg += "File: %s\n" % f
            msg += "\nError message: \n%s\n" % (sys.exc_info()[1])

            with open('error.log', 'a') as logfile:
                logfile.write(msg)

    pmsg = update_msg(success_total, failed_total,  filetotal)
    print(pmsg)

    if filename == 'auto':
        filename = cmt_id

    savemsg = 'Building ahx files:\t'
    for c in get_chan_from_stream(st):
        file_path = "%s/%s/%s" % (savefolder, c, cmt_id)
        savemsg += ' %s' % c
        if os.path.exists('%s/%s.ahx' % (file_path, filename)):
            os.system('cat %s/*-*ahx > %s/%s.NEW.ahx' % (file_path, file_path,
                                                         filename))
        else:
            os.system('cat %s/*-*.ahx > %s/%s.ahx' % (file_path, file_path,
                                                      filename))

        rfiles = glob.glob('%s/???????-*.ahx' % file_path)
        for _f in rfiles:
            os.remove('%s' % _f)
        print(savemsg, end='\r')
        sys.stdout.flush()

        if mean_snr is not None:
            snr_msg = 'Mean SNR: %f' % mean_snr
            print(snr_msg, end='\r')
            sys.stdout.flush()
            ahfile = '%s/%s.ahx' % (file_path, filename)
            st = select_stations(ahfile, min_snr=mean_snr)
            if len(st) != 0:
                ahfile = '%s/%s.ahx' % (file_path, cmt_id)
                print(ahfile)
                st.write("%s" % ahfile, format='ah')

    print()

    return


def data_request(client_name, start=None, end=None, minmag=None, cat=None,
                 inv=None, cat_client_name='globalcmt', net="*", scode="*",
                 channels="*", minlat=None, maxlat=None, minlon=None,
                 maxlon=None, station_minlat=None, station_maxlat=None,
                 station_minlon=None, station_maxlon=None,
                 mindepth=None, maxdepth=None, radialcenterlat=None,
                 radialcenterlon=None, minrad=None, maxrad=None,
                 station_radcenlat=None, station_radcenlon=None,
                 station_minrad=None, station_maxrad=None,
                 azimuth=None, baz=False, record_startt=None, record_endt=None,
                 t_before_first_arrival=1, t_after_first_arrival=9,
                 savefile=False, file_format='SAC',
                 normal_mode_data=False):
    """
    Searches in a given Database for seismic data. Restrictions in terms of
    starttime, endtime, network etc can be made. If data is found it returns a
    stream variable, with the waveforms, an inventory with all station and
    network information and a catalog with the event information.

    :param client_name: Name of desired fdsn client,
                        for a list of all clients see:
                        https://docs.obspy.org/tutorial/code_snippets/retrieving_data_from_datacenters.html
    :type  client_name:  string

    :param start, end: starttime, endtime
    :type : UTCDateTime

    :param minmag: Minimum magnitude of event
    :type  minmag: float

    :param cat_client_name: Name of Event catalog, default is "None", resulting
                            in catalog search, defined by client_name

    :type  cat_client_name: string

    :param net: Network code for which to search data for
    :type  net: string

    :param scode: Station code for which to search data for
    :type  scode: string

    :param channels: Used channels of stations
    :type  channels: string

    :param minlat, maxlat, minlon, maxlon: Coordinate-window of interest
    :type : float

    :param mindepth, maxdepth: depth information of event in km
    :type : float

    :param radialcenterlat, radialcenterlon:
    Centercoordinates of a radialsearch, if radialsearch=True
    :type : float

    :param minrad, maxrad: Minimum and maximum radii for radialsearch
    :type : float

    :param azimuth: Desired range of azimuths of event, station couples in deg
                    as a list [minimum azimuth, maximum azimuth]
    :type  azimuth: list

    :param baz: Desired range of back-azimuths of event, station couples in deg
                as a list [minimum back azimuth, maximum back azimuth]
    :type  baz: list

    :param t_before_first_arrival, t_before_after_arrival:
    Length of the seismograms, startingpoint, minutes before 1st arrival and
    minutes after 1st arrival.
    :type  t_before_first_arrival, t_before_after_arrival: float, int

    :param savefile: if True, Stream, Inventory and Catalog will be saved
                     local, in the current directory.
    :type  savefile: bool

    :param format: File-format of the data, for supported formats see:
    https://docs.obspy.org/packages/autogen/obspy.core.stream.Stream.write.html#obspy.core.stream.Stream.write
    :type  format: string

    returns

    :param: streams, Inventory, Catalog
    :type: list of stream objects, obspy Inventory object, obspy Catalog object



    ### Example 1 ###

    from obspy import UTCDateTime
    from bowpy.util.data_request import data_request

    start = UTCDateTime(2010,1,1,0,0)
    end = UTCDateTime(2010,12,31,0,0)
    minmag = 8
    station = '034A'
    streams, inventory, cat = data_request('IRIS', start, end, minmag,
                                                  net='TA', scode=station)

    st = list_of_stream[0]
    st = st.select(channel='BHZ')
    st.normalize()
    inv = inventory[0]

    st.plot()
    inv.plot()
    cat.plot()

    ### Example 2 ###

    from obspy import UTCDateTime
    from frospy.util.data_request import data_request
    from frospy.preprocessing.data_correction import remove_response

    start = UTCDateTime(2010,1,1,0,0)
    end = UTCDateTime(2010,12,31,0,0)
    minmag = 8
    client = 'IRIS'
    fname = 'my_data'
    inv = read_std_inv()

    streams, inventory, cat = data_request(client, start, end, inv=inv)
    for st in streams:
        st = remove_response(st, sampling_rate=0.1,
                             pre_filt=[0.00015, 0.0002, .5, .51])
        st.sort(['channel'])
        st.sort(['station'])
        print('\n\n Saving file ...')
        st.write(fname + '.ahx', format='ah')

    """
    if not cat and not inv:
        if not start and not end and not minmag:
            print('Neither catalog and inventory specified nor dates')
            return

    stream = Stream()
    streamall = []

    if record_startt is not None and record_endt is not None:
        record_len = record_endt - record_startt
    else:
        record_len = None

    # build in different approach for catalog search, using urllib
    if cat:
        catalog = cat
        client = Client(client_name)
    else:
        if cat_client_name == 'globalcmt':
            catalog = request_gcmt(starttime=start, endtime=end,
                                   minmagnitude=minmag, mindepth=mindepth,
                                   maxdepth=maxdepth, minlatitude=minlat,
                                   maxlatitude=maxlat, minlongitude=minlon,
                                   maxlongitude=maxlon)
            client = Client(client_name)
        else:
            client = Client(client_name)
            try:
                catalog = client.get_events(starttime=start, endtime=end,
                                            minmagnitude=minmag,
                                            mindepth=mindepth,
                                            maxdepth=maxdepth,
                                            latitude=radialcenterlat,
                                            longitude=radialcenterlon,
                                            minradius=minrad,
                                            maxradius=maxrad,
                                            minlatitude=minlat,
                                            maxlatitude=maxlat,
                                            minlongitude=minlon,
                                            maxlongitude=maxlon)

            except FDSNNoDataException:
                print("No events found for given parameters.")
                return

    print("Following events found: \n")
    print(catalog)
    m = TauPyModel(model="ak135")

    for event in catalog:
        if inv:
            origin_t = event.origins[0].time
            inventory = inv

        else:
            print("\n")
            print("########################################")
            print("Looking for available data for event: \n")
            print(event.short_str())
            print("\n")

            origin_t = event.origins[0].time
            station_stime = UTCDateTime(origin_t - 3600*24)
            station_etime = UTCDateTime(origin_t + 3600*24)

            try:
                inventory = client.get_stations(network=net, station=scode,
                                                level="response",
                                                channel=channels,
                                                starttime=station_stime,
                                                endtime=station_etime,
                                                minlatitude=station_minlat,
                                                maxlatitude=station_maxlat,
                                                minlongitude=station_minlon,
                                                maxlongitude=station_maxlon,
                                                latitude=station_radcenlat,
                                                longitude=station_radcenlon,
                                                minradius=station_minrad,
                                                maxradius=station_maxrad)

                msg = "Inventory with %i networks, " +\
                      "containing %i stations found."
                print(msg % (len(inventory),
                             len(inventory.get_contents()['stations'])))
            except FDSNNoDataException:
                print("No Inventory found for given parameters")
                continue

        for net in inventory:

            print("Searching in network: %s" % net.code)
            elat = event.origins[0].latitude
            elon = event.origins[0].longitude
            depth = event.origins[0].depth/1000.

            array_fits = True
            if azimuth or baz:
                cog = center_of_gravity(net)
                slat = cog['latitude']
                slon = cog['longitude']
                epid = locations2degrees(slat, slon, elat, elon)
                arr_time = m.get_travel_times(source_depth_in_km=depth,
                                              distance_in_degree=epid)

                # Checking for first arrival time
                P_arrival_time = arr_time[0]

                Ptime = P_arrival_time.time
                tstart = UTCDateTime(event.origins[0].time + Ptime -
                                     t_before_first_arrival * 60)
                tend = UTCDateTime(event.origins[0].time + Ptime +
                                   t_after_first_arrival * 60)

                center = geometrical_center(net)
                clat = center['latitude']
                clon = center['longitude']
                if azimuth:
                    print("Looking for events in the azimuth range of %f to %f\
                          " % (azimuth[0], azimuth[1]))
                    center_az = gps2dist_azimuth(clat, clon, elat, elon)[1]
                    if center_az > azimuth[1] and center_az < azimuth[0]:
                        print("Geometrical center of Array out of azimuth"
                              + " bounds, \nchecking if single stations fit")
                        array_fits = False

                elif baz:
                    print("Looking for events in the back azimuth " +
                          "range of %f to %f" % (baz[0], baz[1]))
                    center_baz = gps2dist_azimuth(clat, clon, elat, elon)[2]
                    if center_baz > baz[1] and center_baz < baz[0]:
                        print("Geometrical center of Array out of back " +
                              "azimuth bounds, \nchecking if " +
                              "single stations fit")
                        array_fits = False

            # If array fits to azimuth/back azimuth or no azimuth/back azimuth
            # is given
            no_of_stations = 0
            if array_fits:

                for station in net:
                    if record_startt and record_endt:
                            tstart = record_startt
                            tend = record_endt
                    else:
                        epid = locations2degrees(station.latitude,
                                                 station.longitude,
                                                 elat, elon)
                        arr_time = m.get_travel_times(source_depth_in_km=depth,
                                                      distance_in_degree=epid)
                        P_arrival_time = arr_time[0]

                        Ptime = P_arrival_time.time
                        if record_startt:
                            tstart = record_startt
                        else:
                            tstart = UTCDateTime(event.origins[0].time + Ptime
                                                 - t_before_first_arrival * 60)
                        if record_endt:
                            tend = record_endt
                        else:
                            if normal_mode_data:
                                tend = UTCDateTime(event.origins[0].time +
                                                   Ptime + 170 * 60 * 60)
                            else:
                                tend = UTCDateTime(event.origins[0].time +
                                                   Ptime +
                                                   t_after_first_arrival * 60)

                    try:
                        if normal_mode_data:
                            # print(tstart, tend, channels,
                            #       net.code, station.code)
                            st_req = client.get_waveforms(network=net.code,
                                                          station=station.code,
                                                          location='*',
                                                          channel=channels,
                                                          starttime=tstart,
                                                          endtime=tend,
                                                          attach_response=True,
                                                          longestonly=True)
                            # print('after download')
                            # print(st_req)
                        else:
                            st_req = client.get_waveforms(network=net.code,
                                                          station=station.code,
                                                          location='*',
                                                          channel=channels,
                                                          starttime=tstart,
                                                          endtime=tend,
                                                          attach_response=True)
                        no_of_stations += 1
                        msg = "Downloaded data for %i of %i available " +\
                              "stations!"
                        print(msg % (no_of_stations,
                                     net.selected_number_of_stations),
                              end='\r')

                        sys.stdout.flush()
                        if record_len is not None:
                            st_req = remove_short_st(st_req, record_len)
                        # print('After removing')
                        # print(st_req)

                        stream += st_req
                    except Exception:
                        # print("Error: %s\n" % e)
                        pass

                    if savefile == 'station' and len(stream) != 0:
                        stname = str(net.code) + '.' + str(station.code) \
                                 + '.' + str(origin_t).split('.')[0]

                        save_file(stream, origin_t, file_format, stname,
                                  inventory, event)
                        msg = "     Saved data for %i of %i available " +\
                              "stations!"
                        print(msg % (no_of_stations,
                                     net.selected_number_of_stations),
                              end='\r')
                        stream = Stream()

            # If not, checking each station individually.
            else:
                for station in net:
                    epid = locations2degrees(station.latitude,
                                             station.longitude, elat, elon)
                    arr_time = m.get_travel_times(source_depth_in_km=depth,
                                                  distance_in_degree=epid)

                    # Checking for first arrival time
                    P_arrival_time = arr_time[0]

                    Ptime = P_arrival_time.time
                    tstart = UTCDateTime(event.origins[0].time + Ptime -
                                         t_before_first_arrival * 60)
                    tend = UTCDateTime(event.origins[0].time + Ptime +
                                       t_after_first_arrival * 60)

                    fit = False
                    if azimuth:
                        stat_az = gps2dist_azimuth(station.latitude,
                                                   station.longitude,
                                                   elat, elon)[1]
                        if stat_az > azimuth[1] and stat_az < azimuth[0]:
                            fit = True
                    elif baz:
                        stat_baz = gps2dist_azimuth(station.latitude,
                                                    station.longitude,
                                                    elat, elon)[2]
                        if stat_baz > baz[1] and stat_baz < baz[0]:
                            fit = True
                    if fit:
                        try:
                            st_req = client.get_waveforms(network=net.code,
                                                          station=station.code,
                                                          location='*',
                                                          channel=channels,
                                                          startime=tstart,
                                                          endtime=tend,
                                                          attach_response=True)
                            no_of_stations += 1
                            msg = "Downloaded data for %i of %i available " +\
                                  "stations!"
                            print(msg % (no_of_stations,
                                         net.selected_number_of_stations),
                                  end='\r')

                            sys.stdout.flush()
                            if record_len is not None:
                                st_req = remove_short_st(st_req, record_len)
                            stream += st_req
                        except Exception:
                            pass

                    if savefile == 'station' and len(stream) != 0:
                        stname = str(net.code) + '.' + str(station.code) \
                                 + '.' + str(origin_t).split('.')[0]

                        save_file(stream, origin_t, file_format, stname,
                                  inventory, event)
                        msg = "     Saved data for %i of %i available " +\
                              "stations!"
                        print(msg % (no_of_stations,
                                     net.selected_number_of_stations),
                              end='\r')
                        stream = Stream()

            if savefile == 'network' and len(stream) != 0:
                stname = str(origin_t).split('.')[0]

                attach_network_to_traces(stream, inventory)
                attach_coordinates_to_traces(stream, inventory, event)

                save_file(stream, origin_t, file_format, stname,
                          inventory, event)
                print('File Saved: %s' % stname)
                stream = Stream()

            print('\n')

        attach_network_to_traces(stream, inventory)
        attach_coordinates_to_traces(stream, inventory, event)

        if savefile == 'event' and len(stream) != 0:
            stname = str(origin_t).split('.')[0]
            invname = stname + "_inv.xml"
            catname = stname + "_cat.xml"

            save_file(stream, origin_t, file_format, stname,
                      inventory, event)
            inventory.write(invname, format="STATIONXML")
            catalog.write(catname, format="QUAKEML")

            print('File Saved: %s' % stname)

        if not savefile:
            streamall.append(stream)

        stream = Stream()

    return(streamall, inventory, catalog)


def remove_short_st(stream, t_in_s):
    for tr in stream:
        rlen = tr.stats.endtime - tr.stats.starttime
        dt = tr.stats.delta
        if rlen + 10*dt < t_in_s:
            stream.remove(tr)
    return stream


def save_file(stream, origin_t, file_format, stname, inventory=None,
              event=None):

    if file_format == 'ah':
        attach_ah_header(stream, inventory, event)
    stname = stname + "." + file_format
    stream.write(stname, format=file_format)


def create_insta_from_invcat(network, event, database):
    """
    This function creates synthetic data using the given network and
    event information, with the database of instaseis

    :param network: Desired Network, for which the data is generated
    :type  network: obspy.core.inventory.Network

    :param event: Event, for wich the data is generated. The event must have
    stored the moment tensor (e.g. given by glogalcmt.org)
    :type  event: obspy.core.event.Event

    :param database: Link to the database, e.g. the path on your harddrive
    :type  database: str
    """

    db = instaseis.open_db(database)

    tofe = event.origins[0].time
    lat = event.origins[0].latitude
    lon = event.origins[0].longitude
    depth = event.origins[0].depth

    source = instaseis.Source(latitude=lat, longitude=lon, depth_in_m=depth,
                              m_rr=event.MomentTensor.m_rr,
                              m_tt=event.MomentTensor.m_tt,
                              m_pp=event.MomentTensor.m_pp,
                              m_rt=event.MomentTensor.m_rt,
                              m_rp=event.MomentTensor.m_rp,
                              m_tp=event.MomentTensor.m_tp,
                              origin_time=tofe
                              )

    stream = Stream()
    tmp = []
    for station in network:
        rec = instaseis.Receiver(latitude=str(station.latitude),
                                 longitude=str(station.longitude),
                                 network=str(network.code),
                                 station=str(station.code))
        tmp.append(db.get_seismograms(source=source, receiver=rec))

    for x in tmp:
        stream += x

    return stream


def request_gcmt(starttime, endtime, minmagnitude=None, mindepth=None,
                 maxdepth=None, minlatitude=None, maxlatitude=None,
                 minlongitude=None, maxlongitude=None):
    from mechanize import Browser
    import re

    """
    Description
    I am using mechanize. My attempt is just preliminary, for the current
    globalcmt.org site. It is possible to store Moment Tensor information
    in the catalog file.
    """

    # Split numbers and text
    re.compile("([a-zA-Z]+)([0-9]+)")
    br = Browser()
    br.open('http://www.globalcmt.org/CMTsearch.html')
    # Site has just one form
    br.select_form(nr=0)

    br.form['yr'] = str(starttime.year)
    br.form['mo'] = str(starttime.month)
    br.form['day'] = str(starttime.day)
    br.form['oyr'] = str(endtime.year)
    br.form['omo'] = str(endtime.month)
    br.form['oday'] = str(endtime.day)
    br.form['list'] = ['4']
    br.form['itype'] = ['ymd']
    br.form['otype'] = ['ymd']
    # Set to full output
    br.find_control(name="list").items[-1].selected = True

    if minmagnitude:
        br.form['lmw'] = str(minmagnitude)
    if minlatitude:
        br.form['llat'] = str(minlatitude)
    if maxlatitude:
        br.form['ulat'] = str(maxlatitude)
    if minlongitude:
        br.form['llon'] = str(minlongitude)
    if maxlongitude:
        br.form['ulon'] = str(maxlongitude)
    if mindepth:
        br.form['lhd'] = str(mindepth)
    if maxdepth:
        br.form['uhd'] = str(maxdepth)

    req = br.submit()

    data = []
    for line in req:
        data.append(line)

    data_chunked = _chunking_list(keyword='\n', list=data)
    origins = []
    magnitudes = []
    tensor = []

    for i, line in enumerate(data_chunked):
        for element in line:
            if 'Event name:' in element:
                event_infos = data_chunked[i:i+4]
                # Event-id and Year/Month/Day are sitting here
                timing = event_infos[0]
                for l in timing:
                    if 'Event name:' in l:
                        event_id = l.split()[2]
                    if "Region name" in l:
                        region = l.split(':')[1].split('<br>')[0].lstrip()
                        region = region.rstrip()
                    if 'Date' in l:
                        year = int(l.split()[2].split('/')[0])
                        mon = int(l.split()[2].split('/')[1])
                        day = int(l.split()[2].split('/')[2].split('<')[0])

                # hour/min/sec info sits here
                hms = event_infos[1][9].split()
                hour = int(hms[1])
                minute = int(hms[2])
                sec = float(hms[3])

                origin = UTCDateTime(year, mon, day, hour, minute, sec)

                # Magnitudes and tensor are sitting here
                mb = float(event_infos[1][8].split()[7])
                ms = float(event_infos[1][8].split()[8])
                mw = float(event_infos[3][1].split()[2])
                scalar_moment = float(event_infos[3][1].split()[-1])
                latitude = float(event_infos[1][9].split()[4])
                longitude = float(event_infos[1][9].split()[5])
                depth = 1000. * float(event_infos[1][9].split()[6])

                # tensor is sitting here
                ex = event_infos[2][2].split()[4]
                m_rr = float(event_infos[2][4].split()[1] + 'E' + ex)
                m_rr_err = float(event_infos[2][5].split()[1] + 'E' + ex)
                m_tt = float(event_infos[2][4].split()[2] + 'E' + ex)
                m_tt_err = float(event_infos[2][5].split()[2] + 'E' + ex)
                m_pp = float(event_infos[2][4].split()[3] + 'E' + ex)
                m_pp_err = float(event_infos[2][5].split()[3] + 'E' + ex)
                m_rt = float(event_infos[2][4].split()[4] + 'E' + ex)
                m_rt_err = float(event_infos[2][5].split()[4] + 'E' + ex)
                m_rp = float(event_infos[2][4].split()[5] + 'E' + ex)
                m_rp_err = float(event_infos[2][5].split()[5] + 'E' + ex)
                m_tp = float(event_infos[2][4].split()[6] + 'E' + ex)
                m_tp_err = float(event_infos[2][5].split()[6] + 'E' + ex)

                magnitudes.append([("Mw", mw), ("Ms", ms), ("Mb", mb)])
                origins.append((latitude, longitude, depth, origin))
                tensor.append((m_rr, m_tt, m_pp, m_rt, m_rp, m_tp,
                               m_rr_err, m_tt_err, m_pp_err, m_rt_err,
                               m_rp_err, m_tp_err, scalar_moment))

    cat = Catalog()

    for mag, org, ten in zip(magnitudes, origins, tensor):
        # Create event object and append to catalog object.
        event = Event()

        for m in mag:
            # Create magnitude object.
            magnitude = Magnitude()
            magnitude.magnitude_type = m[0]
            magnitude.mag = m[1]
            event.magnitudes.append(magnitude)

        # Write origin object.
        origin = Origin()
        origin.latitude = org[0]
        origin.longitude = org[1]
        origin.depth = org[2]
        origin.time = org[3]
        event.origins.append(origin)

        FM = obspy.core.event.FocalMechanism()
        mt = obspy.core.event.source.MomentTensor()
        mt.tensor = ten[0]
        mt.tensor.m_rr_errors = ten[6]
        mt.tensor.m_tt = ten[1]
        mt.tensor.m_tt_errors = ten[7]
        mt.tensor.m_pp = ten[2]
        mt.tensor.m_pp_errors = ten[8]
        mt.tensor.m_rt = ten[3]
        mt.tensor.m_rt_errors = ten[9]
        mt.tensor.m_rp = ten[4]
        mt.tensor.m_rp_errors = ten[10]
        mt.tensor.m_tp = ten[5]
        mt.tensor.m_tp_errors = ten[11]
        mt.scalar_moment = ten[12]
        FM.moment_tensor = mt
        event.focal_mechanisms.append(FM)

        EDname = EventDescription(text=event_id, type="earthquake name")
        EDregion = EventDescription(text=region, type="region name")
        event.event_descriptions.append(EDname)
        event.event_descriptions.append(EDregion)

        cat.append(event)

    return cat


def _chunking_list(keyword, list):
    """
    taken from
    http://stackoverflow.com/questions/19575702/pythonhow-to-split-
    file-into-chunks-by-the-occurrence-of-the-header-word
    """
    chunks = []
    current_chunk = []

    for line in list:

        if line.startswith(keyword) and current_chunk:
            chunks.append(current_chunk[:])
            current_chunk = []

        current_chunk.append(line)
    chunks.append(current_chunk)

    return chunks


def breqfast_template(inv, event, components, length_in_days, option='file',
                      sender=None):

    if type(event) == str:
        filename = 'breqfast_' + event + '.txt'
        cat_eq = read_std_cat()
        bf_id = event
        for q in cat_eq:
            id = q.resource_id.id.split('/')[1]
            if id == event:
                event = obspy.Catalog()
                event.append(q)
                event = event[0]
                break
        ot = event.origins[0].time

    elif type(event) == obspy.core.event.event.Event:
        ot = event.origins[0].time
        bf_id = str(ot.day) + str(ot.month) + str(ot.year)[-2:] + "EVENT"

    et = ot + 3600 * 24 * length_in_days
    magnitude = event.magnitudes[0].mag
    mtype = str(event.magnitudes[0].magnitude_type)
    lat = event.origins[0].latitude
    lon = event.origins[0].longitude
    depth = event.origins[0].depth / 1000.
    year = ot.year
    month = ot.month
    day = ot.day
    hour = ot.hour
    m = ot.minute
    s = ot.second
    ms = ot.microsecond * 1E-6
    s = s + ms
    bfs = str(s).split('.')

    eyear = et.year
    emonth = et.month
    eday = et.day
    ehour = et.hour
    em = et.minute
    es = et.second
    ems = et.microsecond * 1E-6
    es = es + ems
    ebfs = str(es).split('.')

    chan = components[0]
    cc = len(components)
    for i, _c in enumerate(components):
        if i != 0:
            chan += " " + _c

    if option == 'file':
        with open(filename, 'w') as fh:
            fh.write('.NAME Simon Schneider\n')
            fh.write('.INST Utrecht University\n')
            mail = '.MAIL Department of Earth Sciences, Heidelberglaan 2, '
            mail += 'Utrecht, the Netherlands\n'
            fh.write(mail)
            fh.write('.EMAIL s.a.schneider@uu.nl\n')
            fh.write('.PHONE +31 6 3117 9333\n')
            fh.write('.MEDIA Electronic (aFTP)\n')
            fh.write('.ALTERNATE MEDIA EXABYTE - 2 gigabyte\n')
            fh.write('.ALTERNATE MEDIA EXABYTE - 5 gigabyte\n')
            fh.write('.LABEL %s\n' % bf_id)
            fh.write('.SOURCE ~CMT catalogue~\n')
            h = '.HYPO ~%4i %0.02i %0.02i %0.02i' + \
                '%0.02i %2.2f~%3.3f~%3.3f~%4.2f\n'
            fh.write(h % (year, month, day, hour, m, s, lat, lon, depth))
            fh.write('.MAGNITUDE ~%1.1f~%s~/n' % (magnitude, mtype))
            fh.write('.QUALITY B\n')
            fh.write('.END\n\n')

            if type(inv) == obspy.core.inventory.inventory.Inventory:
                for net in inv:
                    for station in net:
                        req = "%s %s %i %i %i %i %i %i %i %i %i %i %i " + \
                              "%i %i %s --\n"
                        fh.write(req % (str(station.code), str(net.code),
                                 year, month, day, hour,
                                 m, int(bfs[0]), eyear, emonth, eday, ehour,
                                 em, int(ebfs[0]), cc, chan))

            if type(inv) == obspy.core.inventory.network.Network:
                for station in inv:
                    req = "%s %s %i %i %i %i %i %i %i %i %i %i %i " + \
                          "%i %i %s --\n"
                    fh.write(req % (str(station.code), str(inv.code), year,
                             month, day, hour, m, int(bfs[0]), eyear,
                             emonth, eday, ehour, em, int(ebfs[0]), cc,
                             chan))

    elif option == 'send' and sender:
        mail = '.NAME Simon Schneider\n'
        mail += '.INST Utrecht University\n'
        mail += '.MAIL Department of Earth Sciences, Heidelberglaan 2, '
        mail += 'Utrecht, the Netherlands\n'
        mail += '.EMAIL s.a.schneider@uu.nl\n'
        mail += '.PHONE +31 6 3117 9333\n'
        mail += '.MEDIA Electronic (FTP)\n'
        mail += '.ALTERNATE MEDIA EXABYTE - 2 gigabyte\n'
        mail += '.ALTERNATE MEDIA EXABYTE - 5 gigabyte\n'
        mail += '.LABEL %s\n' % bf_id
        mail += '.SOURCE ~CMT catalogue~\n'
        h = '.HYPO ~%4i %0.02i %0.02i %0.02i' + \
            '%0.02i %2.2f~%3.3f~%3.3f~%4.2f\n'
        mail += h % (year, month, day, hour, m, s, lat, lon, depth)
        mail += '.MAGNITUDE ~%1.1f~%s~/n' % (magnitude, mtype)
        mail += '.QUALITY B\n'
        mail += '.END\n\n'

        if type(inv) == obspy.core.inventory.inventory.Inventory:
            for net in inv:
                req = "%s %s %i %i %i %i %i %i %i %i %i %i %i %i %i %s --\n"
                mail += req % (req % (str(station.code), str(net.code),
                               year, month, day, hour, m, int(bfs[0]), eyear,
                               emonth, eday, ehour, em, int(ebfs[0]), cc,
                               chan))

        if type(inv) == obspy.core.inventory.network.Network:
            for station in inv:
                req = "%s %s %i %i %i %i %i %i %i %i %i %i %i %i %i %s --\n"
                mail += req % (str(station.code), str(inv.code), year, month,
                               day, hour, m, int(bfs[0]), eyear, emonth, eday,
                               ehour, em, int(ebfs[0]), cc, chan)

        print("Sending mail...")
        # send_breqfast(mail, sender)
    return


# def send_breqfast(mail, sender):
#     msg = MIMEText(mail)
#     msg['Subject'] = 'request'
#     receiver = 'breq_fast@iris.washington.edu'
#     msg['From'] = sender
#     msg['To'] = receiver
#     s = smtplib.SMTP('localhost')
#     s.sendmail(sender, [receiver], msg.as_string())
#     s.quit()
#     return
