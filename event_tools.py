import os
import h5py
import netCDF4 as nc
import numpy as np
import obspy as op
from obspy.clients.fdsn import Client
from obspy.core.event import Catalog
from obspy import UTCDateTime
from obspy.geodetics.base import locations2degrees, degrees2kilometers
from obspy.taup import TauPyModel
from scipy.signal import butter, filtfilt, detrend

def get_event_list(t1,t2):
    '''
    In
        t1, t2: start and ending timestamps
    Out
        cat : events in USGS AK catalog exceeding GMM threshold
        ttimes : P-wave arrival times from IASPEI91
        otimes : origin times (for first available solution)
    '''
    # Reference point near Homer, AK
    lat0 = 59.4405426
    lon0 = -152.0276765
    # Get catalog events
    #client = Client('USGS')
    #catalog = client.get_events(starttime=t1,endtime=t2,catalog='ak',\
    #        includeallorigins=True,includeallmagnitudes=True)
    client = Client('IRIS')
    catalog = client.get_events(starttime=t1, endtime=t2,\
            includeallorigins=True, includeallmagnitudes=True,\
            latitude=lat0, longitude=lon0, maxradius=10)
    # Parameters for simple GMM
    a = -1
    b = 0.55
    # Set up model for travel-times
    model = TauPyModel(model='iasp91')
    # Remove weak events and get travel-times
    events = []
    otimes = []
    ttimes = []
    for event in catalog:
        lon = event.origins[0]['longitude']
        lat = event.origins[0]['latitude']
        dep = event.origins[0]['depth'] * 1e-3
        dis1 = locations2degrees(lat0,lon0,lat,lon)
        dis = degrees2kilometers(dis1)
        R = np.sqrt(dis**2 + dep**2)
        M = event.magnitudes[0]['mag']
        if (M - 10**(a + b*np.log10(R)) >= 0):
            events.append(event)
            arr = model.get_travel_times(source_depth_in_km=dep,distance_in_degree=dis1)
            ttimes.append(arr[0].time)
            otimes.append(event.origins[0]['time'])
    cat = Catalog(events=events)
    ttimes = np.array(ttimes)
    otimes = np.array(otimes)
    return cat, ttimes, otimes

def get_file_list(fdir,ftemplate,t1,t2):
    '''
    In
        fdir : path to recording directory
        ftemplate : file name format (e.g. 'decimator2_%Y-%m-%d_%H.%M.%S_UTC.h5')
        t1, t2 : start and ending timestamps
    Out
        flist : list of files between t1 and t2 (from file names)
        tlist : list of UTCDateTime from file names in flist
    '''
    # Get all available files
    tmp = os.listdir(fdir)
    flist = np.array([fname for fname in tmp if fname.endswith('h5')])
    tlist = np.array([UTCDateTime.strptime(fname,format=ftemplate) for fname in flist])
    # Select files between t1, t2 (with a few minutes buffer)
    index = np.logical_and((tlist-t1-180)>=0,(tlist-t2+180)<=0)
    flist = flist[index]; tlist = tlist[index]
    # Sort files by time
    index = np.argsort(tlist)
    flist = flist[index]; tlist = tlist[index]
    flist = [os.path.join(fdir,fname) for fname in flist]
    return flist, tlist

'''
# Deprecated, extra space is reserved resulting in unnecessarily large files
def write_event_file(outpath,flist,tlist,stime,nMin):
    #
    #In
    #    outpath : output file path including HDF5 name
    #    flist : sorted list of input files
    #    tlist : sorted UTCDateTime for each file
    #    ttime : modeled P-wave arrival time 
    #    nMin : number of minutes of data to save
    #Out
    #    ...
    #
    etime = stime + nMin*60
    ststamp = stime.timestamp*1e6
    etstamp = etime.timestamp*1e6
    # List of files containing event data
    index = np.searchsorted(tlist,stime,side='left') - 1
    files = flist[index:index+nMin+1]
    if not len(files)==nMin+1:
        print('Files missing for %s' % outpath)
    else:
        # Read and concatenate data
        data = []
        dataTime = []
        dataSample = []
        with h5py.File(files[0],'r') as fp:
            index = np.searchsorted(fp['Acquisition']['Raw[0]']['RawDataTime'][:],ststamp,side='left')
            data.append(fp['Acquisition']['Raw[0]']['RawData'][index:,:])
            dataTime.append(fp['Acquisition']['Raw[0]']['RawDataTime'][index:])
            dataSample.append(fp['Acquisition']['Raw[0]']['RawDataSampleCount'][index:])
        for kk in range(1,nMin):
            with h5py.File(files[kk],'r') as fp:
                data.append(fp['Acquisition']['Raw[0]']['RawData'][:])
                dataTime.append(fp['Acquisition']['Raw[0]']['RawDataTime'][:])
                dataSample.append(fp['Acquisition']['Raw[0]']['RawDataSampleCount'][:])
        with h5py.File(files[-1],'r') as fp:
            index = np.searchsorted(fp['Acquisition']['Raw[0]']['RawDataTime'][:],etstamp,side='left')
            data.append(fp['Acquisition']['Raw[0]']['RawData'][:index,:])
            dataTime.append(fp['Acquisition']['Raw[0]']['RawDataTime'][:index])
            dataSample.append(fp['Acquisition']['Raw[0]']['RawDataSampleCount'][:index])
        data = np.concatenate(data,axis=0)
        dataTime = np.concatenate(dataTime,axis=0)
        dataSample = np.concatenate(dataSample,axis=0)
        # Write output file
        with h5py.File(outpath,'w') as fp_out:
            # Open first file and copy over format
            with h5py.File(files[0],'r') as fp_in:
                fp_in.copy(fp_in['Acquisition'],fp_out['/'],shallow=True)
                del fp_out['Acquisition']['Raw[0]']
                fp_in.copy(fp_in['Acquisition']['Raw[0]'],fp_out['Acquisition'],shallow=True)
                del fp_out['Acquisition']['Raw[0]']['RawData']
                del fp_out['Acquisition']['Raw[0]']['RawDataSampleCount']
                del fp_out['Acquisition']['Raw[0]']['RawDataTime']
                # Add datasets
                fp_out['Acquisition']['Raw[0]'].create_dataset('RawData',data=data)
                fp_out['Acquisition']['Raw[0]'].create_dataset('RawDataTime',data=dataTime)
                fp_out['Acquisition']['Raw[0]'].create_dataset('RawDataSampleCount',data=dataSample)
                # Ammend attributes
                fp_out['Acquisition']['Raw[0]']['RawData'].attrs['Count'] = np.prod(data.shape)
                fp_out['Acquisition']['Raw[0]']['RawData'].attrs['Dimensions'] = 'time, locus'
                fp_out['Acquisition']['Raw[0]']['RawData'].attrs['StartIndex'] = \
                        fp_in['Acquisition']['Raw[0]']['RawData'].attrs['StartIndex']
                StartTime = UTCDateTime(dataTime[0]*1e-6).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                EndTime = UTCDateTime(dataTime[-1]*1e-6).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                fp_out['Acquisition']['Raw[0]']['RawData'].attrs['PartStartTime'] = StartTime
                fp_out['Acquisition']['Raw[0]']['RawData'].attrs['PartEndTime'] = EndTime
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['Count'] = len(dataTime)
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['StartIndex'] = \
                        fp_in['Acquisition']['Raw[0]']['RawDataTime'].attrs['StartIndex']
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['StartTime'] = StartTime
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['EndTime'] = EndTime
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['PartStartTime'] = StartTime
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['PartEndTime'] = EndTime
    return
'''

def write_event_file(outpath,flist,tlist,stime,nMin):
    '''
    In
        outpath : output file path including HDF5 name
        flist : sorted list of input files
        tlist : sorted UTCDateTime for each file
        ttime : modeled P-wave arrival time
        nMin : number of minutes of data to save
    Out
        ...
    '''
    etime = stime + nMin*60
    ststamp = stime.timestamp*1e6
    etstamp = etime.timestamp*1e6
    # List of files containing event data
    index = np.searchsorted(tlist,stime,side='left') - 1
    files = flist[index:index+nMin+1]
    if not len(files)==nMin+1:
        print('Files missing for %s' % outpath)
    else:
        # Read and concatenate data
        data = []
        dataTime = []
        dataSample = []
        with h5py.File(files[0],'r') as fp:
            index = np.searchsorted(fp['Acquisition']['Raw[0]']['RawDataTime'][:],ststamp,side='left')
            data.append(fp['Acquisition']['Raw[0]']['RawData'][index:,:])
            dataTime.append(fp['Acquisition']['Raw[0]']['RawDataTime'][index:])
            dataSample.append(fp['Acquisition']['Raw[0]']['RawDataSampleCount'][index:])
        for kk in range(1,nMin):
            with h5py.File(files[kk],'r') as fp:
                data.append(fp['Acquisition']['Raw[0]']['RawData'][:])
                dataTime.append(fp['Acquisition']['Raw[0]']['RawDataTime'][:])
                dataSample.append(fp['Acquisition']['Raw[0]']['RawDataSampleCount'][:])
        with h5py.File(files[-1],'r') as fp:
            index = np.searchsorted(fp['Acquisition']['Raw[0]']['RawDataTime'][:],etstamp,side='left')
            data.append(fp['Acquisition']['Raw[0]']['RawData'][:index,:])
            dataTime.append(fp['Acquisition']['Raw[0]']['RawDataTime'][:index])
            dataSample.append(fp['Acquisition']['Raw[0]']['RawDataSampleCount'][:index])
        data = np.concatenate(data,axis=0)
        dataTime = np.concatenate(dataTime,axis=0)
        dataSample = np.concatenate(dataSample,axis=0)
        # Write output file
        with h5py.File(outpath,'w') as fp_out:
            # Open first file and copy over format
            with h5py.File(files[0],'r') as fp_in:
                acquisition = fp_out.create_group('Acquisition')
                for attribute in list(fp_in['Acquisition'].attrs):
                    fp_out['Acquisition'].attrs[attribute] = fp_in['Acquisition'].attrs[attribute].copy()
                raw0 = fp_out['Acquisition'].create_group('Raw[0]')
                for attribute in list(fp_in['Acquisition']['Raw[0]'].attrs):
                    fp_out['Acquisition']['Raw[0]'].attrs[attribute] = \
                            fp_in['Acquisition']['Raw[0]'].attrs[attribute].copy()
                # Add datasets
                fp_out['Acquisition']['Raw[0]'].create_dataset('RawData',data=data)
                fp_out['Acquisition']['Raw[0]'].create_dataset('RawDataTime',data=dataTime)
                fp_out['Acquisition']['Raw[0]'].create_dataset('RawDataSampleCount',data=dataSample)
                # Clone dataset attributes
                for attribute in list(fp_in['Acquisition']['Raw[0]']['RawData'].attrs):
                    fp_out['Acquisition']['Raw[0]']['RawData'].attrs[attribute] = \
                            fp_in['Acquisition']['Raw[0]']['RawData'].attrs[attribute].copy()
                for attribute in list(fp_in['Acquisition']['Raw[0]']['RawDataTime'].attrs):
                    fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs[attribute] = \
                            fp_in['Acquisition']['Raw[0]']['RawDataTime'].attrs[attribute].copy()
                for attribute in list(fp_in['Acquisition']['Raw[0]']['RawDataSampleCount'].attrs):
                    fp_out['Acquisition']['Raw[0]']['RawDataSampleCount'].attrs[attribute] = \
                            fp_in['Acquisition']['Raw[0]']['RawDataSampleCount'].attrs[attribute].copy()
                # Ammend attributes
                fp_out['Acquisition']['Raw[0]']['RawData'].attrs['Count'] = np.prod(data.shape)
                fp_out['Acquisition']['Raw[0]']['RawData'].attrs['Dimensions'] = 'time, locus'
                fp_out['Acquisition']['Raw[0]']['RawData'].attrs['StartIndex'] = \
                        fp_in['Acquisition']['Raw[0]']['RawData'].attrs['StartIndex']
                StartTime = UTCDateTime(dataTime[0]*1e-6).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                EndTime = UTCDateTime(dataTime[-1]*1e-6).strftime('%Y-%m-%dT%H:%M:%S.%f+00:00')
                fp_out['Acquisition'].attrs['MeasurementStartTime'] = StartTime
                fp_out['Acquisition']['Raw[0]']['RawData'].attrs['PartStartTime'] = StartTime
                fp_out['Acquisition']['Raw[0]']['RawData'].attrs['PartEndTime'] = EndTime
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['Count'] = len(dataTime)
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['StartIndex'] = \
                        fp_in['Acquisition']['Raw[0]']['RawDataTime'].attrs['StartIndex']
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['StartTime'] = StartTime
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['EndTime'] = EndTime
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['PartStartTime'] = StartTime
                fp_out['Acquisition']['Raw[0]']['RawDataTime'].attrs['PartEndTime'] = EndTime
    return


def process_events(outdir,anames,fdirs,ftemplates,nMin,date):
    '''
    In
        outdir : path to write event data
        anames : names for each array in output
        fdirs : list of paths for raw data from each array
        ftemplates : format string for file naming conventions in each fdir
        nMin : number of minutes to write in output
        date : date (UTCDateTime) to find events
    Out
        ...
    '''
    if not len(fdirs) == len(ftemplates):
        print('Inconsistent fdir and ftemplate inputs')
        raise ValueError
    # Start and end time from 00:00:00 UTC
    t1 = UTCDateTime(date.year,date.month,date.day,0,0,0)
    t2 = t1 + 60*60*24
    # Get list of events
    catalog, ttimes, otimes = get_event_list(t1,t2)
    # Create output directories and write QuakeML
    outpaths = []
    for ii,event in enumerate(catalog):
        outname = otimes[ii].strftime('%Y-%m-%dT%H:%M:%S.%fZ')
        outpath = os.path.join(outdir,outname)
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        try: 
            event.write(os.path.join(outpath,'event.qml'),format='QUAKEML')
        except:
            # USGS catalog output has broken metadata for some events, these are skipped
            os.rmdir(outpath)
            continue
        outpaths.append(outpath)
    # Iterate over arrays
    for fdir,ftemplate,aname in zip(fdirs,ftemplates,anames):
        # Get list of files
        flist, tlist = get_file_list(fdir,ftemplate,t1,t2)
        # Iterate over events
        for outpath,ttime,otime in zip(outpaths,ttimes,otimes):
            if os.path.exists(outpath):
                stime = otime + ttime - 10 # 10 seconds before first arrival
                out = os.path.join(outpath,aname+'.h5')
                write_event_file(out,flist,tlist,stime,nMin)
            else:
                continue
    return outpaths

def assess_data(outdir,anames,fdirs,date,eventpaths):
    '''
    In
        outdir : writing directory for GMT
        anames : list of DAS array names
        fdirs : paths to data
        date : UTCDateTime
        eventpaths : paths to event files
    Out
        ...
    '''
    # Clobber HMS info
    date = UTCDateTime(date.year,date.month,date.day,0,0,0)
    # Open file and write date for title
    fp = open(os.path.join(outdir,'info.txt'),'w')
    fp.write('%s\n' % date.strftime('%B %-d, %Y'))
    # For each array, count files & data volume
    for name,fdir in zip(anames,fdirs):
        flist = [fname for fname in os.listdir(fdir)\
                if fname.startswith(date.strftime('decimator2_%Y-%m-%d'))]
        size = 0
        for fname in flist:
            size += os.path.getsize(os.path.join(fdir,fname))
        size *= 1e-9 # to GB
        nums = len(flist)
        fp.write('%s \t %d \t %d GB\n' % (name,nums,size))
    # Count total events saved
    evnts = len(eventpaths)
    fp.write('%d' % evnts)
    fp.close()
    return

def assess_events(outdir,anames,eventpaths):
    '''
    In
        outdir : writing directory for GMT
        anames : list of DAS arrays
        eventpaths : where the event data is
    Out
        ...
    '''
    # For each event, read QuakeML and write info for GMT
    x = []; y = []; z = []; m = []
    for ii,eventpath in enumerate(eventpaths):
        if not (os.path.exists(os.path.join(eventpath,'%s.h5'%anames[0])) and os.path.exists(os.path.join(eventpath,'%s.h5'%anames[1]))):
            ii-=1
            continue
        try:
            event = op.read_events(os.path.join(eventpath,'event.qml'),format='QUAKEML')[0]
        except:
            ii-=1
            continue
        org = event.origins[0]
        mag = event.magnitudes[0]
        try: 
            evt = event.event_descriptions[0]
        except:
            evt = {'text': 'No description available' }
        x.append(org.longitude); y.append(org.latitude)
        z.append(org.depth*1e-3); m.append(mag.mag)
        with open(os.path.join(outdir,'%d.txt'%ii),'w') as fp:
            fp.write('%s %.2f \t %s \t %.1f km depth \t %s \n' % \
                    (mag.magnitude_type,mag.mag,evt['text'],z[-1],org['time'].strftime('%Y-%m-%d,%H:%M:%S UTC')))
    # Save all events for map plot
    np.savetxt(os.path.join(outdir,'events.xy'),np.column_stack((x,y,z,m)))
    np.savetxt(os.path.join(outdir,'event_list.txt'),np.arange(len(x)),fmt='%s.nc')
    # For each event, create grid file for GMT
    fname1 = anames[0]+'.h5'; fname2 = anames[1]+'.h5'
    for ii,eventpath in enumerate(eventpaths):
        # Interval number of channels to plot (e.g. every 10th channel)
        FACTOR = 10
        # Load data
        if not os.path.exists(os.path.join(eventpath,fname1)):
            ii-=1
            continue
        with h5py.File(os.path.join(eventpath,fname1),'r') as fp:
            data1 = fp['Acquisition']['Raw[0]']['RawData'][:,::FACTOR]
            dx = fp['Acquisition'].attrs['SpatialSamplingInterval']*FACTOR
            fs = fp['Acquisition']['Raw[0]'].attrs['OutputDataRate']
        with h5py.File(os.path.join(eventpath,fname2),'r') as fp:
            data2 = fp['Acquisition']['Raw[0]']['RawData'][:,::FACTOR]
        # Format for plot
        data = np.concatenate((data1[:,::-1],data2),axis=1)
        x = np.arange(-data1.shape[1],data2.shape[1])*dx*1e-3
        t = np.arange(0,data1.shape[0])/fs
        data = data[t<=60,:]
        t = t[t<=60]
        data = data[:,abs(x)<=60]
        x = x[abs(x)<=60]
        # Detrend and filter
        data = detrend(data,axis=0)
        b, a = butter(4,[2,8],btype='bandpass',fs=fs)
        data = filtfilt(b,a,data,axis=0)
        # Normalize for plotting
        data /= 1.5*np.median(abs(data.flatten()))
        # Write to .nc file
        ds = nc.Dataset(os.path.join(outdir,'%d.nc'%ii),'w',format='NETCDF4')
        x_ = ds.createDimension('X',len(x)); xs_ = ds.createVariable('X','f4',('X',))
        y_ = ds.createDimension('Y',len(t)); ys_ = ds.createVariable('Y','f4',('Y',))
        ds_ = ds.createVariable('Data','f4',('Y','X',))
        xs_[:] = x.copy(); ys_[:] = t.copy(); ds_[:,:] = data.copy()
        ds.close()
    return

def generate_report(outdir,anames,fdirs,date,eventpaths):
    # (1) Assess data recorded
    assess_data('/home/efwillia/research/earthquakes/alaska_event_reports/tmp',anames,fdirs,date,eventpaths)

    # (2) Prepare event grids (and metadata for labels)
    assess_events('/home/efwillia/research/earthquakes/alaska_event_reports/tmp',anames,eventpaths)

    # (3) Execute GMT script
    reportname = 'report_%04d-%02d-%02d' % (date.year,date.month,date.day)
    os.system('/bin/bash /home/efwillia/research/earthquakes/alaska_event_reports/make_report.gmt %s' % os.path.join(outdir,reportname))

    # (4) Send GMT plot
    # use 'test_email_list.txt' if debugging
    with open('/home/efwillia/research/earthquakes/alaska_event_reports/email_list.txt','r') as emails:
        TO = ','.join([address.strip() for address in emails.readlines()])
    BODY = '/home/efwillia/research/earthquakes/alaska_event_reports/email.txt'
    ATTCH = '%s.pdf' % os.path.join(outdir,reportname)
    SUBJ = 'Earthquake update for %s' % date.strftime('%B %-d, %Y')
    cmd = 'alpine -passfile ~/.pine-passfile -I ^X,y -subject "%s" -attach "%s" %s < %s' \
            % (SUBJ, ATTCH, TO, BODY)
    os.system(cmd)

    return


