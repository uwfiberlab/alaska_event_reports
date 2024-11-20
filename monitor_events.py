#!/home/efwillia/miniconda3/bin/python

import numpy as np
from obspy import UTCDateTime
from event_tools import process_events, generate_report

# Path to write event data and reports
#outdir1 = '/mnt/disk1/event_data'
outdir1 = '/mnt/disk1/event_data_v2'
outdir2 = '/mnt/disk1/event_reports'

# Actively recording arrays, directories, and file formats
#anames = ['KKFLS','TERRA']
#fdirs = ['/mnt/qnap/KKFL-S_FIberA_25Hz','/mnt/qnap/TERRA_FiberA_25Hz']
anames = ['TERRA','KKFLN']
fdirs = ['/mnt/qnap/TERRA_FiberA_50Hz','/mnt/qnap/KKFL-N_FiberA_50Hz']
ftemplates = ['decimator2_%Y-%m-%d_%H.%M.%S_UTC.h5','decimator2_%Y-%m-%d_%H.%M.%S_UTC.h5']

# Number of minutes to save in each event file
nMin = 2

# Process events from previous day (UTC)
today = UTCDateTime.now()
date = today - 60*60*24
eventpaths = process_events(outdir1,anames,fdirs,ftemplates,nMin,date)

# Generate report
generate_report(outdir2,anames,fdirs,date,eventpaths)


