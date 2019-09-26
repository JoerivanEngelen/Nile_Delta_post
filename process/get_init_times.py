# -*- coding: utf-8 -*-
"""
Created on Thu May 24 13:30:36 2018

@author: engelen
"""

from glob import glob
import numpy as np
import os, sys
import netCDF4 as nc4
from datetime import datetime
import re


#%%

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)] 

def max_date(folder):
    files = glob(os.path.join(folder, "head_*_l1_p001.idf"))
    dates = [int(re.search(r"\d{14}", i).group(0)) for i in files]
    return(max(dates))

def wildcard_to_nr(folder, folder_glob):
    return(int(folder.split(sep = folder_glob.split(sep = "*")[0])[1]))

#%%

##Example arguments
#folder_glob = r"/home/path/results_*"
#inittxt = r"/home/path/init_times_test.txt"
#cont_nr = 1

folder_glob = sys.argv[1]
inittxt = sys.argv[2]

if len(sys.argv)>3:
    cont_nr = int(sys.argv[3])
else:
    cont_nr = None

folders = [f for f in glob(folder_glob) if os.path.isdir(f)]
folders.sort(key=natural_sort_key)

folders_cont = []

if cont_nr is not None:
    init_times = list(np.loadtxt(inittxt, ndmin=1))[:cont_nr]
    for fol in folders:
        nr = wildcard_to_nr(fol, folder_glob)
        if nr >= cont_nr:
           folders_cont.append(fol)
else:
    init_times = [0.]
    folders_cont = folders

for i, f in enumerate(folders_cont):
    try:
        dt = datetime.strptime(str(max_date(f)), '%Y%m%d%H%M%S')
        units = "hours since 2000-01-01 00:00:00.0"
        date = nc4.date2num(dt, units = units, calendar = "gregorian")
        init_time = init_times[-1] + date
        init_times.append(init_time)
    except ValueError:
        pass

np.savetxt(inittxt, np.array(init_times))