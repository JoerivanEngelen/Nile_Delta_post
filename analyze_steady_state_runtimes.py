# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 11:18:44 2019

@author: engelen
"""

import xarray as xr
import pandas as pd
import os
from glob import glob
import re
import numpy as np

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]


#%%
##Model discretization
#dx = 1000
#dy = 1000
#dz = np.array([20.]*21 + [40.]*10 + [50.]*4)
#
##n = 0.10 #porosity
#n = 0.25
#
#volumes = xr.DataArray(dz, coords = [np.arange(1, len(dz)+1)], dims = ["z"]) * dx * dy * n

drop_vars = ["vx", "vy", "vlayer", "constant head", "dcdt", "head", "head dep bounds", "river leakage"]

#Open options
open_opt = dict(decode_times=False, drop_variables = drop_vars, chunks = {"time" : 10, "z" : 35, "y": 258, "x" : 237})

#%%
#def calc_fresh_water

df_paths = pd.read_csv(r"./path_management/conductance/paths2paleo.csv")

steady_paths = df_paths.loc[df_paths["brine"]=="Bot","path_s"]
times=[]
ds_list = []
for i, path in enumerate(steady_paths):
    nc_paths = glob(os.path.join(path, r"results_[0-9][0-9][0-9].nc"))
    nc_paths.sort(key=natural_sort_key)
#    ds_list.append(xr.open_mfdataset(nc_paths, **open_opt))
    times.append(xr.open_dataset(nc_paths[-1], decode_times=False).isel(time=-1).time/365.25/24)
    
    