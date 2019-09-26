# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 11:26:25 2019

@author: engelen
"""

import xarray as xr
import pandas as pd
from glob import glob
import os
import numpy as np
from collections import defaultdict

def get_first_true(df, condition):
    time = df[condition].iloc[0:1].index.values
    if time.size == 0:
        time = df.iloc[-2:-1].index.values
    return(time)

#%%Path management
fw_path = r"./plots/FW_volumes/*-S_fw.nc"
fw_paths = glob(fw_path)

or_path = r"./plots/FW_volumes/*-S_origins.csv"
or_paths = glob(or_path)
#%%Read fresh water volumes
d_fw = {}

open_opt = dict(decode_times=False, 
                drop_variables = ["total_fw_pumpable", "total_onshore"])

for p in fw_paths:
    name = os.path.basename(p).split("_fw.nc")[0]
    d_fw[name] = xr.open_dataset(p, **open_opt)

#%%Differentiate

for name, ds in d_fw.items():
    ds["fw_norm_diff"] = (
            ds["total_fw"]/ds["total_fw"].max()
#            ds["total_fw"]/8734.5725
    ).isel(time=slice(None, -7)).differentiate("time")

#%%time to reach steady state fw_vol

diff = xr.merge(
        [ds["fw_norm_diff"].rename(name) for name, ds in d_fw.items()]
        ).drop(["dx", "dy"]).to_dataframe()
diff = np.log10(np.abs(diff))

time_steady={}

for name in diff.columns:
    time_steady[name]=get_first_true(diff[name], diff[name] < -6)
    
#%%Read origins
colnames = []    

d_or = defaultdict()

for csv in or_paths:
    name = os.path.basename(csv).split("_origins.csv")[0]
    d_or[name] = pd.read_csv(csv, header=0).set_index("time").drop(columns=["dx", "dy"])
    colnames.extend([(name, var) for var in d_or[name].columns])

d_or = pd.concat(d_or, axis=1)

#%%Differentiate
#Use xarray to differentiate, as it automatically differentiates properly
tot_vol = d_or.loc[:, ("C-F-B-S", slice(None))].sum(axis=1).iloc[0]

diff_or = xr.Dataset(d_or/tot_vol).differentiate("time").to_dataframe()
diff_or = np.log10(np.abs(diff_or))

time_steady_or={}

for name in diff_or.columns:
    time_steady_or[name]=get_first_true(diff_or[name], diff_or[name] < -6.25)

#All this stacking, reseting and dropping is to get rid the table in the right format
time_steady_or=pd.DataFrame(time_steady_or).stack().reset_index(level=[0]).drop(columns="level_0")

mx_time_steady_or = time_steady_or[time_steady_or.index=="River"].max(axis=0)
mx_time_steady_or.to_csv(os.path.join(or_path, "..", "time_to_steady.csv"))

#%%