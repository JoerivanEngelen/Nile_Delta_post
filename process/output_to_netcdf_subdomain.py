# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 17:18:53 2018

@author: engelen
"""


from glob import glob
import xarray as xr
from timeit import default_timer as timer
import numpy as np
import os, sys
import netCDF4 as nc4

sys.path.append(os.path.abspath(os.path.join(os.path.abspath(__file__),'..', ".." ,"imodpy", "imod-python")))

from imod import idf, util

#%%
nan_dic = {"conc1" : -1.,
           "conc2" : -1.,
           "conc3" : -1.,
           "conc4" : -1.,
           "conc"  : -1.,
       "constant head" : 0.,
       "dcdt" : 0.,
       "flow front face" : 0.,
       "flow lower face" : 0.,
       "flow right face" : 0.,
       "head" : -9999.,
       "head dep bounds" : 0.,
       "river leakage" : 0.}

#%%
start = timer()

util.Config.np_datetime = False

path_tot = r"g:\3D_Nile_Delta\jengelen.4320447\results_2\head_*_p004*"
path_out = r"g:\3D_Nile_Delta\jengelen.4320447\results_002_p{:03d}.nc"
path_time = r"g:\3D_Nile_Delta\jengelen.4320447\init_times.txt"
res_nr = 2

#path_tot = sys.argv[1]
#path_out = sys.argv[2]
#path_time = sys.argv[3]
#res_nr = sys.argv[4]

#init_time = int(np.loadtxt(path_time)[int(res_nr)-1])

paths = glob(path_tot)

#%%
ds_tot = idf.loadset(path_tot, memmap = False)

nan_dic = {key: nan_dic[key] for key, var in nan_dic.items() if key in ds_tot}

for key in ds_tot.keys():
    if "subdomain" in ds_tot[key].coords:
        subd = ds_tot[key].coords["subdomain"]
        ds_tot[key] = ds_tot[key].drop("subdomain")

path_out = path_out.format(subd.values)

arb_var = key
attrs = ds_tot[arb_var].attrs

ds_tot = xr.Dataset(ds_tot)

ds_tot.attrs = attrs
ds_tot = ds_tot.transpose("time", "layer", "y", "x")
ds_tot = ds_tot.fillna(nan_dic)


ds_tot.time.encoding["units"] = attrs["units"]


ds_tot.to_netcdf(path_out, unlimited_dims = ["time"])

sect4 = timer()
print("Done with Section 4")
print(sect4-start)

#%%
#Explicitly set NetCDF units
testnc = nc4.Dataset(path_out, mode= 'a')
testnc.variables["time"].setncattr("units", attrs["units"])
print(testnc.variables["time"])
testnc.close()

end = timer()
print("Done with Section 5")
print(end-sect4)