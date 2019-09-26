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
def get_q_avg(ds, flow_var, dim, zero_all = None):
    
    if zero_all is None:
    #Find whether out begin or end of dimension equals 0
        zero_all =  [bool((ds[flow_var][{dim : 0}] == 0.).all().values), 
                     bool((ds[flow_var][{dim : -1}] == 0.).all().values)]
        
    if zero_all[0] == zero_all[1]:
        print(flow_var)
        print(dim)
        raise ValueError("zero_all evaluate to the same value: " + str(zero_all))
    
    #All last elements equals 0., preseverve first element
    if zero_all[1]:
        n = 1
        addslice = slice(n, None)
        keepslice = slice(None, n)
    
    #All first elements equals 0., preseverve last element
    elif zero_all[0]:
        n = -1
        addslice = slice(None, n)
        keepslice = slice(n, None)
    
    added = ds[flow_var] + ds[flow_var].shift(**{dim : n})
    added = added[{dim : addslice}]
    
    return(xr.concat([ds[flow_var][{dim : keepslice}], added], dim = dim) * 0.5)
    
def get_vel(ds, flow_var, dim, A, zero_all = None):
    return(get_q_avg(ds, flow_var, dim, zero_all)/A)


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

flow_dic = {"x" : "flow right face", 
            "y" : "flow front face", 
            "layer" : "flow lower face"}

#%%
start = timer()

util.Config.np_datetime = False

##Examples
#path_tot = r"/home/path/head_*"
#path_out = r"/home/path/results_006.nc"
#path_time = r"/home/path/init_times.txt"
#res_nr = 1

path_tot = sys.argv[1]
path_out = sys.argv[2]
path_time = sys.argv[3]
res_nr = sys.argv[4]

init_time = int(np.loadtxt(path_time)[int(res_nr)-1])


paths = glob(path_tot)

#%%
ds_tot = idf.loadset(path_tot, memmap = False)

nan_dic = {key: nan_dic[key] for key, var in nan_dic.items() if key in ds_tot}

sect1 = timer()
print("Done with Section 1")
print(sect1 - start)

#%%
for key in ds_tot.keys():
    if "subdomain" in ds_tot[key].coords:
        ds_tot[key] = ds_tot[key].drop("subdomain")

#ds_tot = {key: ds_tot[key].drop("subdomain") if "subdomain" in ds_tot[key] else ds_tot[key] for key in ds_tot.keys()}

#get attributes of arbitrary variable, otherwise use arb_var = list(ds_tot.keys())[0]
arb_var = key
attrs = ds_tot[arb_var].attrs

ds_tot = xr.Dataset(ds_tot)

ds_tot.attrs = attrs
ds_tot = ds_tot.transpose("time", "layer", "y", "x")

ds_tot.time.encoding["units"] = attrs["units"]


ds_tot.fillna(nan_dic)

sect2 = timer()
print("Done with Section 2")
print(sect2-sect1)



#%%
dx, dy = ds_tot[arb_var].attrs["res"]

reg_top = 20.
thicknesses = np.array([20.]*21 + [40.]*10 + [50.]*4)

indices = ds_tot.layer.values - 1 #-1 because Python indexing and not Fortran indexing

ints = np.insert((reg_top - np.cumsum(thicknesses)), 0, reg_top)
z = ((ints[:-1] + ints[1:])/2)[indices]

A_dic = {"x" : xr.DataArray(thicknesses[indices] * dy, coords = [ds_tot.coords["layer"]], dims = ["layer"]),
         "y" : xr.DataArray(thicknesses[indices] * dx, coords = [ds_tot.coords["layer"]], dims = ["layer"]),
         "layer" : dx*dy}

if "flow front face" in nan_dic: #Check whether fluxes are saved.
    try:
        ds_tot = xr.merge([{"v"+key : get_vel(ds_tot, var, key, A_dic[key]) for key, var in flow_dic.items()}, ds_tot])
        for key, var in flow_dic.items():
            ds_tot = ds_tot.drop(var)
    except ValueError:
        print("Warning could not incorporate velocities")

ds_tot.coords["z"] = ("layer", z)
ds_tot.swap_dims({"layer" : "z"}, inplace = True)

ds_tot = ds_tot.chunk({"time" : 1, "z" : 35, "x" : 237, "y" : 258})
#Offset time with actual initial time
ds_tot["time"] = ds_tot["time"] + init_time

sect3 = timer()
print("Done with Section 3")
print(sect3-sect2)

#%%

ds_tot.to_netcdf(path_out, unlimited_dims = ["time"])

sect4 = timer()
print("Done with Section 4")
print(sect4-sect3)

#%%

#Explicitly set NetCDF units
testnc = nc4.Dataset(path_out, mode= 'a')
testnc.variables["time"].setncattr("units", attrs["units"])
print(testnc.variables["time"])
testnc.close()

end = timer()
print("Done with Section 5")
print(end-sect4)