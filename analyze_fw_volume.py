# -*- coding: utf-8 -*-
"""
Created on Wed Jun 27 16:33:48 2018

@author: engelen
"""
import numpy as np
import xarray as xr
from glob import glob
import re
import matplotlib.pyplot as plt
import os
import pandas as pd
import netCDF4 as nc4
import time
import seaborn as sns
import timeit
from imod import idf
from collections import defaultdict
#%%
def parsecode(geomod, brine, steady, *extra):
    sea_codes = {"Homgn" : "H", "Half_Open" : "hO", "Open" : "O", "Closed" : "C"}
    clayer_codes = {"Homgn" : "", "noCL" : "N", "condCL" : "F", "" : "M"}
    brine_codes = {"Bot" : "B", "Top" : "T"}
    mode_codes = {"Paleo" : "P", "Steady" : "S"}
    
    match = re.search(r"(\w+)_((?<=_)[a-zA-Z]+CL)", geomod)
    if match is None:
        shelf = geomod
        if shelf == "Homgn":
            clayer="Homgn"
        else:
            clayer = ""
    else:
        shelf = match.group(1)
        clayer= match.group(2)
    
    codes = [sea_codes[shelf], clayer_codes[clayer], brine_codes[brine], mode_codes[steady]]
        
    #add final code
    fin_code = r"{}-{}-{}-{}".format(*codes)
    if extra is not None:
        fin_code+=extra[0]
    
    codes.append(fin_code)
    return(codes)

def get_origins(ds, conc_names, nan=-999., keep_attrs = True):
    
    conc_nrs = [int(re.search("conc(\d)", i).group(1)) for i in conc_names]

    origins = xr.concat([ds[conc] for conc in conc_names], dim = pd.Index(conc_nrs, name="origin"))
    origins = origins.argmax(dim = "origin") + conc_nrs[0]       
    origins = origins.where(ds[conc_names[0]]!=nan, other=int(nan))
    
    if keep_attrs:
        origins.attrs = ds.attrs
    
    return(origins, conc_nrs)

def natural_sort_key(s, _nsre=re.compile('([0-9]+)')):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)]

def get_volume_conc(cond, volumes):
    cond = cond.assign_coords(z = np.arange(1, cond["z"].shape[0]+1))   #Change z to layer_nrs, most convenient for multiplying with volumes
    return(cond * volumes * 1e-9)                                                 #Multiply by 1e-9 to convert from m3 to km3    
    
def get_volume_fw(da, volumes):
    fw_where = xr.where((da>-0.1) & (da<1.), 1, 0)                                  #Get cells that have a concentration in between 0. and 1.
    return(get_volume_conc(fw_where, volumes))

def origin_volumes(ds, origins, volumes, conc_nrs):
    for conc_nr in conc_nrs:    
        ds["origin_%d" % conc_nr] = get_volume_conc(xr.where(origins==conc_nr, 1, 0), volumes).sum(dim=["x","y","z"])
    return(ds)

#%%Options for script
skip_or=False
skip_fw=False

steady=False

if steady == True:
    path_col = "path_s"
    pcode_steady = "Steady"
else:
    path_col = "path_h"
    pcode_steady = "Paleo"

#Conductances
cond_sa=True

#%%Path management
coast_path = r"./data/0BP_coastline.IDF"
df_paths = pd.read_csv(r"./path_management/paths2paleo.csv")
df_paths_all = pd.read_csv(r"./path_management/paths2results.csv")
outputfolder = r"./plots/"

if cond_sa==True:
    df_paths = pd.read_csv(r"./path_management/conductance/paths2paleo.csv")
    df_paths_all = pd.read_csv(r"./path_management/conductance/paths2paleo.csv")
    outputfolder=os.path.join(outputfolder, "conductance")

fw_folder = os.path.join(outputfolder, "FW_volumes")
or_folder = os.path.join(outputfolder, "origins")


##Select specific run
#df_paths_all = df_paths_all.loc[(df_paths_all["geomod"] == "Half_Open") & (df_paths_all["brine"] == "Bot")]
#df_paths = df_paths.loc[(df_paths["geomod"] == "Half_Open") & (df_paths["brine"] == "Bot")]

#%%
source_dict={"origin_3" : "Sea", "origin_4" : "River", "origin_5" : "Dunes", "origin_6" : "Brine"}

#%%Model parameters
##Coast
coast = idf.dataarray(coast_path)
coast = coast.fillna(0).astype(np.int32)
onshore = (coast - 1) * -1 #Goofy way of inverting binary matrix, ~ did not seem to work. 

#Model discretization
dx = 1000
dy = 1000
dz = np.array([20.]*21 + [40.]*10 + [50.]*4)

#n = 0.10 #porosity
n = 0.25

volumes = xr.DataArray(dz, coords = [np.arange(1, len(dz)+1)], dims = ["z"]) * dx * dy * n

drop_vars = ["vx", "vy", "vlayer", "constant head", "dcdt", "head", "head dep bounds", "river leakage"]

#Open options
open_opt = dict(decode_times=False, drop_variables = drop_vars, chunks = {"time" : 10, "z" : 35, "y": 258, "x" : 237})

#%%Analyze end state FW
fw = defaultdict(list)

if not skip_fw:
    for i, row in df_paths_all.iterrows():
        if cond_sa == True:
            extra=row["path_a"]
        else:
            exta=None
        
        sea, clayer, _, _, code = parsecode(row["geomod"], row["brine"], row["steady"], extra)
        
        nc_steady = glob(os.path.join(row[path_col], r"results_[0-9][0-9][0-9].nc"))
        nc_steady.sort(key=natural_sort_key)
    
        ds_steady = xr.open_dataset(nc_steady[-1], **open_opt).isel(time=-1)      
        
        fw["total_fw"].append(get_volume_fw(ds_steady["conc1"], volumes).sum(dim=["x","y","z"]).values)
        fw["total_fw_pumpable"].append(get_volume_fw(ds_steady["conc1"].sel(z=slice(None, -300)), volumes).sum(dim=["x","y","z"]).values)
        fw["total_onshore"].append((get_volume_fw(ds_steady["conc1"].sel(z=slice(None, -300)), volumes)*onshore).sum(dim=["x","y","z"]).values)
        fw["codes"].append(code)
        fw["sea"].append(sea)
        fw["clayer"].append(clayer)

    fw = pd.DataFrame(fw, columns=["codes", "sea", "clayer", "total_fw", "total_fw_pumpable", "total_onshore"])
    fw = pd.concat([df_paths_all.loc[:, ["geomod", "brine", "steady"]], fw],axis=1)
    fw.to_csv(os.path.join(fw_folder, "end_volumes.csv"), index=False)

#%%Analyze transience paleo models
if not skip_or:
    corrected_models=[]
#    for i, row in df_paths.loc[(df_paths["geomod"]=="Half_Open") & (df_paths["brine"]=="Top"), :].iterrows():
    for i, row in df_paths.iterrows():
        #skip first one, because we already have this one
        if cond_sa == True:
            extra=row["path_s"]
        else:
            extra=None
        
        code = parsecode(row["geomod"], row["brine"], pcode_steady, extra)[-1]
        
        nc_paths = glob(os.path.join(row[path_col], r"results_[0-9][0-9][0-9].nc"))
        nc_paths.sort(key=natural_sort_key)
    
        ds_last = xr.open_dataset(nc_paths[-1], **open_opt)
    
        ds = xr.concat([xr.open_dataset(nc_path,**open_opt) for nc_path in nc_paths][:-1], dim="time")
            
        #Something went wrong with origins as 1 of the origins concentrations dropped for the last 30 years (seperate run), so concat did not work. 
        missing_vars = set(ds.data_vars.keys()) - set(ds_last.data_vars.keys())
        if len(missing_vars) > 0:
            corrected_models.append(code)
            for var in missing_vars:
                #Hack to get empty array with NANs as there is no nice function for that yet
                ds_last[var] = ds_last[list(ds_last.data_vars.keys())[0]] 
                ds_last[var] = ds_last[var].where(np.nan)
        
        ds = xr.concat([ds, ds_last], dim="time")
        
        ds.attrs["title"] = code
        ds["total_fw"] = get_volume_fw(ds["conc1"], volumes).sum(dim=["x","y","z"])
        ds["total_fw_pumpable"] = get_volume_fw(ds["conc1"].sel(z=slice(None, -300)), volumes).sum(dim=["x","y","z"])
        ds["total_onshore"] = (get_volume_fw(ds["conc1"].sel(z=slice(None, -300)), volumes)*onshore).sum(dim=["x","y","z"])
        
        ds["time"]=ds["time"]/365.25/24
    
    
        #%%Export lines
        path_fw = os.path.join(fw_folder, code+"_fw.nc")
        ds[["total_fw", "total_fw_pumpable", "total_onshore"]].to_netcdf(path_fw, unlimited_dims = ["time"])
    
        #%%
        lines = xr.plot.line(ds["total_fw"], hue="section")
        lines2 = xr.plot.line(ds["total_fw_pumpable"], hue="section")
        lines3 = xr.plot.line(ds["total_onshore"], hue="section")
        
        plt.ylabel("Total fresh water volume (km3)")
        plt.xlabel("Time (years)")
        #plt.ylim(ymin=0, ymax = 5000)
        plt.ylim(ymin=0)
        #ax.legend()
        plt.savefig(os.path.join(fw_folder, code+"_fw.png"))
        plt.close()
    
        #%%
        origin_names = [conc_nr for conc_nr in ds.keys() if re.search("conc[3-9]", conc_nr)]
        if row["brine"] == "Bot":
            origin_names[-2], origin_names[-1] = origin_names[-1], origin_names[-2] 
        
        origins, origin_nrs = get_origins(ds, origin_names)
        attrs = origins.attrs
        
        if row["brine"] == "Top":
            #Correct for error made in sources where 
            origins = origins.where(ds["conc1"] < 37., other=6)
        ds = origin_volumes(ds, origins, volumes, origin_nrs)
        
        df = ds[["origin_%d" % or_nr for or_nr in origin_nrs]].to_dataframe()
        df = df.rename(index = str, columns = source_dict)
        
        df.plot.line(colormap = "viridis")
        df.to_csv(os.path.join(or_folder, code+"_origins.csv"))
        
        plt.ylabel("Total water volume (km3)")
        plt.xlabel("Time (years)")
        plt.savefig(os.path.join(or_folder, code+"_origins.png"))
        plt.close()
    
        #%%
        path_out = os.path.join(or_folder, code+"_origins.nc")
        
        origins["time"] = origins["time"]*365.25*24
        origins = origins.astype(np.float32)
        origins = origins.to_dataset(name = "origins")
        origins = origins.transpose("time", "z", "y", "x")
        
        origins.attrs = attrs
        #origins.time.attrs["units"] = attrs["units"]
        
        origins.to_netcdf(path_out, unlimited_dims = ["time"])
        
        time.sleep(30)
        
        testnc = nc4.Dataset(path_out, mode= 'r')
        testnc.close()
        
        testnc = nc4.Dataset(path_out, mode= 'a')
        testnc.variables["time"].setncattr("units", ds.attrs["units"])
        print(testnc.variables["time"])
        testnc.close()