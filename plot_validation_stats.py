# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 10:21:31 2019

@author: engelen
"""

import geopandas as gpd
import pandas as pd
from glob import glob
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from scipy import stats
import re

def plot_euclids(name_model, nr_nans=None, **data):
    cell_euclid_fig = "Euclid_cells_%s.png"
    for select_type, df in data.items():
        folder="euclid_%s" % select_type
        ax = sns.jointplot(df["TDS (g/l)"], df['cell_dist_'], kind="scatter", alpha = 0.2)
        ax = ax.set_axis_labels("TDS (g/l)", "$\Lambda$ (-)")
        if nr_nans is not None:
            nr_nans_annotate = lambda a, b: nr_nans
            ax = ax.annotate(nr_nans_annotate, template="{stat} = {val}",stat="$N_{nans}$", loc="upper left", fontsize=12) 
        
        plt.ylim(0.0, 0.25)
        plt.tight_layout()
        plt.savefig(os.path.join(outputfolder, folder, cell_euclid_fig % (name_model)))
        plt.close()

def data2dic(name, col = ["name", "code", "sea", "clayer", "brine", "mode", "type", "$\Lambda$ (-)"], **data):
    """Add all data to dic, to quickly initiate a pandas dataframe
    """

    ls = []
    for select_type, arr in data.items():
        df = pd.DataFrame(columns = col)
        df[col[-1]] = arr #first add array so that the name columns and type columns get values for every point
        df[col[-2]] = select_type
        sea, clayer, brine, mode, full_code = parsecode(name)
        df[col[5]] = mode
        df[col[4]] = brine
        df[col[3]] = clayer
        df[col[2]] = sea
        df[col[1]] = full_code
        df[col[0]] = name
        ls.append(df)
    return(ls)

def stats2dic(name, dist_col='cell_dist_', modcol="TDS mod (g/l)", obscol="TDS (g/l)",  **data):
    """Add stats to dic, to quickly initiate a pandas dataframe
    """
    ls = []
    for select_type, arr in data.items():  
        d = {}
        d["model"] = name
        sea, clayer, brine, mode, full_code = parsecode(name)
        d["sea"] = sea
        d["clayer"] = clayer
        d["brine"] = brine
        d["mode"] = mode
        d["code"] = full_code
        d["type"] = select_type
        #cannot do this in list comprehension because they decided it would be more convenient to provide min max in seperate tuple
        d["n"], (d["smin"], d["smax"]), d["sm"], d["sv"], d["sv"], d["sk"] = stats.describe(arr[dist_col])
        d["smed"] = np.median(arr[dist_col])
        d["mae"] = mae(arr[modcol], arr[obscol])
        ls.append(d)
    return(ls)

def parsecode(name):

    sea_codes = {"Homgn" : "O", "Half_Open" : "H", "Open" : "O", "Closed" : "C"}
    clayer_codes = {"Homgn" : "N","noCL" : "N", "condCL" : "F", "" : "M"}
    brine_codes = {"Bot" : "B", "Top" : "T"}
    mode_codes = {"Paleo" : "P", "Steady" : "S"}    
    
    match = re.search(r"(Holocene_(.+?)_Brine_(.+?))_(.+)", name)
    geomod = match.group(2)
    brine  = match.group(3)
    mode = match.group(4)
    
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
    
    codes = [sea_codes[shelf], clayer_codes[clayer], brine_codes[brine], mode_codes[mode]]
    #add final code
    codes.append(r"{}-{}-{}-{}".format(*codes))
    return(codes)

def mae(obs, mod):
    """Calculate mean absolute error
    """
    return(np.mean(np.abs(obs-mod)))

#%%Plotting functions
def plot_barplot(data, file, figsize, legend_loc = "upper right", xlim = None, xcol = "mae"):
    colors = np.array([(0.171875, 0.48046875, 0.7109375),   (0.5703125, 0.7734375, 0.87109375), 
                       (0.99609375, 0.87109375, 0.6015625), (0.83984375, 0.09765625, 0.109375)])
    
    f, ax = plt.subplots(figsize=figsize)
    
    unique_hue = data.loc[:, "type"].unique()
    hue_order = [x for x in ["fresh", "brackish", "saline", "hypersaline"] if x in unique_hue]
    
    sns.barplot(y="code", x= xcol, hue="type", hue_order=hue_order, data=data, palette=colors, orient = "h")

    plt.tight_layout()
    plt.legend(loc=legend_loc)
    plt.savefig(file)
    plt.close()   

def plot_euclid_boxplot(data, file, figsize, legend_loc = "upper right", xlim = None, xcol = "$\Lambda$ (-)"):
    colors = np.array([(0.171875, 0.48046875, 0.7109375),   (0.5703125, 0.7734375, 0.87109375), 
                       (0.99609375, 0.87109375, 0.6015625), (0.83984375, 0.09765625, 0.109375)])
    
    f, ax = plt.subplots(figsize=figsize)
    
    unique_hue = data.loc[:, "type"].unique()
    hue_order = [x for x in ["fresh", "brackish", "saline", "hypersaline"] if x in unique_hue]
    
    sns.boxplot(y="code", x= xcol, hue="type", hue_order=hue_order, data=data, palette=colors)

    ax = correct_colors_boxplot(ax, line_nrs=[0, 1, 2, 3, 5], remove_edge=True)
  
    if xlim is not None:
        plt.xlim(xlim)
    
    sns.despine(trim=True, left=True)
    plt.tight_layout()
    plt.legend(loc=legend_loc)
    plt.savefig(file)
    plt.close()    

def correct_colors_boxplot(ax, line_nrs = range(1,5), remove_edge=False):
    for i,artist in enumerate(ax.artists):
        # Set the linecolor on the artist to the facecolor
        col = artist.get_facecolor()
        if remove_edge:
            artist.set_edgecolor(col)
        artist.set_facecolor(col)
    
        # Each box has 6 associated Line2D objects (to make the whiskers, fliers, etc.)
        # Loop over them here, and use the same colour as above
        #(+0=minimum vertical line, +1=maximum vertical line, +2=minimum horizontal line, 
        #+3=maximum horizontal line, +4=median line, +5=flier)
#        for k in [0, 1, 2, 3, 5]:
#        for k in [5]:
        for k in line_nrs:
            j = i * 6 + k
            line = ax.lines[j]
            line.set_color(col)
            line.set_mfc(col)
            line.set_mec(col)
    return(ax)


def plot_euclid_boxplot_v2(data, file, figsize, legend_loc = "upper right", xlim = None, xcol = "$\Lambda$ (-)"):
    """Seperate groups.
    """
    colors = np.array([(0.171875, 0.48046875, 0.7109375),   (0.5703125, 0.7734375, 0.87109375), 
                       (0.99609375, 0.87109375, 0.6015625), (0.83984375, 0.09765625, 0.109375)])
    
    f, axes = plt.subplots(2,2, figsize=figsize)
    axes = axes.flatten()
    
    unique_hue = data.loc[:, "type"].unique()
    hue_order = [x for x in ["fresh", "brackish", "saline", "hypersaline"] if x in unique_hue]
    
    legend_elements = [Patch(facecolor=colors[i], edgecolor='r',label=hue) for i, hue in enumerate(hue_order)]
    
    for i, hue in enumerate(hue_order):
        ds = data.loc[data["type"]==hue, :]
        code = ds["code"].unique()[0]
        n = len(ds[ds["code"] == code])
        sns.boxplot(ax = axes[i], x="code", y= xcol, data=ds, color = colors[i])
        axes[i] = correct_colors_boxplot(axes[i], line_nrs=[5], remove_edge=False)
        
        if i == 1:
            axes[i].legend(handles=legend_elements, loc=legend_loc)
            
        axes[i].set_xticklabels(axes[i].get_xticklabels(), rotation=45, horizontalalignment="right")
        axes[i].text(0.05, 0.95, "n = {}".format(n), transform=axes[i].transAxes)
        if xlim is not None:
            axes[i].set_ylim(xlim[i])

#        axes[i].axhline(y=0, color = ".80")
    sns.despine(trim=True, bottom=True)
    plt.tight_layout()
    
    name, ext = os.path.splitext(file)
    
    plt.savefig(file, dpi=300)
    plt.savefig(name+".svg")
    plt.savefig(name+".pdf")
    plt.close()

        
#%%Path management
outputfolder = r"./plots/validation/"
TDS_folder = r"./plots/validation/shpfiles/"
br_zone = r"./plots/validation/shpfiles/brackish_zone2.shp"

shp_paths = glob(TDS_folder + r"*.shp")
br_zone = gpd.read_file(br_zone)

obs_model_fig = "Obs_Model_%s.png" 
euclid_fig = "Euclid_%s_%s.png"

#%%Initiate some stuff for data frame
col_stats = ["model", "sea", "clayer", "brine", "mode", "code","type", "n", "smin", "smax", "sm", "sv", "ss", "sk", "smed", "mae"]
ls_stats = []

col_tot = ["name", "code", "sea", "clayer", "brine", "mode", "type", "$\Lambda$"]
ls_tot = []
ls_mae = []
#%%Loop over shapes
for path in shp_paths:
#    name_model = os.path.splitext(os.path.basename(path))[0]
    name_model = os.path.basename(path).split("_long.shp")[0]
    
    obs_long = gpd.read_file(path)
    obs_long["TDS (g/l)"] = obs_long["TDS (mg/l)"]/1000.
    
    obs = obs_long[obs_long["obs"]==1]
    obs["TDS mod (g/l)"] = obs_long.loc[obs_long["obs"]==0, "TDS (g/l)"].reset_index(drop=True)
    
    #%%Select data
    nans = obs['cell_dist_']>998.
    
    _, _, brine, mode, _ = parsecode(name_model)
    
    if brine == "T" and mode == "S":
        obs.loc[obs["TDS (g/l)"] > 35.5, 'cell_dist_'] = np.nan #Ensure cell_dist is nan for hypersaline water under Top Steady.

    nr_nans=np.count_nonzero(nans)
    obs.loc[nans, 'cell_dist_'] = np.nan
    
    brackish    = obs[(obs["TDS (g/l)"] > 1.)   & (obs["TDS (g/l)"] <= 5.)   & ~obs.intersects(br_zone.iloc[2]["geometry"])]
    saline      = obs[(obs["TDS (g/l)"] > 5.)   & (obs["TDS (g/l)"] <= 35.)]
    hypersaline = obs[(obs["Depth (m)"] > 300.) & (obs["TDS (g/l)"] > 35.5)]
    fresh       = obs[obs["TDS (g/l)"] < 1.]
        
    #%%Create data tables
    ls_stats.extend(stats2dic(name_model,
                                    brackish=brackish,
                                    saline=saline,
                                    hypersaline=hypersaline, 
                                    fresh=fresh))

    ls_tot.extend(data2dic(name_model, col = col_tot,
                                    brackish=np.array(brackish['cell_dist_']),
                                    saline=np.array(saline['cell_dist_']),
                                    hypersaline=np.array(hypersaline['cell_dist_']), 
                                    fresh=np.array(fresh['cell_dist_'])))
    


#%%Create dataframes
df_summ = pd.DataFrame(ls_stats, columns=col_stats)
df_tot = pd.concat(ls_tot)

#%%Create boxplot
df_tot = df_tot.rename(columns={"$\Lambda$":"$\Lambda$ (-)"})

gp = ['sea',"clayer", "brine", "type"]

df_p = df_tot.loc[df_tot["mode"]=="P"].set_index(gp)
df_s = df_tot.loc[df_tot["mode"]=="S"].set_index(gp)

df_diff = df_p.copy()
df_diff["$\Lambda$ (-)"]  = df_s["$\Lambda$ (-)"] - df_p["$\Lambda$ (-)"]
diff_colname="$\Delta \; \Lambda$ (-)"
df_diff = df_diff.rename(columns={"$\Lambda$ (-)":diff_colname})
df_diff["code"] = df_diff.loc[:, "code"].replace(to_replace={"-P" : ""}, regex=True)

df_p, df_s, df_diff = [i.reset_index(gp) for i in [df_p, df_s, df_diff]]

plot_euclid_boxplot_v2(df_p, os.path.join(outputfolder, "boxplot_euclid_P_v2.png"), (9, 9), xlim=[(0., 0.25)]*4)
plot_euclid_boxplot_v2(df_diff, os.path.join(outputfolder, "boxplot_euclid_diff_v2.png"), (9,9),
                       xlim=[(-0.105, 0.105)]*3+[(-0.205, 0.205)], xcol = diff_colname)

#plot_barplot(df_summ.loc[(df_summ["mode"]=="P") & (df_summ["type"] != "hypersaline")], os.path.join(outputfolder, "MAE.png"), (6, 12))


#%%
def loc_tup(mode, salt=None, sea=None, clayer=None, brine=None):
    if mode not in ["P", "S"]:
        raise ValueError("")
    return((slice(sea), slice(clayer), slice(brine), slice(salt), mode))

def salt_tup(salt, sea=None, clayer=None, brine=None):
    if salt not in ["brackish", "saline", "fresh", "hypersaline"]:
        raise ValueError("")
    return((slice(sea), slice(clayer), slice(brine), salt))

df_summ = df_summ.sort_values(by=['model', 'sea', 'clayer', "brine"])
mean = df_summ.groupby(['sea',"clayer", "brine", "type", "mode"])["sm"].mean()
median = df_summ.groupby(['sea',"clayer", "brine", "type", "mode"])["smed"].mean()

diff_mean   =   mean.loc[loc_tup("P")] - mean.loc[loc_tup("S")]
diff_median = median.loc[loc_tup("P")] - mean.loc[loc_tup("S")]

print(diff_mean.loc[salt_tup("saline")])
print(diff_mean.loc[salt_tup("brackish")])
print(diff_median.loc[salt_tup("brackish")])
