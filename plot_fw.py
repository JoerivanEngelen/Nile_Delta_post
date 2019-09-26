# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 13:23:35 2019

@author: engelen
"""

import xarray as xr
import seaborn as sns
import os
from glob import glob
import pandas as pd
import matplotlib.pyplot as plt
import SeabornFig2Grid as sfg
import matplotlib.gridspec as gridspec
import numpy as np

#%%Path management
resfolder = r"./plots/FW_volumes/"

validated = ["C-M-B-P", "C-N-T-P","H-M-T-P", "H-F-T-P", "H-N-T-P"]
files = glob(os.path.join(resfolder, r"*.nc"))

opt=dict(axis = 0, sort = False)



#%%Data handling 
#Changes over time
codes = [os.path.basename(f).split("_fw.nc")[0] for f in files]
dfs = [xr.open_dataset(f).to_dataframe() for f in files]

columns = dfs[0].columns

data = pd.concat([pd.concat([df[col] for df in dfs], keys=codes, names = ["model", "time"], 
                            **opt).reset_index().rename(columns = {col : "V (km3)"}) for col in columns], 
                                keys=columns, names = ["var", "old_index"], **opt).reset_index().drop(labels="old_index", axis = 1)

data["kyears BP"] = 32.030-data["time"] / 1000.
data["model"] = data["model"].str.replace("hO", "H")
val=data["model"].isin(validated)

#Sensitivity Analysis
fw_vols = pd.read_csv(os.path.join(resfolder, "end_volumes.csv")).drop(columns=["total_fw_pumpable", "total_onshore"]).rename(columns={"total_fw":"V (km3)"})
#Homgn model got a bit of a sloppy code
fw_vols["sea"].replace({"H" : "O"}, inplace=True)
fw_vols["sea"].replace({"hO" : "H"}, inplace=True)
fw_vols["clayer"].fillna("N", inplace=True)
fw_vols["brine"].replace({"Top" : "T", "Bot" : "B"}, inplace=True)
fw_vols = fw_vols.rename(columns={"brine":"HGw provenance", "sea": "sea connection", "clayer": "clay layer"})

fw_onshore = pd.read_csv(os.path.join(resfolder, "end_volumes.csv"))
fw_onshore = fw_onshore.loc[fw_onshore["codes"].isin(validated), ["codes", "total_fw", "total_fw_pumpable", "total_onshore"]]
fw_onshore.to_csv(os.path.join(resfolder,"table_paper.csv"))

#timeslices
ka = (32000 - np.cumsum(np.array([0, 15000, 3500, 4500, 3000, 3000, 2000, 1000])))/1000

#%%Analyze a bit
fac={}
for i, mod in enumerate(validated):
    select = data.loc[(data["var"] == "total_fw") & (data["model"] == mod), :].reset_index()
    fac[mod] = select.loc[0, "V (km3)"]/select.loc[135, "V (km3)"]
    
print(fac)
#%%Plot
sns.set(style="whitegrid")

#Colors
color_trans = [0.171875, 0.48046875, 0.7109375]
palette = {"Paleo":"xkcd:plum", "Steady":"xkcd:gold"}

#Initialize grid
f = plt.figure(figsize=(8,10))
gs = gridspec.GridSpec(3, 1, height_ratios = [20, 1, 20])
ax = plt.subplot(gs[0])

#Create this one for the title of the bottom plot
with sns.axes_style("white"): #No grid for this axes
    ax_text = plt.subplot(gs[1]) 
    ax_text.text(0.5, -1., "Sensitivity total FGw", horizontalalignment='center', verticalalignment='center', transform=ax_text.transAxes)
    sns.despine(ax=ax_text, left=True, bottom=True) #Remove automatically created spines
    ax_text.set_xticks([], [])                      #Remove automatically created xticks
    ax_text.set_yticks([], [])                      #Remove automatically created yticks
    ax_text.text(0.01, -1., r"$\bf{b)}$", horizontalalignment='center', verticalalignment='center', transform=ax_text.transAxes)

data = data.rename(columns={"kyears BP":"time (ka)"})

#Plot transience
g1 = sns.lineplot(ax = ax, x="time (ka)", y="V (km3)", style="model", data=data.loc[val & (data["var"] == "total_fw"), :], color=color_trans)
g1.axes.invert_xaxis()
g1.axes.grid(b=False, axis="x")
g1.set_ylim(bottom=0., top=g1.get_ylim()[-1])
g1.set_title("Total FGw dynamics")
g1.legend(loc="upper right")
sns.despine(ax=g1, left=True)
g1.axes.text(0.01, 1.04, r"$\bf{a)}$", horizontalalignment='center', verticalalignment='center', transform=g1.axes.transAxes)

#Make lines in legend blue and title bold
leg = g1.get_legend()
texts = leg.get_texts()
texts[0].set_text(r"$\bf{%s}$" % (texts[0].get_text()))
lines=[line.set_color(color_trans) for line in leg.get_lines()]

#Replace ticks and draw vertical dotted lines at timeslice edges
xticklabels = [i.split(".0")[0] for i in list(ka.astype(str))]
g1.set(xticklabels=xticklabels, xticks=ka)
for i, ts in enumerate(ka):
#    g1.axes.axvline(ts, ls=":", color=".60")
    g1.axes.axvline(ts, lw=1, color=".80", zorder=-100)
    if i > 0:
        label_x = (ka[i-1] + ka[i])/2.
        g1.axes.text(label_x, 400, "sp%d" % (i), horizontalalignment="center", fontdict={"weight" : "bold"})



#Create PairGrid
g2 = sns.PairGrid(fw_vols, y_vars=["V (km3)"],
                 x_vars=["HGw provenance", "sea connection", "clay layer"],
                 height=5, aspect=.5, hue="steady", palette = palette)

g2.map(sns.stripplot, alpha="0.6", size=9, jitter=True)
sns.despine(fig=g2.fig, left=True)
g2.add_legend(title="")

handles, labels = g2.axes[0, 0].get_legend_handles_labels()

#Add PairGrid to gridspec
mg2 = sfg.SeabornFig2Grid(g2, f, gs[-1])
plt.legend(handles[::3], labels[::3], title="", loc="upper right")

plt.tight_layout()
f.savefig(os.path.join(resfolder, "_FW_Paleo.png"))
f.savefig(os.path.join(resfolder, "_FW_Paleo_300.png"), dpi=300) #This messes up bottom panel probably due to hacks
f.savefig(os.path.join(resfolder, "_FW_Paleo.svg"))
f.savefig(os.path.join(resfolder, "_FW_Paleo.pdf"))
plt.close()