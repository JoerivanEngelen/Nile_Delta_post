# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 13:19:41 2019

@author: engelen
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from glob import glob
import xarray as xr
import seaborn as sns

def get_files(folder, scenarios, ext):
    return([glob(os.path.join(folder, "*"+scenario+"*"+ext))[0] for scenario in scenarios])

#%%Path management
folder = r"c:\Users\engelen\OneDrive - Stichting Deltares\PhD\Nile_Delta\3D_Model\Results\Plots\conductance\\"

scenarios = ["base", "div2", "prod2"]

or_fol = os.path.join(folder, "origins")
fw_fol = os.path.join(folder, "FW_volumes")

fw_files = get_files(fw_fol, scenarios, ".nc")
or_files = get_files(or_fol, scenarios, ".csv")

#%%Plot fw volumes
fws = [xr.open_dataset(file)["total_fw"] for file in fw_files]

colors = ["xkcd:grey blue", "xkcd:plum", "xkcd:bright orange"]

lines = []
for i, fw in enumerate(fws):
    lines.append(fw.plot(color = colors[i])[0])

plt.legend(lines, scenarios)

plt.savefig(os.path.join(folder, "fw_comparison.png"))
plt.close()

#%%Read and reshuffle dataframe to get accepted by seaborn.
sources = ["HGw", "Sea", "River", "Dunes"]

#Conductance
#scenarios_out = ["1", "0.5", "2"]

#Resistance
scenarios_out = ["1", "2", "0.5"]

ors = [pd.read_csv(file).drop(columns=["dx", "dy"]).set_index("time").rename(columns={"Brine":"HGw"}) for file in or_files]
tot_vols = [df.loc[:, sources].sum(axis=1) for df in ors]
ors = [df.assign(**dict((source, df[source]/tot_vols[i]) for source in sources)) for i, df in enumerate(ors)]
ors = [ori.stack().reset_index().set_axis(["time (ka)", "type", "relative volume (-)"], axis=1, inplace=False) for ori in ors]
ors = [ori.assign(resistance=scenarios_out[i]) for i, ori in enumerate(ors)]
ors = pd.concat(ors, axis=0)

#%%Set times
ors["time (ka)"] = 32.0-ors["time (ka)"] / 1000.

#timeslices
ka = (32000 - np.cumsum(np.array([0, 15000, 3500, 4500, 3000, 3000, 2000, 1000])))/1000

#%%Plot origin volumes

sns.set(style="whitegrid")
f, ax = plt.subplots(figsize=(10,7))

palette = {"HGw" :(0.83984375, 0.09765625, 0.109375), "Sea":(0.9609375, 0.5625, 0.32421875), 
           "River" : (0.171875, 0.48046875, 0.7109375), "Dunes" : (0.99609375, 0.87109375, 0.6015625)}

g1 = sns.lineplot(x="time (ka)", y="relative volume (-)", 
                  hue="type", style="resistance", data=ors, palette = palette, legend="brief")

#Replace ticks and draw vertical dotted lines at timeslice edges
xticklabels = [i.split(".0")[0] for i in list(ka.astype(str))]

g1.axes.set(xticklabels=xticklabels, xticks=ka)
for i, ts in enumerate(ka):
    g1.axes.axvline(ts, lw=1, color=".80", zorder=-100)
    if i > 0:
        label_x = (ka[i-1] + ka[i])/2.
        g1.axes.text(label_x, 0.05, "sp%d" % (i), horizontalalignment="center", fontdict={"weight" : "bold"})

g1.axes.legend(loc = "upper right")

ids = [0,5]
lgtxt = ax.get_legend().get_texts()

[lgtxt[i].set_text(r"$\bf{%s}$" % (lgtxt[i].get_text())) for i in ids]

#g.axes.invert_xaxis()
plt.ylim(0,1)
plt.xlim(32.0, 0)
#plt.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.) #Place legend outside ax
plt.tight_layout()

plt.savefig(os.path.join(folder, "or_comparison.png"))
plt.savefig(os.path.join(folder, "or_comparison.svg"))
plt.savefig(os.path.join(folder, "or_comparison.pdf"))