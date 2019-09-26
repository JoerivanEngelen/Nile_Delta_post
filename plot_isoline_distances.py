# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 16:44:04 2019

@author: engelen
"""

import pandas as pd
import seaborn as sns
from glob import glob
import os
import matplotlib.pyplot as plt
import numpy as np

#%%Estimator
def median_count(arr, lim_count=100):
    if len(arr) < lim_count:
        return(np.nan)
    else:
        return(np.median(arr))

#%%Path Management
plotfolder = r"/plots/"
res_folder = os.path.join(plotfolder, "distance_isolines")

isoval=[3.0, 10., 20., 30.]

sns.set_style("darkgrid", {"axes.facecolor": ".65", "ytick.left": True})

colors = [[0.5703125, 0.7734375, 0.87109375],
          [0.8671875, 0.9375, 0.8125],
          [0.99609375, 0.87109375, 0.6015625],
          [0.9609375, 0.5625, 0.32421875]]


color_pal = dict(i for i in zip(isoval, colors))

#%%Combination plot
fig, axes = plt.subplots(1, 4, figsize=(8, 4))
codes = ["H-T", "hO-F-T", "C-N-T", "hO-N-T",]
files = [os.path.join(res_folder, code+r".csv") for code in codes]
codes = ["O-N-T", "H-F-T", "C-N-T", "H-N-T"]

estimator = np.median
#estimator = median_count

for i, f in enumerate(files):
    dist = pd.read_csv(f)
    
    sns_plot = sns.pointplot(ax = axes[i],  x = "distance (km)", y = "z (m)", data = dist, hue = "TDS (g/l)",
                             orient = "h", ci = None, estimator=estimator, markers=".", palette = color_pal)

    #Make sure pointplot is rendered on top of stripplot.
    plt.setp(sns_plot.lines, zorder=100)
    plt.setp(sns_plot.collections, zorder=100, label="")

    sns_plot = sns.stripplot(ax = axes[i], x = "distance (km)", y = "z (m)", data = dist, hue = "TDS (g/l)", 
                        orient = "h", alpha = 0.05, jitter=True, dodge=True, palette = color_pal, size=2.5)
    
    sns_plot.invert_yaxis()
    sns_plot.set_xlim(-15, 15)
    sns_plot.set_title(codes[i])
    
    set_d = dict(xlabel='$\omega$ (km)')
        
    #Correct yticklabels
    if i == 0:
        nu_labels = sns_plot.get_yticklabels()[::2]
        sns_plot.set_yticks(sns_plot.get_yticks()[::2])
        sns_plot.set_yticklabels(nu_labels)
        set_d["ylabel"]="z (m)"
    else:
        sns_plot.set_yticks([])
        set_d["ylabel"]=""
    
    axes[i].set(**set_d)
    
    if f == files[-1]:
        handles, labels = sns_plot.get_legend_handles_labels()
        l = sns_plot.legend(handles[:len(isoval)], labels[:len(isoval)], title="TDS (g/l)",
                  handletextpad=0, columnspacing=1,
                  loc="lower right", ncol=1, frameon=True)
        for text in l.get_texts():
            text.set_color("white")
        l.get_title().set_color("white")
    else:
        sns_plot.legend_.remove()
    
plt.tight_layout()

plt.savefig(os.path.join(plotfolder, r"distance_isolines", "distance_isolines_three_model.png"), dpi=300)
plt.savefig(os.path.join(plotfolder, r"distance_isolines", "distance_isolines_three_model.svg"))


plt.close()


