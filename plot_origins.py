# -*- coding: utf-8 -*-
"""
Created on Mon Feb  4 13:09:07 2019

@author: engelen
"""

import pandas as pd
import seaborn as sns
import os
from glob import glob
import matplotlib.pyplot as plt
import numpy as np

#%%
class _LinePlotter_minmax(sns.relational._LinePlotter):
    def aggregate(self, vals, grouper, units=None):
        """OVERRIDE bootstrapping ci estimating, we just want the range of the raw data"""
        func = self.estimator
        ci = self.ci
    
        # Group and get the aggregation estimate
        grouped = vals.groupby(grouper, sort=self.sort)
        est = grouped.agg(func)
    
        # Exit early if we don't want a confidence interval
        if ci is None:
            return est.index, est, None
        else:
            _min = grouped.min()
            _max = grouped.max()
            cis = pd.DataFrame(np.c_[_min, _max],
                               index=est.index,
                               columns=["low", "high"]).stack()
        # Unpack the CIs into "wide" format for plotting
        if cis.notnull().any():
            cis = cis.unstack().reindex(est.index)
        else:
            cis = None

        return est.index, est, cis

#%%Path management
folder = r"./plots/origins/"


#%%Data handling
validated = ["CmCBP", "CnCTP", "hOmCTP", "hOfCTP", "hOnCTP"]
brine     = ["Bot",   "Top",   "Top",    "Top",    "Top"]
sources   = ["Sea", "River", "Dunes", "HGw"]

dfs = [pd.read_csv(os.path.join(folder, '{}_origins.csv'.format(mod))).assign(code=mod, brine_source=brine[i]).iloc[slice(0,129)].rename(
        columns={"Brine" : "HGw", "brine_source" : "HGw \; prov"}) for i, mod in enumerate(validated)]
dfs = [df.assign(total=df.loc[:, sources].sum(axis=1)) for df in dfs]
dfs = [df.assign(**dict((source, df[source]/df["total"]) for source in sources)) for df in dfs]
df=pd.concat(dfs, axis=0, sort=False).drop(columns="total")
df_long=pd.concat([df.loc[:, ["time", "code", "HGw \; prov", source]].assign(type=source).rename(columns={source : "fractions"}) for source in sources], axis=0, sort=False)
df_long["kyears BP"] = 32.0-df_long["time"] / 1000.
df_long=df_long.rename(columns={"code" : "model", "fractions":"relative volume (-)", "kyears BP" : "time (ka)"})

#df_long=pd.concat([df_long, df_long.loc[df_long["brine"]=="Bot"]])

#timeslices
ka = (32000 - np.cumsum(np.array([0, 15000, 3500, 4500, 3000, 3000, 2000, 1000])))/1000

#%%Plot
palette = {"HGw" :(0.83984375, 0.09765625, 0.109375), "Sea":(0.9609375, 0.5625, 0.32421875), 
           "River" : (0.171875, 0.48046875, 0.7109375), "Dunes" : (0.99609375, 0.87109375, 0.6015625)}
sns.set(style="whitegrid")

f, ax = plt.subplots(figsize=(10,7))

p = _LinePlotter_minmax(x="time (ka)", y="relative volume (-)", hue="type", 
                                data=df_long, style = "HGw \; prov", dashes=True, style_order=["Top", "Bot"],
                                palette = palette, ci=0, n_boot=0, estimator="median", 
                                err_style="band", legend="brief")

p.plot(ax, dict(lw=3, alpha=0.8))

#Replace ticks and draw vertical dotted lines at timeslice edges
xticklabels = [i.split(".0")[0] for i in list(ka.astype(str))]
ax.set(xticklabels=xticklabels, xticks=ka)
for i, ts in enumerate(ka):
#    g1.axes.axvline(ts, ls=":", color=".60")
    ax.axvline(ts, lw=1, color=".80", zorder=-100)
    if i > 0:
        label_x = (ka[i-1] + ka[i])/2.
        ax.text(label_x, 0.05, "sp%d" % (i), horizontalalignment="center", fontdict={"weight" : "bold"})

plt.legend(loc = "upper right")

ids = [0,5]
lgtxt = ax.get_legend().get_texts()

[lgtxt[i].set_text(r"$\bf{%s}$" % (lgtxt[i].get_text())) for i in ids]

#g.axes.invert_xaxis()
plt.ylim(0,1)
plt.xlim(32.3, 0)
plt.tight_layout()
plt.savefig(os.path.join(folder, "Origins_behavioral_shade.png"), dpi=300)
plt.savefig(os.path.join(folder, "Origins_behavioral_shade.svg"))
plt.savefig(os.path.join(folder, "Origins_behavioral_shade.pdf"))
plt.close()

#%%Plot
f, ax = plt.subplots(figsize=(10,7))
g = sns.lineplot(x="time (ka)", y="relative volume (-)", 
                 hue="type", style="model", data=df_long, palette = palette, lw=5, alpha=0.5)

#g.axes.invert_xaxis()
plt.ylim(0,1)
plt.xlim(32.3, 0)
plt.tight_layout()
plt.savefig(os.path.join(folder, "Origins_behavioral_sep.png"))
plt.close()