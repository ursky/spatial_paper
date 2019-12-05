#!/usr/bin/env python2
import sys
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def load_data(filename):
	data={}
	formation={1:"polygon", 2:"rising polygon", 3:"nodule", 4:"mature nodule"}
	for line in open(filename):
		cut=line.strip().split("\t")
		if cut[0]=="Sites":
			continue
		site = cut[0].split("-")[0]
		if site=="S": site="SG2"
		if site=="NTop": site="SG1-Top"
		if site=="NBot": site="SG1-Bottom"
		transect = cut[0].split("-")[1]
		total=0
		for i in cut[1:]: total+=float(i)
		for i in range(1,5):
			entry = cut[0]+"-"+str(i)
			data[entry]={}
			data[entry]["Formation"] = formation[i]
			data[entry]["Percent coverage"] = 100.0*float(cut[i])/total
			data[entry]["Site"] = site
	return pd.DataFrame.from_dict(data).T


df = load_data("transect_data.tab")
df["Percent coverage"] = df["Percent coverage"].astype(float)
print df

# plotting set-up
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)
plt.rc('font', family='arial')
fig, ax = plt.subplots(1, 1, figsize=(10,6))

#colors = {"SG2": "cyan", "SG1-Top": "magenta", "SG1-Bottom":"gold"}
colors = {"polygon": "cyan", "rising polygon": "magenta", "nodule":"gold", "mature nodule": "red"}

sns.boxplot(hue="Formation", y="Percent coverage", x="Site", data=df, ax=ax, palette=colors, linewidth=1, fliersize=3).set(xlabel="")
sns.swarmplot(hue="Formation", y="Percent coverage", x="Site", data=df, ax=ax, palette=colors, linewidth=0.5, split=True, size=10, label=None, alpha=1).set(xlabel="")

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[:len(handles)/2], labels[:len(labels)/2], ncol=4, loc='lower center', bbox_to_anchor=(0.5, -0.15))
print "what"

plt.tight_layout()
plt.savefig("figure.png", dpi=300)
