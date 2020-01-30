#!/usr/bin/env python2
import sys
import matplotlib.pyplot as plt
plt.style.use('seaborn-whitegrid')
import numpy as np


def load_pcoa(filename):
	f = open(filename)
	line = f.readline()
	data = eval(line.strip())
	return data

def plot_scatter(data, ax):
	explained = data["decomposition"]["percents_explained"]
	samples = data["decomposition"]["sample_ids"]
	coordinates = data["decomposition"]["coordinates"]
	metadata = data["metadata"]
	colors = {}
	for i,sample in enumerate(samples):
		meta = metadata[i]
		if len(meta)==2:
			site=meta[-1]
		else:
			site=meta[-2]
		x = coordinates[i][0]
		y = coordinates[i][1]
		if site=="Top": site="SG1-Top"
		if site=="Bottom": site="SG1-Bottom"
		if colors=={}:
			colors[site]="r"
			ax.scatter(x, y, c=colors[site], edgecolors="k", label=site, s=50, alpha=0.5)
		elif site not in colors:
			colors[site]="b"
			ax.scatter(x, y, c=colors[site], edgecolors="k", label=site, s=50, alpha=0.5)
		else:
			ax.scatter(x, y, c=colors[site], edgecolors="k", s=50, alpha=0.5)
	ax.set_xlabel("PC1 (%d%% explained)" % (int(explained[0])))
	ax.set_ylabel("PC2 (%d%% explained)" % (int(explained[1])))

	# legend ordering
	handles, labels = ax.get_legend_handles_labels()
	if "SG2" in colors:
		labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0], reverse=False))
	else:
		labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0], reverse=True))
	ax.legend(handles, labels, frameon=True)

	ax.grid(ls="--", c="k", alpha=0.2)
		
	


medium_data = load_pcoa("medium.dict")
large_data = load_pcoa("large.dict")


fig, [ax1, ax2] = plt.subplots(1,2, figsize=(8,4))
plot_scatter(large_data, ax1)
plot_scatter(medium_data, ax2)

ax1.annotate("A.", xy=(-0.19, 0.95), xycoords="axes fraction", fontsize=20)
ax2.annotate("B.", xy=(-0.19, 0.95), xycoords="axes fraction", fontsize=20)

plt.tight_layout()
plt.savefig("figure.png", dpi=300)












