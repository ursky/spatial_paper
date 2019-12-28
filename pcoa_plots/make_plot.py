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

		if colors=={}:
			colors[site]="r"
			ax.scatter(x, y, c=colors[site], edgecolors="k", label=site)
		elif site not in colors:
			colors[site]="b"
			ax.scatter(x, y, c=colors[site], edgecolors="k", label=site)
		else:
			ax.scatter(x, y, c=colors[site], edgecolors="k")
	ax.set_xlabel("PC1 (%d%% explained)" % (int(explained[0])))
	ax.set_ylabel("PC2 (%d%% explained)" % (int(explained[1])))
	ax.legend()
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












