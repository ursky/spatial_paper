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

def plot_scatter(data, pc1, pc2, category, ax):
	explained = data["decomposition"]["percents_explained"]
	samples = data["decomposition"]["sample_ids"]
	coordinates = data["decomposition"]["coordinates"]
	metadata = data["metadata"]
	headers = data["metadata_headers"]
	metadata_headers = {}
	for i,meta in enumerate(headers):
		metadata_headers[meta] = i
	color_lib = ["red", "cyan", "magenta"]
	colors = {}
	for i,sample in enumerate(samples):
		meta = metadata[i]
		site=meta[metadata_headers[category]]
		if site=="SG1": site="North"
		if site=="SG2": site="South"
		if site=="Top" and pc1==0: site="North-Top"
		if site=="Bottom" and pc1==0: site="North-Bottom"
		x = coordinates[i][pc1]
		y = coordinates[i][pc2]

		if site not in colors:
			colors[site]=color_lib[len(colors)]
			ax.scatter(x, y, c=colors[site], edgecolors="k", s=50, label=site, alpha=0.9)
		else:
			ax.scatter(x, y, c=colors[site], s=50, edgecolors="k", alpha=0.9)
	ax.set_xlabel("PC%d (%d%% explained)" % (pc1+1, int(explained[pc1])))
	ax.set_ylabel("PC%d (%d%% explained)" % (pc2+1, int(explained[pc2])))

	# legend ordering
        handles, labels = ax.get_legend_handles_labels()
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0], reverse=True))
	print labels
        ax.legend(handles, labels, frameon=True)
	ax.grid(ls="--", c="k", alpha=0.2)	

def plot_continuous(data, pc1, pc2, category, ax):
	explained = data["decomposition"]["percents_explained"]
	samples = data["decomposition"]["sample_ids"]
	coordinates = data["decomposition"]["coordinates"]
	metadata = data["metadata"]
	headers = data["metadata_headers"]
	metadata_headers = {}
	for i,meta in enumerate(headers):
		metadata_headers[meta] = i
	xs=[]; ys=[]; vs=[]
	for i,sample in enumerate(samples):
		meta = metadata[i]
		value = meta[metadata_headers[category]]
		value = float(value)
		x = coordinates[i][pc1]
		y = coordinates[i][pc2]
		xs.append(x)
		ys.append(y)
		vs.append(value)

	sc = ax.scatter(xs, ys, s=50, edgecolors="k", c=vs, cmap="viridis")
	ax.set_xlabel("PC%d (%d%% explained)" % (pc1+1, int(explained[pc1])))
	ax.set_ylabel("PC%d (%d%% explained)" % (pc2+1, int(explained[pc2])))
	ax.legend()
	ax.grid(ls="--", c="k", alpha=0.2)
	fig.colorbar(sc, ax=ax)

small_data = load_pcoa("small.dict")
medium_data = load_pcoa("medium.dict")
large_data = load_pcoa("large.dict")


fig, axs = plt.subplots(2,2, figsize=(8,8))
plot_scatter(large_data, 0, 1, "Location", axs[0,0])
plot_scatter(medium_data, 0, 1, "Location", axs[0,1])
plot_scatter(small_data, 2, 3, "Position", axs[1,0])
plot_continuous(small_data, 2, 3, "Distance to sruface", axs[1,1])


axs[0,0].annotate("A.", xy=(-0.19, 0.95), xycoords="axes fraction", fontsize=20)
axs[0,1].annotate("B.", xy=(-0.19, 0.95), xycoords="axes fraction", fontsize=20)
axs[1,0].annotate("C.", xy=(-0.19, 0.95), xycoords="axes fraction", fontsize=20)
axs[1,1].annotate("D.", xy=(-0.19, 0.95), xycoords="axes fraction", fontsize=20)


plt.tight_layout()
plt.savefig("figure.png", dpi=300)





