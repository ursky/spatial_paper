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
	color_lib = ["red", "magenta", "cyan", "gold", "green", "blue", "pink"]
	shape_lib = ["o", "s", "v", "X", "P"]
	colors = {}
	shapes = {}
	n_colors = 0
	n_shapes = 0
	for i,sample in sorted(enumerate(samples), reverse=True):
		meta = metadata[i]
		site=meta[metadata_headers[category]]
		x = coordinates[i][pc1]
		y = coordinates[i][pc2]

		if site not in colors:
			if n_colors>=len(color_lib):
				n_colors = 0
				n_shapes += 1
			colors[site]=color_lib[n_colors]
			shapes[site]=shape_lib[n_shapes]
			ax.scatter(x, y, c=colors[site], marker=shapes[site], edgecolors="k", s=50, label=site)
			n_colors += 1
		else:
			ax.scatter(x, y, c=colors[site], marker=shapes[site], s=50, edgecolors="k")
	ax.set_xlabel("PC%d (%d%% explained)" % (pc1+1, int(explained[pc1])))
	ax.set_ylabel("PC%d (%d%% explained)" % (pc2+1, int(explained[pc2])))
	
	# legend
	handles, labels = ax.get_legend_handles_labels()
	if category=="Position":
		labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0], reverse=True))
	else:
		labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0], reverse=False))
	if category!="Slice":
		ax.legend(handles, labels, frameon=True)


	ax.grid(ls="--", c="k", alpha=0.2)
		


small_data = load_pcoa("small.dict")

fig, axs = plt.subplots(3,1, figsize=(6,8))
plot_scatter(small_data, 0, 1, "Position", axs[0])
plot_scatter(small_data, 0, 1, "Rock", axs[1])
plot_scatter(small_data, 0, 1, "Slice", axs[2])

axs[0].annotate("A.", xy=(-0.11, 0.95), xycoords="axes fraction", fontsize=20)
axs[1].annotate("B.", xy=(-0.11, 0.95), xycoords="axes fraction", fontsize=20)
axs[2].annotate("C.", xy=(-0.11, 0.95), xycoords="axes fraction", fontsize=20)


plt.tight_layout()
plt.savefig("figure_small_primary.png", dpi=300)












