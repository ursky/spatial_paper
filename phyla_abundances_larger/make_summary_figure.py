#!/usr/bin/env python2
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy import stats


def load_data(filename):
	data = {}
	headers = {}
	toi = set(["D_0__Archaea;D_1__Euryarchaeota", \
		"D_0__Archaea;D_1__Nanohaloarchaeota", \
		"D_0__Bacteria;D_1__Actinobacteria", \
		"D_0__Bacteria;D_1__Bacteroidetes", \
		"D_0__Bacteria;D_1__Cyanobacteria", \
		"D_0__Bacteria;D_1__Proteobacteria", \
		"D_0__Eukarya;D_1__Chlorophyta"])
	for i,line in enumerate(open(filename)):
		cut=line.strip().split(",")
		if i==0:
			head = cut[:]
			for j,col in enumerate(head):
				headers[col]=j
		else:
			sample = cut[0]
			location = cut[headers["Location"]]
			# calculate total sum
			total = 0
			for i,val in enumerate(cut):
				taxon = head[i]
				if "D_1__" in taxon:
					total+=float(val)
			for i,val in enumerate(cut):
				taxon = head[i]
				if "Chloroplast" in taxon:
					taxon = "D_0__Eukarya;D_1__Chlorophyta"
				taxon = ";".join(taxon.split(";")[:2])
				if taxon not in toi:
					continue
				abund = 100.0*float(val)/total
				index = sample+"_"+taxon
				if index not in data:
					data[index] = {}
					data[index]["taxon"]=taxon.split("_")[-1]
					data[index]["location"]=location
					data[index]["sample"]=sample
					data[index]["abundance"]=abund
				else:
					data[index]["abundance"]+=abund
	for k in data:
		if data[k]["abundance"]==0:
			data[k]["abundance"]=0.1
		if data[k]["location"]=="Top":
			data[k]["location"]="North-Top"
		if data[k]["location"]=="Bottom":
			data[k]["location"]="North-Bottom"
		if data[k]["location"]=="SG1":
			data[k]["location"]="North"
		if data[k]["location"]=="SG2":
			data[k]["location"]="South"
	df = pd.DataFrame.from_dict(data).T
	return data, df


def get_max_min_in_dict(dictionary):
        maxs=[]; mins=[]
        for key in dictionary:
                maxs.append(max(dictionary[key]))
                mins.append(min(dictionary[key]))
        return max(maxs), min(mins)


def draw_signifficance_bars(data, labels, ax):
	for x_st,taxon in enumerate(labels):
		groups = {}
		lists = []
		hmax = 0
		for sample in data:
			if taxon not in data[sample]["taxon"]:
				continue
			location = data[sample]["location"]
			abund = data[sample]["abundance"]
			if abund>hmax:
				hmax=abund
			if location not in groups:
				groups[location]=len(lists)
				lists.append([])
			lists[groups[location]].append(abund)

		test=stats.ttest_ind(lists[0], lists[1])
		if test.pvalue > 0.01: continue
		elif test.pvalue > 0.001: m='*'
		elif test.pvalue > 0.0001: m='**'
		else: m='***'
		ax.hlines(y=hmax+5, xmin=x_st-0.25, xmax=x_st+0.25, linewidth=1, color='k')
		ax.text(x_st, hmax+5, m, ha='center', fontsize=12)


def draw_boxplot(filename, ax):
	data, df = load_data(filename)
	taxa = ["Euryarchaeota", "Bacteroidetes", "Cyanobacteria", "Chlorophyta", "Proteobacteria", "Actinobacteria", "Nanohaloarchaeota"]
	df["abundance"] = df["abundance"].astype(float)

	ax = sns.boxplot(x="taxon", y="abundance", data=df, color="w", hue="location", order=taxa, ax=ax)
	ax = sns.swarmplot(x="taxon", y="abundance", data=df, hue="location", split=True, edgecolor="k", order=taxa, ax=ax, alpha=0.6, linewidth=0.5)

	ax.set_ylabel("Relative abundance (%)")
	if "medium" in filename:
		ax.set_xlabel("Phyla")
	else:
		ax.get_xaxis().set_visible(False)

	handles, labels = ax.get_legend_handles_labels()
	l = ax.legend(handles[2:], labels[2:], frameon=True)


	ax.set_xticklabels(taxa, rotation=45)
	#ax.set_ylim(0,100)

	for x in [1,2,3,4,5,6]:
		ax.axvline(x=x-0.5, linewidth=1, color='k', alpha=0.2)

	draw_signifficance_bars(data, taxa, ax)
	ax.set_yscale('log')



fig, axs = plt.subplots(2,1, figsize=(8,8))
sns.set(style="whitegrid")

draw_boxplot("large-3.csv", axs[0])
draw_boxplot("medium-3.csv", axs[1])

axs[0].annotate("A.", xy=(-0.09, 0.95), xycoords="axes fraction", fontsize=20)
axs[1].annotate("B.", xy=(-0.09, 0.95), xycoords="axes fraction", fontsize=20)

axs[0].grid(axis="y", ls="--", alpha=0.2, c="k")
axs[1].grid(axis="y", ls="--", alpha=0.2, c="k")

plt.tight_layout()
plt.savefig("figure.png", dpi=300)
