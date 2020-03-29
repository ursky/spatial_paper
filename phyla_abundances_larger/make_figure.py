#!/usr/bin/env python2
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from collections import OrderedDict
from scipy import stats
import numpy as np
import matplotlib.gridspec as gridspec


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
					if taxon.startswith("D_0"):
						taxon = "D_0__Other;D_1__Other"
					else:
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
			data[k]["location"]="SG1-Top"
		if data[k]["location"]=="Bottom":
			data[k]["location"]="SG1-Bottom"
	df = pd.DataFrame.from_dict(data).T
	return data, df


def reformat_data(data):
	samples=[]
	for k in data:
		taxon = data[k]["taxon"]
		abund = data[k]["abundance"]
		sample = data[k]["sample"]
		location = data[k]["location"]
		name = location+"|"+sample
		if name not in samples:
			samples.append(name)
	samples.sort()
	N=len(samples)
	sample_locs = {}
	for i, sample in enumerate(samples):
		sample_locs[sample]=i

	new_data={}

	for k in data:
		taxon = data[k]["taxon"]
		abund = data[k]["abundance"]
		sample = data[k]["sample"]
		location = data[k]["location"]
		name = location+"|"+sample
		if taxon not in new_data:
			new_data[taxon] = [''] * N
		loc = sample_locs[name]
		new_data[taxon][loc] = abund
	return new_data, samples


def insert_gap(samples, abundances):
	new_abunds = []
	site=""
	for i, abund in enumerate(abundances):
		sample = samples[i]
		new_site = sample.split("|")[0]
		if new_site!=site:
			if site=="":
				site=new_site
			else:
				site=new_site
				new_abunds.append(0)
		
		new_abunds.append(abund)
	return new_abunds


def draw_plot(data, samples, ax):
	N = len(samples) + 1
	ind = np.arange(N)
	width = 1

	base = [0] * N
	plots={}

	taxa = ['Nanohaloarchaeota', 'Actinobacteria', 'Cyanobacteria', 'Chlorophyta', 'Proteobacteria', 'Bacteroidetes', 'Euryarchaeota']	
	for taxon in taxa:
		abundances = data[taxon]
		abundances = insert_gap(samples, abundances)
		plots[taxon] = ax.bar(ind, abundances, width, bottom=base, linewidth=0.5)
		for i, val in enumerate(abundances):
			base[i]+=val

	#ax.set_xlabel('Samples')
	ax.set_ylabel('Relative abundance (%)')
	#ax.set_xticks(ind, samples)
	ax.set_xticks([]) 
	ax.set_xticklabels([])
	ax.set_yticks(np.arange(0, 101, 10))

	# draw legend and labels
	handles = []
	taxa.reverse()
	for taxon in taxa:
		handles.append(plots[taxon][0])
	for i,taxon in enumerate(taxa):
		taxa[i] = "$\it{" + taxon + "}$"
	if "SG1-Top|SG1-TOP-2016-02mt-9b" in samples:
		ax.legend(handles, taxa, bbox_to_anchor=(1, 0.7))
	ax.grid(axis="y", ls="--", c="k", alpha=0.1)
	ax.spines['right'].set_visible(False)
	ax.spines['top'].set_visible(False)
	ax.spines['bottom'].set_visible(False)

	if "SG1-Bottom" in samples[0]:
		ax.text(0, -6, "North-Bottom", fontsize=14)
		ax.text(11, -6, "North-Top", fontsize=14)
	else:
		ax.text(17, -6, "North", fontsize=14)
		ax.text(52, -6, "South", fontsize=14)


########################################################
fig = plt.figure(figsize=(8,8))

gs1 = gridspec.GridSpec(2,10, figure=fig)
#gs1.update(left=0.05, right=0.48, wspace=0.05)
ax1 = plt.subplot(gs1[0, :])
ax2 = plt.subplot(gs1[1, :6])
#ax3 = plt.subplot(gs1[1, 1])


sns.set(style="whitegrid")


data, df = load_data("large-3.csv")
data, samples = reformat_data(data)
draw_plot(data, samples, ax1)

data, df = load_data("medium-3.csv")
data, samples = reformat_data(data)
draw_plot(data, samples, ax2)

ax1.annotate("A.", xy=(-0.09, 0.95), xycoords="axes fraction", fontsize=20)
ax2.annotate("B.", xy=(-0.15, 0.95), xycoords="axes fraction", fontsize=20)

plt.tight_layout()
plt.savefig("figure_full.png", dpi=300)
