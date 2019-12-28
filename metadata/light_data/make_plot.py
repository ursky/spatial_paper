#!/usr/bin/env python
import sys,os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def load_data(dirname):
	data = {}
	for filename in os.listdir(dirname):
		begin = False
		data[filename] = ( [], [] )
		for line in open(dirname+"/"+filename):
			if "Begin" in line:
				begin = True
				continue
			if begin==True:
				cut = line.strip().split()
				if cut[1]=="NaN":
					continue
				if float(cut[0])<500 or float(cut[0])>900:
					continue
				data[filename][0].append(float(cut[0]))
				data[filename][1].append(float(cut[1]))
	return data

def load_solar_data(filename):
	x=[]
	y=[]
	for line in open(filename):
		if line.startswith("Wvl"):
			continue
		cut = line.strip().split()
		if float(cut[0])<500 or float(cut[0])>900:
			continue
		x.append(float(cut[0]))
		y.append(float(cut[2]))
	for i, val in enumerate(y):
		y[i] = val/max(y)
	return (x, y)


def rename_labels(data):
	out = {}
	for k in data:
		if k=="ha04_a30_transmission.txt": new_name="Halite 1, middle"
		elif k=="ha04_a50_transmission.txt": new_name="Halite 1, top"
		elif k=="ha05_c16_transmission.txt": new_name="Halite 2, middle"
		elif k=="ha05_c30_transmission.txt": new_name="Halite 2, top"
		elif k=="ha06_a35_transmission.txt": new_name="Halite 3, middle"
		elif k=="ha06_a68_transmission.txt": new_name="Halite 3, top"
		else: new_name=k
		out[new_name] = data[k]
	return out
		

def standardize_data(data):
	means = {}
	for i in ["1","2","3"]:
		sample = "Halite "+i+", top"
		mean = np.mean(data[sample][1])
		means[i] = mean
	mean = np.mean(means.values())
	for i in ["1","2","3"]:
		adjustment = means[i]/mean
		sample1 = "Halite "+i+", top"
		sample2 = "Halite "+i+", middle"
		for j,absorbance in enumerate(data[sample1][1]):
			data[sample1][1][j] = absorbance/adjustment
		for j,absorbance in enumerate(data[sample2][1]):
			data[sample2][1][j] = absorbance/adjustment
	return data

def plot_data(data):
	font = {'family': 'arial', 'weight': 'normal', 'size': 12}
	plt.rc('font', **font)
	fig, ax = plt.subplots(figsize=(10, 7))
	
	for k in sorted(data, reverse=False):
		if "Sun" in k:
			c="k"
			ls="-"
		elif "Halite 1" in k:
			c="r"
		elif "Halite 2" in k:
			c="b"
		else:
			c="c"
		if "iddle" in k:
			ls = "--"
		else:
			ls = "-"
		ax.plot(data[k][0], data[k][1], label=k, color=c, ls=ls)

	ax.grid(which='both', color='k', linestyle='--', alpha=0.1)
	ax.set_yscale('log')
	ax.set_xlim(500, 900)
	ax.set_ylim(0.0000001, 1.2)
	ax.legend()
	ax.set_xlabel("Wavelength (nm)")
	ax.set_ylabel("Transmission")
	plt.tight_layout()
	plt.savefig("figure.png")



data = load_data("20180410")
data = rename_labels(data)
data = standardize_data(data)
data["Sunlight"] = load_solar_data("sunlight.txt")
plot_data(data)
