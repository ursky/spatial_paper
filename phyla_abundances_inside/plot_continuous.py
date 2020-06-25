#!/usr/bin/env python
import sys
import numpy as np
from scipy import stats
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods


def load_distances(filename):
	distances={}
	for line in open(filename):
		cut=line.strip().split("\t")
		if cut[0]=="Sample":
			continue
		distances[cut[0]] = float(cut[-1])
	return distances


def import_asv_abund(file_name):
	data={}
	for i,line in enumerate(open(file_name)):
		cut=line.strip().split("\t")
		if i==0:
			continue
		elif i==1:
			head = cut[:]
		else:
			otu=cut[0]
			data[otu]={}
			for i,val in enumerate(cut):
				if i==0:
					continue
				sample = head[i]
				data[otu][sample] = float(val)
	return data


def standardize_table(data):
	sums={}
	for otu in data:
		for sample in data[otu]:
			if sample not in sums:
				sums[sample]=0
			val = data[otu][sample]
			sums[sample]+=val
	for otu in data:
		for sample in data[otu]:
			val = data[otu][sample]
			new_val = val * 100.0 / sums[sample]
			data[otu][sample] = new_val
	return data


def import_taxonomy(file_name):
	data={}
	for i,line in enumerate(open(file_name)):
		if i==0:
			continue
		else:
			cut=line.strip().split("\t")
			full_taxon = cut[1].split(";")
			taxon=[]
			for t in full_taxon:
				if t.split("__")[-1]=="":
					break
				taxon.append(t)
			for i in range(len(taxon)):
				taxon[i] = taxon[i].split("__")[-1]
			data[cut[0]] = taxon
	return data


def collapse_taxa ( abund, taxa, level=2 ):
	data={}
	for otu in taxa:
		taxon = taxa[otu]
		if len(taxon)<level:
			continue
		taxon = ";".join(taxon[:level])
		if "Plancto" in taxon:
			continue
		if taxon not in data:
			data[taxon]={}
			for sample in abund[otu]:
				data[taxon][sample]=0
		for sample in abund[otu]:
			data[taxon][sample]+=abund[otu][sample]
	return data
		

def standardize_to_average( table, distances, only_center=False ):
	data={}
	halites=["H2", "H3", "H4", "H5", "H6", "H7"]
	slices=["A", "B", "C"]
	positions=["T", "M", "B"]

	for taxon in table:
		# get average abundance in middle
		tot=0; ct=0
		for sample in table[taxon]:
			tot+=table[taxon][sample]
			ct+=1
		total_average = tot/ct
		if total_average < 0.1: 
			continue
		
		
		data[taxon]=([], [], [])
		for position in positions:
			for halite in halites:
				for sli in slices:
					if only_center==True:
						if halite=="H2":
							if sli!="A":
								continue
						elif halite=="H3" or halite=="H4" or halite=="H5":
							if sli!="B":
								continue
						else:
							continue
					sample = "-".join([halite, sli, position])+"_ASD"
					if sample not in table[taxon]:
						continue
					if sample=="H5-A-M_ASD":
						continue
					if "-".join([halite, sli, "T"])+"_ASD" not in table[taxon] or "-".join([halite, sli, "M"])+"_ASD" not in table[taxon] or "-".join([halite, sli, "B"])+"_ASD" not in table[taxon]:
						continue

					abund = table[taxon][sample]
					ave_slice_abund = compute_ave_abund_in_slice(table[taxon], halite, sli)

					if ave_slice_abund==0:
						continue
					#new_abund = total_average * abund/ave_slice_abund
					new_abund = abund/ave_slice_abund
					data[taxon][1].append(new_abund)
					distance = distances["_".join(sample.split("_")[:-1])]
					data[taxon][0].append(distance)
					data[taxon][2].append(position)
	return data


def compute_ave_abund_in_slice (data, halite, sli):
	sample_slice = halite+"-"+sli
	tot=0
	ct=0
	for sample in data:
		if sample.startswith(sample_slice):
			tot+=data[sample]
			ct+=1
	return tot/ct


def print_data ( data ):
	positions=["T","M","B"]
	for taxon in data:
		print ""
		print taxon
		#print len(data[taxon]["T"]), len(data[taxon]["M"]), len(data[taxon]["B"])
		print "Top\tMiddle\tBottom"
		for i in range(100):
			line=[]
			for position in positions:
				if i<len(data[taxon][position]):
					line.append(str(data[taxon][position][i]))
				else:
					line.append("")	
			if line[0]=="" and line[1]=="" and line[2]=="": 
				break
			print "\t".join(line)
					

def stats_tests( data ):
	for taxon in data:
		list1=data[taxon][0]
		list2=data[taxon][1]
		ratios=[]
		for i in range(len(list1)):
			if list2[i]==0:
				ratios.append(1)
			else:
				ratios.append(list1[i]/list2[i])
		average=np.mean(list1+list2)
		ratio=np.mean(ratios)
		if average==0:
			continue
		test = stats.ttest_rel(list1, list2)
		if test.pvalue>1: 
			continue
		print "\t".join([taxon, str(average), str(ratio), str(test.pvalue)])


def import_cell_counts(filename):
	data={}
	positions={"top":"T", "bot":"B", "mid":"M"}
	for line in open(filename):
		cut=line.strip().split("\t")
		name = cut[0].split("-")
		name[2]=positions[name[2]]
		name = "-".join(name)
		data[name+"_ASD"]=float(cut[1])
	return data
	

def standardize_to_cell_counts(taxa_table, cell_counts):
	out_table={}
	for taxon in taxa_table:
		out_table[taxon]={}
		for sample in taxa_table[taxon]:
			if sample not in cell_counts:
				continue
			else:
				out_table[taxon][sample] = taxa_table[taxon][sample] * cell_counts[sample]

	return out_table


def convert_to_df(std_data):
	taxa = ["Archaea;Euryarchaeota", "Bacteria;Bacteroidetes", "Bacteria;Cyanobacteria", "Bacteria;Actinobacteria", "Archaea;Nanohaloarchaeota", "Bacteria;Proteobacteria"]
	positions={"T":"Top", "B":"Bot", "M":"Mid"}

	df = pd.DataFrame(columns=["taxon", "position", "taxon_position", "value"])
	dic = {}
	for taxon in taxa:
		tax = taxon.split(";")[1]
		for position in std_data[taxon]:
			dic[tax+"_"+positions[position]]=[]
			for i,val in enumerate(std_data[taxon][position]):
				ID = "_".join([tax, position, str(i)])

				df = df.append({"taxon":tax, 
					"position":positions[position], 
					"taxon_position":tax+"_"+positions[position], 
					"value":val}, ignore_index=True)
				dic[tax+"_"+positions[position]].append(val)
	return df, dic
				

def get_max_min_in_dict(dictionary):
        maxs=[]; mins=[]
        for key in dictionary:
                maxs.append(max(dictionary[key]))
                mins.append(min(dictionary[key]))
        return max(maxs), min(mins)



def make_plot(data):
	fig, axs = plt.subplots(2,3, figsize=(8,6), sharey=True)
	for i,taxon in enumerate(data):
		if i<=2:
			ax = axs[0][i]
		else:
			ax = axs[1][i-3]
		xs = data[taxon][0]
		ys = data[taxon][1]
		positions = data[taxon][2]
		colors = []
		for position in positions:
			if position=="T": colors.append("red")
			if position=="M": colors.append("magenta")
			if position=="B": colors.append("cyan")
		ax.scatter(xs, ys, alpha=0.6, c=colors, edgecolor='k')
	
		# line of best fit
		grid = np.r_[0.5:5:512j]
		k0 = smooth.NonParamRegression(xs, ys, method=npr_methods.SpatialAverage())
		k0.fit()
		ax.plot(grid, k0(grid), "k", linewidth=2)

		test = stats.pearsonr(data[taxon][0], data[taxon][1])
		print taxon, test

		ax.set_title(taxon.split(";")[-1], fontsize=12)
		ax.grid(ls="--")
		ax.set_xticks([0,1,2,3,4,5])
	#	if i!=0:
	#		ax.set_yticklabels([]) 
	
	fig.add_subplot(111, frameon=False)
	plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
	plt.xlabel("Distance to surface (cm)")
	plt.ylabel("Standardized abundance of taxon")

	plt.tight_layout()
	plt.savefig("figure.png", dpi=300)





####################  MAIN  ####################

distances = load_distances("distance_metadata.txt")
asv_table = import_asv_abund("asv_table.txt")
cell_counts = import_cell_counts("cell_counts.txt")
std_asv_table = standardize_table(asv_table)
taxonomy = import_taxonomy("taxonomy.txt")
taxa_table = collapse_taxa ( std_asv_table, taxonomy, level=2 )

# optional standardize to average cell counts in each replicate
# taxa_table = standardize_to_cell_counts(taxa_table, cell_counts)

std_data = standardize_to_average(taxa_table, distances, only_center=False)
make_plot(std_data)







