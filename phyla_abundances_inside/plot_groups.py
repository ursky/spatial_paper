#!/usr/bin/env python
import sys
import numpy as np
from scipy import stats
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


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
			if "Chloroplast" in cut[1]:
				cut[1] = "D_0__Eukarya;D_1__Chlorophyta"
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
		if taxon not in data:
			data[taxon]={}
			for sample in abund[otu]:
				data[taxon][sample]=0
		for sample in abund[otu]:
			data[taxon][sample]+=abund[otu][sample]
	return data
		

def standardize_to_average( table, only_center=False ):
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
		
		
		data[taxon]={}
		for position in positions:
			data[taxon][position]=[]
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
					data[taxon][position].append(new_abund)
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


def analyze_standardized_counts (data):
	for taxon in data:
		print "\t".join([taxon, str(np.mean(data[taxon]["T"])), str(np.mean(data[taxon]["M"])), str(np.mean(data[taxon]["B"]))])
		test = stats.f_oneway(data[taxon]["T"], data[taxon]["M"], data[taxon]["B"])
		print test
		print ""


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


def draw_signifficance_bars(df, labels, ax):
	max_val,min_val = get_max_min_in_dict(df)
	inc = (max_val-min_val)*0.05

	taxon=""
	for x_st,s1 in enumerate(labels):
		for x_fi,s2 in enumerate(labels):
			if s1>=s2: continue
			taxon1=s1.split("_")[0]
			taxon2=s2.split("_")[0]
			if taxon1!=taxon2: continue

			if taxon!=taxon1:
				taxon=taxon1
				h=(max(df[s1]+df[s2]))*1.05

			test=stats.ttest_ind(df[s1], df[s2])
			if test.pvalue > 0.01: continue
			elif test.pvalue > 0.001: m='*'
			elif test.pvalue > 0.0001: m='**'
			else: m='***'
			ax.hlines(y=h, xmin=x_st, xmax=x_fi, linewidth=1, color='k')
			ax.text(float(x_fi+x_st)/2, h+0, m, ha='center', fontsize=12)
			h+=inc


def make_boxplot(df, dic):
	sns.set(rc={'figure.figsize':(8,5)})
	sns.set(style="whitegrid")

	ax = sns.boxplot(x="taxon_position", y="value", data=df, color="w")
	ax = sns.swarmplot(x="taxon_position", y="value", data=df, hue="position")
	#ax.set_yscale("log")

	labels = [item.get_text() for item in ax.get_xticklabels()]
	draw_signifficance_bars(dic, labels, ax)

	# simplify x labels
	for i,label in enumerate(labels):
		if "Mid" in label:
			labels[i]=label.split("_")[0]
		else:
			labels[i]=""
	ax.set_xticklabels(labels, rotation=30)
	for x in [3,6,9,12,15]:
		ax.axvline(x=x-0.5, linewidth=1, color='k', alpha=0.2)

	ax.set_xlabel("Taxon")
	ax.set_ylabel("Relative abundance standardized to slice average")
	
	plt.tight_layout()
	plt.savefig("figure_gouped.png")
	#plt.show()








asv_table = import_asv_abund("asv_table.txt")
cell_counts = import_cell_counts("cell_counts.txt")
std_asv_table = standardize_table(asv_table)
taxonomy = import_taxonomy("taxonomy.txt")
taxa_table = collapse_taxa ( std_asv_table, taxonomy, level=2 )
print taxa_table["Eukarya;Chlorophyta"]

# optional standardize to average cell counts in each replicate
#taxa_table = standardize_to_cell_counts(taxa_table, cell_counts)

std_data = standardize_to_average(taxa_table, only_center=False)
analyze_standardized_counts(std_data)
#print_data(std_data)

df, dic = convert_to_df(std_data)

make_boxplot(df, dic)







