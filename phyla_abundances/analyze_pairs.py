#!/usr/bin/env python
# ./analyze.py Top Bottom 3
import sys
import numpy as np
from scipy import stats


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
		if taxon not in data:
			data[taxon]={}
			for sample in abund[otu]:
				data[taxon][sample]=0
		for sample in abund[otu]:
			data[taxon][sample]+=abund[otu][sample]
	return data
		

def get_paired_abundances( table, pos1, pos2 ):
	data={}
	halites=["H2", "H3", "H4", "H5", "H6", "H7"]
	slices=["A", "B", "C"]
	pos1=pos1[0]
	pos2=pos2[0]

	for taxon in table:
		data[taxon]=[[],[]]
		for halite in halites:
			for sli in slices:
				sample1 = "-".join([halite, sli, pos1])+"_ASD"
				sample2 = "-".join([halite, sli, pos2])+"_ASD"
				if sample1 not in table[taxon] or sample2 not in table[taxon]:
					continue
				abund1 = table[taxon][sample1]
				abund2 = table[taxon][sample2]
				data[taxon][0].append(abund1)
				data[taxon][1].append(abund2)
	return data
					

def stats_tests( data, pos1, pos2 ):
	print "\t".join(["Taxon", "Abundance", pos1+":"+pos2+" ratio", "pvalue"])
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


position1=sys.argv[1]
position2=sys.argv[2]

abund = import_asv_abund("asv_table.txt")
std_abund = standardize_table(abund)
taxonomy = import_taxonomy("taxonomy.txt")
taxa_table = collapse_taxa ( std_abund, taxonomy, level=int(sys.argv[3]) )

paired_data = get_paired_abundances(taxa_table, position1, position2)
stats_tests(paired_data, position1, position2)













