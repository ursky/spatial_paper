#!/usr/bin/env python2.7
print "loading libs..."
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats
import random

print "loading metadata..."
sites={}
site1=set(); site2=set()
for line in open("large_metadata.txt"):
	if line[0]=="#":
		continue
	cut=line.strip().split("\t")
	sample=cut[0].split("_")[0]
	site=cut[4]
	sites[sample]=site
	if site=="SG1":
		site1.add(sample)
	else:
		site2.add(sample)



print "loading abundance data..."
df=pd.read_csv("large_otu_table.txt", sep='\t', index_col=0, header=1)


# store away the taxonomy info
taxa = df["taxonomy"]
df.drop(["taxonomy"], inplace=True, axis=1)


# rename columns
sg1=1
sg2=1
for col in df:
	if col in site1:
		new_name = "SG1_"+str(sg1)
		sg1+=1
	elif col in site2:
		new_name = "SG2_"+str(sg2)
		sg2+=1
	df[new_name]=df[col]
	df.drop([col], inplace=True, axis=1)

df = df.reindex(sorted(df.columns), axis=1)




print "standardizing data..."
# standardize columns by total sum in each column
df = df.div(df.sum(axis=0), axis=1)
df=1000000*df

# remove rows with low total
df = df[(df.T != 0).any()]
df['TOTAL'] = df.sum(axis=1)
df = df[df['TOTAL'] > 1000]
df.drop(["TOTAL"], inplace=True, axis=1)

# log standardize:
df+=0.01; df=np.log(df)

# standardize rows by maximum value in each row
#df = df.div(df.max(axis=1), axis=0)



print "computing distance statistics..."
correlations=[]
for sample1 in df.columns.values:
	for sample2 in df.columns.values:
		if sample1==sample2:
			continue
		if sample1.startswith("SG1") and sample2.startswith("SG2"):
			l1 = df[sample1].values
			l2 = df[sample2].values
			corr, pval = stats.spearmanr(l1, l2)
			#print corr
			correlations.append(corr)
print "mean correlation between between SG1 and SG2 OTUs is ", np.mean(correlations)
		



print "calculating taxon enrichment"
if len(sys.argv)>1:
        toi=sys.argv[1]
else:
        toi=None

all_taxa=[]
select_taxa=[]
target_ct=0
for index in df.index.values:
	if toi==None:
		continue
        taxon_all = taxa[index].split(';')
        if len(taxon_all)<2:
                taxon="unknown"
        else:
                taxon = taxa[index].split(';')[1].split("_")[-1]

	set1 = []
	set2 = []
	for sample in df.columns.values:
		site = sample.split("_")[0]
		val = df.loc[index, sample]
		if val<0:
			continue
		if site=="SG1":
			set1.append(val)
		else:
			set2.append(val)
	if len(set1)==0:
		overrepresentation=100
	else:
		overrepresentation=len(set2)/len(set1)
	if overrepresentation>3 and len(set2)>=5:
		if toi in taxon:
			target_ct+=1
		select_taxa.append(taxon)
	all_taxa.append(taxon)

if toi!=None:
	passed_ct=0
	for i in range(10000):
		ct = 0
		rand_subset = random.sample(all_taxa, k=len(select_taxa))
		for taxon in rand_subset:
			if toi in taxon:
				ct+=1
		if ct>=target_ct:
			passed_ct+=1
	print passed_ct*1.0/10000

	





print "adding color code..."
taxa_colors = {"Bacteroidetes":"r", "Cyanobacteria":"g", "Euryarchaeota":"y", "Nanohaloarchaeota":"b", "Proteobacteria":"m"}

col_colors=[]
for sample in df.columns.values:
	if "SG1" in sample: col_colors.append('r')
	elif "SG2" in sample: col_colors.append('b')
	else: col_colors.append('w')
row_colors=[]

if len(sys.argv)>1:
	toi=sys.argv[1]
else:
	toi=None

for index in df.index.values:
	taxon_all = taxa[index].split(';')
	if len(taxon_all)<2:
		taxon="unknown"
	else:
		taxon = taxa[index].split(';')[1].split("_")[-1]
	if toi!=None:
		if toi not in taxon:
			df.drop([index], inplace=True, axis=0)
			continue


	if taxon in taxa_colors:
		row_colors.append(taxa_colors[taxon])
	else:
		row_colors.append("w")



print "plotting..."
sns.set(font_scale=1)
#g = sns.clustermap(df, figsize=(8,8), col_colors=col_colors, row_colors=row_colors, col_cluster=False, yticklabels=False, cmap="magma")
g = sns.clustermap(df, figsize=(8,8), col_colors=col_colors, col_cluster=False, yticklabels=False, cmap="magma")

# adjust axis labels
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)



plt.savefig("figure.png", bbox_inches='tight', dpi=300)


