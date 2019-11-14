#!/usr/bin/env python2.7
print "loading libs..."
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

print "loading metadata..."
samples={}
for line in open("metadata.txt"):
	if line[0]=="#":
		continue
	cut=line.strip().split("\t")
	sample=cut[0].split("_")[0]
	halite=cut[1]
	samples[sample]=halite



print "loading abundance data..."
df=pd.read_csv("small_otu_table.txt", sep='\t', index_col=0, header=1)

# store away the taxonomy info
taxa = df["taxonomy"]
df.drop(["taxonomy"], inplace=True, axis=1)

# rename columns
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
intra_correlations=[]
inter_correlations=[]
for sample1 in df.columns.values:
	for sample2 in df.columns.values:
		if sample1==sample2:
			continue
		halite1=samples[sample1]
		halite2=samples[sample2]
		l1 = df[sample1].values
		l2 = df[sample2].values
		corr, pval = stats.spearmanr(l1, l2)
		if halite1==halite2:
			intra_correlations.append(corr)
		else:
			inter_correlations.append(corr)
			#print corr
print "Inter-nodule correlation is", np.mean(inter_correlations)
print "Intra-nodule correlation is", np.mean(intra_correlations)



print "adding color code..."
taxa_colors = {"Bacteroidetes":"r", "Cyanobacteria":"g", "Euryarchaeota":"y", "Nanohaloarchaeota":"b", "Proteobacteria":"m"}
halite_colors = {"H2":"y", "H3":"r", "H4":"g", "H5":"b", "H6":"c", "H7":"m"}

col_colors=[]
for sample in df.columns.values:
	halite = sample.split("-")[0]
	if halite in halite_colors:
		col_colors.append(halite_colors[halite])
	else:
		col_colors.append('w')
row_colors=[]

# see if needed to filter by taxa
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
g = sns.clustermap(df, figsize=(8,8), col_colors=col_colors, row_colors=row_colors, col_cluster=False, yticklabels=False, cmap="magma")

# adjust axis labels
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)



plt.savefig("figure.png", bbox_inches='tight', dpi=300)


