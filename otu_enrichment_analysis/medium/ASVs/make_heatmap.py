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
for line in open("medium_metadata.txt"):
	if line[0]=="#":
		continue
	cut=line.strip().split("\t")
	sample=cut[0].split("_")[0]
	site=cut[1]
	sites[sample]=site
	if site=="Bottom":
		site1.add(sample)
	else:
		site2.add(sample)



print "loading abundance data..."
df=pd.read_csv("medium_otu_table.txt", sep='\t', index_col=0, header=1)

# rename columns
for col in df:
	sample=col.split("_")[0]
	new_name = sample
	if col!=sample:
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



print "adding color code..."

col_colors=[]
for sample in df.columns.values:
	site = sample.split("_")[0]
	if site=="SG1": col_colors.append('r')
	elif site=="SG2": col_colors.append('b')
	else: col_colors.append('w')

print "plotting..."
sns.set(font_scale=1)
g = sns.clustermap(df, figsize=(8,8), col_colors=col_colors, col_cluster=False, yticklabels=False, cmap="magma")

# adjust axis labels
plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
plt.setp(g.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)



plt.savefig("figure.png", bbox_inches='tight', dpi=300)


