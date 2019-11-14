#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import random
from scipy import stats

df_dic={"labels":[],"values":[]}


# load large scale metrics
large_site={}
for line in open("metadata_large.txt"):
	if line[0]=="#":
		continue
	cut=line.split("\t")
	sample=cut[0].split("_")[0]
	site=cut[4]
	large_site[sample]=site

for i,line in enumerate(open("dissimilarity_large.txt")):
	cut=line.rstrip().split("\t")
	if i==0:
		head=cut
		continue
	for j,val in enumerate(cut):
		if j==0:
			sample1=val
			continue
		sample2=head[j]
		if sample1<=sample2:
			continue
		val=float(val)
		site1=large_site[sample1]
		site2=large_site[sample2]
		if site1!=site2:
			df_dic["labels"].append("large")
			df_dic["values"].append(val)
		
	
# load medium scale metrics
medium_site={}
for line in open("metadata_medium.txt"):
	if line[0]=="#":
		continue
	cut=line.split("\t")
	sample=cut[0].split("_")[0]
	site=cut[1]
	medium_site[sample]=site

for i,line in enumerate(open("dissimilarity_medium.txt")):
	cut=line.rstrip().split("\t")
	if i==0:
		head=cut
		continue
	for j,val in enumerate(cut):
		if j==0:
			sample1=val
			continue
		sample2=head[j]
		if sample1<=sample2:
			continue
		val=float(val)
		site1=medium_site[sample1]
		site2=medium_site[sample2]
		if site1!=site2:
			df_dic["labels"].append("medium")
			df_dic["values"].append(val)	


# load small scale metrics
small_site={}
for line in open("metadata_small.txt"):
	if line[0]=="#":
		continue
	cut=line.split("\t")
	sample=cut[0].split("_")[0]
	site=cut[1]
	small_site[sample]=site

for i,line in enumerate(open("dissimilarity_small.txt")):
	cut=line.rstrip().split("\t")
	if i==0:
		head=cut
		continue
	for j,val in enumerate(cut):
		if j==0:
			sample1=val
			continue
		sample2=head[j]
		if sample1<=sample2:
			continue
		val=float(val)
		site1=small_site[sample1]
		site2=small_site[sample2]
		slice1=sample1.split("-")[1]
		slice2=sample2.split("-")[1]
		pos1=sample1.split("-")[2]
		pos2=sample2.split("-")[2]
		if site1!=site2:
			df_dic["labels"].append("inter-halite")
			df_dic["values"].append(val)
		else:
			df_dic["labels"].append("intra-halite")
			df_dic["values"].append(val)
		#	if slice1==slice2:
		#		df_dic["labels"].append("intra-slice")
		#		df_dic["values"].append(val)
		#	else:
		#		df_dic["labels"].append("inter-slice")
		#		df_dic["values"].append(val)
			if pos1==pos2:
				df_dic["labels"].append("same-position")
				df_dic["values"].append(val)
		#	else:
		#		df_dic["labels"].append("inter-position")
		#		df_dic["values"].append(val)

data={}
labels=df_dic["labels"]
values=df_dic["values"]
for i,val in enumerate(values):
	label=labels[i]
	if label not in data:
		data[label]=[]
	data[label].append(val)

df = pd.DataFrame.from_dict(df_dic)

ax = sns.boxplot(x="labels", y="values", data=df, color="gray")
#ax = sns.violinplot(x="labels", y="values", data=df)
#ax = sns.swarmplot(x="labels", y="values", data=df)

#add signifficance bars
labels = [item.get_text() for item in ax.get_xticklabels()]
print labels
h=0.7
for x_st,s1 in enumerate(labels):
	for x_fi,s2 in enumerate(labels):
		if s1>=s2: continue
		if abs(x_st-x_fi)>1: continue
		test=stats.ttest_ind(data[s1], data[s2])
		if test.pvalue > 0.05: continue
		elif test.pvalue > 0.01: m='*'
		elif test.pvalue > 0.001: m='**'
		else: m='***'
		ax.hlines(y=h, xmin=x_st, xmax=x_fi, linewidth=1, color='k')
		ax.text(float(x_fi+x_st)/2, h+0, m, ha='center', fontsize=12)
		h-=0.05





plt.ylabel("Bray-Curtis Dissimilarity")
plt.xlabel("Spatial Scale")
plt.xticks(rotation=45)

plt.tight_layout()
plt.savefig("figure_bray_curtis.png", dpi=300)
