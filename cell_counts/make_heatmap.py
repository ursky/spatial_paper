#!/usr/bin/env python
import sys
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

positions={"T":"Top", "M":"Middle", "B":"Bottom"}
data={}
data_df={}
data_df["Top"]={}
data_df["Middle"]={}
data_df["Bottom"]={}


def count_to_density(count):
	eFOVs = {100:0.031, 40:0.203, 10:3.198}
	eFOV = eFOVs[40]
	filter_area = 226.98
	counts_per_filter = count * 226.98 / eFOV
	counts_per_gram = counts_per_filter * (6.5/2) / 0.5
	return counts_per_gram / 1000000


for line in open("cell_counts.txt"):
	cut = line.strip().split("\t")
	if cut[0]=="Membrane-FOV":
		continue
	count = int(cut[1])
	halite = cut[0].split("-")[0]
	sli = cut[0].split("-")[1][0]
	position = cut[0].split("-")[1][-1]
	position = positions[position]
	membrane = cut[0].split("-")[2][0]
	fov = cut[0].split("-")[2][2]
	
	name = "-".join([halite, sli])
	if name not in data_df[position]:
		data_df[position][name] = []
	data_df[position][name].append( count_to_density(count) )




# remove outliers: counts that are more than 1 stdev away from the mean in a given biological replicate
for position in data_df:
	for sli in data_df[position]:
		counts = data_df[position][sli]
		mean = np.mean(counts)
		std = np.std(counts)*1
		new_counts = []
		for ct in counts:
			if ct>mean-std and ct<mean+std:
				new_counts.append(ct)
		data_df[position][sli] = np.mean(new_counts)


df = pd.DataFrame.from_dict(data_df)
print df.T
df.T.to_csv("cell_count_table.tab", sep="\t")

df = df.div(df.max(axis=1), axis=0)
sns.clustermap(df.T, figsize=(6,5), cmap="magma_r", row_cluster=False)


plt.savefig("heatmap.png", dpi=300, bbox_inches='tight')
#plt.show()



















