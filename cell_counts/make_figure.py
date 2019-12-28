#!/usr/bin/env python
import sys
import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

positions={"T":"top", "M":"mid", "B":"bot"}
data={}
data_df={}
data_df["halite"]={}
data_df["slice"]={}
data_df["position"]={}
data_df["replicate"]={}
data_df["membrane"]={}
data_df["fov"]={}
data_df["count"]={}


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
	
	name = "-".join([halite, sli, position, membrane, fov])
	data_df["halite"][name] = halite
	data_df["slice"][name] = "-".join([halite, sli])
	data_df["replicate"][name] = "-".join([halite, sli, position])
	data_df["position"][name] = position
	data_df["membrane"][name] = "-".join([halite, sli, position, membrane])
	data_df["fov"][name] = fov
	data_df["count"][name] = count_to_density(count)


	if halite not in data:
		data[halite]={}
	if sli not in data[halite]:
		data[halite][sli]={}
	if position not in data[halite][sli]:
		data[halite][sli][position]={}
	if membrane not in data[halite][sli][position]:
		data[halite][sli][position][membrane]=[]
	data[halite][sli][position][membrane].append(count)


df = pd.DataFrame.from_dict(data_df)

# remove outliers: counts that are more than 1 stdev away from the mean in a given biological replicate
for halite in data:
	for sli in data[halite]:
		for position in data[halite][sli]:
			counts=[]
			replicate="-".join([halite, sli, position])
			for membrane in data[halite][sli][position]:
				counts+=data[halite][sli][position][membrane]
			mean = np.mean(counts)
			std = np.std(counts)*1
			print replicate +"\t"+ str(mean)
			for membrane in data[halite][sli][position]:
				for i,ct in enumerate(data[halite][sli][position][membrane]):
					name = "-".join([halite, sli, position, membrane, str(i+1)])
					if ct<mean-std or ct>mean+std:
						df = df.drop([name])
						


if len(sys.argv)>1:
	category = sys.argv[1]
else:
	category = "slice"

sns.set(rc={'figure.figsize':(8,6)})
sns.set(style="whitegrid")


ax = sns.boxplot(x=category, y="count", data=df, color="w")
ax = sns.swarmplot(x=category, y="count", data=df, hue="position", palette=sns.xkcd_palette(["cyan","magenta","gold"]))

ax.set_ylabel("Million cells per gramm of halite")
ax.set_xlabel("Slice")


plt.tight_layout()
plt.savefig("figure.png", dpi=300)
#plt.show()



















