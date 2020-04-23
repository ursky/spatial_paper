#!/usr/bin/env python
import sys
from scipy import stats


def paired_ttest(data, A, B):
	a=[]; b=[]
	for key in data:
		if A not in data[key] or B not in data[key]:
			continue
		a.append(data[key][A])
		b.append(data[key][B])
	stat, pval = stats.ttest_rel(a, b)
	print "%s vs %s\t%s\t%s" %(A, B, stat, pval)
	return stat, pval

def anova_test(data, A="T", B="M", C="B"):
	a=[]; b=[]; c=[]
	for key in data:
		if A not in data[key] or B not in data[key] or C not in data[key]:
			continue
		a.append(data[key][A])
		b.append(data[key][B])
		c.append(data[key][C])
	stat, pval = stats.f_oneway(a,b,c)
	print "Anova test: %s, pval=%s" %(stat, pval)
	return stat, pval


####################################

slices = {}

for line in open(sys.argv[1]):
	cut = line.strip().split("\t")
	if len(cut)!=2:
		continue
	val = float(cut[1])
	sample = cut[0].split("_")[0]
	sli = "-".join(sample.split("-")[:-1])
	pos = sample.split("-")[-1]
	if sli not in slices:
		slices[sli]={}
	slices[sli][pos] = val

#anova_test(slices)
paired_ttest(slices, "T", "M")
paired_ttest(slices, "B", "M")
paired_ttest(slices, "T", "B")
