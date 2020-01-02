#!/usr/bin/env python2
import sys
from scipy import stats
import numpy as np

data = {}
for line in open("alpha-diversity.tsv"):
	if "faith" in line:
		continue
	cut = line.strip().split("\t")
	sample = cut[0].split("_")[0]
	faith = float(cut[1])
	sli = "-".join(sample.split("-")[:-1])
	position = sample.split("-")[-1]

	if sli not in data:
		data[sli]={}
	data[sli][position]=faith


def test_difference(s1, s2):
	site1=[]
	site2=[]

	for sli in data:
		if s1 in data[sli] and s2 in data[sli]:
			site1.append(data[sli][s1])
			site2.append(data[sli][s2])
	print "\n%s vs %s paired test:" % (s1, s2)
	print stats.ttest_rel(site1, site2)

	diffs=[]
	for i, A in enumerate(site1):
		diffs.append( A-site2[i] )
	print "The average diversity difference between %s and %s is %f" % (s1, s2, np.mean(diffs))


test_difference("T", "M")
test_difference("B", "M")
test_difference("T", "B")
		


