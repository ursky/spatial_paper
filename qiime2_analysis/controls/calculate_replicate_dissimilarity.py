#!/usr/bin/env python2
import sys
import numpy as np

def get_sample(name):
	sample = name.strip().split("_")[0]
	return name, "-".join(sample.split("-")[1:3])

diss = []
for line in open (sys.argv[1]):
	cut = line.split("\t")
	if cut[0]=="":
		header = cut
		continue
	id1, sample1 = get_sample(cut[0])
	for i in range(1, len(cut)):
		val = float(cut[i])
		id2, sample2 = get_sample(header[i])
		id1=id1.strip()
		id2=id2.strip()
		if "MedCtrl" not in id1 or "MedCtrl" not in id2:
			continue
		if id1>=id2:
			continue
		if sample1==sample2:
			diss.append(val)
			print id1, id2, val
print diss
print np.mean(diss)
print np.std(diss)
	
