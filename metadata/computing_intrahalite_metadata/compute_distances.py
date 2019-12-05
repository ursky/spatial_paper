#!/usr/bin/env python
import sys
import numpy as np

# load vertical positions
vertical_position={}
for line in open("vertical_component.txt"):
	if line[0]=="#": 
		continue
	cut=line.strip().split("\t")

	halite=cut[0]
	slice=cut[1]
	position=cut[2]
	distance=float(cut[3])
	
	name = "Halite"+halite+"_"+slice+"_"+position
	vertical_position[name]=distance


# calculate vertical distances
vertical_distance={}
for sample1 in vertical_position:
	for sample2 in vertical_position:
		if "Total" in sample1 or "Total" in sample2:
			continue
		cut1=sample1.split("_")
		cut2=sample2.split("_")
		if sample1==sample2:
			v_distance=0
		elif cut1[0]!=cut2[0]:
			v_distance=0
		else:
			total_hight1 = vertical_position["_".join(cut1[:-1]+["Total"])]
			total_hight2 = vertical_position["_".join(cut2[:-1]+["Total"])]
			distance1 = vertical_position[sample1]
			distance2 = vertical_position[sample2]
			distance_to_bottom1 = total_hight1 - distance1
			distance_to_bottom2 = total_hight2 - distance2
			v_distance = max(distance_to_bottom1-distance_to_bottom2, distance_to_bottom2-distance_to_bottom1)
		if sample1 not in vertical_distance:
			vertical_distance[sample1] = {}
		vertical_distance[sample1][sample2] = v_distance


# load horisontal positions
slice_distances={}
for line in open("horisontal_component.txt"):
	if line[0]=="#":
		continue
	cut=line.strip().split("\t")
	halite="Halite"+cut[0].split()[-1]
	order=cut[1]
	left_wing=float(cut[2])
	right_wing=float(cut[3])
	left_slice=float(cut[4])
	right_slice=float(cut[5])

	# compute distances between 3 possible pairs of slices
	slice1=halite+"_"+order[0]
	slice2=halite+"_"+order[1]
	slice3=halite+"_"+order[2]

	slice_distances[slice1]={}
	slice_distances[slice2]={}
	slice_distances[slice3]={}

	distance=left_slice
	slice_distances[slice1][slice2]=distance
	slice_distances[slice2][slice1]=distance

	distance=right_slice
	slice_distances[slice2][slice3]=distance
	slice_distances[slice3][slice2]=distance

	distance=left_slice+right_slice
	slice_distances[slice1][slice3]=distance
	slice_distances[slice3][slice1]=distance



# calculate 2d distances
distances={}
for sample1 in vertical_distance:
	for sample2 in vertical_distance:
		cut1=sample1.split("_")
		cut2=sample2.split("_")
		halite1=cut1[0]
		halite2=cut2[0]
		slice1=cut1[1]
		slice2=cut2[1]
	
		v_distance=vertical_distance[sample1][sample2]
		if halite1==halite2:
			if slice1==slice2:
				h_distance=0
			else:
				h_distance=slice_distances[halite1+"_"+slice1][halite2+"_"+slice2]
		else:
			h_distance=100

		distance = np.sqrt(v_distance*v_distance + h_distance*h_distance)
		if sample1 not in distances:
			distances[sample1]={}
		distances[sample1][sample2]=distance


renamed={}
for k1 in distances:
	cut1=k1.split("_")
	cut1[0]="H"+cut1[0][-1]
	cut1[2]=cut1[2][0]+"_ASD"
	k1_re = "-".join(cut1)
	renamed[k1_re]={}
	for k2 in distances:
		cut2=k2.split("_")
		cut2[0]="H"+cut2[0][-1]
		cut2[2]=cut2[2][0]+"_ASD"
		k2_re = "-".join(cut2)
		renamed[k1_re][k2_re]=distances[k1][k2]
distances=renamed



if len(sys.argv)>1:
	halite=sys.argv[1]
else:
	halite="any"



cut=[""]
for k1 in sorted(distances):
	if halite=="any" or k1.split("-")[0]==halite:
		cut.append(k1)
print "\t".join(cut)


for k1 in sorted(distances):
	if halite=="any" or k1.split("-")[0]==halite:
		line=k1
		for k2 in sorted(distances):
			if halite=="any" or k2.split("-")[0]==halite:
				line=line+"\t"+str(distances[k1][k2])[:5]
		print line




















