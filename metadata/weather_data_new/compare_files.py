#!/usr/bin/env python 
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats


def convert_time(time):
	cut = time.split()
	day_cut=cut[0].split("/")
	time_cut = cut[1].split(":")
	hour = int(time_cut[0])
	if int(time_cut[1]) > 30: 
		hour+=1
		if hour==13 and cut[2]=="AM":
			hour=1
		elif hour==13 and cut[2]=="PM":
			cut[2]="AM"
		elif hour==12 and cut[2]=="PM":
			cut[2]="AM"
	if cut[2]=="PM" and hour!=12:
		hour+=12
	if cut[2]=="AM" and hour==12:
		hour=0
	if hour<10:
		hour="0"+str(hour)
	else:
		hour=str(hour)
	if hour=="25":
		hour="01"
		day_cut[1] = str(int(day_cut[1])+1)
	day = "-".join([   day_cut[2], day_cut[0], day_cut[1]  ])
	return day, hour


def time_to_int(time):
	out=[]
	ct=0
	for point in time:
		cut = point.split("_")
		date_cut = cut[0].split("-")
		year=int(date_cut[0]) *1000000
		month=int(date_cut[1]) *10000
		day=int(date_cut[2]) *100
		hour=int(cut[1]) *1
		int_time = year+month+day+hour
		#out.append(int_time)
		out.append(ct)
		ct+=1
	return out


def load_data(filename):
	data={}
	farenheit=False
	for i,line in enumerate(open(filename)):
		if i<2: continue
		cut=line.strip().split(',')
		day,hour = convert_time(cut[1])
		temp = float(cut[2])
		if temp>50:
			farenheit=True
		if farenheit==True:
			temp = (temp-32) * (5.0/9.0)
		humi = float(cut[3])

		if day not in data:
			data[day]={}
		if hour in data[day]:
			temp = (temp+data[day][hour][0])/2
			humi = (humi+data[day][hour][1])/2
		data[day][hour]=(temp, humi)	
	return data
	

def pair_data(set1, set2, position):
	out1=[]
	out2=[]
	time_out=[]
	for i,day in enumerate(sorted(set1)):
		if day not in set2:
			continue
		#if i<3 or i>len(set1)-2:
		#	continue
		for hour in sorted(set1[day]):
			if hour not in set2[day]:
				continue
			point1 = set1[day][hour][position]
			point2 = set2[day][hour][position]
			time = day+"_"+hour
			out1.append(point1)
			out2.append(point2)
			time_out.append(time)
	if time_out==[]:
		print "No matching times found. Exiting..."
		quit()
	return time_out, out1, out2

def pair_dates(set1, set2):
	out={}
	for day in sorted(set1):
		if day not in set2:
			continue
		for hour in sorted(set1[day]):
			if hour not in set2[day]:
				continue
			point1 = set1[day][hour]
			point2 = set2[day][hour]
			if day not in out:
				out[day]={}
			out[day][hour]=[point1, point2]
	return out


def add_diagonal(ax):
	lims = [
	    np.min([ax.get_xlim(), ax.get_ylim()]),  # min of both axes
	    np.max([ax.get_xlim(), ax.get_ylim()]),  # max of both axes
	]
	ax.plot(lims, lims, 'b-', zorder=0)
	ax.set_aspect('equal')
	ax.set_xlim(lims)
	ax.set_ylim(lims)


def draw_data_vs_data(ax, data1, data2, title, alpha, time_points):
	print "plotting data vs data ("+title+")"
	time=[]
	for point in time_points:
		time.append(int(point.split("_")[1]))
	sc = ax.scatter(data1, data2, alpha=alpha, s=10, c=time, cmap="hsv")
	ax.set_xlabel(sys.argv[1].split(".")[0].split("/")[-1])
	ax.set_ylabel(sys.argv[2].split("/")[-1].split(".")[0])
	ax.set_title(title)
	ax.grid(linestyle='--', alpha=0.5)
	fig.colorbar(sc, ax=ax)
	add_diagonal(ax)


def draw_data_vs_time(ax, time, data1, data2, title, alpha):
	print "plotting data vs time ("+title+")"
	time = time_to_int(time)
	ax.scatter(time, data1, alpha=alpha, s=10, label=sys.argv[1].split(".")[0])
	ax.scatter(time, data2, alpha=alpha, s=10, label=sys.argv[2].split(".")[0])
	ax.legend()
	ax.set_ylabel(title)
	ax.set_title(title)
	ax.set_xlabel("Time (hours)")
	ax.grid(linestyle='--', alpha=0.5)



def draw_ratio_vs_time(ax, time, data1, data2, title, alpha):
	print "plotting ratio vs time ("+title+")"
	time = time_to_int(time)
	ratios=[]
	for i,point1 in enumerate(data1):
		point2 = data2[i]
		ratio=point2/point1
		ratios.append(ratio)
	ax.scatter(time, ratios, alpha=alpha, s=10, c='k')
	ax.axhline(y=1, color='b', linestyle='-')
	ax.set_xlabel("Time (hours)")
	samp1=sys.argv[1].split(".")[0].split("/")[-1]
	samp2=sys.argv[2].split(".")[0].split("/")[-1]
	ax.set_ylabel("Ratio ("+samp2+"/"+samp1+")")
	ax.set_title(title)
	ax.grid(linestyle='--', alpha=0.5)
	

def adjacent_values(vals, q1, q3):
	upper_adjacent_value = q3 + (q3 - q1) * 1.5
	upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

	lower_adjacent_value = q1 - (q3 - q1) * 1.5
	lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
	return lower_adjacent_value, upper_adjacent_value


def draw_deliquescence_time(ax, dates):
	set1_percent=[]
	set2_percent=[]
	set1_hours=[]
	set2_hours=[]
	for day in sorted(dates):
		set1=0; set2=0
		for hour in dates[day]:
			if dates[day][hour][0][1]>=75:
				set1+=1
			if dates[day][hour][1][1]>=75:
				set2+=1
		set1_percent.append(set1*100.0/24)
		set2_percent.append(set2*100.0/24)
		set1_hours.append(set1)
		set2_hours.append(set2)

	data = [set1_hours, set2_hours]
	ax.violinplot(data, showmeans=False, showmedians=False, showextrema=False)
	quartile1, medians, quartile3 = np.percentile(data, [25, 50, 75], axis=1)
	whiskers = np.array([ adjacent_values(sorted_array, q1, q3) for sorted_array, q1, q3 in zip(data, quartile1, quartile3)])
	whiskersMin, whiskersMax = whiskers[:, 0], whiskers[:, 1]
	inds = np.arange(1, len(medians) + 1)
	ax.scatter(inds, medians, marker='o', color='white', s=30, zorder=3)
	ax.vlines(inds, quartile1, quartile3, color='k', linestyle='-', lw=5)
	ax.vlines(inds, whiskersMin, whiskersMax, color='k', linestyle='-', lw=1)

	labels = [sys.argv[1].split(".")[0].split("/")[-1], sys.argv[2].split(".")[0].split("/")[-1]]

	ax.get_xaxis().set_tick_params(direction='out')
	ax.xaxis.set_ticks_position('bottom')
	ax.set_xticks(np.arange(1, len(labels) + 1))
	ax.set_xticklabels(labels)
	ax.set_xlim(0.25, len(labels) + 0.75)
	ax.set_ylabel("Time under deliquescence (h)")
	ax.set_ylim(0,24)
	ax.grid(linestyle='--', axis="y", alpha=0.5)


	


# MAIN

# loading data
data1 = load_data(sys.argv[1])
data2 = load_data(sys.argv[2])
time, temp_1, temp_2 = pair_data(data1, data2, 0)
time, humi_1, humi_2 = pair_data(data1, data2, 1)
dates = pair_dates(data1, data2)


# statistics
temp_test = stats.ttest_rel(temp_1, temp_2)
humi_test = stats.ttest_rel(humi_1, humi_2)
print "Temperature comparison test results: "+str(temp_test)
print "Humidity comparison test results: "+str(humi_test)


# plotting set-up
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)
plt.rc('font', family='arial')
fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, figsize=(10,10))
if len(sys.argv)>3:
	alpha=float(sys.argv[3])
else:
	alpha=0.5



# seperate plots
draw_data_vs_data(ax1, temp_1, temp_2, "Temperature", alpha, time)
ax1.set_xticks(np.arange(10, 40, 5))
ax1.set_yticks(np.arange(10, 40, 5))

draw_data_vs_data(ax2, humi_1, humi_2, "Humidity", alpha, time)
ax2.set_xticks(np.arange(20, 100, 10))
ax2.set_yticks(np.arange(20, 100, 10))

draw_data_vs_time(ax3, time, temp_1, temp_2, "Temperature", alpha)
draw_deliquescence_time(ax3, dates)
draw_data_vs_time(ax4, time, humi_1, humi_2, "Humidity", alpha)

draw_ratio_vs_time(ax5, time, temp_1, temp_2, "Temperature", alpha)
draw_ratio_vs_time(ax6, time, humi_1, humi_2, "Humidity", alpha)





#plt.tight_layout()
plt.savefig("figure.png", dpi=300)
plt.show()




