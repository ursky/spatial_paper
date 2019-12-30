#!/usr/bin/env python 
import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats
from matplotlib.patches import Patch
import pyqt_fit.nonparam_regression as smooth
from pyqt_fit import npr_methods


def convert_time(time):
	cut = time.split()
	day_cut=cut[0].split("/")
	time_cut = cut[1].split(":")
	hour = float(time_cut[0])
	mins = float(time_cut[1])
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
	return day, float(hour)+mins/60


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

		if day.split('-')[1]!="04" and day.split('-')[1]!="03":
			continue
		if day not in data:
			data[day]={}
		if hour in data[day]:
			temp = (temp+data[day][hour][0])/2
			humi = (humi+data[day][hour][1])/2
		data[day][hour]=(temp, humi)	
	return data
	


def get_average_day(data):
	hours=[]
	temps=[]
	humis=[]
	for day in data:
		for hour in sorted(data[day]):
			hours.append(float(hour))
			temps.append(data[day][hour][0])
			humis.append(data[day][hour][1])
	return sorted(hours), [x for _,x in sorted(zip(hours,temps))], [x for _,x in sorted(zip(hours,humis))]


def plot_data (xs, ys, color, ax):
	ax.scatter(xs, ys, alpha=0.01, c=color)
	if color=="gold": color="orange"
	if color=="viotet": color="magenta"
	
	# line of best fit
	grid = np.r_[0:24:512j]
	k0 = smooth.NonParamRegression(xs, ys, method=npr_methods.SpatialAverage())
	k0.fit()
	ax.plot(grid, k0(grid), color, linewidth=2)
	ax.set_xticks([0,4,8,12,16,20,24])
	ax.set_xlim(0,24)

def fit(xs, ys):
	est = smooth.NonParamRegression(xs, ys, method=npr_methods.LocalPolynomialKernel(q=2))
	est.fit()
	return est

def plot_temp(data, color, ax):
	hours, temps, humis = get_average_day(data)
	plot_data(hours, temps, color, ax)
	

def plot_humi(data, color, ax):
	hours, temps, humis = get_average_day(data)
	plot_data(hours, humis, color, ax)
	ax.set_ylim(19,100)



data = {}
data["top"] = load_data("SG1_In-Top.csv")
data["mid"] = load_data("SG1_In-Mid.csv")
data["bot"] = load_data("SG1_In-Bot.csv")
data["out"] = load_data("SG1_Top.csv")

# plotting set-up
font = {'family': 'arial', 'weight': 'normal', 'size': 12}
plt.rc('font', **font)
plt.rc('font', family='arial')
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6,6))
colors=["gold", "red", "violet", "cyan"]

for i,position in enumerate(["out", "top", "mid", "bot"]):
	print "plotting", position, "temperature..."
	plot_temp(data[position], colors[i], ax1)
	print "plotting", position, "humidity..."
	plot_humi(data[position], colors[i], ax2)


ax1.set_xlabel("Time after midnight (h)")
ax2.set_xlabel("Time after midnight (h)")
ax1.set_ylabel("Temperature ($^\circ$C)")
ax2.set_ylabel("Relative humidity (%)")
ax1.grid(ls="--", c="k", alpha=0.2)
ax2.grid(ls="--", c="k", alpha=0.2)


#making legend
ax = fig.add_axes([0.0, 0.0, 1, 0.05])
ax.axis("off")
legend_elements = [Patch(facecolor=colors[0], edgecolor='k', label='Air', linewidth=1),
        Patch(facecolor=colors[1], edgecolor='k', label='Top', linewidth=1),
        Patch(facecolor=colors[2], edgecolor='k', label='Middle', linewidth=1),
        Patch(facecolor=colors[3], edgecolor='k', label='Bottom', linewidth=1)]
ax.legend(handles=legend_elements, loc="lower center", framealpha=1, frameon=True, facecolor='w', ncol=4, columnspacing=1, handlelength=1, prop={'size': 16})

ax1.annotate("A.", xy=(-0.13, 1.02), xycoords="axes fraction", fontsize=20)
ax2.annotate("B.", xy=(-0.13, 1.02), xycoords="axes fraction", fontsize=20)

plt.tight_layout(rect=[0, 0.1, 1, 0.97])
plt.savefig("average_day.png", dpi=300)




