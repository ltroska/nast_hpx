#!/usr/bin/env python
from __future__ import print_function
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import sys
import re
import csv
import glob
import numpy as np
import math

dpi = 600

def cmpit(a):
	print(a)
	pattern = '.*_([0-9]*)x[0-9]*(.*)\..*'
	offset = 1 if 'SOR' in re.search(pattern, a).group(2) else 0
	return int(re.search(pattern, a).group(1)) + offset

colors = ['r', 'b', 'g', 'k', 'm']
markers = ['o', 's', '^', 'x', '*']
linestyles = [[1,0], [16,8], [4, 2], [16, 4, 2, 4, 2, 4, 2, 4, 2, 4], [16, 4, 2, 4, 2, 4]]

def inccolor(this_size, last_size, color_num, marker_num, linestyle_num):
	if last_size is not None:	
		if this_size != last_size:
				color_num += 1
				marker_num = 0
		else:
				marker_num += 1

	last_size = this_size
	return this_size,  last_size, color_num, marker_num, linestyle_num

def incmarker(this_size, last_size, color_num, marker_num, linestyle_num):
	if last_size is not None:	
		if last_size is not None:
			if this_size != last_size:
				marker_num += 1
				color_num = 0
			else:
				color_num += 1

	last_size = this_size
	return this_size,  last_size, color_num, marker_num, linestyle_num

def inccolorandmarker(this_size, last_size, color_num, marker_num, linestyle_num):
	if last_size is not None:	
		if this_size != last_size:
				color_num += 1
				marker_num += 1
				linestyle_num = 0
		else:
				linestyle_num += 1

	last_size = this_size
	return this_size,  last_size, color_num, marker_num, linestyle_num

def incmarkerandstyle(this_size, last_size, color_num, marker_num, linestyle_num):
	if last_size is not None:	

		marker_num += 1

		if marker_num >= len(markers):
			marker_num = 0
			linestyle_num += 1
	
		if linestyle_num >= len(linestyles):
			color_num += 1

	last_size = this_size
	return this_size,  last_size, color_num, marker_num, linestyle_num

def inccolorandstyle(this_size, last_size, color_num, marker_num, linestyle_num):
	if last_size is not None:	
		if this_size != last_size:
			color_num += 1
			linestyle_num += 1
			if linestyle_num >= len(linestyles):
				linestyle_num = 0
			marker_num = 0

		else:
			marker_num += 1
	
	last_size = this_size
	return this_size,  last_size, color_num, marker_num, linestyle_num

incfunc = inccolorandmarker

if len(sys.argv) < 4:
	print("Usage:", sys.argv[0] , "<outfile> <title> [input grid]... ", file=sys.stderr)
	sys.exit(0)

adjust = 10*1000

outfile_path = sys.argv[1]
title = sys.argv[2]

infile_paths = []
labels = []
all_labels = []


if len(sys.argv) == 4:
	if sys.argv[3].endswith('.dat'):
		infile_paths.append(sys.argv[3])
	else:
		infile_paths = sorted(glob.glob(sys.argv[3] + ('*.dat' if sys.argv[3].endswith('/') else '/*.dat')), key=cmpit)

else:
	for idx in range(3, len(sys.argv)):
		infile_paths.append(sys.argv[idx])
		

	infile_paths = sorted(infile_paths, key=cmpit)

num_files = len(infile_paths)
num_plots = num_files/2

pattern = '.*_([0-9]*x[0-9]*).*\..*'
idx = 0
while idx < num_files:
	labels.append(re.search(pattern, infile_paths[idx]).group(1))
	idx += 2

pattern = '.*_([0-9]*x[0-9]*.*)\..*'

for infile in infile_paths:
	all_labels.append(re.search(pattern, infile).group(1))

print(infile_paths)
print(labels)
print(all_labels)
processors = []
avg_runtimes = []
all_jac_avg_runtimes = []
all_sor_avg_runtimes = []
all_jac_avg_idle_rates = []
all_sor_avg_idle_rates = []

num = 0
idx = 0
while idx < num_files:
	processors.append([])
	jac_avg_runtimes = []
	sor_avg_runtimes = []
	jac_avg_idle_rate = []
	sor_avg_idle_rate = []

	with open(infile_paths[idx]) as jacobi_file:
		print('processing', infile_paths[idx])
		for line in jacobi_file:
			splitline = line.split()
			
			processors[num].append(int(splitline[0]))				
			jac_avg_runtimes.append(float(splitline[1])/adjust)
			jac_avg_idle_rate.append(float(splitline[4])/100)
	
	with open(infile_paths[idx+1]) as sor_file:
		print('processing', infile_paths[idx+1])
		for line in sor_file:
			splitline = line.split()
			
			sor_avg_runtimes.append(float(splitline[1])/adjust)
			sor_avg_idle_rate.append(float(splitline[4])/100)

	avg_runtimes.append([sor/jac for jac, sor in zip(jac_avg_runtimes, sor_avg_runtimes)])
	all_jac_avg_runtimes.append(jac_avg_runtimes)
	all_sor_avg_runtimes.append(sor_avg_runtimes)
	all_jac_avg_idle_rates.append(jac_avg_idle_rate)
	all_sor_avg_idle_rates.append(sor_avg_idle_rate)
	idx += 2
	num += 1
	print('next')

f = plt.gcf()
ax = plt.gca()

incfunc = inccolorandmarker
color_num = 0
marker_num = 0
linestyle_num = 0
last_size = None
for idx, times in enumerate(avg_runtimes):
	this_size = labels[idx].split('x')[0]
	this_size, last_size, color_num, marker_num, linestyle_num = incfunc(this_size, last_size, color_num, marker_num, linestyle_num)
	
	plt.plot(processors[idx], times, dashes=linestyles[linestyle_num], marker=markers[marker_num], color=colors[color_num], label=labels[idx])

val_max = max([max(time) for time in avg_runtimes])
val_min = min([min(time) for time in avg_runtimes])

val_max = math.ceil(val_max * 10.0) / 10
val_min = math.floor(val_min * 10.0) / 10

plt.xlabel("Number of Processors")
plt.ylabel(r"$\frac{t_{SOR}}{t_{Jacobi}}$", rotation=0, fontsize=26)
plt.ylim([val_min, val_max])
plt.yticks(np.arange(val_min, val_max + 0.1, 0.1))
ax = plt.gca()
#ax.tick_params(length=0)
ax.hlines(1, 0, 600)
minor_locator = AutoMinorLocator(4)
ax.xaxis.set_minor_locator(minor_locator)

ax.grid(True, which='both')
handles, labels = ax.get_legend_handles_labels()
lgd = plt.legend(handles[::-1], labels[::-1], loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=len(labels), title='Grid Sizes', handlelength=3)
plt.savefig(outfile_path + '_ratio.png', bbox_inches='tight', format='png', dpi=dpi)

fig = plt.figure()
ax = plt.gca()
fig.set_figwidth(15)

incfunc = inccolorandmarker
color_num = 0
marker_num = 0
linestyle_num = 0
last_size = None
for idx in range(len(all_jac_avg_idle_rates)):
	label = all_labels[2*idx]
	this_size = label.split('x')[0]
	this_size, last_size, color_num, marker_num, linestyle_num = incfunc(this_size, last_size, color_num, marker_num, linestyle_num)
	
	plt.plot(processors[idx], all_jac_avg_idle_rates[idx], dashes=linestyles[linestyle_num], marker=markers[marker_num], color=colors[color_num], label=label)

	label = all_labels[1 + 2*idx]
	this_size = label.split('x')[0]
	this_size, last_size, color_num, marker_num, linestyle_num = incfunc(this_size, last_size, color_num, marker_num, linestyle_num)
	
	plt.plot(processors[idx], all_sor_avg_idle_rates[idx], dashes=linestyles[linestyle_num], marker=markers[marker_num], color=colors[color_num], label=label)


plt.xlabel("Number of Processors")
plt.ylabel("Idle Rate (%)")
plt.yticks(np.arange(0, 110, 10))
ax = plt.gca()
#ax.tick_params(length=0)
minor_locator = AutoMinorLocator(4)
ax.xaxis.set_minor_locator(minor_locator)

ax.grid(True, which='both')
handles, labelsd = ax.get_legend_handles_labels()
lgd = plt.legend(handles[::-1], labelsd[::-1], loc='upper center', bbox_to_anchor=(0.5, -0.1), ncol=len(labels), title='Grid Sizes', handlelength=3)
plt.savefig(outfile_path + '_idle_rate.png', bbox_inches='tight', format='png', dpi=dpi)

fig = plt.figure()
ax = plt.gca()

indices = np.arange(len(processors[-1]))
width = 0.9/2

jac_bar = ax.bar(indices, all_jac_avg_runtimes[-1], width, color='r')
sor_bar = ax.bar(indices+width, all_sor_avg_runtimes[-1], width, color='b')

ax.set_ylabel('Average Runtime Per Iteration (s)')
ax.set_xlabel('Number Of Processors')
ax.set_xticks(indices+width)
ax.set_xticklabels(processors[-1])
ax.legend((jac_bar, sor_bar), ('Jacobi', 'SOR'), title='Grid: ' + labels[-1])

plt.savefig(outfile_path + '_bar.png', bbox_inches='tight', format='png', dpi=dpi)

