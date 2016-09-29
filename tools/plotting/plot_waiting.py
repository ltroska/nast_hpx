#!/usr/bin/env python
from __future__ import print_function
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import sys
import re
import csv
import glob
import numpy as np

dpi = 600

wait_factor = 1

def cmpit(a):
	pattern = '.*_([0-9]*)x[0-9]* (.*) (.*) waiting\..*'
	m = re.search(pattern, a)
	return int(m.group(1)) + (10 if m.group(2).startswith('S') else 0) + (1 if 'without' not in m.group(3) else 0)

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

outfile_path = sys.argv[1]
title = sys.argv[2]

infile_paths = []
labels = []
pattern = '.*_(.*) .* .*\..*'

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

idx = 0
while idx < num_files:
	labels.append(re.search(pattern, infile_paths[idx]).group(1))
	idx += 2

print(infile_paths)
print(labels)

processors = []
avg_runtimes = []

num = 0
idx = 0
while idx < num_files:
	processors.append([])
	nowait_avg_runtimes = []
	wait_avg_runtimes = []

	with open(infile_paths[idx]) as nowait_file:
		print('processing', infile_paths[idx])
		for line in nowait_file:
			splitline = line.split()
			
			processors[num].append(int(splitline[0]))				
			nowait_avg_runtimes.append(float(splitline[1]))

	with open(infile_paths[idx+1]) as wait_file:
		print('processing', infile_paths[idx+1])
		for line in wait_file:
			splitline = line.split()
			
			wait_avg_runtimes.append(wait_factor*float(splitline[1]))

	avg_runtimes.append([wait/nowait for wait, nowait in zip(wait_avg_runtimes, nowait_avg_runtimes)])
	idx += 2
	num += 1

f = plt.gcf()
ax = plt.gca()
f.set_figheight(15)
f.set_figwidth(15)

color_num = 0
marker_num = 0
linestyle_num = 0
last_size = None
for idx, times in enumerate(avg_runtimes):
	this_size = labels[idx].split('x')[0]
	this_size, last_size, color_num, marker_num, linestyle_num = incfunc(this_size, last_size, color_num, marker_num, linestyle_num)
	
	plt.plot(processors[idx], times, dashes=linestyles[linestyle_num], marker=markers[marker_num], color=colors[color_num], label=labels[idx])

plt.xlabel("Number of Processors")
plt.ylabel(r"$\frac{t_w}{t_n}$", rotation=0, fontsize=26)
plt.yticks(np.arange(0.4, 1.6, 0.1))
ax = plt.gca()
ax.hlines(1, 0, 600)
minor_locator = AutoMinorLocator(3)
ax.xaxis.set_minor_locator(minor_locator)
ax.yaxis.set_minor_locator(minor_locator)
ax.tick_params(axis='x',which='minor',bottom='off')
ax.grid(True, which='both')
plt.legend(loc='best', title='Grid Sizes')
plt.savefig(outfile_path + '.png', bbox_inches='tight', format='png', dpi=dpi)

