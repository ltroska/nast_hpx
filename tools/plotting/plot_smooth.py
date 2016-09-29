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

reverse = False

colors = ['r', 'b', 'g', 'k', 'm', 'c']
markers = ['o', 's', '^', 'x', '*', 'p']
#linestyles = ['-', '--', '-.', ':']
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
	global reverse
	if last_size is not None:	
		if this_size != last_size:
				color_num += 1
				if color_num >= len(colors):
					color_num = 0
					reverse = not reverse
				
				marker_num += 1 if not reverse else -1
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
	
		if linestyle_num >= len(linestyles) or this_size != last_size:
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

if len(sys.argv) < 5:
	print("Usage:", sys.argv[0] , "<outfile> <title> <weak/strong> [input grid]... ", file=sys.stderr)
	sys.exit(0)

outfile_path = sys.argv[1]
title = sys.argv[2]
strong = (sys.argv[3] == "strong" or sys.argv[3] == "1")

infile_paths = []
labels = []
pattern = '.*/(.*)\..*'

def cmpit(a):
	return int(a.split()[0].split('/')[-1])
	

if len(sys.argv) == 5:
	if sys.argv[4].endswith('.dat'):
		infile_paths.append(sys.argv[4])
		labels.append(re.search(pattern, sys.argv[4]).group(1))
	else:
		infile_paths = sorted(glob.glob(sys.argv[4] + ('*.dat' if sys.argv[4].endswith('/') else '/*.dat')), key=cmpit)
		
else:
	for idx in range(4, len(sys.argv)):
		infile_paths.append(sys.argv[idx])
		

if any('step' in infile_path for infile_path in infile_paths):
	colors.insert(0, 'm')
	markers.insert(0, '*')

for path in infile_paths:
	labels.append(re.search(pattern, path).group(1))

print(labels)

num_sizes = len(set([label.split('x')[0] for label in labels]))

processors = []
avg_runtimes = []
min_runtimes = []
max_runtimes = []

avg_idle = []
min_idle = []
max_idle = []

num_files = len(labels)

has_serial = [True] * len(labels)

for idx, file_path in enumerate(infile_paths):
	with open(file_path, 'rb') as infile:
		lines = infile.readlines()
		processors.append([])
		avg_runtimes.append([])
		min_runtimes.append([])
		max_runtimes.append([])

		avg_idle.append([])
		min_idle.append([])
		max_idle.append([])

		for line in lines:
			splitline = line.split()
			
			processors[idx].append(int(splitline[0]))				
			avg_runtimes[idx].append(float(splitline[1]))
			min_runtimes[idx].append(float(splitline[2]))
			max_runtimes[idx].append(float(splitline[3]))

			avg_idle[idx].append(float(splitline[4]) / 100)
			min_idle[idx].append(float(splitline[5]) / 100)
			max_idle[idx].append(float(splitline[6]) / 100)

	print(processors[idx])
	print(avg_runtimes[idx])

	if processors[idx][0] != 1:
		#avg_runtimes[idx].insert(0, avg_runtimes[idx][0] * processors[idx][0])
		#min_runtimes[idx].insert(0, min_runtimes[idx][0] * processors[idx][0])
		#max_runtimes[idx].insert(0, max_runtimes[idx][0] * processors[idx][0])
		#avg_idle[idx].insert(0, 0)
		#min_idle[idx].insert(0, 0)
		#max_idle[idx].insert(0, 0)
		#processors[idx].insert(0, 1)
		has_serial[idx] = False

if strong:
	runtime_errors = []
	speedups = []
	max_speedups = []
	min_speedups = []
	speedup_errors = []
	efficiencies = []
	max_efficiencies = []
	min_efficiencies = []
	efficiency_errors = []

	for i, times in enumerate(max_runtimes):
		t0 = times[0]
		max_speedups.append([t0/time for time in times])
		max_efficiencies.append([speed/processors[i][idx] for idx, speed in enumerate(max_speedups[i])])

	for i, times in enumerate(min_runtimes):
		t0 = times[0]
		min_speedups.append([t0/time for time in times])
		min_efficiencies.append([speed/processors[i][idx] for idx, speed in enumerate(min_speedups[i])])

	for i, times in enumerate(avg_runtimes):
		t0 = times[0]
		speedups.append([t0/time for time in times])
		speedup_errors.append([[speedups[i][j] - min_speedups[i][j] for j in range(len(speedups[i]))], [max_speedups[i][j] - speedups[i][j] for j in range(len(speedups[i]))]])

		efficiencies.append([speed/processors[i][idx] for idx, speed in enumerate(speedups[i])])
		efficiency_errors.append([[efficiencies[i][j] - min_efficiencies[i][j] for j in range(len(efficiencies[i]))], [max_efficiencies[i][j] - efficiencies[i][j] for j in range(len(efficiencies[i]))]])

	print(processors[0])
	print(avg_runtimes[0])

	#ideal = [t[0]/p for p in processors] 

	plt.suptitle(title)

	f = plt.gcf()
	f.set_figheight(15)
	f.set_figwidth(15)
	
	if False:
		ax = plt.subplot(4, 1, 1)
		ax.set_xscale("log", nonposx='clip')
		ax.set_yscale("log", nonposy='clip')

		color_num = 0
		marker_num = 0
		linestyle_num = 0
		last_size = None
		for idx, times in enumerate(avg_runtimes):
			this_size = labels[idx].split()[2][:-1]
			print(this_size)
			this_size, last_size, color_num, marker_num, linestyle_num = incfunc(this_size, last_size, color_num, marker_num, linestyle_num)
			
			plt.loglog(processors[idx], times, dashes=linestyles[linestyle_num], marker=markers[marker_num], color=colors[color_num], label=labels[idx])
			#print(processors[idx], times)
		#plt.loglog(p, ideal, 'b-o', label="ideal scaling")
		plt.title("Runtime")
		plt.xlabel("Number of Processors")
		plt.ylabel("Runtime (s)")
		ax = plt.gca()
		ax.set_aspect('auto')
		plt.grid(True, which='major')
		minor_locator = AutoMinorLocator(2)
		ax.xaxis.set_minor_locator(minor_locator)
		#ax.set_xscale("log", nonposx='clip')
		plt.subplot(4, 1, 2)

		color_num = 0
		marker_num = 0
		linestyle_num = 0
		last_size = None

		for idx, speedup in enumerate(speedups):
			this_size = labels[idx].split()[2][:-1]
			print(this_size)
			this_size, last_size, color_num, marker_num, linestyle_num = incfunc(this_size, last_size, color_num, marker_num, linestyle_num)

			if has_serial[idx]:
				plt.plot(processors[idx], speedup, dashes=linestyles[linestyle_num], marker=markers[marker_num], color=colors[color_num], label=labels[idx])

		plt.title("Speedup")
		plt.xlabel("Number of Processors")
		plt.ylabel("Speedup")
		ax = plt.gca()
		ax.set_aspect('auto')
		plt.grid(True, which='both')
		minor_locator = AutoMinorLocator(2)
		ax.xaxis.set_minor_locator(minor_locator)
		#ax.set_xscale("log", nonposx='clip')

	plt.subplot(4, 1, 3)

	color_num = 0
	marker_num = 0
	linestyle_num = 0
	last_size = None

	for idx, efficiency in enumerate(efficiencies):
		this_size = labels[idx].split()[0]
		print(this_size)
		this_size, last_size, color_num, marker_num, linestyle_num = incfunc(this_size, last_size, color_num, marker_num, linestyle_num)

		if has_serial[idx]:
			plt.plot(processors[idx], efficiency, dashes=linestyles[linestyle_num], marker=markers[marker_num], color=colors[color_num], label=labels[idx])

	plt.title("Efficiency")
	plt.xlabel("Number of Processors")
	plt.ylabel("Efficiency")
	plt.yticks(np.arange(0, 1.1, 0.1))
	ax = plt.gca()
	ax.set_aspect('auto')
	plt.grid(True, which='both')
	minor_locator = AutoMinorLocator(2)
	ax.xaxis.set_minor_locator(minor_locator)
	#ax.set_xscale("log", nonposx='clip')
	plt.subplot(4, 1, 4)
	color_num = 0
	marker_num = 0
	linestyle_num = 0
	last_size = None

	plot_data = []
	for idx, idle_rate in enumerate(avg_idle):
		this_size = labels[idx].split()[0]
		this_size, last_size, color_num, marker_num, linestyle_num = incfunc(this_size, last_size, color_num, marker_num, linestyle_num)

		plt.plot(processors[idx], idle_rate, marker=markers[marker_num], dashes=linestyles[linestyle_num], color=colors[color_num], label=labels[idx])

	plt.title("Average Idle Rate")
	plt.xlabel("Number of Processors")
	plt.ylabel("Idle Rate (%)")
	plt.yticks(np.arange(0, 110, 10))
	ax = plt.gca()
	ax.set_aspect('auto')
	plt.grid(True, which='both')
	minor_locator = AutoMinorLocator(2)
	ax.xaxis.set_minor_locator(minor_locator)
	#ax.set_xscale("log", nonposx='clip')


	lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=5, handlelength=3)
	#lgd = plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

	plt.tight_layout()
	plt.subplots_adjust(top=0.9)
	plt.subplots_adjust(hspace=.35)

	plt.savefig(outfile_path + ".png", bbox_extra_artists=(lgd,), bbox_inches='tight', format='png', dpi=dpi)
else:
	ideal = [t[0] for pv in p] 
	plt.suptitle(title)

	f = plt.gcf()
	f.set_figheight(15)
	f.set_figwidth(15)

	plt.loglog(p, t, 'r-o', label="achieved runtimes")
	plt.loglog(p, ideal, 'b-o', label="ideal scaling")
	plt.title("Runtime")
	plt.xlabel("Number of Processors")
	plt.ylabel("Runtime (s)")
	ax = plt.gca()
	ax.set_ylim([10,100])
	ax.grid(True, which='both')
	plt.legend()
	plt.savefig(outfile_path + '.svg', bbox_inches='tight', format='svg', dpi=dpi)


