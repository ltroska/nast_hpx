#!/usr/bin/env python
from __future__ import print_function
import matplotlib.pyplot as plt
import sys
import re

import csv

if len(sys.argv) < 5:
	print("Usage:", sys.argv[0] , "<outfile> <title> <weak/strong> [input grid]... ", file=sys.stderr)
	sys.exit(0)

outfile_path = sys.argv[1]
title = sys.argv[2]
strong = (sys.argv[3] == "strong" or sys.argv[3] == "1")

infile_paths = []
labels = []

pattern = '.*_(.*)\..*'
for idx in range(4, len(sys.argv)):
	infile_paths.append(sys.argv[idx])
	labels.append(re.search(pattern, sys.argv[idx]).group(1))

processors = []
avg_runtimes = []
min_runtimes = []
max_runtimes = []

avg_idle = []
min_idle = []
max_idle = []

num_files = len(labels)

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

	if processors[idx][0] != 1:
		avg_runtimes[idx].insert(0, avg_runtimes[idx][0] * processors[idx][0])
		min_runtimes[idx].insert(0, min_runtimes[idx][0] * processors[idx][0])
		max_runtimes[idx].insert(0, max_runtimes[idx][0] * processors[idx][0])
		processors[idx].insert(0, 1)



styles = ['r-o', 'b-x', 'g-*', 'k-+', 'm-h', 'c-v']

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

	ax = plt.subplot(4, 1, 1)
	ax.set_xscale("log", nonposx='clip')
	ax.set_yscale("log", nonposy='clip')

	for idx, times in enumerate(avg_runtimes):
		print(idx)
		plt.loglog(processors[idx], times, styles[idx], label=labels[idx])
		#print(processors[idx], times)
	#plt.loglog(p, ideal, 'b-o', label="ideal scaling")
	plt.title("Runtime")
	plt.xlabel("Number of Processors")
	plt.ylabel("Runtime (s)")
	ax = plt.gca()
	#ax.set_ylim([10,100])
	ax.grid(True, which='both')

	plt.subplot(4, 1, 2)

	for idx, speedup in enumerate(speedups):
		plt.errorbar(processors[idx], speedup, yerr=speedup_errors[idx], fmt=styles[idx], label=labels[idx])

	plt.title("Speedup")
	plt.xlabel("Number of Processors")
	plt.ylabel("Speedup")
	ax.set_aspect('auto')
	plt.grid()

	plt.subplot(4, 1, 3)

	for idx, efficiency in enumerate(efficiencies):
		plt.errorbar(processors[idx], efficiency, yerr=efficiency_errors[idx], fmt=styles[idx], label=labels[idx])

	plt.title("Efficiency")
	plt.xlabel("Number of Processors")
	plt.ylabel("Efficiency")
	ax.set_aspect('auto')
	plt.grid()

	plt.subplot(4, 1, 4)

	for idx, idle_rate in enumerate(avg_idle):
		plt.plot(processors[idx], idle_rate, styles[idx], label=labels[idx])

	plt.title("Average Idle Rate")
	plt.xlabel("Number of Processors")
	plt.ylabel("Idle Rate (%)")
	ax.set_aspect('auto')
	plt.grid()
	lgd = plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=len(avg_runtimes))

	plt.tight_layout()
	plt.subplots_adjust(top=0.9)

	plt.savefig(outfile_path + ".png", bbox_extra_artists=(lgd,), bbox_inches='tight')
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
	plt.savefig(outfile_path + ".png")


