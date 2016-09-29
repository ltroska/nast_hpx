#!/usr/bin/env python
from __future__ import print_function
import matplotlib.pyplot as plt
import sys

import csv

if len(sys.argv) != 4:
	print("Usage:", sys.argv[0] , "<input grid> <title> <weak/strong>", file=sys.stderr)
	sys.exit(0)

infile_path = sys.argv[1]
title = sys.argv[2]
strong = (sys.argv[3] == "strong" or sys.argv[3] == "1")
outfile_path = infile_path[:-4]

p = []
t = []


with open(infile_path, 'rb') as csvfile:
	gridreader = csv.reader(csvfile, delimiter=' ', quotechar='"')
	
	for row in gridreader:
		if len(row) > 1:
			p.append(float(row[0]))
			t.append(float(row[1]))

print(p)
print(t)

if strong:
	ideal = [t[0]/pv for pv in p] 
	print(ideal)

	plt.suptitle(title)

	f = plt.gcf()
	f.set_figheight(15)
	f.set_figwidth(15)

	plt.subplot(3, 1, 1)
	plt.loglog(p, t, 'r-o', label="achieved runtimes")
	#plt.loglog(p, ideal, 'b-o', label="ideal scaling")
	plt.title("Runtime")
	plt.xlabel("#processors")
	plt.ylabel("runtime (s)")
	ax = plt.gca()
	#ax.set_ylim([10,100])
	ax.grid(True, which='both')
	plt.legend()

	plt.subplot(3, 1, 2)

	t0 = t[0]

	speedup = [t0/time for time in t]

	plt.plot(p, speedup, 'r-o')
	plt.title("Speedup")
	plt.xlabel("#processors")
	plt.ylabel("speedup")
	ax.set_aspect('auto')
	plt.grid()


	plt.subplot(3, 1, 3)

	efficiency = [speed/p[pv] for pv, speed in enumerate(speedup)]

	plt.plot(p, efficiency, 'r-o')
	plt.title("Efficiency")
	plt.xlabel("#processors")
	plt.ylabel("efficiency")
	ax.set_aspect('auto')
	plt.grid()

	plt.tight_layout()
	plt.savefig(outfile_path + ".png")
else:
	ideal = [t[0] for pv in p] 
	plt.suptitle(title)

	f = plt.gcf()
	f.set_figheight(15)
	f.set_figwidth(15)

	plt.loglog(p, t, 'r-o', label="achieved runtimes")
	plt.loglog(p, ideal, 'b-o', label="ideal scaling")
	plt.title("Runtime")
	plt.xlabel("#processors")
	plt.ylabel("runtime (ms)")
	ax = plt.gca()
	ax.set_ylim([10,100])
	ax.grid(True, which='both')
	plt.legend()
	plt.savefig(outfile_path + ".png")


