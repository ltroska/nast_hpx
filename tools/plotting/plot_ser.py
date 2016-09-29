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

N = 3

dpi = 100

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


grid_sizes = []
regular = []
mpi = []
hpx = []

with open(infile_paths[0], 'rb') as infile:
	for line in infile:
		splitline = line.split()
		grid_sizes.append(splitline[0])
		reg = float(splitline[1])
		regular.append(1)
		hpx.append(float(splitline[2])/reg)
		mpi.append(float(splitline[3])/reg)

if N == 3:
	fig = plt.figure()
	ax = plt.gca()

	indices = np.arange(len(grid_sizes))
	width = 0.9/3

	reg_bar = ax.bar(indices, regular, width, color='r')
	hpx_bar = ax.bar(indices+width, hpx, width, color='b')
	mpi_bar = ax.bar(indices+2*width, mpi, width, color='g')

	ax.set_ylabel('Runtime (normalized to runtime of regular C++)')
	ax.set_xlabel('Grid Sizes')
	ax.set_xticks(indices+1.5*width)
	ax.set_xticklabels(grid_sizes)
	ax.legend((reg_bar, hpx_bar, mpi_bar), ('regular C++', 'HPX', 'MPI'), title='Code Version', loc='center left', bbox_to_anchor=(1, 0.5))

	plt.savefig(outfile_path + '_bar.png', bbox_inches='tight', format='png', dpi=dpi)
elif N == 2:
	fig = plt.figure()
	ax = plt.gca()

	indices = np.arange(len(grid_sizes))
	width = 0.9/2

	reg_bar = ax.bar(indices, regular, width, color='r')
	hpx_bar = ax.bar(indices+width, hpx, width, color='b')

	ax.set_ylabel('Runtime (normalized to runtime of regular C++)')
	ax.set_xlabel('Grid Sizes')
	ax.set_xticks(indices+width)
	ax.set_xticklabels(grid_sizes)
	ax.legend((reg_bar, hpx_bar), ('regular C++', 'HPX'), title='Code Version', loc='center left', bbox_to_anchor=(1, 0.5))

	plt.savefig(outfile_path + '_bar.png', bbox_inches='tight', format='png', dpi=dpi)

