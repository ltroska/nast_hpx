#!/usr/bin/env python
from __future__ import print_function
import sys

import csv

if len(sys.argv) != 3:
	print("Usage:", sys.argv[0] , "<input grid> <output grid>", file=sys.stderr)
	sys.exit(0)

infile_path = sys.argv[1]

in_grid = []

with open(infile_path, 'rb') as csvfile:
	gridreader = csv.reader(csvfile, delimiter=',', quotechar='"')
	
	for row in gridreader:
		in_grid.append(row)

outfile_path = sys.argv[2]
with open(outfile_path, 'wb') as outfile:
	outfile.write( str(len(in_grid)) + " " + str(len(in_grid[0])) + " 1\n")

	for i, row in enumerate(in_grid):
		outfile.write(" ".join(row))

		if i != len(in_grid) - 1:
			outfile.write(" ")
	
