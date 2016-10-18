#!/usr/bin/env python
from __future__ import print_function
import sys
import csv
from grid_3d import grid

if len(sys.argv) != 6:
	print("Usage:\t", sys.argv[0] , "<#repetitions> <height> <distance> <input grid> <output grid>", file=sys.stderr)
	sys.exit(0)


num_repetitions = int(sys.argv[1])
height = int(sys.argv[2])
distance = int(sys.argv[3])

infile_path = sys.argv[4]
outfile_path = sys.argv[5]

in_grid = []

with open(infile_path, 'rb') as csvfile:
	gridreader = csv.reader(csvfile, delimiter=',', quotechar='"')
	
	for row in gridreader:
		in_grid.append(row)


i_max = len(in_grid[0]) - 2
j_max = len(in_grid) - 2
k_max = 2 * distance + (num_repetitions - 1) * distance + height * num_repetitions

b = grid(i_max, j_max, k_max, all_fluid=False)

k = 1 + distance

for _ in range(num_repetitions):
	for _ in range(height):
		b.add_xy_plane(k, in_grid)
		k += 1

	k += distance

print("Writing grid of size", i_max, j_max, k_max, "to", outfile_path, ".")
b.write_to(outfile_path)

