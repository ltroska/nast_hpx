#!/usr/bin/env python
from __future__ import print_function
import sys
from grid_3d import grid

if len(sys.argv) < 6:
	print("Usage:\t", sys.argv[0] , "<input grid> <output grid> <i_max> <j_max> <k_max> [x_length] [y_length] [z_length] [ignore_before] [ignore_after]\n", file=sys.stderr)
	sys.exit(0)

i_max = int(sys.argv[3])
j_max = int(sys.argv[4])
k_max = int(sys.argv[5])

x_length = sys.argv[6] if len(sys.argv) > 6 else 1
y_length = sys.argv[7] if len(sys.argv) > 7 else 1
z_length = sys.argv[8] if len(sys.argv) > 8 else 1
ignore_before = int(sys.argv[9]) if len(sys.argv) > 9 else 0
ignore_after = int(sys.argv[10]) if len(sys.argv) > 10 else 0

infile_path = sys.argv[1]
outfile_path = sys.argv[2]

b = grid(i_max, j_max, k_max + 2, x_length, y_length, z_length, all_fluid=False)

with open(infile_path) as infile:
	for k in range(k_max + 2):
		for _ in range(ignore_before):
			infile.readline()

		add = [infile.readline().split() for _ in range(j_max + 1, -1, - 1)]
		b.add_xy_plane(k + 1, add, invert=False)
	
		for _ in range(ignore_after):
			infile.readline()

b.refine().write_to(outfile_path)
