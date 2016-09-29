#!/usr/bin/env python
from __future__ import print_function
import sys
from grid_3d import grid

if len(sys.argv) != 6:
	print("Usage:\t", sys.argv[0] , "<i_max> <j_max> <k_max> <type> <output grid>\n\t", "where <type> is one of: driven, step", file=sys.stderr)
	sys.exit(0)

i_max = int(sys.argv[1])
j_max = int(sys.argv[2])
k_max = int(sys.argv[3])

type_string = sys.argv[4]
outfile_path = sys.argv[5]

b = grid(i_max, j_max, k_max)

if type_string.lower() == "driven":
	print("Creating 'Driven Cavity' grid at", outfile_path, ".")
elif type_string.lower() == "step":
	print("Creating 'Flow Over A Step' grid at", outfile_path, ".")
	b.insert_rectangle(0, 1 + i_max/4, 0, j_max + 1, 0, 1 + k_max/2)
else:
	print("Type not supported!", file=sys.stderr)
	sys.exit(0)

								
b.write_to(outfile_path)
