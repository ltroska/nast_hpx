#!/usr/bin/env python
from __future__ import print_function
from grid import *
import sys
import random

if len(sys.argv) < 9:
  print("Usage:", sys.argv[0] , "<output grid> <#cols> <#rows> <width> <lower bound for length of rect> <upper bound for length of rect> <lower bound for distance> <upper bound for distance> [deadend = 0.1] [max recursion=0('infinite')]", file=sys.stderr)
  sys.exit(0)

outfile_path = sys.argv[1]
num_cols = int(sys.argv[2])
num_rows = int(sys.argv[3])
width = int(sys.argv[4])
lengthlow = int(sys.argv[5])
lengthhigh = int(sys.argv[6])
distancelow = int(sys.argv[7])
distancehigh = int(sys.argv[8])

dead_perc = 10
if len(sys.argv) == 10:
	dead_perc = int(sys.argv[9])

max_n = sys.maxint
if len(sys.argv) == 11:
	max_n = int(sys.argv[10])

def make_branch(grid, pos, width, distance, n):
	if n <= max_n and not (pos[0] > num_cols - 2 or pos[1] > num_rows - 2 or pos[0] < 1 or pos[1] < 1):
		lengthrand = random.randint(lengthlow, lengthhigh)
		distrand = random.randint(distancelow, distancehigh)
		grid.add_shape(rectangle(pos, [pos[0] + lengthrand, pos[1] + width - 1]))
		if random.randint(0, 100) >= dead_perc:
			grid.add_shape(zigzag([pos[0] + 1 + lengthrand - width, pos[1] + width], width, distrand))
			make_branch(grid, [pos[0] + 1 + lengthrand + distrand, pos[1] + distrand + width], width, distrand, n + 1)
		if random.randint(0, 100) >= dead_perc:
			grid.add_shape(zigzag([pos[0] + 1 + lengthrand - width, pos[1] - 1], width, -distrand))
			make_branch(grid, [pos[0] + 1 + lengthrand + distrand, pos[1] - 1 - distrand - width], width, distrand, n + 1)

g = grid(num_cols, num_rows, inverted=True)

start_y = 1 + (num_rows - 2) / 2 - width/2

make_branch(g, [1, start_y], width, random.randint(distancelow, distancehigh), 1)

g.apply_shapes()
g.write_to(outfile_path)

#for split in range(num_splits + 1):
	
