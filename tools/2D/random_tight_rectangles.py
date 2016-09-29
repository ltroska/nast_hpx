#!/usr/bin/env python
from __future__ import print_function
from grid import *
import sys
import random
import itertools

if len(sys.argv) < 8:
  print("Usage:", sys.argv[0] , "<output grid> <#cols> <#rows> <#iterations> <minimal side length> <maximal side length> <distance> [include_boundary=0]", file=sys.stderr)
  sys.exit(0)

outfile_path = sys.argv[1]
num_cols = int(sys.argv[2])
num_rows = int(sys.argv[3])
max_iter = int(sys.argv[4])
min_side_length = int(sys.argv[5])
max_side_length = int(sys.argv[6])
distance = int(sys.argv[7])
  
include_boundary = False
if len(sys.argv) >= 9:
  include_boundary = (int(sys.argv[8]) == 1)
 
g = grid(num_cols, num_rows)

x_range = range(1, num_cols - 2) if include_boundary else range(0, num_cols)
y_range = range(1, num_rows - 2) if include_boundary else range(0, num_rows)

indices = itertools.product(x_range, y_range)

iters = 0
while iters < max_iter:  
  for x, y in indices:
    length1 = random.randint(min_side_length, max_side_length)
    length2 = random.randint(min_side_length, max_side_length)

    if g.add_shape(rectangle([x, y], [x + length1, y + length2]), distance, include_boundary):
      iters += 1
      print(iters)
      
      if iters >= max_iter:
	break
  
  iters += 1
  
g.apply_shapes()
g.write_to(outfile_path)