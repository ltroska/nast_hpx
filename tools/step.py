#!/usr/bin/env python
from __future__ import print_function
from grid import *
import sys
import random
import itertools

if len(sys.argv) != 9:
  print("Usage:", sys.argv[0] , "<output grid> <#cols> <#rows> <bottom step width> <bottom step height> <distance to obstacle> <obstacle width> <obstacle height>", file=sys.stderr)
  sys.exit(0)

outfile_path = sys.argv[1]
num_cols = int(sys.argv[2])
num_rows = int(sys.argv[3])
bs_width = int(sys.argv[4])
bs_height = int(sys.argv[5])
distance = int(sys.argv[6])
obs_width = int(sys.argv[7])
obs_height = int(sys.argv[8])
  
g = grid(num_cols, num_rows)

g.add_shape(rectangle([1, 1], [1 + bs_width, 1 + bs_height]))
g.add_shape(rectangle([1 + bs_width + distance, num_rows - 1 - obs_height], [1 + bs_width + distance + obs_width, num_rows - 2]))
  
g.apply_shapes()
g.write_to(outfile_path)
