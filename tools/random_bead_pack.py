#!/usr/bin/env python
from __future__ import print_function
from grid import *
import sys
import random

if len(sys.argv) < 7:
  print("Usage:", sys.argv[0] , "<output grid> <#cols> <#rows> <#beads> <minimal radius> <maximal radius> [distance=0] [include_boundary=0] [maximal number of failures=2*num_beads]", file=sys.stderr)
  sys.exit(0)

outfile_path = sys.argv[1]
num_cols = int(sys.argv[2])
num_rows = int(sys.argv[3])
num_beads = int(sys.argv[4])
min_radius = int(sys.argv[5])
max_radius = int(sys.argv[6])

distance = 0
if len(sys.argv) >= 8:
  distance = int(sys.argv[7])
  
include_boundary = False
if len(sys.argv) >= 9:
  include_boundary = (int(sys.argv[8]) == 1)
  
max_fails = 2 * num_beads
if len(sys.argv) >= 10:
  max_fails = int(sys.argv[9])
  
g = grid(num_cols, num_rows)

count = 0
failed = 0
resets = 0
while count < num_beads:
  center_x = random.randint(1, num_cols - 2)
  center_y = random.randint(1, num_rows - 2)
  
  radius = random.randint(min_radius, max_radius)
  
  success = g.add_shape(circle([center_x, center_y], radius), distance, include_boundary) 
  count += success
  print("trying to add", center_x, center_y, radius, "now at", count)
  failed += 1 - success
  
  if failed > max_fails:
    print("Resetting...")
    resets += 1
    count = 0
    failed = 0
    g.reset()

g.apply_shapes()
g.write_to(outfile_path)
print("Wrote grid of size", num_cols, "x", num_rows, "after", resets, "resets with", count, "insertions and", failed, "failures!")
