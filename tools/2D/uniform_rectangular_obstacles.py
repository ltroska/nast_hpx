#!/usr/bin/env python
from __future__ import print_function
from grid import *
import sys
import random

if len(sys.argv) < 7:
  print("Usage:", sys.argv[0] , "<output grid> <#cols> <#rows> <x length of rectangles> <y length of rectangles> <distance> [offset=0] [center=0] [include_boundary=0]", file=sys.stderr)
  sys.exit(0)

outfile_path = sys.argv[1]
num_cols = int(sys.argv[2])
num_rows = int(sys.argv[3])
length_x = int(sys.argv[4])
length_y = int(sys.argv[5])
distance = int(sys.argv[6])
  
offset = 0
if len(sys.argv) >= 8:
  offset = int(sys.argv[7])
 
center = False
if len(sys.argv) >= 9:
  center = (int(sys.argv[8]) == 1)
  
include_boundary = False
if len(sys.argv) >= 10:
  include_boundary = (int(sys.argv[9]) == 1)
  
g = grid(num_cols, num_rows)
  
num_rectangles_x = (num_cols - 2) / (distance + length_x)
num_rectangles_y = (num_rows - 2) / (distance + length_y)

num_gaps_x = num_rectangles_x - 1
num_gaps_y = num_rectangles_y - 1

start_x = 1 + offset
start_y = 1

end_y = start_y + length_y * (num_rectangles_y + 1) - 1 + distance * num_rectangles_y

offset_y = (num_rows - 2 - end_y) / 2

if center:
  start_y += offset_y

for num_y in range(0, num_rectangles_y + 1):
  y = start_y + num_y * (distance + length_y)
  
  for num_x in range(0, num_rectangles_x + 1):
    x = start_x + num_x * (distance + length_x)
    
    
    g.add_shape(rectangle([x, y], [x + length_x - 1, y + length_y - 1]), distance, include_boundary)
  
g.apply_shapes()
g.write_to(outfile_path)
print("Wrote grid of size", num_cols, "x", num_rows, "with", num_rectangles_x, "x", num_rectangles_y, "rectangles!")