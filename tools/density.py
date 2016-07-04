#!/usr/bin/env python
from __future__ import print_function
import csv
import sys

if len(sys.argv) != 2:
  print("Usage:", sys.argv[0], "<input grid>", file=sys.stderr)
  sys.exit(0)

infile_path = sys.argv[1]

m = sum(1 for line in open(infile_path))
n = 0

surface = 0
volume = 0
boundary = 0


with open(infile_path, 'rb') as csvfile:
  in_grid = csv.reader(csvfile, delimiter=',', quotechar='"')

  for row in in_grid:
    n = len(row)
    for col in row:
			if (int(col) & (1 << 5)) != 0:
				boundary += 1
			elif (int(col) & (1 << 4)) == 0 and col != "0":
				surface += 1
				volume += 1
			elif col == "0":
				volume += 1

	
total = m*n
 
print("This file has the following densities:")
print("Boundary:", 100 * boundary / float(total), "% (", boundary, "/", total, ")")
print("Obstacle Surface:", 100 * surface / float(total), "% (", surface, "/", total, ")")
print("Obstacle Volume:", 100 * volume / float(total), "% (", volume, "/", total, ")")
print("Total:", 100 * (boundary + volume) / float(total), "% (", boundary + volume, "/", total, ")")


