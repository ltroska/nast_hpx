#!/usr/bin/env python
from __future__ import print_function
from PIL import Image
import csv
import sys

if len(sys.argv) != 3:
	print("Usage:", sys.argv[0], "<input grid> <output image>", file=sys.stderr)
	sys.exit(0)

infile_path = sys.argv[1]
outfile_path = sys.argv[2]

m = sum(1 for line in open(infile_path))
n = 0

img = None
pixels = None

with open(infile_path, 'rb') as csvfile:
  in_grid = csv.reader(csvfile, delimiter=',', quotechar='"')

  for j, row in enumerate(in_grid):
  if j == 0:
    n = len(row)	  
    img = Image.new('RGB', (n, m), "white")
    pixels = img.load()

  for i, col in enumerate(row):
    if col == '1':
      pixels[i,j] = (0, 0, 0)

img.save(outfile_path + ".png")
