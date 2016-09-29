#!/usr/bin/env python
from __future__ import print_function
from PIL import Image
import numpy as np
import csv
import sys

if len(sys.argv) != 4:
	print("Usage:", sys.argv[0] , "<input image> <output grid> <add boundary (1)/no boundary (0)>", file=sys.stderr)
	sys.exit(0)

add_boundary = int(sys.argv[3]) == 1
infile = sys.argv[1]
outfile = sys.argv[2]

col = Image.open(infile)
pixels = col.transpose(Image.FLIP_LEFT_RIGHT).transpose(Image.ROTATE_90).load()


n, m = col.size

with open(outfile, 'wb') as csvfile:
	gridwriter = csv.writer(csvfile, delimiter=',', quotechar='"')
	
	if add_boundary:
		gridwriter.writerow([1]*(n+2))

	for i in range(m):
		row = []
		if add_boundary:
			row.append(1)
		for j in range(n):
			row.append(1 if pixels[i,j] != (255, 255, 255, 255) else 0)

		if add_boundary:
			row.append(1)
		gridwriter.writerow(row)

	if add_boundary:
		gridwriter.writerow([1]*(n+2))

col.save("bw.png")

