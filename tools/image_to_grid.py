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
gray = col.convert('L')
arr = np.array(gray)


m = len(arr)
n = len(arr[0])

with open(outfile, 'wb') as csvfile:
	gridwriter = csv.writer(csvfile, delimiter=',', quotechar='"')
	
	if add_boundary:
		gridwriter.writerow([1]*(n+2))

	for i in range(m):
		row = []
		if add_boundary:
			row.append(1)
		for j in range(n):
			row.append(0 if arr[i][j] > 128 else 1)

		if add_boundary:
			row.append(1)
		gridwriter.writerow(row)		

	if add_boundary:
		gridwriter.writerow([1]*(n+2))

#gray.point(lambda x: 0 if x<128 else 255, '1').save("result_bw.png")
