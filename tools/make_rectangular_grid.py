#!/usr/bin/env python
from __future__ import print_function
import sys

import csv

if len(sys.argv) != 4:
	print("Usage:", sys.argv[0] , "n m <output grid>", file=sys.stderr)
	sys.exit(0)

n = int(sys.argv[1])
m = int(sys.argv[2])
outfile_path = sys.argv[3]

with open(outfile_path, 'wb') as csvfile:
	gridwriter = csv.writer(csvfile, delimiter=',', quotechar='"')

	for j in range(m):
		row = []
		if j == 0 or j == m - 1:
			row = [1]*n

		else:
			row.append(1)
			row += [0]*(n-2)
			row.append(1)

		gridwriter.writerow(row)


