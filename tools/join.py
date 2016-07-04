#!/usr/bin/env python
from __future__ import print_function
import sys

if len(sys.argv) < 3:
	print("Usage:", sys.argv[0] , "<outfile> [files]...", file=sys.stderr)
  	sys.exit(0)

outfile_path = sys.argv[1]
file_paths = []
for idx in range(2, len(sys.argv)):
	file_paths.append(sys.argv[idx])


processors = []
outrows = []

for idx, file_path in enumerate(file_paths):
	with open(file_path, 'rb') as infile:
		lines = infile.readlines()
		for idy, line in enumerate(lines):
			splitline = line.split()
			
			if idx == 0:
				processors.append(splitline[0])
				outrows.append([])
				
			outrows[idy].append(splitline[1])


with open(outfile_path, 'wb') as outfile:
	for idx, row in enumerate(outrows):
		outfile.write(processors[idx] + " " + " ".join(row) + "\n")
