from __future__ import print_function
import sys

import csv

if len(sys.argv) != 3:
	print("Usage:", sys.argv[0] , "<input grid> <output grid>", file=sys.stderr)
	sys.exit(0)

infile_path = sys.argv[1]

in_grid = []

with open(infile_path, 'rb') as csvfile:
	gridreader = csv.reader(csvfile, delimiter=',', quotechar='"')
	
	for row in gridreader:
		in_grid.append(row)

out_grid = []

for row in range(len(in_grid)):
	out_row = []

	for col in range(len(in_grid[0])):
		cell = 1 - int(in_grid[row][col])

		type_int = 16*cell
		
		#north
		if row != 0 and in_grid[row-1][col] != "1":
			type_int += 1

		#south
		if row != len(in_grid) - 1 and in_grid[row+1][col] != "1":
			type_int += 2

		#west
		if col != 0 and in_grid[row][col-1] != "1":
			type_int += 4
		
		#east
		if col != len(in_grid[0]) - 1 and in_grid[row][col+1] != "1":
			type_int += 8

		out_row.append(type_int)
	
	out_grid.append(out_row)		
		

outfile_path = sys.argv[2]

with open(outfile_path, 'wb') as csvfile:
	gridwriter = csv.writer(csvfile, delimiter=',', quotechar='"')

	for row in out_grid:
		gridwriter.writerow(row)
