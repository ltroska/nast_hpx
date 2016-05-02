#!/usr/bin/env python
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

		if row == 0 and col != 0 and col != len(in_grid[0]) - 1:
			type_int = 0 + 1 + 4 + 8

		elif row == len(in_grid) - 1 and col != 0 and col != len(in_grid[0]) - 1:
			type_int = 0 + 2 + 4 + 8

		elif col == 0 and row != 0 and row != len(in_grid) - 1:
			type_int = 0 + 1 + 2 + 4

		elif col == len(in_grid[0]) - 1 and row != 0 and row != len(in_grid) - 1:
			type_int = 0 + 1 + 2 + 8

		else:
			count = 0
			#north
			if row != 0 and in_grid[row-1][col] != "1":
				type_int += 1
				count += 1 - cell
			 
			#south
			if row != len(in_grid) - 1 and in_grid[row+1][col] != "1":
				type_int += 2
				count += 1 - cell
			#west
			if col != 0 and in_grid[row][col-1] != "1":
				type_int += 4
				count += 1 - cell
			#east
			if col != len(in_grid[0]) - 1 and in_grid[row][col+1] != "1":
				type_int += 8
				count += 1 - cell
			
			if count > 2:
				print("FAILURE (illegal grid): obstacle cell is neighbors with more than two fluid cells at index i =",col, ", j =",row, " (origin top left)!", file=sys.stderr)
				sys.exit(0)

		out_row.append(type_int)
	
	out_grid.append(out_row)		
		

outfile_path = sys.argv[2]

with open(outfile_path, 'wb') as csvfile:
	gridwriter = csv.writer(csvfile, delimiter=',', quotechar='"')

	for row in out_grid:
		gridwriter.writerow(row)
