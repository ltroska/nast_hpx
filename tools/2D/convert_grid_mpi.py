#!/usr/bin/env python
from __future__ import print_function
import sys

import csv

if len(sys.argv) != 3:
	print("Usage:", sys.argv[0] , "<input grid> <output grid>", file=sys.stderr)
	sys.exit(0)

infile_path = sys.argv[1]

in_grid = []

imax = 0
jmax = 0

with open(infile_path, 'rb') as csvfile:
	gridreader = csv.reader(csvfile, delimiter=',', quotechar='"')
	
	for row in gridreader:
		in_grid.append(row)
		jmax += 1
		imax = len(row)

outfile_path = sys.argv[2]

jmax -= 1

with open(outfile_path, 'wb') as csvfile:
	gridwriter = csv.writer(csvfile, delimiter=',', quotechar='"')

	for row in range(len(in_grid)):
		out_row = [jmax]
		jmax -= 1

		for col in range(len(in_grid[0])):

			cell = 1 - int(in_grid[row][col])
			#is_fluid bit
			type_int = 16*cell	
		
			count = 0
			#has_fluid_X bits
			#north
			if row != 0 and in_grid[row-1][col] != "1":
				type_int += 1
				count += 1
			#south
			if row != len(in_grid) - 1 and in_grid[row+1][col] != "1":
				type_int += 2
				count += 1
			#west
			if col != 0 and in_grid[row][col-1] != "1":
				type_int += 4
				count += 1
			#east
			if col != len(in_grid[0]) - 1 and in_grid[row][col+1] != "1":
				type_int += 8
				count += 1

			other_count = 0
			#counting fluid neighbors
			#north west
			if row != 0 and col != 0 and in_grid[row-1][col-1] != "1":
				other_count += 1		 
			#north east
			if row != 0 and col != len(in_grid[0]) - 1 and in_grid[row-1][col+1] != "1":
				other_count += 1
			#south west
			if row != len(in_grid) - 1 and col != 0 and in_grid[row+1][col-1] != "1":
				other_count += 1
			#south east
			if row != len(in_grid) - 1 and col != len(in_grid[0]) - 1 and in_grid[row+1][col+1] != "1":
				other_count += 1
			
			if count > 2 and cell == 0:
				print("FAILURE (illegal grid): obstacle cell is neighbors with more than two fluid cells at index i =",col, ", j =",row, " (origin top left)!", file=sys.stderr)
				#sys.exit(0)
			
			if count + other_count == 0:
				type_int = 0

			out_row.append(type_int)
		
		gridwriter.writerow(out_row)

	lastrow = [None]
	for i in range(imax):
		lastrow.append(i)

	gridwriter.writerow(lastrow)		
			


