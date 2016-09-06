#!/usr/bin/env python
from __future__ import print_function
import sys
import numpy as np

import csv

class grid:	
	def __init__(self, i_max, j_max, k_max):
		self.__data = [[[1 for k in range(k_max + 2)] for j in range(j_max + 2)] for i in range(i_max + 2)]
		self.__i_max = i_max
		self.__j_max = j_max
		self.__k_max = k_max
		
	def __getitem__(self, tup):
		i, j, k = tup
		return self.__data[k][j][i]
  
	def __setitem__(self, tup, val):
		i, j, k = tup
		self.__data[k][j][i] = val
		
	def __str__(self):
		output = ""
		for k, plane in enumerate(self.__data):
			for j, row in enumerate(plane):
			  output += str(row)
			  if j != len(plane) - 1:
				output += "\n"
				
			output += "\n"
			output += "\n"
		
		return output
		
	def write_to(self, file_path):
		with open(file_path, 'wb') as csvfile:
			gridwriter = csv.writer(csvfile, delimiter=',', quotechar='"')

			for k, plane in enumerate(self.__data):
				for j, row in enumerate(plane):
					outrow = [0] * len(row)
										
					for i, cell in enumerate(row):
						outvalue = 0
						
						if cell == 1:
							outvalue = 2**6
						else:
							outvalue = 2**7	
						
						if i != 0 and self[i - 1, j, k] == 1:
							outvalue += 2**0
						if i != self.__i_max + 1 and self[i + 1, j, k] == 1:
							outvalue += 2**1
						if j != 0 and self[i, j - 1, k] == 1:
							outvalue += 2**4
						if j != self.__j_max + 1 and self[i, j + 1, k] == 1:
							outvalue += 2**5
						if k != 0 and self[i, j, k - 1] == 1:
							outvalue += 2**2
						if k != self.__k_max + 1 and self[i, j, k + 1] == 1:
							outvalue += 2**3
							
						outrow[i] = outvalue				  
				  
					gridwriter.writerow(outrow)
		

if len(sys.argv) != 5:
	print("Usage:", sys.argv[0] , "i_max j_max k_max <output grid>", file=sys.stderr)
	sys.exit(0)

i_max = int(sys.argv[1])
j_max = int(sys.argv[2])
k_max = int(sys.argv[3])

outfile_path = sys.argv[4]

b = grid(i_max, j_max, k_max)

for k in range(k_max + 2):
	for j in range(j_max + 2):
		for i in range(i_max + 2):
			if i == 0 or j == 0 or k == 0 or i == i_max + 1 or j == j_max + 1 or k == k_max + 1:
				b[i, j, k] = 0
								
b.write_to(outfile_path)
				

				
