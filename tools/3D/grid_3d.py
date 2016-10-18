#!/usr/bin/env python
import numpy as np
import itertools

import csv

HAS_FLUID_LEFT = 2**0
HAS_FLUID_RIGHT = 2**1
HAS_FLUID_BOTTOM = 2**2
HAS_FLUID_TOP = 2**3
HAS_FLUID_FRONT = 2**4
HAS_FLUID_BACK = 2**5
IS_FLUID = 2**6
IS_OBSTACLE = 2**7
IS_BOUNDARY = 2**8

def reverse_enumerate(iterable):
    """
    Enumerate over an iterable in reverse order while retaining proper indexes
    """
    return itertools.izip(reversed(xrange(len(iterable))), reversed(iterable))


class grid:	
	def __init__(self, i_max, j_max, k_max, x_length = 1, y_length = 1, z_length = 1, all_fluid=True):
		self.__data = [[[1 if all_fluid else 0 for i in range(i_max + 2)] for j in range(j_max + 2)] for k in range(k_max + 2)]
		self.__i_max = i_max
		self.__j_max = j_max
		self.__k_max = k_max
		self.__x_length = x_length
		self.__y_length = y_length
		self.__z_length = z_length

		for k in range(self.__k_max + 2):
			for j in range(self.__j_max + 2):
				for i in range(self.__i_max + 2):
					if i == 0 or j == 0 or k == 0 or i == self.__i_max + 1 or j == self.__j_max + 1 or k == self.__k_max + 1:
						self[i, j, k] = 0
		
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

	def insert_rectangle(self, i_begin, i_end, j_begin, j_end, k_begin, k_end):
		for k in range(k_begin, k_end + 1):
			for j in range(j_begin, j_end + 1):
				for i in range(i_begin, i_end + 1):
					self[i, j, k] = 0
		
	def write_to(self, file_path):
		with open(file_path, 'wb') as csvfile:
			csvfile.write(str(self.__i_max)+'\n')
			csvfile.write(str(self.__j_max)+'\n')
			csvfile.write(str(self.__k_max)+'\n')
			csvfile.write(str(self.__x_length)+'\n')
			csvfile.write(str(self.__y_length)+'\n')
			csvfile.write(str(self.__z_length)+'\n')
			gridwriter = csv.writer(csvfile, delimiter=',', quotechar='"')

			for k, plane in enumerate(self.__data):
				for j, row in enumerate(plane):
					outrow = [0] * len(row)
							
					for i, cell in enumerate(row):
						outvalue = 0
						
						if cell == 1:
							outvalue = IS_FLUID
						elif i == 0 or i == self.__i_max + 1 or j == 0 or j == self.__j_max + 1 or k == 0 or k == self.__k_max + 1:
							outvalue = IS_BOUNDARY + IS_OBSTACLE
						else:
							outvalue = IS_OBSTACLE
						
						if i != 0 and self[i - 1, j, k] == 1:
							outvalue += HAS_FLUID_LEFT
						if i != self.__i_max + 1 and self[i + 1, j, k] == 1:
							outvalue += HAS_FLUID_RIGHT
						if j != 0 and self[i, j - 1, k] == 1:
							outvalue += HAS_FLUID_FRONT
						if j != self.__j_max + 1 and self[i, j + 1, k] == 1:
							outvalue += HAS_FLUID_BACK
						if k != 0 and self[i, j, k - 1] == 1:
							outvalue += HAS_FLUID_BOTTOM
						if k != self.__k_max + 1 and self[i, j, k + 1] == 1:
							outvalue += HAS_FLUID_TOP
							
						outrow[i] = outvalue				  
				  
					gridwriter.writerow(outrow)

	def add_xy_plane(self, k, plane, invert=True):
		add = [[1 - int(x) if invert else int(x) for x in row] for row in plane]
		self.__data[k] = add			
				
