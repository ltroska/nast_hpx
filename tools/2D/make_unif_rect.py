#!/usr/bin/env python
from subprocess import Popen
import sys
import random

init_size = 3360
init_rect = 200
max_iter = 4 

num_beads = 60

processes = []


#Rectangles
size = init_size
rect = init_rect
for _ in range(max_iter):
	str_size = str(size)
	str_rect = str(rect)
	str_rect_10 = str(rect/10)
	cmd = "python uniform_rectangular_obstacles.py uniform_rect_" + str_size + "_" + str_size + ".csv " + str_size + " " + str_size + " " + str_rect + " " + str_rect + " " + str_rect + " 0 1 0"

	processes.append([Popen(cmd, shell=True)])
	cmd = "python uniform_rectangular_obstacles.py uniform_rect_small_" + str_size + "_" + str_size + ".csv " + str_size + " " + str_size + " " + str_rect + " " + str_rect + " " + str_rect_10 + " 0 1 0"

	processes.append([Popen(cmd, shell=True)])

	size *= 2
	rect *= 2

for l in processes:
	for p in l:
		p.wait()

processes = []

size = init_size
for _ in range(max_iter):
	str_size = str(size)
	str_size_4 = str(size/4)
	cmd = "python convert_grid.py uniform_rect_small_" + str_size + "_" + str_size + ".csv  ../examples/geometry/uniform_rect_small_" + str_size + "_" + str_size + ".csv  && rm uniform_rect_small_" + str_size + "_" + str_size + ".csv"

	processes.append([Popen(cmd, shell=True)])
	size *= 2

for l in processes:
	for p in l:
		p.wait()


