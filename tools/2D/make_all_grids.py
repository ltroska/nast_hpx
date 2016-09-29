#!/usr/bin/env python
from subprocess import Popen
import sys
import random

init_size = 10
max_iter = 5 

num_beads = 60

processes = []

#DRIVEN
size = init_size
for _ in range(max_iter):
	str_size = str(size)
	commands = ["python make_rectangular_grid.py " + str_size + " " + str_size + " driven_cavity_" + str_size + "_" + str_size + ".csv"]
	
	processes.append([Popen(cmd, shell=True) for cmd in commands])


	size *= 2

for l in processes:
	for p in l:
		p.wait()

processes = []

size = init_size
for _ in range(max_iter):
	str_size = str(size)
	commands = ["python convert_grid.py driven_cavity_" + str_size + "_" + str_size + ".csv ../examples/geometry/driven_cavity_" + str_size + "_" + str_size + ".csv && rm driven_cavity_" + str_size + "_" + str_size + ".csv"]
	
	processes.append([Popen(cmd, shell=True) for cmd in commands])
	size *= 2

for l in processes:
	for p in l:
		p.wait()

#STEP
size = init_size
for _ in range(max_iter):
	if size >= 20:
		str_size = str(size)
		str_size_8 = str(size/8)
		str_size_4 = str(size/4)
		commands = ["python step.py step_" + str_size + "_" + str_size_4 + ".csv " + str_size + " " + str_size_4 + " " + str_size_4 + " " + str_size_8 + " " + str_size_8 + " " + str_size_4 + " " + str_size_8]
	
		processes.append([Popen(cmd, shell=True) for cmd in commands])


	size *= 2

for l in processes:
	for p in l:
		p.wait()

processes = []

size = init_size
for _ in range(max_iter):
	if size >= 5040:
		str_size = str(size)
		str_size_4 = str(size/4)
		commands = ["python convert_grid.py step_" + str_size + "_" + str_size_4 + ".csv ../examples/geometry/step_" + str_size + "_" + str_size_4 + ".csv && rm step_" + str_size + "_" + str_size_4 + ".csv"]
	
		processes.append([Popen(cmd, shell=True) for cmd in commands])
	size *= 2

for l in processes:
	for p in l:
		p.wait()

#Rectangles
size = init_size
for _ in range(max_iter):
	str_size = str(size)
	str_size_8 = str(size/8)
	str_size_40 = str(size/40)
	commands = ["python uniform_rectangular_obstacles.py uniform_rect_" + str_size + "_" + str_size + ".csv " + str_size + " " + str_size + " " + str_size_8 + " " + str_size_8 + " " + str_size_40 + " 0 1 0"]

	processes.append([Popen(cmd, shell=True) for cmd in commands])


	size *= 2

for l in processes:
	for p in l:
		p.wait()

processes = []

size = init_size
for _ in range(max_iter):
	str_size = str(size)
	str_size_4 = str(size/4)
	commands = ["python convert_grid.py uniform_rect_" + str_size + "_" + str_size + ".csv  ../examples/geometry/uniform_rect_" + str_size + "_" + str_size + ".csv  && rm uniform_rect_" + str_size + "_" + str_size + ".csv"]

	processes.append([Popen(cmd, shell=True) for cmd in commands])
	size *= 2

for l in processes:
	for p in l:
		p.wait()

#Porous
size /= 2
str_size = str(size)
str_size_8 = str(size/8)
str_size_40 = str(size/40)

cmd = ["python random_bead_pack.py bead_pack_" + str_size + "_" + str_size + ".csv " + str_size + " " + str_size + " " + str(num_beads) + " 4 " + str_size_8 + " " + str_size_40 + " 0 " + str(sys.maxint)]


p = Popen(cmd, shell=True)
p.wait()

cmd = ["python convert_grid.py bead_pack_" + str_size + "_" + str_size + ".csv  ../examples/geometry/bead_pack_" + str_size + "_" + str_size + ".csv  && rm bead_pack_" + str_size + "_" + str_size + ".csv"]

p = Popen(cmd, shell=True)
p.wait()

