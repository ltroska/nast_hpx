#!/usr/bin/env python
from __future__ import print_function
import sys
import re
import csv

def get_time(output):
	pattern = '\/threads\{locality#0\/total\}\/time\/overall,[0-9\.]*,([0-9\.]*),\[s\],[0-9\.]*,\[ns\]\n'
	m = re.search(pattern, output)
	return m.group(1)

def get_idle_rate(output):
	pattern = '\/threads\{locality#0\/total\}\/idle-rate,.,[0-9\.]*,\[s\],([0-9]*),\[0.01%\]\n'
	m = re.search(pattern, output)
	return m.group(1) 

def get_rate(output):
	pattern = 'Rate \(MB\/s\):\t([0-9.]*?)\n'
	m = re.search(pattern, output)
	return m.group(1)

if len(sys.argv) < 4:
	print("Usage:", sys.argv[0] , "<task name> <id low> <id high> [nodes low] [proc low]", file=sys.stderr)
	sys.exit(0)

infile_path = sys.argv[1]
id_low = int(sys.argv[2])
id_high = int(sys.argv[3])
nodes_low = 1
procs_low = 1

outfile_path = infile_path

if len(sys.argv) > 4:
	nodes_low = int(sys.argv[4])

if len(sys.argv) > 5:
	procs_low = int(sys.argv[5])


nodes = nodes_low
procs = procs_low

i = 0
with open(infile_path + '_data.dat', 'wb') as outfile:	
	for id in range(id_low, id_high + 1):
		with open(infile_path + '.o' + str(id), 'rb') as infile:
			output = infile.read()
			outfile.write(str(nodes * procs) + ' ' + get_time(output) + '\n')

		if procs < 16:
			procs *= 2
		elif i == 0:
			nodes = 2
			i += 1 
		else:
			i += 1
			nodes = i ** 2
			


