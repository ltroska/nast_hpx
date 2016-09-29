#!/usr/bin/env python
from __future__ import print_function
import sys
import random
import csv
import subprocess

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

if len(sys.argv) != 3:
	print("Usage:", sys.argv[0] , "<input grid> <output grid>", file=sys.stderr)
	sys.exit(0)

infile_path = sys.argv[1]
outfile_path = sys.argv[2]

num_rows = file_len(infile_path)
num_cols = 0

with open(infile_path, 'rb') as csvfile:
	gridreader = csv.reader(csvfile, delimiter=',', quotechar='"')

	with open(outfile_path, 'wb') as csvfile:
		gridwriter = csv.writer(csvfile, delimiter=',', quotechar='"')
	
		for i, row in enumerate(gridreader):
			if i == 0 or i == num_rows - 1:
				num_cols = len(row)
				outrow = [1] * 2 * len(row)
				gridwriter.writerow(outrow)			
			else:
				outrow = [1]
				outrow.append(row[1])
				for j in range(1, len(row) - 1):
					outrow.extend([row[j]]*2)
				outrow.append(row[len(row) - 2])
				outrow.append(1)

				gridwriter.writerow(outrow)
				gridwriter.writerow(outrow)
				if i == 1 or i == num_rows - 2:
					gridwriter.writerow(outrow)
				
			
