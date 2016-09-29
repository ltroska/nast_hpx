from __future__ import print_function
import sys
import csv
import subprocess

def file_len(fname):
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE, 
                                              stderr=subprocess.PIPE)
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])

if len(sys.argv) != 5:
	print("Usage:", sys.argv[0] , "<input grid> <output grid> <x_length> <y_length>", file=sys.stderr)
	sys.exit(0)

infile_path = sys.argv[1]
outfile_path = sys.argv[2]
x_length = sys.argv[3]
y_length = sys.argv[4]

num_rows = file_len(infile_path)
num_cols = 0

with open(infile_path, 'rb') as csvfile:
	gridreader = csv.reader(csvfile, delimiter=',', quotechar='"')
	with open(outfile_path, 'wb') as outfile:

		for i, row in enumerate(gridreader):
			outrow = []
			if i == 0:
				num_cols = len(row)
				outfile.write('size = ' + str(num_rows) + ' ' + str(num_cols) + '\n')
				outfile.write('length = ' + x_length + ' ' + y_length + '\n')
				outfile.write('pressure = 0\n')
				outfile.write('geometry = free\n')
				outrow = ['#'] * len(row)
			else:
				if row[0] != row[1]:
					outrow = ['H']
				else:
					outrow = ['#']

				for j in range(1, len(row) - 1):
					outrow.append('#' if row[j] == '1' else ' ')

				outrow.append('#')

			outfile.write(''.join(outrow) + '\n')
				
			


