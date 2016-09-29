import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.ticker import AutoMinorLocator
import sys
import glob
import re
import math

dpi = 600
directory = 'heatmaps/'

if len(sys.argv) < 3:
	print 'Usage', sys.argv[0], '<geometry file> [data]...'
	sys.exit(1)

def cmpit(a):
	pattern = '.*_([0-9]*)x[0-9]*\(([0-9]*)x[0-9]*\)\..*'
	return int(re.search(pattern, a).group(1) + re.search(pattern,a).group(2))

def read_data(infile):
	processors = []
	idle_rates = []

	for line in infile:
		splitline = line.split()
		nloc = max(int(splitline[0])/16, 1)
		if nloc == 2:
			nloc_x = 2
			nloc_y = 1
		else:
			nloc_x = nloc_y = int(math.sqrt(nloc))
		processors.append([nloc_x, nloc_y])
		
		offset = 7
		idle_rate = []
		for idx in range(nloc):
			idle_rate.append(float(splitline[offset + idx * 3]) / 100)
		
		idle_rates.append(np.array(idle_rate).reshape(nloc_y, nloc_x))
	return processors, idle_rates

geometry_file = sys.argv[1]
geometry_img = mpimg.imread(geometry_file)
aspect_ratio = float(len(geometry_img)) / len(geometry_img[0])

infile_paths = []
labels = []
pattern = '.*\/(.*)\..*'

if len(sys.argv) == 3:
	if sys.argv[2].endswith('.dat'):
		infile_paths.append(sys.argv[2])
	else:
		infile_paths = sorted(glob.glob(sys.argv[2] + ('*.dat' if sys.argv[2].endswith('/') else '/*.dat')), key=cmpit)
else:
	for idx in range(2, len(sys.argv)):
		infile_paths.append(sys.argv[idx])
		
	infile_paths = sorted(infile_paths, key=cmpit)

for path in infile_paths:
	labels.append(re.search(pattern, path).group(1))

for infile_path, label in zip(infile_paths, labels):
	with open(infile_path, 'rb') as infile:	
		print 'processing', infile_path
		processors, idle_rates = read_data(infile)

		for proc, idle_data in zip(processors, idle_rates):
			nloc_x, nloc_y = proc

			print 'processors', nloc_x, nloc_y

			fig, ax = plt.subplots()

			#imshow portion
			im = plt.imshow(idle_data, interpolation='none', aspect=aspect_ratio, cmap='RdBu_r')
			plt.imshow(geometry_img, alpha=0.75, extent=[-0.5, nloc_x-0.5, -0.5, nloc_y-0.5], aspect=aspect_ratio)

			#text portion
			x_ind_array = np.arange(0, nloc_x)
			y_ind_array = np.arange(0, nloc_y)
			x, y = np.meshgrid(x_ind_array, y_ind_array)

			for x_val, y_val in zip(x.flatten(), y.flatten()):
			    c = '%.2f' % idle_data[y_val, x_val]
			    ax.text(x_val, y_val, c, va='center', ha='center', backgroundcolor='w')

			#set tick marks for grid
			ax.set_xticks(np.arange(0, nloc_x))
			ax.set_yticks(np.arange(0, nloc_y))
			ax.set_xticklabels(np.arange(1, nloc_x + 1))
			ax.set_yticklabels(np.arange(1, nloc_x + 1))
			ax.set_xlim(-0.5, nloc_x-0.5)
			ax.set_ylim(-0.5, nloc_y-0.5)
			minor_locator = AutoMinorLocator(2)
			if nloc_x >= 2:
				ax.xaxis.set_minor_locator(minor_locator)
			if nloc_y >= 2:
				ax.yaxis.set_minor_locator(minor_locator)
			ax.tick_params(length=0)
			#ax.grid(which='minor', linestyle='-')
			ax.tick_params(axis=u'both', which=u'both',length=0)
			if aspect_ratio < 1:
				cbar = plt.colorbar(im, orientation='horizontal')
				cbar.ax.set_xlabel('Idle Rate (%)')
			else:
				cbar = plt.colorbar(im)
				cbar.ax.set_ylabel('Idle Rate (%)')

			cbar.set_alpha(1)
			cbar.draw_all()

			plt.savefig(directory + label + '_' + str(proc[0]) + 'x' + str(proc[1]) + '_heatmap.png', bbox_inches='tight', format='png', dpi=dpi)
