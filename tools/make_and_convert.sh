#!/bin/bash

if [ $# -ne 3 ]
	then
		echo "Usage:" $0 "name x y"
		exit
fi

python make_rectangular_grid.py $2 $3 ../examples/geometry_orig/$1
python convert_grid.py ../examples/geometry_orig/$1 ../examples/geometry/$1
