#!/bin/bash
case $# in
0,1)  	echo "Usage: runner.sh i0 ifinal  // runs from i=i0 to i<ifinal";
	exit 1 ;;
2)
	i0=$1;
	ifinal=$2
	for ((i=${i0};i<=${ifinal};i=i+1)) do
		b3d default_${i} ${i} ${i};
	done
esac
