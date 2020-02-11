#!/bin/bash
case $# in
	0|1)
		echo "Usage: runner_b3d.sh ifirst ifinal  // runs from i=ifirst to <=ifinal";
  	exit 1 ;;
	2)
		firsti=$1
		lasti=$2
		for ((ii=${firsti};ii<=${lasti};ii++))
		do
			echo "____________ b3d for run number " ${ii} ______________;
			./b3d ${ii};
		done
esac
