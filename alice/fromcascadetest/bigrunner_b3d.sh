#!/bin/bash
case $# in
	0)
		echo "Usage: bigrunner_b3d.sh iproc0 // runs from idefault*1000 to idefault*1000 +999";
  	exit 1 ;;
	1)
		iproc0=$1
		nproc=12
		iprocf=`expr ${iproc0} + ${nproc}`
		make b3d_fromcascade
		for((i=iproc0;i<iprocf;i++))
		do
			`./runner_b3d.sh ${i} > logfiles/b3d_${i}.txt &` ;
		done
		
esac
