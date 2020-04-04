#!/bin/bash
NEVENTS=1000
case $# in
	0)
		echo "Usage: runner_balance.sh iproc0 // runs from idefault*1000 to idefault*1000 +999";
  	exit 1 ;;
	1)
		iproc0=$1
		nproc=24
		iprocf=`expr ${iproc0} + ${nproc}`
		make balance_fromcascade
		do
			`./runner_balance.sh ${i} > logfiles/balance_${i}.txt &` ;
		done
		
esac