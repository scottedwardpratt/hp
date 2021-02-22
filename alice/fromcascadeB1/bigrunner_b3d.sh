#! /bin/bash
case $# in
	0|1)
		echo "Usage: bigrunner_b3d.sh idefault0 NEVENTS  // runs from i=ifirst to <=ifinal";
  	exit 1 ;;
	2)
	NPROC=24
	idefault0=$1
	NEVENTS=$2
	idefaultf=`expr ${idefault0} + ${NPROC}`
	for ((idefault=${idefault0};idefault<${idefaultf};idefault+=1))
	do
		firsti=`expr ${idefault} \* ${NEVENTS}`;
		lasti=`expr ${firsti} + ${NEVENTS} - 1`;
		rm -f logfiles/b3d_${idefault}.txt;
		`./runner_b3d.sh ${idefault} ${firsti} ${lasti} > logfiles/b3d_${idefault}.txt &` ;
	done
esac