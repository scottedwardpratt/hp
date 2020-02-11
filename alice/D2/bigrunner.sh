#!/bin/bash
nproc=24
make b3d;
rm logfiles/crap*.dat
nruns=30;
for ((i=0;i<${nproc};i=i+1)) do
	event0=`expr ${i} \\* ${nruns}`
	eventf=`expr ${event0} + ${nruns} - 1`
	echo event0=${event0}
	echo eventf=${eventf}
	runner.sh ${event0} ${eventf} > logfiles/crap_${event0}_${eventf}.dat &
done

