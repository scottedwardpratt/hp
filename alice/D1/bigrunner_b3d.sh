#!/bin/bash
make b3d
nproc=24
nruns=30
for ((i=0;i<${nproc};i+=1))
do
	firsti=`expr ${i} \* ${nruns}`;
	lasti=`expr ${firsti} + ${nruns}`;
	echo firsti=${firsti}, lasti=${lasti};
	`./runner_b3d.sh ${firsti} ${lasti} > logfiles/b3d_${firsti}_${lasti}.txt &` ;
done
