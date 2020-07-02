#!/bin/bash
make b3d_fromhydro
nproc=12
nruns=67
for ((i=0;i<${nproc};i+=1))
do
	firsti=`expr ${i} \* ${nruns}`;
	lasti=`expr ${firsti} + ${nruns}`;
	rm -f logfiles/b3d_${firsti}_${lasti}.txt;
	echo starting runs with firsti=${firsti}, lasti=${lasti};
	`./runner_b3d.sh ${firsti} ${lasti} > logfiles/b3d_${firsti}_${lasti}.txt &` ;
done
