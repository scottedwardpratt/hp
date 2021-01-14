#!/bin/bash
make hydro2uds
nproc=24
nruns=30
for ((i=0;i<${nproc};i+=1))
do
	firsti=`expr ${i} \* ${nruns}`;
	lasti=`expr ${firsti} + ${nruns} - 1`;
	echo firsti=${firsti}, lasti=${lasti};
	`./runner_both.sh ${firsti} ${lasti} > logfiles/both_${firsti}_${lasti}.txt &` ;
done
