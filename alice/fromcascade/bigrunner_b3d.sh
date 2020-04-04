#!/bin/bash
make b3d_fromcascade
nproc=24
iproc0=0
iprocf=`expr ${iproc0} + ${nproc}`
for ((i=iproc0;i<procf;iproc+=1))
do
	`./runner_b3d.sh ${i} > logfiles/b3d_${i}.txt &` ;
done
