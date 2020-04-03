#!/bin/bash
make b3d_fromcascade
nproc=24
for ((i=0;i<${nproc};i+=1))
do
	`./runner_b3d.sh ${i} > logfiles/b3d_${i}.txt &` ;
done
