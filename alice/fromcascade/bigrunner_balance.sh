#!/bin/bash
make balance_fromcascade
nproc=24
for ((i=0;i<${nproc};i+=1))
do
	`./runner_balance.sh ${i} > logfiles/balance_${i}.txt &` ;
done
