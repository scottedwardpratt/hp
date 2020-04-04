#!/bin/bash
make balance_fromcascade
nproc=24
iproc0=0
iprocf=`expr ${iproc0} + ${nproc}`
for ((i=iproc0;i<procf;iproc+=1))
do
	`./runner_balance.sh ${i} > logfiles/balance_${i}.txt &` ;
done