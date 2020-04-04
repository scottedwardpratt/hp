#!/bin/bash
case $# in
	0)
		echo "Usage: runner_balance.sh idefault // runs from idefault*1000 to idefault*1000 +999";
  	exit 1 ;;
	1)
		idefault=$1
		NEVENTS=1000
		firsti=`expr ${idefault} \* ${NEVENTS}`
		lasti=`expr ${firsti} + ${NEVENTS}`
		echo idefault=${idefault}, firsti=${firsti}, lasti=${lasti}
		./balance_fromcascade default_${idefault} ${firsti} ${lasti}
		
esac
