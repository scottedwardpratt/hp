#!/bin/bash
NEVENTS=250
case $# in
	0)
		echo "Usage: runner_balance.sh idefault // runs from idefault*1000 to idefault*1000 +999";
  	exit 1 ;;
	1)
		idefault=$1
		firsti=`expr ${idefault} \* ${NEVENTS}`
		lasti=`expr ${firsti} + ${NEVENTS}`
		echo idefault=${idefault}, firsti=${firsti}, lasti=${lasti};
		./b3d_fromcascade default_${idefault} ${firsti} ${lasti};
		./balance_fromcascade default_${idefault} ${firsti} ${lasti}
esac
