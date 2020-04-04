#!/bin/bash
NEVENTS=4000
case $# in
	0)
		echo "Usage: runner_b3d.sh idefault // runs from idefault*100 to idefault*100+99";
  	exit 1 ;;
	1)
		idefault=$1
		firsti=`expr ${idefault} \* ${NEVENTS}`
		lasti=`expr ${firsti} + ${NEVENTS}`
		echo idefault=${idefault}, firsti=${firsti}, lasti=${lasti}
		./b3d_fromcascade default_${idefault} ${firsti} ${lasti}
		
esac
