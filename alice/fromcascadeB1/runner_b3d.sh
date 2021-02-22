#!/bin/bash
case $# in
	0|1|2)
		echo "Usage: idefault first lasti  // runs from i=ifirst to <=ifinal";
  	exit 1 ;;
	*)
		idefault=$1
		firsti=$2
		lasti=$3
		../bin/b3d default_${idefault} ${firsti} ${lasti};
esac
