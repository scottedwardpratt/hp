include ../makefile_defs.mk

HBLIBS = -lgsl -lgslcblas -lhydrobalance -larmadillo
HBDIR = ../install
BASELIBS = -lcoralutils -lgsl -lgslcblas -lrhicutils
INCLUDELIBDIRS = -I${MADAI_INSTALLDIR}/include/madai -I${MADAI_GSLPATH}/include -I${MADAI_EIGEN_HOME} -L${MADAI_GSLPATH}/lib -L${MADAI_INSTALLDIR}/lib/madai


balhydrolibs :
	make -C ../software install

hydro2uds : hydro2uds.cc balhydrolibs
	${CPP} -o hydro2uds ${OPT} -I ${HBDIR}/include/hydrobalance -I ${GSLPATH}/include hydro2uds.cc -L ${HBDIR}/lib/hydrobalance -L ${MADAI_GSLPATH}/lib ${HBLIBS}

corallibs :
	make -C ${MADAI_HOME}/rhic/coral install

rhicutilslibs :
	make -C ${MADAI_HOME}/rhic/rhicutils install

b3dlibs : corallibs rhicutilslibs
	make -C ${MADAI_HOME}/rhic/b3d2 install

b3d : b3dmain.cc b3dlibs
	${CPP} -o b3d b3dmain.cc ${OPT} ${INCLUDELIBDIRS} -lb3d ${BASELIBS}
