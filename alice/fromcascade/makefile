include ../makefile_defs.mk

BASELIBS = -lmadai_common -lhydrobalance -lgsl -lgslcblas  -lb3d
INCLUDELIBDIRS = -I${MADAI_INSTALLDIR}/include -I${MADAI_GSLPATH}/include -I${MADAI_EIGEN_HOME} -I${MADAI_INSTALLDIR}/include/hydrobalance -L${MADAI_GSLPATH}/lib -L${MADAI_INSTALLDIR}/lib

balhydrolibs :
	make -C ../software install

hydro2uds : hydro2udsmain.cc balhydrolibs
	${MADAI_CPP} -o hydro2uds ${MADAI_CFLAGS} hydro2udsmain.cc ${INCLUDELIBDIRS} -lhydrobalance ${BASELIBS}

b3d : b3dmain.cc balhydrolibs
	${MADAI_CPP} -o b3d ${MADAI_CFLAGS} b3dmain.cc ${INCLUDELIBDIRS} -lb3d  ${BASELIBS}
