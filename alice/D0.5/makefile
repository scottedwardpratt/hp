include ../makefile_defs.mk
BASELIBS = -lcoralutils -lgsl -lgslcblas -lrhicutils -larmadillo
INCLUDELIBDIRS = -I${MADAI_INSTALLDIR}/include/madai -I${MADAI_GSLPATH}/include -I${MADAI_EIGEN_HOME} -L${MADAI_GSLPATH}/lib -L${MADAI_INSTALLDIR}/lib/madai

corallibs :
	make -C ${MADAI_HOME}/rhic/coral install

rhicutilslibs :
	make -C ${MADAI_HOME}/rhic/rhicutils install

b3dlibs : corallibs rhicutilslibs
	make -C ${MADAI_HOME}/rhic/b3d2 install

b3d : b3dmain.cc b3dlibs
	${CPP} -o b3d b3dmain.cc ${OPT} ${INCLUDELIBDIRS} -lb3d ${BASELIBS}


