include ../../makefile_defs.mk

#e.g. /Users/scottpratt/local
INC=-I include -I ${GSLPATH}/include -I ${MADAI_EIGEN_HOME} -I ../commonutils/include -I ../b3d2/include

hydrobalance : commonutils b3d2libs lib/libhydrobalance.a

commonutils : 
	make -C ../commonutils install

install : ${MADAI_INSTALLDIR}/lib/libhydrobalance.a ${MADAI_INSTALLDIR}/include/hydro2uds.h ${MADAI_INSTALLDIR}/include/charge.h ${MADAI_INSTALLDIR}/include/eos.h

clean :
	rm -f build/*.o lib/*.a

uninstall :
	rm -f ${MADAI_INSTALLDIR}/lib/libhydrobalance.a ${MADAI_INSTALLDIR}/include/hydro2uds.h ${MADAI_INSTALLDIR}/include/eos.h ${MADAI_INSTALLDIR}/include/charge.h

HFILES_HYDROBALANCE = include/hydro2uds.h include/eos.h include/charge.h ../b3d2/include/hyper.h ../b3d2/include/resonances.h

OBJFILES =  build/hbio.o build/hydro2uds.o build/eos.o build/mesh.o build/hbio.o build/hyperfind.o build/hbcharge.o
	
######################

lib/libhydrobalance.a : commonutils ${OBJFILES}
	rm -f lib/libhydrobalance.a;\
	ar -ru lib/libhydrobalance.a ${OBJFILES};

build/hydro2uds.o : src/hydro2uds.cc ${HFILES_HYDROBALANCE}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/hydro2uds.o src/hydro2uds.cc

build/hbio.o : src/hbio.cc ${HFILES_HYDROBALANCE}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/hbio.o src/hbio.cc

build/eos.o : src/eos.cc ${HFILES_HYDROBALANCE}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/eos.o src/eos.cc

build/mesh.o : src/mesh.cc ${HFILES_HYDROBALANCE}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/mesh.o src/mesh.cc

build/hyperfind.o : src/hyperfind.cc ${HFILES_HYDROBALANCE}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/hyperfind.o src/hyperfind.cc
	
build/hbcharge.o : src/hbcharge.cc ${HFILES_HYDROBALANCE}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/hbcharge.o src/hbcharge.cc
	
############

${MADAI_INSTALLDIR}/lib/libhydrobalance.a : lib/libhydrobalance.a
	cp -f lib/libhydrobalance.a ${MADAI_INSTALLDIR}/lib/

${MADAI_INSTALLDIR}/include/hydro2uds.h : include/hydro2uds.h
	cp -f include/hydro2uds.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/eos.h : include/eos.h
	cp -f include/eos.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/charge.h : include/charge.h
	cp -f include/charge.h ${MADAI_INSTALLDIR}/include/
