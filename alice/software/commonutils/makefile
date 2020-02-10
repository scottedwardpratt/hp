include ../../makefile_defs.mk
#e.g. /Users/scottepratt/local
INC=-I include -I ${MADAI_GSLPATH}/include -I ${MADAI_EIGEN_HOME} -I ../b3d2/include -I ../hydrobalance/include -I ../coral/include

INSTALL_HFILES_COMMONUTILS = ${MADAI_INSTALLDIR}/include/arrays.h ${MADAI_INSTALLDIR}/include/commondefs.h ${MADAI_INSTALLDIR}/include/constants.h ${MADAI_INSTALLDIR}/include/gslmatrix.h ${MADAI_INSTALLDIR}/include/misc.h ${MADAI_INSTALLDIR}/include/misc.h ${MADAI_INSTALLDIR}/include/parametermap.h ${MADAI_INSTALLDIR}/include/randy.h ${MADAI_INSTALLDIR}/include/sf.h ${MADAI_INSTALLDIR}/include/qualifier.h

HFILES_COMMONUTILS = include/arrays.h include/commondefs.h include/constants.h include/gslmatrix.h include/misc.h include/parametermap.h include/randy.h include/sf.h include/qualifier.h

OBJFILES =  build/misc.o build/parametermap.o build/gslmatrix.o build/ylm.o build/bess.o build/chcalc.o build/coulwave.o build/mclist.o build/ylmarray.o build/arraycalc.o build/mcprlist.o build/charray.o build/threed.o build/randy.o build/qualifier.o

madai_common : lib/libmadai_common.a bin/bmaker

install : installdirs ${MADAI_INSTALLDIR}/lib/libmadai_common.a ${INSTALL_HFILES_COMMONUTILS} ${MADAI_INSTALLDIR}/bin/bmaker

installdirs :
	mkdir -p ${MADAI_INSTALLDIR}/include;\
	mkdir -p ${MADAI_INSTALLDIR}/lib;\
	mkdir -p ${MADAI_INSTALLDIR}/bin;\
	mkdir -p ${MADAI_INSTALLDIR}/progdata/resinfo;

${MADAI_INSTALLDIR}/lib/libmadai_common.a : lib/libmadai_common.a
	/bin/cp -f lib/libmadai_common.a ${MADAI_INSTALLDIR}/lib/

lib/libmadai_common.a : ${OBJFILES}
	rm -f lib/libmadai_common.a;\
	ar -ru lib/libmadai_common.a ${OBJFILES};

clean :
	rm -f build/*.o  lib/*.a bin/*

uninstall :
	/bin/rm -f ${MADAI_INSTALLDIR}/lib/libmadai_common.a 
	/bin/rm -f ${MADAI_INSTALLDIR}/include/commondefs.h ${MADAI_INSTALLDIR}/include/constants.h ${MADAI_INSTALLDIR}/include/gslmatrix.h ${MADAI_INSTALLDIR}/include/misc.h ${MADAI_INSTALLDIR}/include/randy.h ${MADAI_INSTALLDIR}/include/sf.h ${MADAI_INSTALLDIR}/include/parametermap.h ${MADAI_INSTALLDIR}/include/qualifier.h ${MADAI_INSTALLDIR}/include/arrays.h 
	/bin/rm -f ${MADAI_INSTALLDIR}/bin/bmaker

########### BUILD OBJECT FILES ###########

build/misc.o : ${HFILES_COMMONUTILS} src/Misc/misc.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/misc.o src/Misc/misc.cc
build/randy.o : ${HFILES_COMMONUTILS} src/NumMath/randy.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/randy.o src/NumMath/randy.cc
build/parametermap.o : ${HFILES_COMMONUTILS} src/ParameterMap/parametermap.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/parametermap.o src/ParameterMap/parametermap.cc
build/gslmatrix.o : ${HFILES_COMMONUTILS} src/GSLMatrix/gslmatrix.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/gslmatrix.o src/GSLMatrix/gslmatrix.cc
build/ylm.o : ${HFILES_COMMONUTILS} src/SpecialFunctions/Ylm/ylm.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/ylm.o src/SpecialFunctions/Ylm/ylm.cc
build/bess.o : ${HFILES_COMMONUTILS} src/SpecialFunctions/Bessel/bess.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/bess.o src/SpecialFunctions/Bessel/bess.cc
build/chcalc.o : ${HFILES_COMMONUTILS} src/SpecialFunctions/CartesianHarmonics/chcalc.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/chcalc.o src/SpecialFunctions/CartesianHarmonics/chcalc.cc
build/coulwave.o : ${HFILES_COMMONUTILS} src/SpecialFunctions/CoulWave/coulwave.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/coulwave.o src/SpecialFunctions/CoulWave/coulwave.cc
build/mclist.o : ${HFILES_COMMONUTILS} src/Arrays/mclist.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/mclist.o src/Arrays/mclist.cc
build/ylmarray.o : ${HFILES_COMMONUTILS} src/Arrays/ylmarray.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/ylmarray.o src/Arrays/ylmarray.cc
build/arraycalc.o : ${HFILES_COMMONUTILS} src/Arrays/arraycalc.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/arraycalc.o src/Arrays/arraycalc.cc
build/mcprlist.o : ${HFILES_COMMONUTILS} src/Arrays/mcprlist.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/mcprlist.o src/Arrays/mcprlist.cc
build/charray.o : ${HFILES_COMMONUTILS} src/Arrays/charray.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/charray.o src/Arrays/charray.cc
build/threed.o : ${HFILES_COMMONUTILS} src/Arrays/threed.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/threed.o src/Arrays/threed.cc
build/qualifier.o : ${HFILES_COMMONUTILS} src/Qualifiers/qualifier.cc
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC} -o build/qualifier.o src/Qualifiers/qualifier.cc
		
############ INSTALL HEADER FILES ############

${MADAI_INSTALLDIR}/include/arrays.h : include/arrays.h
	/bin/cp -f include/arrays.h ${MADAI_INSTALLDIR}/include/arrays.h
${MADAI_INSTALLDIR}/include/commondefs.h : include/commondefs.h
	/bin/cp -f include/commondefs.h ${MADAI_INSTALLDIR}/include/commondefs.h
${MADAI_INSTALLDIR}/include/constants.h : include/constants.h
	/bin/cp -f include/constants.h ${MADAI_INSTALLDIR}/include/constants.h
${MADAI_INSTALLDIR}/include/gslmatrix.h : include/gslmatrix.h
	/bin/cp -f include/gslmatrix.h ${MADAI_INSTALLDIR}/include/gslmatrix.h
${MADAI_INSTALLDIR}/include/misc.h : include/misc.h
	/bin/cp -f include/misc.h ${MADAI_INSTALLDIR}/include/misc.h
${MADAI_INSTALLDIR}/include/parametermap.h : include/parametermap.h
	/bin/cp -f include/parametermap.h ${MADAI_INSTALLDIR}/include/parametermap.h
${MADAI_INSTALLDIR}/include/randy.h : include/randy.h
	/bin/cp -f include/randy.h ${MADAI_INSTALLDIR}/include/randy.h
${MADAI_INSTALLDIR}/include/sf.h : include/sf.h
	/bin/cp -f include/sf.h ${MADAI_INSTALLDIR}/include/sf.h
${MADAI_INSTALLDIR}/include/qualifier.h : include/qualifier.h
	/bin/cp -f include/qualifier.h ${MADAI_INSTALLDIR}/include/qualifier.h

############### EXECUTABLES #################

${MADAI_INSTALLDIR}/bin/bmaker : bin/bmaker
	/bin/cp -f bin/bmaker ${MADAI_INSTALLDIR}/bin/
bin/bmaker : src/BMaker/bmaker.cc
	${MADAI_CPP} ${MADAI_CFLAGS} -o bin/bmaker src/BMaker/bmaker.cc ${INC}
