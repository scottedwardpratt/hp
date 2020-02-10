include ../../makefile_defs.mk

# enter make all, make b3d, make install, make clean or make uninstall
### the following variables are set in the include statement above
#e.g. /opt/local
#e.g. /Users/scottpratt/local
GSLPATH=${MADAI_GSLPATH}
#e.g. /opt/local
#########################################################################
#HOPEFULLY, everything below need not be touched

b3d: commonutils b3d_dirs lib/libb3d.a ${HFILES_B3D2}

commonutils : 
	make -C ../commonutils install
	
INSTALL_HFILES_B3D2 = ${MADAI_INSTALLDIR}/include/b3d.h  ${MADAI_INSTALLDIR}/include/sampler.h ${MADAI_INSTALLDIR}/include/resonances.h  ${MADAI_INSTALLDIR}/include/resonances.h ${MADAI_INSTALLDIR}/include/part.h  ${MADAI_INSTALLDIR}/include/inelastic.h ${MADAI_INSTALLDIR}/include/balancearrays.h  ${MADAI_INSTALLDIR}/include/inelastic.h ${MADAI_INSTALLDIR}/include/acceptance.h  ${MADAI_INSTALLDIR}/include/inelastic.h ${MADAI_INSTALLDIR}/include/hyper.h  ${MADAI_INSTALLDIR}/include/inelastic.h ${MADAI_INSTALLDIR}/include/mutinfo.h ${MADAI_INSTALLDIR}/include/seinfo.h  ${MADAI_INSTALLDIR}/include/cell.h ${MADAI_INSTALLDIR}/include/regen.h ${MADAI_INSTALLDIR}/include/resonances.h  ${MADAI_INSTALLDIR}/progdata/resinfo/resonances_pdg.dat  ${MADAI_INSTALLDIR}/progdata/resinfo/resonances_pdg_weak.dat  ${MADAI_INSTALLDIR}/progdata/resinfo/decays_pdg.dat  ${MADAI_INSTALLDIR}/progdata/resinfo/decays_pdg_weak.dat  ${MADAI_INSTALLDIR}/progdata/resinfo/resonances_nodecays.dat  ${MADAI_INSTALLDIR}/progdata/resinfo/pionsonly.dat  ${MADAI_INSTALLDIR}/progdata/resinfo/decays_none.dat

install : installdirs ${MADAI_INSTALLDIR}/lib/libb3d.a ${INSTALL_HFILES_B3D2}

installdirs :
	mkdir -p ${MADAI_INSTALLDIR}/include;\
	mkdir -p ${MADAI_INSTALLDIR}/lib;\
	mkdir -p ${MADAI_INSTALLDIR}/bin;\
	mkdir -p ${MADAI_INSTALLDIR}/progdata/resinfo;\

clean :
	rm -f build/*.o lib/*.a

INC=-I include -I ${GSLPATH}/include -I ${MADAI_EIGEN_HOME} -I ../commonutils/include -I ../coral/include -I ../hydrobalance/include

uninstall :
	rm -f ${INSTALL_HFILES_B3D2}
	rm -f ${MADAI_INSTALLDIR}/lib/libb3d.a

HFILES_B3D2 = include/b3d.h\
include/sampler.h\
include/resonances.h\
include/inelastic.h\
include/part.h\
include/mutinfo.h\
include/seinfo.h\
include/regen.h\
include/cell.h\
include/balancearrays.h\
include/acceptance.h\
include/hyper.h

OBJFILES = build/b3d.o\
build/b3dsubs.o\
build/b3dio.o\
build/addaction.o\
build/cell.o\
build/collide.o\
build/scatter.o\
build/annihilate.o\
build/decay.o\
build/action.o\
build/action_perform.o\
build/action_perform_activate.o\
build/action_perform_collide.o\
build/action_perform_decay.o\
build/action_perform_regen.o\
build/action_perform_exitcell.o\
build/findcollision.o\
build/sampler.o\
build/mutinfo.o\
build/seinfo.o\
build/regen.o\
build/volume_element.o\
build/part.o\
build/inelastic.o\
build/resonances.o\
build/anal_spectra.o\
build/anal_hbt.o\
build/anal_v2.o\
build/consolidate.o\
build/balancearrays.o\
build/bfnumer.o\
build/acceptance.o\
build/acceptance_cheap.o\
build/acceptance_STAR.o\
build/acceptance_ALICE.o\
build/charge.o\
build/makepart.o\
build/hyper.o

#######################

b3d_dirs :
	mkdir -p build;\
	mkdir -p include;\
	mkdir -p lib;\
	mkdir -p ${MADAI_INSTALLDIR}/include;\
	mkdir -p ${MADAI_INSTALLDIR}/lib;\
	mkdir -p ${MADAI_INSTALLDIR}/progdata
	mkdir -p ${MADAI_INSTALLDIR}/progdata/resinfo

#######################

lib/libb3d.a : ${HFILES_B3D2} ${OBJFILES}
	rm -f lib/libb3d.a;\
	ar -ru lib/libb3d.a ${OBJFILES}

build/b3d.o : src/b3d.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/b3d.o src/b3d.cc

build/b3dsubs.o : src/b3dsubs.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/b3dsubs.o src/b3dsubs.cc
	
build/b3dio.o : src/b3dio.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/b3dio.o src/b3dio.cc

build/cell.o : src/cell.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/cell.o src/cell.cc

build/collide.o : src/collide.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/collide.o src/collide.cc

build/addaction.o : src/addaction.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/addaction.o src/addaction.cc

build/scatter.o : src/scatter.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/scatter.o src/scatter.cc

build/annihilate.o : src/annihilate.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/annihilate.o src/annihilate.cc

build/decay.o : src/decay.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/decay.o src/decay.cc

build/action.o : src/action.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/action.o src/action.cc

build/action_perform.o : src/action_perform.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/action_perform.o src/action_perform.cc

build/action_perform_activate.o : src/action_perform_activate.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/action_perform_activate.o src/action_perform_activate.cc

build/action_perform_collide.o : src/action_perform_collide.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/action_perform_collide.o src/action_perform_collide.cc

build/action_perform_decay.o : src/action_perform_decay.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/action_perform_decay.o src/action_perform_decay.cc

build/action_perform_regen.o : src/action_perform_regen.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/action_perform_regen.o src/action_perform_regen.cc

build/action_perform_exitcell.o : src/action_perform_exitcell.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/action_perform_exitcell.o src/action_perform_exitcell.cc

build/findcollision.o : src/findcollision.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/findcollision.o src/findcollision.cc

build/sampler.o : src/sampler.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/sampler.o src/sampler.cc

build/mutinfo.o : src/mutinfo.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/mutinfo.o src/mutinfo.cc

build/seinfo.o : src/seinfo.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/seinfo.o src/seinfo.cc
	
build/regen.o : src/regen.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/regen.o src/regen.cc

build/volume_element.o : src/volume_element.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/volume_element.o src/volume_element.cc

build/part.o : src/part.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/part.o src/part.cc

build/resonances.o : 	src/resonances.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/resonances.o src/resonances.cc

build/inelastic.o : src/inelastic.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/inelastic.o src/inelastic.cc
	
build/anal_spectra.o : src/anal_spectra.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/anal_spectra.o src/anal_spectra.cc

build/acceptance.o : src/acceptance.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/acceptance.o src/acceptance.cc

build/acceptance_cheap.o : src/acceptance_cheap.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/acceptance_cheap.o src/acceptance_cheap.cc

build/acceptance_STAR.o : src/acceptance_STAR.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/acceptance_STAR.o src/acceptance_STAR.cc

build/acceptance_ALICE.o : src/acceptance_ALICE.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/acceptance_ALICE.o src/acceptance_ALICE.cc

build/charge.o : src/charge.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/charge.o src/charge.cc

build/makepart.o : src/makepart.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/makepart.o src/makepart.cc

build/hyper.o : src/hyper.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/hyper.o src/hyper.cc
	
build/anal_hbt.o : src/anal_hbt.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/anal_hbt.o src/anal_hbt.cc

build/anal_v2.o : src/anal_v2.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/anal_v2.o src/anal_v2.cc

build/consolidate.o : src/consolidate.cc ${B3D_HFIES}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/consolidate.o src/consolidate.cc
		
build/balancearrays.o : src/balancearrays.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/balancearrays.o src/balancearrays.cc

build/bfnumer.o : src/bfnumer.cc ${HFILES_B3D2}
	${MADAI_CPP} -c ${MADAI_CFLAGS} ${INC}  -o build/bfnumer.o src/bfnumer.cc

############

${MADAI_INSTALLDIR}/lib/libb3d.a : lib/libb3d.a
	cp -f lib/libb3d.a ${MADAI_INSTALLDIR}/lib

${MADAI_INSTALLDIR}/include/b3d.h : include/b3d.h
	cp -f include/b3d.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/sampler.h : include/sampler.h
	cp -f include/sampler.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/resonances.h : include/resonances.h
	cp -f include/resonances.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/part.h : include/part.h
	cp -f include/part.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/inelastic.h : include/inelastic.h
	cp -f include/inelastic.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/balancearrays.h : include/balancearrays.h
	cp -f include/balancearrays.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/acceptance.h : include/acceptance.h
	cp -f include/acceptance.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/hyper.h : include/hyper.h
	cp -f include/hyper.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/mutinfo.h : include/mutinfo.h
	cp -f include/mutinfo.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/seinfo.h : include/seinfo.h
	cp -f include/seinfo.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/regen.h : include/regen.h
	cp -f include/regen.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/include/cell.h : include/cell.h
	cp -f include/cell.h ${MADAI_INSTALLDIR}/include/

${MADAI_INSTALLDIR}/progdata/resinfo/resonances_pdg_new.dat : resinfo/resonances_pdg_new.dat
	cp -f resinfo/resonances_pdg_new.dat ${MADAI_INSTALLDIR}/progdata/resinfo/

${MADAI_INSTALLDIR}/progdata/resinfo/resonances_pdg.dat : resinfo/resonances_pdg.dat
	cp -f resinfo/resonances_pdg.dat ${MADAI_INSTALLDIR}/progdata/resinfo/

${MADAI_INSTALLDIR}/progdata/resinfo/resonances_nodecays.dat : resinfo/resonances_nodecays.dat
	cp -f resinfo/resonances_nodecays.dat ${MADAI_INSTALLDIR}/progdata/resinfo/

${MADAI_INSTALLDIR}/progdata/resinfo/resonances_pdg_weak.dat : resinfo/resonances_pdg_weak.dat
	cp -f resinfo/resonances_pdg_weak.dat ${MADAI_INSTALLDIR}/progdata/resinfo/

${MADAI_INSTALLDIR}/progdata/resinfo/pionsonly.dat : resinfo/pionsonly.dat
	cp -f resinfo/pionsonly.dat ${MADAI_INSTALLDIR}/progdata/resinfo/

${MADAI_INSTALLDIR}/progdata/resinfo/decays_pdg.dat : resinfo/decays_pdg.dat
	cp -f resinfo/decays_pdg.dat ${MADAI_INSTALLDIR}/progdata/resinfo/

${MADAI_INSTALLDIR}/progdata/resinfo/decays_pdg_weak.dat : resinfo/decays_pdg_weak.dat
	cp -f resinfo/decays_pdg_weak.dat ${MADAI_INSTALLDIR}/progdata/resinfo/

${MADAI_INSTALLDIR}/progdata/resinfo/decays_none.dat : resinfo/decays_none.dat
	cp -f resinfo/decays_none.dat ${MADAI_INSTALLDIR}/progdata/resinfo/
