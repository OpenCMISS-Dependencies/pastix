#
# File: NewMake
#
# A brand new makefile which should deal better with dependencies
#
# usage : make [debug] target
#   targets :
#      help           : print this help
#      all            : build PaStiX library
#      drivers        : build matrix drivers library
#      examples       : build examples, require drivers and all
#      murge_examples : build murge examples
#      pypastix       : Build Python interface to PaStiX, using Cython
#      murge_up       : synchronize murge directory (dev only)
#      clean          : remove all binaries and objects directories
#      cleanall       : remove all binaries, objects and dependencies directories
#
#
include config.in
include all_macros.mk


PWD           := $(shell pwd)
PREFIX        ?= ${PWD}/../install
BUILD_INCDIR   = ${PREFIX}
BUILD_BINDIR   = ${PREFIX}
BUILD_LIBDIR   = ${PREFIX}
BUILD_TSTDIR   = ${PREFIX}/test
INCLUDEDIR    ?= ${PREFIX}/include
LIBDIR        ?= ${PREFIX}/lib
BINDIR        ?= ${PREFIX}/bin
PYTHON_PREFIX ?= ${PREFIX}
VERSIONTYPE    = ${VERSIONINT}${VERSIONPRC}${VERSIONFLT}
VERSIONNAME    = ${VERSIONBIT}${VERSIONMPI}${VERSIONSMP}${VERSIONBUB}${VERSIONTYPE}${VERSIONORD}_${HOSTARCH}

MODULES := common order symbol fax perf blend kass sopalin utils matrix_drivers wrapper sparse-matrix

#
#We don't need to use dependencies when we're making these targets
#
NODEPS:=clean tags cleanall help info man

# DBL_FLG 0:simple 1:double
DBL_FLG = 0
ifeq (1, $(words $(findstring -DPREC_DOUBLE, $(CCTYPESFLT))))
DBL_FLG = 1
endif
# CPLX_FLG 0:real 1:complex
CPLX_FLG = 0
ifeq (1, $(words $(findstring -DTYPE_COMPLEX, $(CCTYPESFLT))))
CPLX_FLG = 1
endif

# INT_FLG 0:int 1:int64_t 2:int32_t 3:long
INT_FLG = 0
ifeq (1, $(words $(findstring -DFORCE_INT64, $(CCTYPES))))
INT_FLG = 1
endif
ifeq (1, $(words $(findstring -DINTSIZE64, $(CCTYPES))))
INT_FLG = 1
endif
ifeq (1, $(words $(findstring -DINTSSIZE64, $(CCTYPES))))
INT_FLG = 1
endif
ifeq (1, $(words $(findstring -DFORCE_INT32, $(CCTYPES))))
INT_FLG = 2
endif
ifeq (1, $(words $(findstring -DINTSIZE32, $(CCTYPES))))
INT_FLG = 2
endif
ifeq (1, $(words $(findstring -DFORCE_LONG, $(CCTYPES))))
INT_FLG = 3
endif
ifeq (1, $(words $(findstring -DLONG, $(CCTYPES))))
INT_FLG = 2
endif

# look for include files in
#   each of the modules
CCINC = $(patsubst %,-I./%/src,$(MODULES)) -Imurge/include

# SVN Version number
VERSION = $(shell ./myversion.sh)
VERSION := '"${VERSION}"'

COMMON_FLAGS := $(COMMON_FLAGS) -DVERSION=${VERSION} -DX_ARCH$(HOSTARCH) -DDOF_CONSTANT
CCINC:=$(COMMON_FLAGS) $(CCINC)

CCOPT=$(CCFOPT) $(COMMON_FLAGS)
FOPT=$(CCFOPT)

# Debug mode
CUDA_MEMCHECK=
ifeq (1, $(words $(findstring debug, $(MAKECMDGOALS))))
   DEBUG_FLAGS   = -DCOMMON_DEBUG -DSYMBOL_DEBUG        \
        -DDOF_DEBUG -DSOLVER_DEBUG -DORDER_DEBUG        \
        -DFAX_DEBUG -DSCOTCH_DEBUG -DDEBUG_BLEND        \
        -DBUBBLE_DEBUG -DPASTIX_DEBUG
   CCOPT=$(CCFDEB) $(DEBUG_FLAGS)
   FOPT=$(CCFDEB)
   CUDA_MEMCHECK=cuda-memcheck
endif



# extra libraries if required
#LIBS :=

# each module will add to this
OBJ_LIB     :=
OBJ_DRIVERS :=
OBJ_MURGE   :=

# All rules uses CCPRG/MPCCPRG
CCPRG := $(CCPROG)
MPCCPRG := $(MPCCPROG)

# Directory to store built files
DEPTREE := $(patsubst %,%/dep,$(MODULES))
TREE := $(patsubst %,%/obj/$(HOSTARCH),$(MODULES))      \
        $(patsubst %,%/bin/$(HOSTARCH),$(MODULES))      \
        $(BUILD_LIBDIR) ${BUILD_BINDIR} ${BUILD_INCDIR} \
        ${BUILD_TSTDIR}

# General rule
all:    murge						\
	$(patsubst %,%/.,$(TREE))                       \
        $(patsubst %,%/.,$(DEPTREE))                    \
        $(BUILD_LIBDIR)/libpastix$(LIB)                 \
        $(BUILD_LIBDIR)/libpastix_murge$(LIB)           \
        $(BUILD_INCDIR)/murge.h                         \
        $(BUILD_INCDIR)/smurge.h                        \
        $(BUILD_INCDIR)/dmurge.h                        \
        $(BUILD_INCDIR)/cmurge.h                        \
        $(BUILD_INCDIR)/zmurge.h                        \
        $(BUILD_INCDIR)/MURGE_GetLocalElementList.c     \
        $(BUILD_INCDIR)/pastix.h                        \
        $(BUILD_INCDIR)/murge.inc                       \
        $(BUILD_INCDIR)/pastix_fortran.h                \
        $(BUILD_INCDIR)/pastix_nompi.h                  \
        $(BUILD_BINDIR)/pastix-conf

.PHONY: murge_up

.svn/entries:
	mkdir .svn
	touch $@

MURGE_TAG=7addb68bbf78d3dbfe33c53b8e9ff1cef81c9a74
murge_up: murge
	-(cd murge; git branch | grep \* | grep feature > /dev/null ; if [ $$? -eq 1 ]; then git pull; git checkout ${MURGE_TAG}; fi; )

murge:
	git clone git://scm.gforge.inria.fr/murge/murge.git $@
	cd murge; git checkout ${MURGE_TAG}
	-(cd murge; git branch | grep \* | grep feature > /dev/null ; if [ $$? -eq 1 ]; then git pull; git checkout ${MURGE_TAG}; fi; )

ifeq (1, ${SHARED})
all:    $(BUILD_LIBDIR)/libpastix${LIB_SO} \
        $(BUILD_LIBDIR)/libpastix_murge${LIB_SO}
all:    $(BUILD_LIBDIR)/libpastix${LIB_SO}.${SO_VERSION} \
        $(BUILD_LIBDIR)/libpastix_murge${LIB_SO}.${SO_VERSION}
endif

tests: examples murge

%${LIB_SO}.${SO_VERSION}: %${LIB_SO}
	(cd $(shell dirname $@) && ln -sf $(shell basename $<) $(shell basename $@))

test:   $(patsubst %,%/.,$(TREE))			\
	$(patsubst %,%/.,$(DEPTREE))			\

#Provide some help
help info man:
	@echo "usage : make [debug] target"
	@echo "  targets :"
	@echo "     help           : print this help"
	@echo "     all            : build PaStiX library"
	@echo "     drivers        : build matrix drivers library"
	@echo "     examples       : build examples, require drivers and all"
	@echo "     murge_examples : build murge examples"
	@echo "     pypastix       : Build Python interface to PaStiX, using Cython"
	@echo "     murge_up       : synchronize murge directory (dev only)"
	@echo "     clean          : remove all binaries and objects directories"
	@echo "     cleanall       : remove all binaries, objects and dependencies directories"

ifeq (0, $(words $(filter-out debug, $(MAKECMDGOALS))))
debug:  all
else
debug:
	@echo "debug mode"
endif

#Cleaning rule
clean:
	-rm -rf $(TREE)
	-rm -rf example/bin

cleanall: clean
	-rm -rf $(TREE)
	-rm -rf example/bin
	-rm -rf $(DEPTREE)

%/.:
	mkdir -p $@

tree:
	-mkdir -p $(patsubst %,%/dep,$(MODULES))
	-mkdir -p $(patsubst %,%/bin/$(HOSTARCH),$(MODULES))
	-mkdir -p $(patsubst %,%/obj/$(HOSTARCH),$(MODULES))
	-mkdir -p ${BUILD_BINDIR}
	-mkdir -p ${BUILD_INCDIR}
	-mkdir -p ${BUILD_LIBDIR}
	-mkdir -p ${BUILD_TSTDIR}

python: $(BUILD_LIBDIR)/_pastix.so
	PYTHONPATH=$(MPI4PY_LIBDIR):$(BUILD_LIBDIR) python example/src/pastix_python.py

#Making examples
examples: all drivers install
	$(MKPROG) -C example/src

examples_s: all drivers install
	$(MKPROG) -C example/src examples_s

examples_d: all drivers install
	$(MKPROG) -C example/src examples_d

examples_c: all drivers install
	$(MKPROG) -C example/src examples_c

examples_z: all drivers install
	$(MKPROG) -C example/src  examples_z

examples_common: all drivers install
	$(MKPROG) -C example/src  examples_common


murge/Makefile.inc: config.in
	cp config.in murge/Makefile.inc
	echo "">> murge/Makefile.inc
	echo "LDSOLVER = -L$(LIBDIR) -lpastix_murge -lpastix" >> murge/Makefile.inc
	echo "CFLAGS  := $(CFLAGS) -I$(INCLUDEDIR)" >> murge/Makefile.inc

murge/pymurge/src/%.pxd: wrapper/src/pypastix/src/%.pxd
	cp $< $@

murge/pymurge/src/%.pxi: wrapper/src/pypastix/src/%.pxi
	cp $< $@

murge/pymurge/src/MURGE.pxd:
	echo include \"_pastix.pxd\" >  $@
	echo include \"_murge.pxd\"  >> $@
	echo include \"PASTIX.pxd\"  >> $@

pymurge:murge/Makefile.inc		\
	pypastix
	$(MKPROG) -C murge pymurge

murge_examples: install murge/Makefile.inc

	cp config.in murge/Makefile.inc
	$(MKPROG) -C murge/examples
ifeq (1, $(words ${CTAGSPROG}))
ifeq (0, $(shell ${CTAGSPROG} -x --c-kinds=fp include/murge.h > /dev/null 2>&1 && echo $$?))
sopalin/src/smurge_pastix.h sopalin/src/dmurge_pastix.h sopalin/src/cmurge_pastix.h sopalin/src/zmurge_pastix.h: sopalin/src/murge_pastix.h

murge-arithmetics : murge/include/murge.h murge/Makefile.inc
	$(MKPROG) -C murge arithmetics
murge/include/smurge.h: murge/include/murge.h murge/Makefile.inc
	$(MKPROG) -C murge arithmetics
murge/include/dmurge.h: murge/include/murge.h murge/Makefile.inc
	$(MKPROG) -C murge arithmetics
murge/include/cmurge.h: murge/include/murge.h murge/Makefile.inc
	$(MKPROG) -C murge arithmetics
murge/include/zmurge.h: murge/include/murge.h murge/Makefile.inc
	$(MKPROG) -C murge arithmetics
endif
endif

drivers: $(BUILD_LIBDIR)/libmatrix_driver$(LIB)
drivers: $(BUILD_INCDIR)/get_options.h
drivers: $(BUILD_INCDIR)/read_matrix.h
ifeq (1, ${SHARED})
drivers: $(BUILD_LIBDIR)/libmatrix_driver${LIB_SO}
drivers: $(BUILD_LIBDIR)/libmatrix_driver${LIB_SO}.${SO_VERSION}
endif

INSTALLDIRS=${INCLUDEDIR} ${BINDIR} ${LIBDIR}
${INSTALLDIRS} :
	mkdir -p $@

install: all drivers
install: $(foreach dir, ${INSTALLDIRS}, $(filter-out $(wildcard $(dir)), $(dir)))
install: ${BINDIR}/pastix-conf
install: ${LIBDIR}/libpastix${LIB}
ifeq (1, ${SHARED})
install: ${LIBDIR}/libpastix${LIB_SO}
install: ${LIBDIR}/libmatrix_driver${LIB_SO}
install: ${LIBDIR}/libpastix_murge${LIB_SO}
install: ${LIBDIR}/libpastix${LIB_SO}.${SO_VERSION}
install: ${LIBDIR}/libmatrix_driver${LIB_SO}.${SO_VERSION}
install: ${LIBDIR}/libpastix_murge${LIB_SO}.${SO_VERSION}
endif
install: ${LIBDIR}/libmatrix_driver${LIB}
install: ${LIBDIR}/libpastix_murge${LIB}

install: $(INCLUDEDIR)/pastix.h
install: $(INCLUDEDIR)/cscd_utils.h
install: $(INCLUDEDIR)/murge.h
install: $(INCLUDEDIR)/MURGE_GetLocalElementList.c
install: $(INCLUDEDIR)/murge.inc
install: $(INCLUDEDIR)/pastix_fortran.h
install: $(INCLUDEDIR)/pastix_nompi.h
install: $(INCLUDEDIR)/read_matrix.h
install: $(INCLUDEDIR)/get_options.h

${LIBDIR}/% :  $(filter-out $(wildcard ${LIBDIR}), ${LIBDIR})
${LIBDIR}/%${LIB}: ${BUILD_LIBDIR}/%${LIB} $(filter-out $(wildcard ${LIBDIR}), ${LIBDIR})
	cp $< $@

${INCLUDEDIR}/%: ${BUILD_INCDIR}/% $(filter-out $(wildcard ${INCLUDEDIR}), ${INCLUDEDIR})
	cp $< $@

$(BUILD_INCDIR)/get_options.h: matrix_drivers/src/get_options.h
	cp matrix_drivers/src/get_options.h  $(BUILD_INCDIR)/get_options.h

$(BUILD_INCDIR)/read_matrix.h: matrix_drivers/src/read_matrix.h
	cp matrix_drivers/src/read_matrix.h  $(BUILD_INCDIR)/read_matrix.h

TAGS tags: */src/*.c */src/*.h
	find . -regex ".*\.[ch]"  | grep -v "redefine_functions.h" | grep -v "oldsrc" | xargs etags


# include the description for
#   each module
include $(patsubst %,%/src/module.mk,$(MODULES))

$(BUILD_LIBDIR)/libpastix$(LIB): $(OBJ_LIB)
	$(AR) $(ARFLAGS) $(BUILD_LIBDIR)/libpastix${VERSIONNAME}${LIB} $(^)
	-$(RANLIB) $(BUILD_LIBDIR)/libpastix${VERSIONNAME}${LIB}
	(cd $(BUILD_LIBDIR) && ln -sf libpastix${VERSIONNAME}${LIB} libpastix$(LIB))

ifeq (1, ${SHARED})
DYLIB_OPT=
ifeq (i686_mac, ${HOSTARCH})
DYLIB_OPT=${LDFLAGS} -L__LIBDIR__ -lpastix
endif
%/libmatrix_driver${LIB_SO}: $(OBJ_DRIVERS)
	$(MPCCPROG) $(subst __SO_NAME__,$@, ${SHARED_FLAGS}) -o $(dir $@)libmatrix_driver${VERSIONNAME}${LIB_SO} $(^) $(subst __LIBDIR__,$(dir $@),${DYLIB_OPT})
	(cd $(dir $@) && ln -sf libmatrix_driver${VERSIONNAME}${LIB_SO} libmatrix_driver${LIB_SO})

%/libpastix${LIB_SO}: $(OBJ_LIB)
	$(MPCCPROG)  $(subst __SO_NAME__,$@, ${SHARED_FLAGS}) -o $(dir $@)libpastix${VERSIONNAME}${LIB_SO} $(^) $(filter-out -lpastix, $(subst __LIBDIR__,$(dir $@),${DYLIB_OPT}))
	(cd $(dir $@) && ln -sf libpastix${VERSIONNAME}${LIB_SO} libpastix${LIB_SO})

%/libpastix_murge${LIB_SO}: $(OBJ_MURGE) %/libpastix${LIB}
	$(MPCCPROG)  $(subst __SO_NAME__,$@, ${SHARED_FLAGS}) -o $(dir $@)libpastix_murge${VERSIONNAME}${LIB_SO} $(^) $(subst __LIBDIR__,$(dir $@),${DYLIB_OPT})
	(cd $(dir $@) && ln -sf libpastix_murge${VERSIONNAME}${LIB_SO} libpastix_murge${LIB_SO})

endif


$(BUILD_LIBDIR)/libmatrix_driver$(LIB): $(OBJ_DRIVERS)
	$(AR) $(ARFLAGS) $(BUILD_LIBDIR)/libmatrix_driver${VERSIONNAME}${LIB} $(^)
	-$(RANLIB) $(BUILD_LIBDIR)/libmatrix_driver${VERSIONNAME}${LIB}
	(cd $(BUILD_LIBDIR) && ln -sf libmatrix_driver${VERSIONNAME}${LIB} libmatrix_driver$(LIB))

$(BUILD_LIBDIR)/libpastix_murge$(LIB): $(OBJ_MURGE)
	$(AR) $(ARFLAGS) $(BUILD_LIBDIR)/libpastix_murge${VERSIONNAME}${LIB} $(^)
	-$(RANLIB) $(BUILD_LIBDIR)/libpastix_murge${VERSIONNAME}${LIB}
	(cd $(BUILD_LIBDIR) && ln -sf libpastix_murge${VERSIONNAME}${LIB} libpastix_murge$(LIB))


$(BUILD_BINDIR)/pastix-conf:	config.in 		\
			utils/src/pastix-conf.sh	\
			$(filter-out $(wildcard $(BUILD_BINDIR)), $(BUILD_BINDIR))
	cp utils/src/pastix-conf.sh $@${VERSIONNAME}.sh
	cat $@${VERSIONNAME}.sh						   | \
	   sed -e s%\_LIBS\_%"-L$(BUILD_LIBDIR) -lpastix ${EXTRALIB}"%	   | \
	   sed -e s%\_INC\_%"-I$(BUILD_INCDIR)"%			   | \
	   sed -e s%\_CC\_%"${MPCCPROG}"%				   | \
	   sed -e s%\_CCOPT\_%"${CCOPT}"%				   | \
	   sed -e s%\_CXX\_%"${MPCXXPROG}"%				   | \
	   sed -e s%\_CXXOPT\_%"${CXXOPT}"%				   | \
	   sed -e s%\_FC\_%"${MCFPROG} ${CF90CCPOPT}"%			   | \
	   sed -e s%\_FCOPT\_%"${FCOPT}"%				   | \
	   sed -e s%\_CL\_%"${MPCCPROG}"%				   | \
	   sed -e s%\_FL\_%"${MCFPROG} "%				   | \
	   sed -e s%\_BLAS\_%"${BLASLIB}"%				   | \
	   sed -e s%\_VERSION\_%"`./myversion.sh`"% | \
	   sed -e s%\_OPTS\_%"${CCTYPES} ${CCTYPESFLT} ${CCPASTIX}"%         \
	   > $@${VERSIONNAME}.sh.tmp
	mv $@${VERSIONNAME}.sh.tmp $@${VERSIONNAME}.sh
	chmod +x $@${VERSIONNAME}.sh
	(cd $(BUILD_BINDIR) && ln -fs $@${VERSIONNAME}.sh pastix-conf)

${BINDIR}/pastix-conf:	config.in 			\
			utils/src/pastix-conf.sh	\
			$(filter-out $(wildcard $(BINDIR)), $(BINDIR))
	cp utils/src/pastix-conf.sh $@
	echo "MM${MPCXXPROG}MM"
	cat $@								   | \
	   sed -e s%\_LIBS\_%"-L$(LIBDIR) -lpastix ${EXTRALIB}"%	   | \
	   sed -e s%\_INC\_%"-I$(INCLUDEDIR)"%				   | \
	   sed -e s%\_CC\_%"${MPCCPROG}"%				   | \
	   sed -e s%\_CXX\_%"${MPCXXPROG}"%				   | \
	   sed -e s%\_CCOPT\_%"${CCOPT}"%				   | \
	   sed -e s%\_FC\_%"${MCFPROG} ${CF90CCPOPT}"%			   | \
	   sed -e s%\_FCOPT\_%"${FCOPT}"%				   | \
	   sed -e s%\_CL\_%"${MPCCPROG}"%				   | \
	   sed -e s%\_FL\_%"${MCFPROG} "%				   | \
	   sed -e s%\_BLAS\_%"${BLASLIB}"%				   | \
	   sed -e s%\_VERSION\_%"`cd sopalin/src && ../../myversion.sh `"% | \
	   sed -e s%\_OPTS\_%"${CCTYPES} ${CCTYPESFLT} ${CCPASTIX}"% > $@.tmp
	mv $@.tmp $@
	chmod +x $@


$(BUILD_INCDIR)/pastix_nompi.h: common/src/nompi.h
	cp  $< $@

ifeq (1, ${CPLX_FLG})
ifeq (1, ${DBL_FLG})
PASTIX_FLOAT_T=double complex
else
PASTIX_FLOAT_T=float complex
endif
else
ifeq (1, ${DBL_FLG})
PASTIX_FLOAT_T=double
MURGE_REAL=double
else
PASTIX_FLOAT_T=float
MURGE_REAL=float
endif
endif
MURGE_COEF=${PASTIX_FLOAT_T}

ifeq (1, ${CPLX_FLG})
PRINT_DEFINES="\#ifndef TYPE_COMPLEX\n\#define TYPE_COMPLEX\n\#endif\n"
endif
ifeq (1, ${DBL_FLG})
PRINT_DEFINES:=${PRINT_DEFINES}"\#ifndef PREC_DOUBLE\n\#define PREC_DOUBLE\n\#endif\n"
endif

PASTIX_INT_T=int32_t
MURGE_INTS=int32_t
MURGE_INTL=int32_t

ifeq (1, ${INT_FLG})
PRINT_DEFINES:=${PRINT_DEFINES}"\#ifndef INTSSIZE64\n\#define INTSSIZE64\n\#endif\n"
PASTIX_INT_T=int64_t
MURGE_INTS=int64_t
MURGE_INTL=int64_t
endif
ifeq (2, ${INT_FLG})
PRINT_DEFINES:=${PRINT_DEFINES}"\#ifndef INTSIZE32\n\#define INTSIZE32\n\#endif\n"
PASTIX_INT_T=int32_t
MURGE_INTS=int32_t
MURGE_INTL=int32_t
endif
ifeq (3, ${INT_FLG})
PRINT_DEFINES:=${PRINT_DEFINES}"\#ifndef LONG\n\#define LONG\n\#endif\n"
PASTIX_INT_T=int32_t
MURGE_INTS=int32_t
MURGE_INTL=int64_t
endif


${BUILD_INCDIR}/MURGE_GetLocalElementList.c : sopalin/src/MURGE_GetLocalElementList.c
	cp $< $@

$(BUILD_INCDIR)/murge.h: murge/include/murge.h sopalin/src/murge_pastix.h
	echo -e ${PRINT_DEFINES} > $@
	cat  murge/include/murge.h >> $@
	cat sopalin/src/murge_pastix.h >> $@

$(BUILD_INCDIR)/smurge.h: murge/include/smurge.h sopalin/src/smurge_pastix.h
	echo -e ${PRINT_DEFINES} > $@
	cat  murge/include/smurge.h >> $@
	cat sopalin/src/smurge_pastix.h >> $@

$(BUILD_INCDIR)/dmurge.h: murge/include/dmurge.h sopalin/src/dmurge_pastix.h
	echo -e ${PRINT_DEFINES} > $@
	cat  murge/include/dmurge.h >> $@
	cat sopalin/src/dmurge_pastix.h >> $@

$(BUILD_INCDIR)/cmurge.h: murge/include/cmurge.h sopalin/src/cmurge_pastix.h
	echo -e ${PRINT_DEFINES} > $@
	cat  murge/include/cmurge.h >> $@
	cat sopalin/src/cmurge_pastix.h >> $@

$(BUILD_INCDIR)/zmurge.h: murge/include/zmurge.h sopalin/src/zmurge_pastix.h
	echo -e ${PRINT_DEFINES} > $@
	cat  murge/include/zmurge.h >> $@
	cat sopalin/src/zmurge_pastix.h >> $@

$(BUILD_INCDIR)/pastix${VERSIONTYPE}.h:	utils/bin/${HOSTARCH}/genheader$(EXE) \
					sopalin/src/murge_pastix.inc	\
					murge/include/murge.h		\
					murge/scripts/genfort.pl	\
					utils/src/insert-fortran-inc.sh \
					sopalin/src/pastix.h

	utils/bin/${HOSTARCH}/genheader$(EXE) 			\
		$(BUILD_INCDIR)/pastix${VERSIONTYPE}.h 		\
		$(BUILD_INCDIR)/pastix_fortran${VERSIONTYPE}.h	\
		$(BUILD_INCDIR)/murge${VERSIONTYPE}.inc		\
		sopalin/src/murge_pastix.inc			\
		murge/include/murge.h				\
		murge/scripts/genfort.pl			\
		utils/src/insert-fortran-inc.sh C
	cat sopalin/src/pastix.h >> $(BUILD_INCDIR)/pastix${VERSIONTYPE}.h
	cat common/src/api.h     >> $(BUILD_INCDIR)/pastix${VERSIONTYPE}.h

$(BUILD_INCDIR)/pastix.h: $(BUILD_INCDIR)/pastix${VERSIONTYPE}.h
	(cd $(BUILD_INCDIR) && ln -sf pastix${VERSIONTYPE}.h           pastix.h)

$(BUILD_INCDIR)/pastix_fortran${VERSIONTYPE}.h:				\
				utils/bin/${HOSTARCH}/genheader$(EXE) 	\
				murge/include/murge.h	      		\
				murge/scripts/genfort.pl      		\
				sopalin/src/pastix.h	      		\
				common/src/api.h	      		\
				common/src/pastix_fortran.h   		\
				common/src/pastix_fortran.inc 		\
				sopalin/src/murge_pastix.inc  		\
				utils/src/insert-fortran-inc.sh
	utils/bin/${HOSTARCH}/genheader$(EXE) 			\
		$(BUILD_INCDIR)/pastix${VERSIONTYPE}.h 		\
		$(BUILD_INCDIR)/pastix_fortran${VERSIONTYPE}.h	\
		$(BUILD_INCDIR)/murge${VERSIONTYPE}.inc		\
		sopalin/src/murge_pastix.inc			\
		murge/include/murge.h				\
		murge/scripts/genfort.pl			\
		utils/src/insert-fortran-inc.sh Fortran
	cat common/src/pastix_fortran.h >> \
	   $(BUILD_INCDIR)/pastix_fortran${VERSIONTYPE}.h
$(BUILD_INCDIR)/pastix_fortran.h: $(BUILD_INCDIR)/pastix_fortran${VERSIONTYPE}.h
	(cd $(BUILD_INCDIR) && ln -sf pastix_fortran${VERSIONTYPE}.h   pastix_fortran.h)

$(BUILD_INCDIR)/murge${VERSIONTYPE}.inc: 					\
				utils/bin/${HOSTARCH}/genheader$(EXE) 	\
				murge/include/murge.h			\
				murge/scripts/genfort.pl		\
				sopalin/src/pastix.h			\
				common/src/api.h			\
				common/src/pastix_fortran.h		\
				common/src/pastix_fortran.inc		\
				sopalin/src/murge_pastix.inc		\
				utils/src/insert-fortran-inc.sh
	utils/bin/${HOSTARCH}/genheader$(EXE) 			\
		$(BUILD_INCDIR)/pastix${VERSIONTYPE}.h 		\
		$(BUILD_INCDIR)/pastix_fortran${VERSIONTYPE}.h	\
		$(BUILD_INCDIR)/murge${VERSIONTYPE}.inc		\
		sopalin/src/murge_pastix.inc			\
		murge/include/murge.h				\
		murge/scripts/genfort.pl			\
		utils/src/insert-fortran-inc.sh Murge
	cat common/src/pastix_fortran.inc >> \
	    $@

$(BUILD_INCDIR)/murge.inc: $(BUILD_INCDIR)/murge${VERSIONTYPE}.inc
	(cd $(BUILD_INCDIR) && ln -sf murge${VERSIONTYPE}.inc murge.inc)


$(BUILD_INCDIR)/cscd_utils${VERSIONTYPE}.h: 				\
				utils/bin/${HOSTARCH}/genheader$(EXE) 	\
				sopalin/src/murge_pastix.inc		\
				murge/include/murge.h			\
				murge/scripts/genfort.pl		\
				utils/src/insert-fortran-inc.sh 	\
				sopalin/src/cscd_utils.h

	utils/bin/${HOSTARCH}/genheader$(EXE) 			\
		$(BUILD_INCDIR)/cscd_utils${VERSIONTYPE}.h 	\
		$(BUILD_INCDIR)/pastix_fortran${VERSIONTYPE}.h	\
		$(BUILD_INCDIR)/murge${VERSIONTYPE}.inc		\
		sopalin/src/murge_pastix.inc			\
		murge/include/murge.h				\
		murge/scripts/genfort.pl			\
		utils/src/insert-fortran-inc.sh C
	cat sopalin/src/cscd_utils.h >> $@

$(BUILD_INCDIR)/cscd_utils.h: $(BUILD_INCDIR)/cscd_utils${VERSIONTYPE}.h
	(cd $(BUILD_INCDIR) && ln -sf cscd_utils${VERSIONTYPE}.h   cscd_utils.h)

#Don't create dependencies when we're cleaning, for instance
# Find a way to do with "make clean all"
ifeq (0, $(words $(foreach goal, ${MAKECMDGOALS}, $(findstring ${goal}, ${NODEPS}))))
  # include the C include
  #   dependencies
  -include $(patsubst %.o,%.d,$(DEP_LIB))
  -include $(patsubst %.o,%.d,$(DEP_DRIVERS))
  -include $(patsubst %.o,%.d,$(DEP_MURGE))
endif

MAKEDEPEND=@./depend.sh "$(CPP)" $(shell echo $@ |sed -e 's@/obj/$(HOSTARCH)/@/dep/@' -e 's/\.o/\.d/g') $< obj/$(HOSTARCH) $(CCOPT) $(CCINC)
