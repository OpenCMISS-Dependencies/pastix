MODULE=wrapper

MY_SRC:=
MY_SRC_MULT_ARCH :=
MY_SRC_FACT_TYPE :=

MY_SRC_MURGE 		:= 
MY_SRC_MURGE_MULT_ARCH  := 
include all_modules.mk
include all_rules.mk

$(OJTDIR) : 
	mkdir -p $@

$(OJTDIR)/pastix.i : $(SRCDIR)/pastix.i $(BUILD_INCDIR)/pastix.h $(filter-out $(wildcard $(OJTDIR)), $(OJTDIR))

$(OJTDIR)/%.i : $(SRCDIR)/%.i
	cat $< > $@.tmp2
	cat $(BUILD_INCDIR)/pastix.h >> $@.tmp2

ifeq (0, $(DBL_FLG))
ifeq (0, $(CPLX_FLG))
	sed -e 's/pastix_float_t/float         /g' $@.tmp2 > $@.tmp
else
	sed -e 's/pastix_float_t/float complex /g' $@.tmp2 > $@.tmp
endif
else
ifeq (0, $(CPLX_FLG))
	sed -e 's/pastix_float_t/double         /g' $@.tmp2 > $@.tmp
else
	sed -e 's/pastix_float_t/double complex /g' $@.tmp2 > $@.tmp
endif
endif

ifeq (1, $(INT_FLG))
	sed -e 's/pastix_int_t/int64_t     /g' $@.tmp > $@
else
ifeq (2, $(INT_FLG))
	sed -e 's/pastix_int_t/int32_t     /g' $@.tmp > $@
else
ifeq (3, $(INT_FLG))
	sed -e 's/pastix_int_t/long        /g' $@.tmp > $@
else
	sed -e 's/pastix_int_t/int         /g' $@.tmp > $@
endif
endif
endif
	rm $@.tmp
	rm $@.tmp2

$(OJTDIR)/%.c :$(filter-out $(wildcard $(OJTDIR)), $(OJTDIR))
$(OJTDIR)/%.o :$(filter-out $(wildcard $(OJTDIR)), $(OJTDIR))

$(OJTDIR)/pastix_wrap_python.c : $(OJTDIR)/pastix.i $(BUILD_INCDIR)/pastix.h 
	swig -python -I$(MPI4PY_INC)/mpi4py -I$(BUILD_INCDIR) -o $@ $<
	mv $(patsubst %.i,%.py, $<) $(LIBDIR)

$(OJTDIR)/pastix_wrap_octave.c : $(OJTDIR)/pastix.i $(BUILD_INCDIR)/pastix.h 
	swig -octave -I$(MPI4PY_INC)/mpi4py -I$(BUILD_INCDIR)  -o $@ $<

$(OJTDIR)/pastix_wrap_python.o : $(OJTDIR)/pastix_wrap_python.c 
	$(MPCCPRG) -I$(MPI4PY_INC) $(shell ${BUILD_BINDIR}/pastix-conf --incs) -c -fpic $< -I${PYTHON_INC} -o $@


$(LIBDIR)/_pastix.so : $(BUILD_LIBDIR)/libpastix$(LIB) ${BUILD_BINDIR}/pastix-conf 
$(LIBDIR)/_pastix.so : $(filter-out $(wildcard $(LIBDIR)), $(LIBDIR))
$(LIBDIR)/_pastix.so : $(OJTDIR)/pastix_wrap_python.o
	$(MPCCPRG) -shared $< $(shell ${BUILD_BINDIR}/pastix-conf --libs) $(shell ${BUILD_BINDIR}/pastix-conf --blas)  -o $@

ifeq (1, $(words $(findstring pxd, ${MAKECMDGOALS})))
pxd: wrapper/src/pypastix/src/_pastix.pxd 	\
	wrapper/src/pypastix/src/_murge.pxd	\
	wrapper/src/pypastix/src/_smurge.pxd	\
	wrapper/src/pypastix/src/_dmurge.pxd	\
	wrapper/src/pypastix/src/_cmurge.pxd	\
	wrapper/src/pypastix/src/_zmurge.pxd	\
	wrapper/src/pypastix/src/PASTIX.pxd	\
	wrapper/src/pypastix/src/_cscd_utils.pxd\
	wrapper/src/pypastix/src/pastix_types.pxi

wrapper/src/pypastix/src/pastix_types.pxi:
	echo 'cdef extern from "pastix.h":' > $@
	echo '    ctypedef ${PASTIX_INT_T} pastix_int_t' >> $@
	echo 'ctypedef np.${PASTIX_INT_T} np_pastix_int_t' >> $@
	echo 'ctypedef np.${PASTIX_FLOAT_T}_t  np_pastix_float_t' >> $@
	echo 'ctypedef np.complex64_t   np_pastix_float_t' >> $@
	echo 'ctypedef np.complex128_t  np_pastix_float_t' >> $@
	echo "" >> $@
	echo 'cdef extern from "murge.h":' >> $@
	echo '    ctypedef ${MURGE_COEF} COEF' >> $@
	echo '    ctypedef ${MURGE_REAL} REAL' >> $@
	echo '    ctypedef int INTS' >> $@
	echo '    ctypedef int INTL' >> $@
	echo 'ctypedef np.${MURGE_INTS} np_INTS' >> $@
	echo 'ctypedef np.${MURGE_INTL} np_INTL' >> $@
	echo 'ctypedef np.${MURGE_COEF}_t  np_COEF' >> $@
	echo 'ctypedef np.${MURGE_REAL}_t  np_REAL' >> $@

wrapper/src/pypastix/src/_murge.pxd :$(INCLUDEDIR)/murge.h
	echo  "#include <complex.h>" > /tmp/$(notdir $<)
	cat $< >> /tmp/$(notdir $<)
	python wrapper/src/pypastix/genpxd.py wrapper/src/pypastix/src /tmp/$(notdir $<)
	sed -e 's/int \([^ ]*comm\)/MPI_Comm \1/g' $@ > $@.tmp
	mv $@.tmp $@
	echo "    int MURGE_SetCommunicator(int id, MPI_Comm comm)" >> $@

wrapper/src/pypastix/src/_smurge.pxd :$(INCLUDEDIR)/smurge.h
	echo  "#include <complex.h>" > /tmp/$(notdir $<)
	cat $< >> /tmp/$(notdir $<)
	python wrapper/src/pypastix/genpxd.py wrapper/src/pypastix/src /tmp/$(notdir $<)
	sed -e 's/int \([^ ]*comm\)/MPI_Comm \1/g' $@ > $@.tmp
	mv $@.tmp $@
	echo "    int SMURGE_SetCommunicator(int id, MPI_Comm comm)" >> $@

wrapper/src/pypastix/src/_dmurge.pxd :$(INCLUDEDIR)/dmurge.h
	echo  "#include <complex.h>" > /tmp/$(notdir $<)
	cat $< >> /tmp/$(notdir $<)
	python wrapper/src/pypastix/genpxd.py wrapper/src/pypastix/src /tmp/$(notdir $<)
	sed -e 's/int \([^ ]*comm\)/MPI_Comm \1/g' $@ > $@.tmp
	mv $@.tmp $@
	echo "    int DMURGE_SetCommunicator(int id, MPI_Comm comm)" >> $@

wrapper/src/pypastix/src/_cmurge.pxd :$(INCLUDEDIR)/cmurge.h
	echo  "#include <complex.h>" > /tmp/$(notdir $<)
	cat $< >> /tmp/$(notdir $<)
	python wrapper/src/pypastix/genpxd.py wrapper/src/pypastix/src /tmp/$(notdir $<)
	sed -e 's/int \([^ ]*comm\)/MPI_Comm \1/g' $@ > $@.tmp
	mv $@.tmp $@
	echo "    int CMURGE_SetCommunicator(int id, MPI_Comm comm)" >> $@

wrapper/src/pypastix/src/_zmurge.pxd :$(INCLUDEDIR)/zmurge.h
	echo  "#include <complex.h>" > /tmp/$(notdir $<)
	cat $< >> /tmp/$(notdir $<)
	python wrapper/src/pypastix/genpxd.py wrapper/src/pypastix/src /tmp/$(notdir $<)
	sed -e 's/int \([^ ]*comm\)/MPI_Comm \1/g' $@ > $@.tmp
	mv $@.tmp $@
	echo "    int ZMURGE_SetCommunicator(int id, MPI_Comm comm)" >> $@

wrapper/src/pypastix/src/_%.pxd :$(INCLUDEDIR)/%.h
	echo  "#include <complex.h>" > /tmp/$(notdir $<)
	cat $< >> /tmp/$(notdir $<)
	python wrapper/src/pypastix/genpxd.py wrapper/src/pypastix/src /tmp/$(notdir $<)
	sed -e 's/int \([^ ]*comm\)/MPI_Comm \1/g' $@ > $@.tmp
	mv $@.tmp $@


wrapper/src/pypastix/src/PASTIX.pxd: wrapper/src/pypastix/src/template.py
wrapper/src/pypastix/src/PASTIX.pxd: wrapper/src/pypastix/parseapi.py
wrapper/src/pypastix/src/PASTIX.pxd: common/src/api.h
	python wrapper/src/pypastix/parseapi.py $< wrapper/src/pypastix/src/template.py > $@

endif

MYLDSOLVER = -L$(LIBDIR) -lpastix_murge -lpastix
MYLIBRARIES=$(patsubst -l%,%,$(filter -l%,  ${MYLDSOLVER} ${LDFLAGS}))
MYLIBRARYDIR=$(patsubst -L%,%,$(filter -L%, ${MYLDSOLVER} ${LDFLAGS}))
MYLIBRARYDIR:=$(shell echo ${MYLIBRARYDIR} | tr ' ' ':')
MYINCLUDEDIR=$(patsubst -I%,%,$(filter -I%,${CFLAGS}))
MYINCLUDEDIR:=$(shell echo ${INCLUDEDIR} ${MYINCLUDEDIR} | tr ' ' ':')

wrapper/src/pypastix/setup.cfg: config.in
	echo [build_ext]                  >  $@
	echo libraries = ${MYLIBRARIES}     >> $@
	echo library_dirs = ${MYLIBRARYDIR} >> $@
	echo include_dirs = ${MYINCLUDEDIR} >> $@

pypastix: wrapper/src/pypastix/build

.PHONY:wrapper/src/pypastix/build
wrapper/src/pypastix/build : wrapper/src/pypastix/src/pypastix.pyx
wrapper/src/pypastix/build : wrapper/src/pypastix/src/_pastix.pxd
wrapper/src/pypastix/build : wrapper/src/pypastix/src/_murge.pxd
wrapper/src/pypastix/build : wrapper/src/pypastix/src/_cscd_utils.pxd
wrapper/src/pypastix/build : wrapper/src/pypastix/src/PASTIX.pxd
#wrapper/src/pypastix/build : wrapper/src/pypastix/src/MURGE.pxi
#wrapper/src/pypastix/build : wrapper/src/pypastix/src/MURGE.pyx
#wrapper/src/pypastix/build : wrapper/src/pypastix/src/Matrix.pyx
wrapper/src/pypastix/build : wrapper/src/pypastix/src/pastix_types.pxi
wrapper/src/pypastix/build : wrapper/src/pypastix/setup.cfg
	@echo "Copying murge pyx files into PaStiX"
	@for file in MURGE.pxi MURGE.pyx Matrix.pyx;\
	do \
	  if [ -e murge/pymurge/src/$$file ]; \
          then \
	    cp murge/pymurge/src/$$file wrapper/src/pypastix/src/;\
	  else \
	    touch wrapper/src/pypastix/src/$$file;\
	  fi;\
	done;
	CC="${CC}" make -C wrapper/src/pypastix PASTIX_PREFIX=${PYTHON_PREFIX}
	touch $@
	#make -C wrapper/src/pypastix test
