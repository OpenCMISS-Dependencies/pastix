PASTIX_BIN ?= ../../../install
PASTIX_LIB ?= ../../../install
PASTIX_INC ?= ../../../install
PASTIX_CONF = ${PASTIX_BIN}/pastix-conf
CFLAGS       = $(shell ${PASTIX_CONF} --ccopts) $(shell ${PASTIX_CONF} --opts)
FFLAGS       = $(shell ${PASTIX_CONF} --fcopts)
CXXFLAGS     = $(shell ${PASTIX_CONF} --cxxopts) $(shell ${PASTIX_CONF} --opts)
FC          := $(shell ${PASTIX_CONF} --fc) $(FFLAGS)
FL          := $(shell ${PASTIX_CONF} --fl) $(FFLAGS)
CC          := $(shell ${PASTIX_CONF} --cc) $(CFLAGS)
CXX         := $(shell ${PASTIX_CONF} --cxx) $(CXXFLAGS)
CL          := $(shell ${PASTIX_CONF} --cl) $(CFLAGS)
INC         := $(shell ${PASTIX_CONF} --incs)
OPTS        := $(shell ${PASTIX_CONF} --opts)

# To activate openmp in examples
ifeq (1, $(words $(findstring -DMURGE_THREADSAFE, ${OPTS})))
  OMP_FLAG = -openmp
endif

EXTRALIB    := -lmatrix_driver $(shell ${PASTIX_CONF} --libs)
EXTRALIB_MURGE = $(shell ${PASTIX_CONF} --libs_murge)

DST         =  ../bin
OBJ         = ../obj
###################################################################
#                              BLAS                               #
###################################################################

# Choose Blas library (Only 1)
# Do not forget to set BLAS_HOME if it is not in your environnement
# BLAS_HOME=/path/to/blas
#---- BLAS from config.in ----
BLAS = BLAS from pastix-conf
EXTRALIB         := ${EXTRALIB} $(shell ${PASTIX_CONF} --blas)
EXTRALIB_MURGE   := ${EXTRALIB_MURGE} $(shell ${PASTIX_CONF} --blas)
#----  Blas    ----
#BLAS       = REFBLAS
#EXTRALIB  := ${EXTRALIB} -lblas
#---- Gotoblas ----
#BLAS       = GOTOBLAS
#EXTRALIB  := ${EXTRALIB} -L${BLAS_HOME} -lgoto
#----  MKL     ----
#BLAS       = MKL
#EXTRALIB  := ${EXTRALIB} -L${BLAS_HOME}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core
#----  Acml    ----
#BLAS       = ACML
#EXTRALIB  := ${EXTRALIB} -L${BLAS_HOME}/lib -lacml

EXAMPLES_NOTYPE  =                      ${DST}/simple
EXAMPLES_NOTYPE := ${EXAMPLES_NOTYPE}   ${DST}/cppsimple
EXAMPLES_NOTYPE := ${EXAMPLES_NOTYPE}   ${DST}/simple_param
EXAMPLES_NOTYPE := ${EXAMPLES_NOTYPE}   ${DST}/murge
ifeq (0, $(words $(findstring -DFORCE_NOMPI, ${OPTS})))
  EXAMPLES_NOTYPE := ${EXAMPLES_NOTYPE} ${DST}/cppsimple_dist
  EXAMPLES_NOTYPE := ${EXAMPLES_NOTYPE} ${DST}/simple_dist
  EXAMPLES_NOTYPE := ${EXAMPLES_NOTYPE} ${DST}/simple_dist_redist
endif
EXAMPLES_NOTYPE := ${EXAMPLES_NOTYPE} ${DST}/plot_memory_usage

EXAMPLES_S = $(patsubst ${DST}/%, ${DST}/s%, ${EXAMPLES_NOTYPE})
EXAMPLES_D = $(patsubst ${DST}/%, ${DST}/d%, ${EXAMPLES_NOTYPE})
EXAMPLES_C = $(patsubst ${DST}/%, ${DST}/c%, ${EXAMPLES_NOTYPE})
EXAMPLES_Z = $(patsubst ${DST}/%, ${DST}/z%, ${EXAMPLES_NOTYPE})
.NOTPARALLEL : promptblas ${DST}/%
all :       promptblas                              \
            ${EXAMPLES_NOTYPE}                      \
            ${EXAMPLES_S}                           \
            ${EXAMPLES_D}                           \
            ${EXAMPLES_C}                           \
            ${EXAMPLES_Z}                           \
            ${DST}/do_not_redispatch_rhs            \
            ${DST}/step-by-step                     \
            ${DST}/reentrant                        \
            ${DST}/schur                            \
            ${DST}/schur2                           \
            ${DST}/refinement                       \
            ${DST}/isolate_zeros

ifeq (1, $(words $(findstring -DMULT_SMX, ${OPTS})))
all                     : ${DST}/mult-rhs
ifeq (1, $(words $(findstring -DDISTRIBUTED, ${OPTS})))
all                     : ${DST}/mult-rhs-dist
endif
endif
ifeq (0, $(words $(findstring -DFORCE_NOMPI, ${OPTS})))
all                     :       ${DST}/fsimple                          \
                                ${DST}/fstep-by-step                    \
                                ${DST}/step-by-step_dist                \
                                ${DST}/multi-comm                       \
                                ${DST}/multi-comm-step

all                     :       ${DST}/fmurge                           \
                                ${DST}/murge                            \
                                ${DST}/murge_2Dlap                      \
                                ${DST}/murge_sequence                   \
                                ${DST}/murge-product
endif

examples_notypes        :       promptblas                              \
                                ${EXAMPLES_NOTYPE}                      \
                                ${DST}/step-by-step                     \
                                ${DST}/step-by-step_dist                \
                                ${DST}/reentrant                        \
                                ${DST}/multi-comm                       \
                                ${DST}/multi-comm-step                  \
                                ${DST}/schur                            \
                                ${DST}/schur2                           \
                                ${DST}/refinement                       \
                                ${DST}/isolate_zeros
ifeq (0, $(words $(findstring -DFORCE_NOMPI, ${OPTS})))
examples_notypes        :       ${DST}/fsimple                          \
                                ${DST}/fstep-by-step                    \
                                ${DST}/fmurge                           \
                                ${DST}/murge
endif

examples_s              :       promptblas                              \
                                echo ${EXAMPLES_S}

examples_d              :       promptblas                              \
                                ${EXAMPLES_D}

examples_c              :       promptblas                              \
                                ${EXAMPLES_C}

examples_z              :       promptblas                              \
                                ${EXAMPLES_Z}

DEPEND =        $(filter-out $(wildcard ${DST}), ${DST}) \
                ${PASTIX_LIB}/libpastix.a                \
                ${PASTIX_LIB}/libmatrix_driver.a         \
                ${PASTIX_INC}/pastix.h                   \
                ${PASTIX_INC}/get_options.h              \
                ${PASTIX_INC}/read_matrix.h              \
                ${PASTIX_INC}/cscd_utils.h               \
                utils.h sutils.h cutils.h                \
                dutils.h zutils.h                        \
                mem_trace.h ${DST}/mem_trace.o

EXTRALIB       :=  ${EXTRALIB}       ${DST}/mem_trace.o
EXTRALIB_MURGE :=  ${EXTRALIB_MURGE} ${DST}/mem_trace.o

DEPEND_O = $(filter-out $(filter %.o, ${DEPEND}), ${DEPEND})
${DST}/%.o : %.c ${DEPEND_O}
	${CC} -c $< -o $@ ${INC}

${DST}/% : %.c ${DEPEND}
	${CL} $< -o $@ ${INC} ${EXTRALIB}

${DST}/plot_memory_usage : plot_memory_usage.c mem_trace.c ${DEPEND}
	${CL} $< -o $@ ${INC} ${EXTRALIB}

%.C : %.c ${DEPEND}
	${CL} $< -o $@ -C -E ${INC}

DEPEND_CPP =    $(filter-out $(wildcard ${DST}), ${DST}) \
                ${PASTIX_LIB}/libpastix.a                \
                ${PASTIX_LIB}/libmatrix_driver.a         \
                ${PASTIX_INC}/pastix.h                   \
                ${PASTIX_INC}/get_options.h              \
                ${PASTIX_INC}/read_matrix.h              \
                ${PASTIX_INC}/cscd_utils.h               \
                utils.hpp sutils.hpp cutils.hpp          \
                dutils.hpp zutils.hpp                    \
                mem_trace.h ${DST}/mem_trace.o

${DST}/% : %.cpp ${DEPEND_CPP}
	${CXX} $< -o $@ ${INC} ${EXTRALIB}
%.CPP : %.cpp ${DEPEND_CPP}
	${CXX} $< -o $@ -C -E ${INC}

DEPEND_MURGE =  $(filter-out $(wildcard ${DST}), ${DST}) \
                ${PASTIX_LIB}/libpastix.a                \
                ${PASTIX_LIB}/libpastix_murge.a          \
                ${PASTIX_INC}/pastix.h                   \
                ${PASTIX_INC}/murge.h                    \
                mem_trace.h mem_trace.c

${DST}/murge_2Dlap : ${PASTIX_INC}/MURGE_GetLocalElementList.c
${DST}/murge : murge.c ${DEPEND_MURGE}
	${CL} ${OMP_FLAG} $< -o $@ ${INC} ${EXTRALIB_MURGE}
${DST}/smurge : smurge.c ${DEPEND_MURGE}
	${CL} ${OMP_FLAG} $< -o $@ ${INC} ${EXTRALIB_MURGE}
${DST}/dmurge : dmurge.c ${DEPEND_MURGE}
	${CL} ${OMP_FLAG} $< -o $@ ${INC} ${EXTRALIB_MURGE}
${DST}/cmurge : cmurge.c ${DEPEND_MURGE}
	${CL} ${OMP_FLAG} $< -o $@ ${INC} ${EXTRALIB_MURGE} -DTYPE_COMPLEX
${DST}/zmurge : zmurge.c ${DEPEND_MURGE}
	${CL} ${OMP_FLAG} $< -o $@ ${INC} ${EXTRALIB_MURGE} -DTYPE_COMPLEX

${DST}/murge_2Dlap : murge_2Dlap.c ${DEPEND_MURGE}
	${CL} ${OMP_FLAG} $< -o $@ ${INC} ${EXTRALIB_MURGE}
${DST}/murge_sequence :  murge_sequence.c ${DEPEND_MURGE}
	${CL} ${OMP_FLAG} $< -o $@ ${INC} ${EXTRALIB_MURGE}
${DST}/murge-product :  murge-product.c ${DEPEND_MURGE}
	${CL} ${OMP_FLAG} $< -o $@ ${INC} ${EXTRALIB_MURGE}

FORTRAN_DEPEND = $(filter-out $(wildcard ${DST}), ${DST})       \
                 ${OBJ}/utils.o                                 \
                 ${PASTIX_INC}/pastix_fortran.h                 \
                 ${PASTIX_LIB}/libmatrix_driver.a               \
                 ${PASTIX_LIB}/libpastix.a

${DST}/% : %.F90 ${FORTRAN_DEPEND}
	${FL} ${INC} $< ${OBJ}/utils.o -o $@ ${EXTRALIB}

MURGE_FORTRAN_DEPEND =  $(filter-out $(wildcard ${DST}), ${DST})\
                        ${PASTIX_INC}/murge.inc                 \
                        ${PASTIX_LIB}/libpastix.a               \
                        ${PASTIX_LIB}/libpastix_murge.a

${DST}/fmurge : ${OBJ}/fmurge.F90 ${MURGE_FORTRAN_DEPEND}
	${FL} ${OMP_FLAG} ${INC} ${OBJ}/fmurge.F90 -o $@  ${EXTRALIB_MURGE}

${OBJ}/fmurge.F90 : fmurge.F90 ${DEPEND_MURGE}
	sh replaceCOEF.sh $< $@ `${PASTIX_CONF} --opts`

${OBJ}/utils.o : utils.F90 $(filter-out $(wildcard ${OBJ}), ${OBJ})\
                 ${PASTIX_INC}/pastix_fortran.h
	${FC} ${INC} $< -c -o $@

${DST} :
	mkdir -p $@

${OBJ} :
	mkdir -p $@

clean :
	-rm -f ../bin/*
	-rm -f ../obj/*
	-rm -f *.mod
	-rm -f ${EXAMPLES_S} ${EXAMPLES_D} ${EXAMPLES_C} ${EXAMPLES_Z}
	-rm -f $(patsubst ${DST}/%, %.c, ${EXAMPLES_S} ${EXAMPLES_D} \
               ${EXAMPLES_C} ${EXAMPLES_Z})


PtoS =          sed -e 's/PaStiX::pastix_float_t/float/g'               \
            -e 's/pastix_float_t/float/g'                               \
            -e 's/MPI_PASTIX_FLOAT/MPI_FLOAT/g'                         \
            -e 's/\"utils.h\"/\"sutils.h\"/g'                           \
            -e 's/\"utils.hpp\"/\"sutils.hpp\"/g'                       \
            -e 's/ read_matrix(/ s_read_matrix(/g'                      \
            -e 's/ dread_matrix(/ s_dread_matrix(/g'                    \
            -e 's/ pastix_checkMatrix(/ s_pastix_checkMatrix(/g'        \
            -e 's/ pastix_getLocalNodeNbr(/ s_pastix_getLocalNodeNbr(/g'\
            -e 's/ pastix_getLocalNodeLst(/ s_pastix_getLocalNodeLst(/g'\
            -e 's/ cscd_redispatch(/ s_cscd_redispatch(/g'              \
            -e 's/ pastix(/ s_pastix(/g'                                \
            -e 's/PRINT_RHS(/PRINT_RHS_REAL(/g'                         \
            -e 's/CONJ(/CONJ_REAL(/g'                                   \
            -e 's/ dpastix(/ s_dpastix(/g'                              \
            -e 's/\([ (]\)MURGE_\([^ ]*\)*(/\1SMURGE_\2(/g'             \
            -e 's/murge.h/smurge.h/g'                                   \
            -e 's/:read_matrix(/:s_read_matrix(/g'                      \
            -e 's/:dread_matrix(/:s_dread_matrix(/g'                    \
            -e 's/:pastix_checkMatrix(/:s_pastix_checkMatrix(/g'        \
            -e 's/:pastix_getLocalNodeNbr(/:s_pastix_getLocalNodeNbr(/g'\
            -e 's/:pastix_getLocalNodeLst(/:s_pastix_getLocalNodeLst(/g'\
            -e 's/:cscd_redispatch(/:s_cscd_redispatch(/g'              \
            -e 's/:pastix(/:s_pastix(/g'                                \
            -e 's/:dpastix(/:s_dpastix(/g' $< > $@

PtoD =          sed -e 's/PaStiX::pastix_float_t/double/g'              \
            -e 's/pastix_float_t/double/g'                              \
            -e 's/MPI_PASTIX_FLOAT/MPI_DOUBLE/g'                        \
            -e 's/\"utils.h\"/\"dutils.h\"/g'                           \
            -e 's/\"utils.hpp\"/\"dutils.hpp\"/g'                       \
            -e 's/ read_matrix(/ d_read_matrix(/g'                      \
            -e 's/ dread_matrix(/ d_dread_matrix(/g'                    \
            -e 's/ pastix_checkMatrix(/ d_pastix_checkMatrix(/g'        \
            -e 's/ pastix_getLocalNodeNbr(/ d_pastix_getLocalNodeNbr(/g'\
            -e 's/ pastix_getLocalNodeLst(/ d_pastix_getLocalNodeLst(/g'\
            -e 's/ cscd_redispatch(/ d_cscd_redispatch(/g'              \
            -e 's/ pastix(/ d_pastix(/g'                                \
            -e 's/PRINT_RHS(/PRINT_RHS_REAL(/g'                         \
            -e 's/CONJ(/CONJ_REAL(/g'                                   \
            -e 's/ dpastix(/ d_dpastix(/g'                              \
            -e 's/\([ (]\)MURGE_\([^ ]*\)*(/\1DMURGE_\2(/g'             \
            -e 's/murge.h/dmurge.h/g'                                   \
            -e 's/:read_matrix(/:d_read_matrix(/g'                      \
            -e 's/:dread_matrix(/:d_dread_matrix(/g'                    \
            -e 's/:pastix_checkMatrix(/:d_pastix_checkMatrix(/g'        \
            -e 's/:pastix_getLocalNodeNbr(/:d_pastix_getLocalNodeNbr(/g'\
            -e 's/:pastix_getLocalNodeLst(/:d_pastix_getLocalNodeLst(/g'\
            -e 's/:cscd_redispatch(/:d_cscd_redispatch(/g'              \
            -e 's/:pastix(/:d_pastix(/g'                                \
            -e 's/:dpastix(/:d_dpastix(/g' $< > $@

PtoC =  sed -e 's/PaStiX::pastix_float_t/std::complex<float>/g'         \
            -e 's/pastix_float_t/float complex/g'                       \
            -e 's/MPI_PASTIX_FLOAT/MPI_COMPLEX/g'                       \
            -e 's/\"utils.h\"/\"cutils.h\"/g'                           \
            -e 's/\"utils.hpp\"/\"cutils.hpp\"/g'                       \
            -e 's/ read_matrix(/ c_read_matrix(/g'                      \
            -e 's/ dread_matrix(/ c_dread_matrix(/g'                    \
            -e 's/ pastix_checkMatrix(/ c_pastix_checkMatrix(/g'        \
            -e 's/ pastix_getLocalNodeNbr(/ c_pastix_getLocalNodeNbr(/g'\
            -e 's/ pastix_getLocalNodeLst(/ c_pastix_getLocalNodeLst(/g'\
            -e 's/ cscd_redispatch(/ c_cscd_redispatch(/g'              \
            -e 's/ pastix(/ c_pastix(/g'                                \
            -e 's/PRINT_RHS(/PRINT_RHS_CPLX(/g'                         \
            -e 's/CONJ(/CONJ_CPLX(/g'                                   \
            -e 's/ dpastix(/ c_dpastix(/g'                              \
            -e 's/\([ (]\)MURGE_\([^ ]*\)*(/\1CMURGE_\2(/g'             \
            -e 's/murge.h/cmurge.h/g'                                   \
            -e 's/:read_matrix(/:c_read_matrix(/g'                      \
            -e 's/:dread_matrix(/:c_dread_matrix(/g'                    \
            -e 's/:pastix_checkMatrix(/:c_pastix_checkMatrix(/g'        \
            -e 's/:pastix_getLocalNodeNbr(/:c_pastix_getLocalNodeNbr(/g'\
            -e 's/:pastix_getLocalNodeLst(/:c_pastix_getLocalNodeLst(/g'\
            -e 's/:cscd_redispatch(/:c_cscd_redispatch(/g'              \
            -e 's/:pastix(/:c_pastix(/g'                                \
            -e 's/:dpastix(/:c_dpastix(/g' $< > $@

PtoZ =  sed -e 's/PaStiX::pastix_float_t/std::complex<double>/g'        \
            -e 's/pastix_float_t/double complex/g'                      \
            -e 's/MPI_PASTIX_FLOAT/MPI_DOUBLE_COMPLEX/g'                \
            -e 's/\"utils.h\"/\"zutils.h\"/g'                           \
            -e 's/\"utils.hpp\"/\"zutils.hpp\"/g'                       \
            -e 's/ read_matrix(/ z_read_matrix(/g'                      \
            -e 's/ dread_matrix(/ z_dread_matrix(/g'                    \
            -e 's/ pastix_checkMatrix(/ z_pastix_checkMatrix(/g'        \
            -e 's/ pastix_getLocalNodeNbr(/ z_pastix_getLocalNodeNbr(/g'\
            -e 's/ pastix_getLocalNodeLst(/ z_pastix_getLocalNodeLst(/g'\
            -e 's/ cscd_redispatch(/ z_cscd_redispatch(/g'              \
            -e 's/ pastix(/ z_pastix(/g'                                \
            -e 's/PRINT_RHS(/PRINT_RHS_CPLX(/g'                         \
            -e 's/CONJ(/CONJ_DCPLX(/g'                                  \
            -e 's/ dpastix(/z_dpastix(/g'                               \
            -e 's/murge.h/zmurge.h/g'                                   \
            -e 's/\([ (]\)MURGE_\([^ ]*\)*(/\1ZMURGE_\2(/g'             \
            -e 's/:read_matrix(/:z_read_matrix(/g'                      \
            -e 's/:dread_matrix(/:z_dread_matrix(/g'                    \
            -e 's/:pastix_checkMatrix(/:z_pastix_checkMatrix(/g'        \
            -e 's/:pastix_getLocalNodeNbr(/:z_pastix_getLocalNodeNbr(/g'\
            -e 's/:pastix_getLocalNodeLst(/:z_pastix_getLocalNodeLst(/g'\
            -e 's/:cscd_redispatch(/:z_cscd_redispatch(/g'              \
            -e 's/:pastix(/:z_pastix(/g'                                \
            -e 's/:dpastix(/:z_dpastix(/g' $< > $@

sutils.h : utils.h
	${PtoS}
sutils.hpp : utils.hpp
	${PtoS}
s%.c	: %.c
	${PtoS}
s%.cpp	: %.cpp
	${PtoS}

dutils.h : utils.h
	${PtoD}
dutils.hpp : utils.hpp
	${PtoD}
d%.c	: %.c
	${PtoD}
d%.cpp	: %.cpp
	${PtoD}

cutils.h : utils.h
	${PtoC}
cutils.hpp : utils.hpp
	${PtoC}
c%.c	: %.c
	${PtoC}
c%.cpp	: %.cpp
	${PtoC}

zutils.h : utils.h
	${PtoZ}
zutils.hpp : utils.hpp
	${PtoZ}
z%.c	: %.c
	${PtoZ}
z%.cpp	: %.cpp
	${PtoZ}

promptblas :
	@echo "************************************************************************************"
	@echo "* Your are using ${BLAS}, you can edit this in examples/src/makefile *"
	@echo "************************************************************************************"
	@sleep 2
