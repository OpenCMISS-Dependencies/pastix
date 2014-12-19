# Makefile.inc must define integer type INTS and INTL by defining 
# INTSIZE64 (INTL 64 bits, INTS 32 bits), INTSSIZE64 (All 64 bits), 
# INTSIZE32 (All 32 bits) or none (all int), these defines have to be 
# stored into CFLAGS environnement variable
#
# Makefile.inc also defines FLOAT and REAL types by defining
# PREC_SIMPLE or PREC_DOUBLE for simple or double precision
# TYPE_REAL or TYPE_COMPLEX for real or complex (only for FLOAT), 
# these defines have to be stored into CFLAGS environnement variable
#
# Makefile.inc must define : 
#   MAKE     - Makefile application
#   CC       - C compiler
#   CFLAGS   - C compilation flags
#   FC       - Fortran compiler
#   FFLAGS   - Fortran compilation flag
#   LDFLAGS  - Linking flags
#   LDSOLVER - Solver linking path

include Makefile.inc

CFLAGS := $(CFLAGS) -I./include/

example	: Makefile.inc
	(cd examples && $(MAKE))

fortran	:
	./scripts/genfort.pl include/murge.h > include/murge.inc

bin	:
	mkdir bin

clean   :
	(cd examples && make clean)

ifeq (0, $(shell ${CTAGS} -x --c-kinds=fp include/murge.h > /dev/null 2>&1 && echo $$?))
arithmetics: include/zmurge.h include/dmurge.h \
	include/cmurge.h include/smurge.h

include/murge_unprotected.h: include/murge.h
	echo "#ifdef MURGE_H"                  >  $@
	echo "#  define MURGE_WAS_INCLUDED"    >> $@
	echo "#  undef  MURGE_H"               >> $@
	echo "#endif /* MURGE_H */"            >> $@
	echo ""                                >> $@
	cat $<                                 >> $@
	echo ""                                >> $@
	echo "#undef MURGE_H"                  >> $@
	echo "#ifdef MURGE_WAS_INCLUDED"       >> $@
	echo "#  define MURGE_H"               >> $@
	echo "#  undef MURGE_WAS_INCLUDED"     >> $@
	echo "#endif /* MURGE_WAS_INCLUDED */" >> $@

include/zmurge.h: include/murge_unprotected.h Makefile.inc
	echo "#ifndef ZMURGE_H" > $@
	echo "#define ZMURGE_H" >> $@
	echo "" >> $@
	${CTAGS} -x --c-kinds=fp include/murge.h |awk '{printf "#define %-30s\t Z%-30s\n", $$1, $$1}' >> $@
	echo "" >> $@
	echo "#ifndef PREC_DOUBLE" >> $@
	echo "#  define PREC_DOUBLE" >> $@
	echo "#endif /* PREC_DOUBLE */" >> $@
	echo "#ifndef TYPE_COMPLEX" >> $@
	echo "#  define TYPE_COMPLEX" >> $@
	echo "#endif /* TYPE_COMPLEX */" >> $@
	echo "" >> $@
	cat include/murge_unprotected.h >> $@
	echo "" >> $@
	echo "#ifndef OVERWRITE_MURGE_CALLS" >> $@
	${CTAGS} -x --c-kinds=fp include/murge.h |awk '{printf "#  undef %-30s\n", $$1}' >> $@
	echo "#endif /* OVERWRITE_MURGE_CALLS */" >> $@
	echo "#endif /* ZMURGE_H */" >> $@

include/dmurge.h: include/murge_unprotected.h Makefile.inc
	echo "#ifndef DMURGE_H" > $@
	echo "#define DMURGE_H" >> $@
	echo "" >> $@
	${CTAGS} -x --c-kinds=fp include/murge.h |awk '{printf "#define %-30s\t D%-30s\n", $$1, $$1}' >> $@
	echo "" >> $@
	echo "#ifndef PREC_DOUBLE" >> $@
	echo "#  define PREC_DOUBLE" >> $@
	echo "#endif /* PREC_DOUBLE */" >> $@
	echo "#ifdef TYPE_COMPLEX" >> $@
	echo "#  undef TYPE_COMPLEX" >> $@
	echo "#endif /* TYPE_COMPLEX */" >> $@
	echo "" >> $@
	cat include/murge_unprotected.h >> $@
	echo "" >> $@
	echo "#ifndef OVERWRITE_MURGE_CALLS" >> $@
	${CTAGS} -x --c-kinds=fp include/murge.h |awk '{printf "#  undef %-30s\n", $$1}' >> $@
	echo "#endif /* OVERWRITE_MURGE_CALLS */" >> $@
	echo "#endif /* DMURGE_H */" >> $@

include/cmurge.h: include/murge_unprotected.h Makefile.inc
	echo "#ifndef CMURGE_H" > $@
	echo "#define CMURGE_H" >> $@
	echo "" >> $@
	${CTAGS} -x --c-kinds=fp include/murge.h |awk '{printf "#define %-30s\t C%-30s\n", $$1, $$1}' >> $@
	echo "" >> $@
	echo "#ifdef PREC_DOUBLE" >> $@
	echo "#  undef PREC_DOUBLE" >> $@
	echo "#endif /* PREC_DOUBLE */" >> $@
	echo "#ifndef TYPE_COMPLEX" >> $@
	echo "#  define TYPE_COMPLEX" >> $@
	echo "#endif /* TYPE_COMPLEX */" >> $@
	echo "" >> $@
	cat include/murge_unprotected.h >> $@
	echo "" >> $@
	echo "#ifndef OVERWRITE_MURGE_CALLS" >> $@
	${CTAGS} -x --c-kinds=fp include/murge.h |awk '{printf "#  undef %-30s\n", $$1}' >> $@
	echo "#endif /* OVERWRITE_MURGE_CALLS */" >> $@
	echo "#endif /* CMURGE_H */" >> $@

include/smurge.h: include/murge_unprotected.h Makefile.inc
	echo "#ifndef SMURGE_H" > $@
	echo "#define SMURGE_H" >> $@
	echo "" >> $@
	${CTAGS} -x --c-kinds=fp include/murge.h |awk '{printf "#define %-30s\t S%-30s\n", $$1, $$1}' >> $@
	echo "" >> $@
	echo "#ifdef PREC_DOUBLE" >> $@
	echo "#  undef PREC_DOUBLE" >> $@
	echo "#endif /* PREC_DOUBLE */" >> $@
	echo "#ifdef TYPE_COMPLEX" >> $@
	echo "#  undef TYPE_COMPLEX" >> $@
	echo "#endif /* TYPE_COMPLEX */" >> $@
	echo "" >> $@
	cat include/murge_unprotected.h >> $@
	echo "" >> $@
	echo "#ifndef OVERWRITE_MURGE_CALLS" >> $@
	${CTAGS} -x --c-kinds=fp include/murge.h |awk '{printf "#  undef %-30s\n", $$1}' >> $@
	echo "#endif /* OVERWRITE_MURGE_CALLS */" >> $@
	echo "#endif /* SMURGE_H */" >> $@
else
arithmetics:
	echo "To build [sdcz]murge.h you need a version of ctags supporting '-x --c-kinds=fp'"
endif
