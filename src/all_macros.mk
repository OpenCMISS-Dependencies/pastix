##
##  Macros.
##

SHELL		= bash

VPATH		= $(OJTDIR):$(DSTDIR)

EXE		= $(EXEEXT)
LIB		= $(LIBEXT)
LIB_SO          = $(SOEXT)
OBJ		= $(OBJEXT)

AR		= $(ARPROG)
CCINC           = -I. -I$(DSTDIR) -I$(INSDIR)
CAT		= cat
CP		= cp
MKDIR		= mkdir
MV              = mv
RANLIB		= ranlib
RM		= rm -f
TOUCH		= touch
TAIL		= tail
HEAD		= head
VERSION         = `../../myversion.sh`
SO_VERSION      = 5.2

COMMONPASTIX_H	= 	$(INSDIR)/common_pastix.h	\
			$(INSDIR)/api.h			\
			$(INSDIR)/debug.h		\
			$(INSDIR)/errors.h		\
			$(INSDIR)/redefine_functions.h

to_S =  "\#define WITH_TYPE_PREFIX\n"	\
	"\#ifdef TYPE_COMPLEX\n"	\
	"\#  undef TYPE_COMPLEX\n"	\
	"\#endif\n"			\
	"\#ifdef PREC_DOUBLE\n"	\
	"\#  undef PREC_DOUBLE\n"	\
	"\#endif"			\

to_D =  "\#define WITH_TYPE_PREFIX\n"	\
	"\#ifdef TYPE_COMPLEX\n"	\
	"\#  undef TYPE_COMPLEX\n"	\
	"\#endif\n"			\
	"\#ifndef PREC_DOUBLE\n"	\
	"\#  define PREC_DOUBLE\n"	\
	"\#endif"			\

to_C =  "\#define WITH_TYPE_PREFIX\n"	\
	"\#ifndef TYPE_COMPLEX\n"	\
	"\#  define TYPE_COMPLEX\n"	\
	"\#endif\n"			\
	"\#ifdef PREC_DOUBLE\n"	\
	"\#  undef PREC_DOUBLE\n"	\
	"\#endif"			\

to_Z =  "\#define WITH_TYPE_PREFIX\n"	\
	"\#ifndef TYPE_COMPLEX\n"	\
	"\#  define TYPE_COMPLEX\n"	\
	"\#endif\n"			\
	"\#ifndef PREC_DOUBLE\n"	\
	"\#  define PREC_DOUBLE\n"	\
	"\#endif"			\

to_GE = "\#define SOPALIN_LU\n"
to_SY = "\#undef CHOL_SOPALIN\n"
to_HE = "\#undef CHOL_SOPALIN\n\#define HERMITIAN\n"
