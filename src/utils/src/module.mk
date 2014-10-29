MODULE:=utils
MY_SRC :=
MY_SRC_MULT_ARCH := 

include all_modules.mk
include all_rules.mk

$(DSTDIR)/genheader$(EXE)	:	$(SRCDIR)/genheader.c $(DSTDIR)/. config.in
	$(CCPRG) $(CCINC) $(CCTYPES) $(CCTYPESFLT) $(LKOPT) $(<) -o $(@) $(LKLIB)

PO_FLAGS = $(CCINC) $(CCTYPES) $(LKOPT) -DVERSION=${VERSION}	\
	-DX_ARCH$(HOSTARCH) -I../../sopalin/src -I../../common/src		\
	-I../../symbol/src  -I../../solver.src -I../../blend/src 

$(DSTDIR) :
	mkdir -p $@

$(DSTDIR)/print_options$(EXE)	:	$(SRCDIR)/print_options.c	\
					sopalin/src/sopalin_option.c 	\
					sopalin/src/sopalin_option.h	\
					$(filter-out $(wildcard $(DSTDIR)), $(DSTDIR))
	$(MPCCPRG) $(PO_FLAGS) $(CCTYPESFLT) $(<) sopalin/src/sopalin_option.c -o $(@) $(LKLIB)

$(DSTDIR)/print_options_s$(EXE)	:	$(SRCDIR)/print_options.c	\
					sopalin/src/sopalin_option.c 	\
					sopalin/src/sopalin_option.h	\
					$(filter-out $(wildcard $(DSTDIR)), $(DSTDIR))
	$(MPCCPRG) $(PO_FLAGS) $(<) sopalin/src/sopalin_option.c -o $(@) $(LKLIB)

$(DSTDIR)/print_options_d$(EXE)	:	$(SRCDIR)/print_options.c	\
					sopalin/src/sopalin_option.c 	\
					sopalin/src/sopalin_option.h	\
					$(filter-out $(wildcard $(DSTDIR)), $(DSTDIR))
	$(MPCCPRG) $(PO_FLAGS) -DPREC_DOUBLE $(<) sopalin/src/sopalin_option.c -o $(@) $(LKLIB)

$(DSTDIR)/print_options_c$(EXE)	:	$(SRCDIR)/print_options.c	\
					sopalin/src/sopalin_option.c 	\
					sopalin/src/sopalin_option.h	\
					$(filter-out $(wildcard $(DSTDIR)), $(DSTDIR))
	$(MPCCPRG) $(PO_FLAGS) -DTYPE_COMPLEX $(<) sopalin/src/sopalin_option.c -o $(@) $(LKLIB)

$(DSTDIR)/print_options_z$(EXE)	:	$(SRCDIR)/print_options.c	\
					sopalin/src/sopalin_option.c 	\
					sopalin/src/sopalin_option.h	\
					$(filter-out $(wildcard $(DSTDIR)), $(DSTDIR))
	$(MPCCPRG) $(PO_FLAGS) -DTYPE_COMPLEX -DPREC_DOUBLE $(<) sopalin/src/sopalin_option.c -o $(@) $(LKLIB)
