HOSTARCH    = _HOSTARCH_
VERSIONBIT  = _VERSIONBIT_
EXEEXT      =
OBJEXT      = .o
LIBEXT      = .a
CCPROG      = _CCPROG_
CFPROG      = _F77PROG_
CF90PROG    = _F77PROG_ 
MCFPROG     = _MCFPROG_
CF90CCPOPT  = _F90FLAGS_
# Compilation options for optimization (make expor)
CCFOPT      = _CCFOPT_
# Compilation options for debug (make | make debug)
CCFDEB      = _CCFDEB_

_ORDPATH_
_HWLPATH_

LKFOPT      =
MKPROG      = make 
MPCCPROG    = _MPCCPROG_
ARFLAGS     = _ARFLAGS_
ARPROG      = _ARPROG_
EXTRALIB    = _LDFLAGS_

_INT_
_PRC_
_FLT_
_VMPI_
_SMP_
_ORD_
_DYN_
CCTYPES     = _CCTYPES_
CCTYPESFLT  = _CCTYPESFLT_
CCPASTIX    = _CCFLAGS_


###################################################################
#                          DO NOT TOUCH                           #
###################################################################

FOPT      := $(CCFOPT)
FDEB      := $(CCFDEB)
CCHEAD    := $(CCPROG) $(CCTYPES) $(CCFOPT)
CCFOPT    := $(CCFOPT) $(CCTYPES) $(CCPASTIX)
CCFDEB    := $(CCFDEB) $(CCTYPES) $(CCPASTIX)


###################################################################
#                        MURGE COMPATIBILITY                      #
###################################################################

MAKE     = $(MKPROG)
CC       = $(MPCCPROG)
CFLAGS   = $(CCTYPESFLT) $(CCFOPT) 
FC       = $(MCFPROG) 
FFLAGS   = $(CCFOPT)
LDFLAGS  = $(EXTRALIB) $(BLASLIB)



