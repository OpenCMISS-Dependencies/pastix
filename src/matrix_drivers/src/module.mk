MODULE:=matrix_drivers
MY_SRC := skitf.f		\
	  iohb.c		\
	  mmio.c		\
          common_drivers.c	\
	  get_options.c

MY_SRC_MULT_ARCH := read_matrix.c	\
		    rsaread.c		\
		    hbread.c		\
		    mmread.c		\
		    mmdread.c		\
		    petscread.c		\
		    cccread.c		\
		    olafread.c		\
		    chbread.c		\
		    cscdread.c		\
		    peerread.c		\
		    threefilesread.c	\
		    fdupread.c		\
		    laplacian.c

OJTDIR:=$(MODULE)/obj/$(HOSTARCH)
DSTDIR:=$(MODULE)/bin/$(HOSTARCH)
SRCDIR:=$(MODULE)/src
DEPDIR:=$(MODULE)/dep

MY_SRC 		 := $(patsubst %,$(SRCDIR)/%,$(MY_SRC))
MY_SRC_MULT_ARCH := $(patsubst %,$(SRCDIR)/%,$(MY_SRC_MULT_ARCH))

MY_OBJ :=  $(patsubst %.c,%.o,$(filter %.c,$(MY_SRC))) 		\
	$(patsubst %.f,%.o,$(filter %.f,$(MY_SRC))) 		\
	$(patsubst %.c,%.o,$(filter %.c,$(MY_SRC_MULT_ARCH))) 	\
	$(patsubst %.c,%_s.o,$(filter %.c,$(MY_SRC_MULT_ARCH))) \
	$(patsubst %.c,%_d.o,$(filter %.c,$(MY_SRC_MULT_ARCH))) \
	$(patsubst %.c,%_c.o,$(filter %.c,$(MY_SRC_MULT_ARCH))) \
	$(patsubst %.c,%_z.o,$(filter %.c,$(MY_SRC_MULT_ARCH)))

MY_OBJ := $(patsubst $(SRCDIR)/%,$(OJTDIR)/%,$(MY_OBJ))
MY_OBJ := ${MY_OBJ} ${OJTDIR}/api_str_to_int${OBJ}

OBJ_DRIVERS := $(OBJ_DRIVERS) $(MY_OBJ)

DEP_DRIVERS := $(DEP_DRIVERS) $(patsubst $(OJTDIR)/%,$(DEPDIR)/%,$(MY_OBJ))

$(DEPDIR):
	mkdir -p $(DEPDIR)

##
##  General inference rules.
##
$(OBJ_DRIVERS) : $(BUILD_INCDIR)/pastix.h $(BUILD_INCDIR)/cscd_utils.h

$(OJTDIR)/%_s$(OBJ)	:	$(SRCDIR)/%.c config.in
	$(MAKEDEPEND)$(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR) -DFORGET_TYPE;
	$(MPCCPROG) $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR) -c $(<) -o $(@) -DFORGET_TYPE

$(OJTDIR)/%_d$(OBJ)	:	$(SRCDIR)/%.c config.in
	$(MAKEDEPEND) $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR) -DPREC_DOUBLE -DPREC_DOUBLE -DFORGET_TYPE;
	$(MPCCPROG) $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR) -c $(<) -o $(@) -DPREC_DOUBLE -DPREC_DOUBLE -DFORGET_TYPE

$(OJTDIR)/%_c$(OBJ)	:	$(SRCDIR)/%.c config.in
	$(MAKEDEPEND) $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR)  -DTYPE_COMPLEX -DTYPE_COMPLEX -DFORGET_TYPE;
	$(MPCCPROG) $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR) -c $(<) -o $(@) -DTYPE_COMPLEX -DTYPE_COMPLEX -DFORGET_TYPE

$(OJTDIR)/%_z$(OBJ)	:	$(SRCDIR)/%.c config.in
	$(MAKEDEPEND) $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR) -DPREC_DOUBLE -DPREC_DOUBLE -DTYPE_COMPLEX -DTYPE_COMPLEX -DFORGET_TYPE;
	$(MPCCPROG) $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR) -c $(<) -o $(@) -DPREC_DOUBLE -DPREC_DOUBLE -DTYPE_COMPLEX -DTYPE_COMPLEX -DFORGET_TYPE

$(OJTDIR)/%$(OBJ)	:	$(SRCDIR)/%.c config.in
	$(MAKEDEPEND) $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR)  $(CCTYPESFLT);
	$(MPCCPROG) $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR) $(CCTYPESFLT) -c $(<) -o $(@)

$(OJTDIR)/%$(OBJ)	:       $(SRCDIR)/%.f config.in
	echo "$@ $(patsubst $(DEPDIR)/%,$(OJTDIR)/%,$(patsubst %.d,%$(OBJ),$@)) : $<" > $@
	$(CFPROG) -c $< -o $@

$(OJTDIR)/api_str_to_int.c : common/src/api.h
	echo "#include \"api_str_to_int.h\"" > $@
	echo "#include \"string.h\"" >> $@
	echo "int api_str_to_int( char * string, int * value) {" >> $@
	grep "[ID]PARM_[A-Z_]*[ ]*=[ ]*[0-9]*" common/src/api.h | sed -e 's/\([ID]PARM_[A-Z_]*\)[ ]*=[ ]*\([0-9]*\).*/  if(!strcmp("\1", string)) { *value = \2\; return 0\;}/' >> $@
	grep "API_[A-Z_]*[ ]*=[ ]*[0-9]*" common/src/api.h | sed -e 's/\(API_[A-Z_]*\)[ ]*=[ ]*\([0-9]*\).*/  if(!strcmp("\1", string)) { *value = \2\; return 0\;}/' >> $@
	echo "  return 1;" >> $@
	echo "}" >> $@


${OJTDIR}/%${OBJ} : ${OJTDIR}/%.c config.in
	$(MAKEDEPEND) $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR) -Imatrix_drivers/src
	$(MPCCPROG)   $(CCOPT) $(COMMON_FLAGS) -I$(BUILD_INCDIR) -Imatrix_drivers/src -c $< -o $@
