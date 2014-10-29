OJTDIR:=$(MODULE)/obj/$(HOSTARCH)
DSTDIR:=$(MODULE)/bin/$(HOSTARCH)
SRCDIR:=$(MODULE)/src
DEPDIR:=$(MODULE)/dep

MY_SRC 		 := $(patsubst %,$(SRCDIR)/%,$(MY_SRC))
MY_SRC_MULT_ARCH := $(patsubst %,$(SRCDIR)/%,$(MY_SRC_MULT_ARCH))

MY_OBJ :=  $(patsubst %.c,%.o,$(filter %.c,$(MY_SRC))) 		\
	$(patsubst %.c,%.o,$(filter %.c,$(MY_SRC_MULT_ARCH))) 	\
	$(patsubst %.c,%_s.o,$(filter %.c,$(MY_SRC_MULT_ARCH))) \
	$(patsubst %.c,%_d.o,$(filter %.c,$(MY_SRC_MULT_ARCH))) \
	$(patsubst %.c,%_c.o,$(filter %.c,$(MY_SRC_MULT_ARCH))) \
	$(patsubst %.c,%_z.o,$(filter %.c,$(MY_SRC_MULT_ARCH)))

MY_OBJ := $(patsubst $(SRCDIR)/%,$(OJTDIR)/%,$(MY_OBJ))

OBJ_LIB := $(OBJ_LIB) $(MY_OBJ)

DEP_LIB := $(DEP_LIB) $(patsubst $(OJTDIR)/%.o,$(DEPDIR)/%.d,$(MY_OBJ))

$(DEPDIR): 
	mkdir -p $@