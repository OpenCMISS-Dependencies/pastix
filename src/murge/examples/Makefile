include ../Makefile.inc
CFLAGS := $(CFLAGS) -I../include 

all : examples
include ../Makefile.rules

examples : 	../Makefile.inc 	\
	$(BINDIR)/murge-dist$(EXE) 	\
	$(BINDIR)/murge-global$(EXE) 	\
	$(BINDIR)/murge-mm$(EXE)	\
	$(BINDIR)/murge-save-load$(EXE)	\
	$(BINDIR)/multi-matrix-and-rhs$(EXE)

$(BINDIR)/murge-mm$(EXE) : $(OBJDIR)/murge-mm$(OBJ) $(OBJDIR)/mmio$(OBJ) $(BINDIR)
	$(CC) $(CFLAGS) $(OBJDIR)/murge-mm$(OBJ) $(OBJDIR)/mmio$(OBJ) -o $@ $(LDSOLVER) $(LDFLAGS) 

