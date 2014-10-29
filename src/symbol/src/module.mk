MODULE:=symbol

MY_SRC = dof.c				\
	dof_io.c			\
	symbol.c			\
	symbol_base.c			\
	symbol_check.c			\
	symbol_cost.c			\
	symbol_draw.c			\
	symbol_io.c			\
	symbol_keep.c			\
	symbol_levf.c			\
	symbol_nonzeros.c		\
	symbol_tree.c



MY_SRC_MULT_ARCH = 

include all_modules.mk
include all_rules.mk

$(OJTDIR)/dof$(OBJ) :	$(SRCDIR)/dof.c config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC);
	$(MPCCPRG) $(CCOPT) $(CCINC) -c $< -o $@  

