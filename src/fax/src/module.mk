MODULE:=fax

MY_SRC = symbol_compact.c		\
	symbol_costi.c			\
	symbol_fax.c			\
	symbol_fax_graph.c		\
	symbol_faxi.c			\
	symbol_faxi_graph.c

MY_SRC_MULT_ARCH = 


include all_modules.mk
include all_rules.mk

$(OJTDIR)/symbol_faxi_graph$(OBJ):	$(SRCDIR)/symbol_faxi_graph.c config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC);
	$(MPCCPRG) $(CCOPT) $(CCINC) -c $< -o $@  
$(OJTDIR)/symbol_fax_graph$(OBJ):	$(SRCDIR)/symbol_fax_graph.c config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC);
	$(MPCCPRG) $(CCOPT) $(CCINC) -c $< -o $@



