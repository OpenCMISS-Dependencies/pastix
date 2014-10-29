MODULE:=kass

MY_SRC = kass.c				\
	compact_graph.c			\
	amalgamate.c			\
	ifax.c				\
	sparRow.c			\
	SF_Direct.c			\
	SF_level.c			\
	find_supernodes.c		\
	KSupernodes.c


MY_SRC_MULT_ARCH = sort_row.c

include all_modules.mk
include all_rules.mk

$(OJTDIR)/kass$(OBJ)		: $(SRCDIR)/kass.c config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC);
	$(MPCCPRG) $(CCOPT) $(CCINC) -c $(<) -o $(@)

$(OJTDIR)/amalgamate$(OBJ)	: $(SRCDIR)/amalgamate.c config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC);
	$(MPCCPRG) $(CCOPT) $(CCINC) -c $(<) -o $(@)
