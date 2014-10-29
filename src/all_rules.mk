##
##  General inference rules.
##

## OBJ : SRC

${OJTDIR}/%_s.cu : ${SRCDIR}/%.cu
	@echo -e ${to_S}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@
${OJTDIR}/%_s.c  : ${SRCDIR}/%.c
	@echo -e ${to_S}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@

${OJTDIR}/%_d.cu : ${SRCDIR}/%.cu
	@echo -e ${to_D}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@
${OJTDIR}/%_d.c  : ${SRCDIR}/%.c
	@echo -e ${to_D}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@

${OJTDIR}/%_c.cu : ${SRCDIR}/%.cu
	@echo -e ${to_C}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@
${OJTDIR}/%_c.c  : ${SRCDIR}/%.c
	@echo -e ${to_C}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@

${OJTDIR}/%_z.c  : ${SRCDIR}/%.c
	@echo -e ${to_Z}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@
${OJTDIR}/%_z.cu : ${SRCDIR}/%.cu
	@echo -e ${to_Z}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@


${OJTDIR}/%${OBJ} :	$(SRCDIR)/%.c config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC) $(CCTYPESFLT);
	${CCPRG}      $(CCOPT) $(CCINC) $(CCTYPESFLT) -c $(<) -o $(@)
${OJTDIR}/%${OBJ} :	$(OJTDIR)/%.c config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC) $(CCTYPESFLT);
	$(CCPRG)      $(CCOPT) $(CCINC) $(CCTYPESFLT) -c $(<) -o $(@)

## EXE : SRC
$(DSTDIR)/%$(EXE)	:	$(SRCDIR)/%.c config.in
	$(CCPRG) $(CCOPT) $(CCINC) $(LKOPT) $(<) -o $(@) $(LKLIB)
