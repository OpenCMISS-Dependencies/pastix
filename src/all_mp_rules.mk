##
##  General inference rules.
##

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





${OJTDIR}/%_s.cu : ${OJTDIR}/%.cu
	@echo -e ${to_S}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@
${OJTDIR}/%_s.c  : ${OJTDIR}/%.c
	@echo -e ${to_S}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@

${OJTDIR}/%_d.cu : ${OJTDIR}/%.cu
	@echo -e ${to_D}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@
${OJTDIR}/%_d.c  : ${OJTDIR}/%.c
	@echo -e ${to_D}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@

${OJTDIR}/%_c.cu : ${OJTDIR}/%.cu
	@echo -e ${to_C}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@
${OJTDIR}/%_c.c  : ${OJTDIR}/%.c
	@echo -e ${to_C}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@

${OJTDIR}/%_z.c  : ${OJTDIR}/%.c
	@echo -e ${to_Z}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@
${OJTDIR}/%_z.cu : ${OJTDIR}/%.cu
	@echo -e ${to_Z}                     > $@;
	@echo "#include \""$(notdir $<)"\"" >> $@



${OJTDIR}/%_ge.c : ${SRCDIR}/%.c
	@echo -e ${to_GE}                    > $@
	@echo "#include \""$(notdir $<)"\"" >> $@

${OJTDIR}/%_sy.c : ${SRCDIR}/%.c
	@echo -e ${to_SY}                    > $@
	@echo "#include \""$(notdir $<)"\"" >> $@

${OJTDIR}/%_he.c : ${SRCDIR}/%.c
	@echo -e ${to_HE}                    > $@
	@echo "#include \""$(notdir $<)"\"" >> $@

$(OJTDIR)/%$(OBJ) :	$(SRCDIR)/%.c config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC) $(CCTYPESFLT) -DCHOL_SOPALIN;
	$(MPCCPROG)   $(CCOPT) $(CCINC) $(CCTYPESFLT) -DCHOL_SOPALIN -c $(<) -o $(@)

$(OJTDIR)/%$(OBJ) :	$(OJTDIR)/%.c config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC) $(CCTYPESFLT) -DCHOL_SOPALIN;
	$(MPCCPROG)   $(CCOPT) $(CCINC) $(CCTYPESFLT) -DCHOL_SOPALIN -c $(<) -o $(@)

$(DSTDIR)/%$(EXE) :	$(SRCDIR)/%.c config.in
	$(MPCCPROG) $(CCOPT) $(CCINC) $(LKOPT) $(<) -o $(@) $(LKLIB)

ifeq (1, $(words $(findstring -DWITH_STARPU, ${CCOPT})))
ifeq (0, $(words $(findstring -DFORCE_NO_CUDA, ${CCOPT})))
$(OJTDIR)/%$(OBJ) : $(SRCDIR)/%.cu config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC) $(CCTYPESFLT);
	nvcc ${NVCCOPT} $(CCINC) $(CCTYPESFLT) -c $< -o $@

$(OJTDIR)/%$(OBJ) : $(OJTDIR)/%.cu config.in
	$(MAKEDEPEND) $(CCOPT) $(CCINC) $(CCTYPESFLT);
	nvcc ${NVCCOPT} $(CCINC) $(CCTYPESFLT) -c $< -o $@
endif
endif