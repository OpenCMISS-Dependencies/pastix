MODULE:=blend
MY_SRC =assemblyGener.c 		\
	blend.c 			\
	blend_symbol_cost.c		\
	blendctrl.c			\
	bulles.c			\
	cost.c	 			\
	costfunc.c 			\
	distribPart.c			\
	elimin.c 			\
	eliminfunc.c 			\
	extendVector.c 			\
	extrastruct.c 			\
	fanboth2.c			\
	param_blend.c			\
	partbuild.c			\
	queue.c				\
	simu.c				\
	smart_cblk_split.c		\
	solverMatrixGen.c		\
	solverRealloc.c 		\
	solver_check.c			\
	solver_io.c			\
	splitfunc.c			\
	splitpart.c			\
	splitpartlocal.c		\
	symbolrand.c			\
        task.c				\
	write_ps.c              \
	smart_cblk_split.c		\
	blend_distributeOnGPU.c

MY_SRC_MULT_ARCH =

OJTDIR:=blend/obj/$(HOSTARCH)
SRCDIR:=blend/src
DEPDIR:=blend/dep

include all_modules.mk
include all_rules.mk
