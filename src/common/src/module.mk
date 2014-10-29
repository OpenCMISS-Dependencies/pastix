MODULE:=common
MY_SRC := common_integer.c 	\
	common_error.c	 	\
	common_memory.c 	\
	trace.c

MY_SRC_MULT_ARCH := common.c

include all_modules.mk
include all_rules.mk
