MODULE:= order

MY_SRC = order.c			\
	order_base.c			\
	order_check.c			\
	order_io.c

MY_SRC_MULT_ARCH = 

include all_modules.mk
include all_rules.mk
