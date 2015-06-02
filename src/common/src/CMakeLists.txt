SET(COMMON_SRC
    common_integer.c
	common_error.c
	common_memory.c
	trace.c)
SET(COMMON_SRC_MA common.c)

if (MINGW)
    list(APPEND COMMON_SRC mingw_fix.c)
endif()