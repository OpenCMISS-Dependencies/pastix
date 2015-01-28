SET(PERF_EXECS
reg_axpy
reg_copy
reg_geam
reg_gemm
reg_gemv
reg_pof
reg_ppf
reg_scal
reg_trsm
reg_trsv
)
add_library(num_recipes num_recipes/covsrt.c num_recipes/gaussj.c num_recipes/lfit.c num_recipes/nrutil.c)
foreach (PERF_EXEC ${PERF_EXECS})
    add_executable(${PERF_EXEC} ${PERF_EXEC}.c)
    target_include_directories(${PERF_EXEC} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/num_recipes>)
    target_link_libraries(${PERF_EXEC} num_recipes)
    if (UNIX)
        target_link_libraries(${PERF_EXEC} m)
    endif()
endforeach()