# Build examples.
if (MPI_FOUND)
    include_directories(${MPI_CXX_INCLUDE_PATH})
    include_directories(${MPI_C_INCLUDE_PATH})
    include_directories(${MPI_Fortran_INCLUDE_PATH})
endif()

macro(ADDEXEC name src)
    #message(STATUS "adding executable ${name} with incdir: ${PASTIX_INCLUDE_DIR_PREC}")
    add_executable(${name} ${src})
    target_link_libraries(${name} pastix pastix-matrix-driver)
    target_include_directories(${name} PRIVATE
        ${PASTIX_INCLUDE_DIR_PREC}
    )    
endmacro()

macro(ALLEXECS precision)
    ADDEXEC(${precision}simple simple.c)
    ADDEXEC(${precision}simple_param simple_param.c)
    ADDEXEC(${precision}plot_memory_usage "plot_memory_usage.c;mem_trace.c")
    ADDEXEC(${precision}cppsimple cppsimple.cpp)
    if (MPI_FOUND)
        ADDEXEC(${precision}cppsimple_dist cppsimple_dist.cpp)
        ADDEXEC(${precision}simple_dist simple_dist.c)
        ADDEXEC(${precision}simple_dist_redist simple_dist_redist.c)
    endif()
endmacro()

# NoType examples still need an include dir (any)
SETINCDIR(${NEUTRAL_PRECISION})
ALLEXECS(notype)
foreach(prec s d c z)
    if(BUILD_PRECISION MATCHES [${prec}])
        execute_process(COMMAND sed -f sed_${prec}.in utils.h
            OUTPUT_FILE ${prec}utils.h
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
        SETINCDIR(${prec})
        ALLEXECS(${prec})
    endif()
endforeach()

ADDEXEC(schur schur.c)
ADDEXEC(schur2 schur2.c)
ADDEXEC(reentrant reentrant.c)
ADDEXEC(refinement refinement.c)
ADDEXEC(step-by-step step-by-step.c)
ADDEXEC(do_not_redispatch_rhs do_not_redispatch_rhs.c)
ADDEXEC(isolate_zeros isolate_zeros.c)
ADDEXEC(mult-rhs mult-rhs.c)
ADDEXEC(mult-rhs-dist mult-rhs-dist.c)
if (MPI_FOUND)
    ADDEXEC(step-by-step_dist step-by-step_dist.c)
    ADDEXEC(multi-comm multi-comm.c) 
    ADDEXEC(multi-comm-step multi-comm-step.c)
    # Fortran examples
    add_library(futils OBJECT utils.F90)
    add_dependencies(futils pastix-matrix-driver)
    target_include_directories(futils PRIVATE
        ${PASTIX_INCLUDE_DIR_PREC}
    ) 
    macro(ADDEXECF name src)
        add_executable(${name} ${src} $<TARGET_OBJECTS:futils>)
        target_link_libraries(${name} pastix pastix-matrix-driver)
        target_include_directories(${name} PRIVATE
            ${PASTIX_INCLUDE_DIR_PREC}
        )    
    endmacro()
    ADDEXECF(fstep-by-step fstep-by-step.F90)
    ADDEXECF(fsimple fsimple.F90)
endif()