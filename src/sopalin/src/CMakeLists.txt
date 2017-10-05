SET(SOPALIN_SRC bordi.c
	 sopalin_thread.c
	 compute_context_nbr.c
)
SET(SRC_FAC
    sopalin3d.c		
	starpu_submit_tasks.c	
	csc_intern_compute.c	
	raff_functions.c
	starpu_updo.c)
SET(SOPALIN_SRC_MA
    coefinit.c
 	csc_intern_build.c	
	#csc_intern_compute.c	
	csc_intern_io.c		
	csc_intern_solve.c	
	csc_intern_updown.c	
	csc_utils.c		
	cscd_utils.c		
	cscd_utils_fortran.c    
	debug_dump.c		
	ooc.c			
	pastix.c		
	pastix_fortran.c	
	#sopalin3d.c
	sopalin_init.c		
	sopalin_option.c	
	sparse_gemm_cpu.c	
	tools.c
	${SRC_FAC}
)
if (USE_STARPU)
    LIST(APPEND SOPALIN_SRC starpu_pastix_sched_policy.c)
    if (USE_CUDA)
        LIST(APPEND SOPALIN_SRC_MA sparse_gemm.cu geadd_cuda.cu getra_cuda.cu)
    endif()
    if (USE_MAGMABLAS)
        LIST(APPEND SOPALIN_SRC_MA 
            zgetrf_stapiv_gpu.c
            zpotrf_stapiv_gpu.c
            zhetrf_stapiv_gpu.c
            zsytrf_stapiv_gpu.c)
    endif()
endif()

# add tests
if (BUILD_TESTS)
    macro(ADDTEST basename)
        add_executable(${basename} ${basename}.c)
        target_link_libraries(${basename} pastix ${NEUTRAL_PRECISION}matrix_driver)
        target_include_directories(${basename} PRIVATE
            ${PASTIX_INCLUDE_DIR_PREC}
            ../../common/src
        )
        target_compile_definitions(${basename} PRIVATE ${COMMON_DEFINITIONS})
        add_test(${basename} ${basename} ${args})
        foreach(prec s d c z)
            if(BUILD_PRECISION MATCHES [${prec}])
                SETINCDIR(${prec})
                add_executable(${prec}${basename} ${basename}.c)
                target_link_libraries(${prec}${basename} pastix ${prec}matrix_driver)
                target_include_directories(${prec}${basename} PRIVATE
                    ${PASTIX_INCLUDE_DIR_PREC}
                    ../../common/src
                )
                target_compile_definitions(${prec}${basename} PRIVATE ${COMMON_DEFINITIONS})
                add_test(${prec}${basename} ${prec}${basename} ${args}) 
            endif()
        endforeach()
    endmacro()
    macro(ADDTEST_MPI basename np)
        add_executable(${basename} ${basename}.c)
        target_link_libraries(${basename} pastix ${NEUTRAL_PRECISION}matrix_driver mpi-c)
        target_compile_definitions(${basename} PRIVATE NODESIZE=8 ${COMMON_DEFINITIONS})
        target_include_directories(${basename} PRIVATE
            ${PASTIX_INCLUDE_DIR_PREC}
            ${MPI_C_INCLUDE_PATH}
            ../../common/src
        ) 
        add_test(${basename} 
           ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${np} ${CMAKE_CURRENT_BINARY_DIR}/${basename})
    endmacro()
    #ADDTEST(test_gemm) seems to hang up
    ADDTEST(test_axpy)
    ADDTEST(test_copy)
    ADDTEST(test_scal)
    ADDTEST(test_gemv)
    ADDTEST(test_trsv)
    if (USE_CUDA)
        ADDTEST(test_sparse_gemm "1000;100;100")
        ADDTEST(test_geadd "1000;10000")
        ADDTEST(test_getra "1000;10000")
        ADDTEST(test_getrf "1000;10000")
    endif()
    if (MPI_FOUND)
        SET(NODESIZE2 16)
        ADDTEST_MPI(pingpong_1 2)
        ADDTEST_MPI(pingpong_1e ${NODESIZE2})
        ADDTEST_MPI(pingpong_2 4)
        ADDTEST_MPI(pingpong_2e ${NODESIZE2})
        ADDTEST_MPI(pingpong_4 8)
        ADDTEST_MPI(pingpong_4e ${NODESIZE2})
        ADDTEST_MPI(pingpong_8 16)
        ADDTEST_MPI(pingpong_8e ${NODESIZE2})
        ADDTEST_MPI(pingpong_n1 2)
        ADDTEST_MPI(pingpong_n1e ${NODESIZE2})
        ADDTEST_MPI(pingpong_n2 4)
        ADDTEST_MPI(pingpong_n2e ${NODESIZE2})
        ADDTEST_MPI(pingpong_n4 8)
        ADDTEST_MPI(pingpong_n4e ${NODESIZE2})
        ADDTEST_MPI(pingpong_n8 16)
        ADDTEST_MPI(pingpong_n8e ${NODESIZE2})
    endif()
endif()
