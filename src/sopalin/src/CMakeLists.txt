SET(SOPALIN_SRC bordi.c
	 sopalin_thread.c
	 compute_context_nbr.c
	 coefinit.c
)
SET(SOPALIN_SRC_MULTIARCH	 		
 	csc_intern_build.c	
	csc_intern_compute.c	
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
	sopalin3d.c		
	sopalin_init.c		
	sopalin_option.c	
	sparse_gemm_cpu.c	
	tools.c
	sopalin3d.c		
	starpu_submit_tasks.c	
	csc_intern_compute.c	
	raff_functions.c	
	starpu_updo.c
)
if (USE_STARPU)
    LIST(APPEND SOPALIN_SRC starpu_pastix_sched_policy.c)
    if (USE_CUDA)
        LIST(APPEND SOPALIN_SRC_MULTIARCH sparse_gemm.cu geadd_cuda.cu getra_cuda.cu)
    endif()
    if (USE_MAGMABLAS)
        LIST(APPEND SOPALIN_SRC_MULTIARCH 
            zgetrf_stapiv_gpu.c
            zpotrf_stapiv_gpu.c
            zhetrf_stapiv_gpu.c
            zsytrf_stapiv_gpu.c)
    endif()
endif()

SET(SOPALIN_SRC ${SOPALIN_SRC} ${SOPALIN_SRC_MULTIARCH})