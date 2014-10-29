  !File: fstep-by-step.F90
  !
  ! Simple example in Fortran calling PaStiX in step-by-step mode.
  ! Runs first step then 2 factorizations with 2 solves each.
  !
  
#include "pastix_fortran.h"
  
  ! Program: step_by_step_f
  ! 
  ! Parses the program options 
  ! Reads the corresponding matrix
  ! check the matrix
  ! Runs Pastix on it
  !
  ! parameters:
  ! 
  ! Gets parameters from command line (see <get_option> in <utils.f90>)
  !
  program step_by_step_f

    use utils
    implicit none

    include 'mpif.h'

    pastix_data_ptr_t                         :: pastix_data ! Structure to keep information in PaStiX (0 for first call)
    integer                                   :: pastix_comm ! MPI communicator used by pastix
    pastix_int_t                              :: n           ! Number of columns in the matrix
    pastix_int_t   ,dimension(:), allocatable :: ia          ! Index of first element of each column in ja and avals
    pastix_int_t   ,dimension(:), allocatable :: ja          ! Row of each element
    pastix_float_t ,dimension(:), allocatable :: avals       ! Value of each element
    pastix_int_t   ,dimension(:), allocatable :: perm        ! permutation tabular
    pastix_int_t   ,dimension(:), allocatable :: invp        ! reverse permutation tabular
    pastix_float_t ,dimension(:), allocatable :: rhs         ! Right hand side
    pastix_float_t ,dimension(:), allocatable :: rhssaved    ! Copy of Right hand side 
    pastix_int_t                              :: nrhs        ! right hand side number (only one possible)
    pastix_int_t                              :: iparm(IPARM_SIZE) ! Integer parameters
    double precision                          :: dparm(DPARM_SIZE) ! Floating poin parameters
    Integer                                   :: driver_num  ! Driver number
    Character(len=64)                         :: filename    ! Path to the matrix
    Character(len=4)                          :: type        ! type of the matrix
    Character(len=4)                          :: rhstype     ! type of the right-hand-side member
    pastix_int_t                              :: nbthread    ! Number of threads in PaStiX
    pastix_int_t                              :: verbose     ! Verbose level
    Integer                                   :: ierr        ! Error retrun value
    Integer                                   :: required    ! MPI thread level required  
    Integer                                   :: provided    ! MPI thread level provided 
    Integer                                   :: StatInfo    ! Info returned by MPI
    Integer                                   :: NbrOfFact   ! Number of factorizations
    Integer                                   :: NbrOfSolv   ! Number of rhs for solve
    Integer                                   :: rank        ! MPI rank
    Integer                                   :: i,j


    NbrOfFact = 2
    NbrOfSolv = 2

    !
    ! initiate MPI communication
    !
    pastix_comm = MPI_COMM_WORLD
    required=MPI_THREAD_MULTIPLE
    call MPI_Init_thread(required, provided, StatInfo)
    call MPI_Comm_rank  (pastix_comm, rank, StatInfo);
    !
    ! Get options ftom command line
    !
    Call get_option(driver_num, filename, nbthread, verbose, n)

    !
    ! reads the matrix
    !
    call read_matrix(driver_num, filename, &
         n, ia, ja, avals, rhs,            &
         Type, rhstype, pastix_comm, ierr)

    !
    ! First PaStiX call to initiate parameters
    !
    pastix_data = 0
    nrhs        = 1
    iparm(IPARM_MODIFY_PARAMETER) = API_NO
    iparm(IPARM_START_TASK)       = API_TASK_INIT
    iparm(IPARM_END_TASK)         = API_TASK_INIT

    call pastix_fortran(pastix_data ,pastix_comm, &
         n,ia,ja,avals,perm,invp,rhs,nrhs,iparm,dparm)


    !
    ! Customize some parameters
    !
    iparm(IPARM_THREAD_NBR) = nbthread  
    iparm(IPARM_VERBOSE)    = verbose

    if (type(2:2) == 'S') then

       iparm(IPARM_SYM)           = API_SYM_YES
       iparm(IPARM_FACTORIZATION) = API_FACT_LDLT

    else
       iparm(IPARM_SYM)           = API_SYM_NO
       iparm(IPARM_FACTORIZATION) = API_FACT_LU
    End if
    iparm(IPARM_MATRIX_VERIFICATION) = API_YES
    iparm(IPARM_RHS_MAKING)          = API_RHS_1

    allocate(perm(n))
    allocate(invp(n))
    !
    ! Call PaStiX first steps (Scotch - Fax - Blend
    !
    iparm(IPARM_START_TASK)       = API_TASK_ORDERING
    iparm(IPARM_END_TASK)         = API_TASK_ANALYSE

    call pastix_fortran(pastix_data ,pastix_comm, &
         n,ia,ja,avals,perm,invp,rhs,nrhs,iparm,dparm)

    
    allocate(rhssaved(n))
    rhssaved = rhs
    
    Do i = 1, NbrOfFact
    
       !
       ! Call PaStiX factorization
       !
       iparm(IPARM_START_TASK)       = API_TASK_NUMFACT
       iparm(IPARM_END_TASK)         = API_TASK_NUMFACT
       If (rank == 0) print *, "      > Factorisation number",i,"<"
       call pastix_fortran(pastix_data ,pastix_comm, &
            n,ia,ja,avals,perm,invp,rhs,nrhs,iparm,dparm)
       
       Do j = 1, NbrOfSolv
          !
          ! Call PaStiX updown and refinement
          !
          iparm(IPARM_START_TASK)       = API_TASK_SOLVE
          iparm(IPARM_END_TASK)         = API_TASK_REFINE
          ! rhs has been changed to solution by previous solve call
          rhs = rhssaved
          If (rank == 0) print *, "      >> Solve step number",i," <<"
          call pastix_fortran(pastix_data ,pastix_comm, &
               n,ia,ja,avals,perm,invp,rhs,nrhs,iparm,dparm)
       End Do
    End Do
    !
    ! Call PaStiX clean
    !
    iparm(IPARM_START_TASK)       = API_TASK_CLEAN
    iparm(IPARM_END_TASK)         = API_TASK_CLEAN
    
    call pastix_fortran(pastix_data ,pastix_comm, &
         n,ia,ja,avals,perm,invp,rhs,nrhs,iparm,dparm)
    
    deallocate(ia)
    deallocate(ja)
    deallocate(avals)
    deallocate(perm)
    deallocate(invp)
    deallocate(rhssaved)
    deallocate(rhs)
    
    ! Clean up MPI
    call MPI_FINALIZE(IERR)   


  end program step_by_step_f
