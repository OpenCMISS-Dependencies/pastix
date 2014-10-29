!
! File: fmurge.F90
!
! Example that generate A Laplacian with multiple
! degrees of freedom and solves it.
!
! That example only works with the -DDISTRIBUTED option.
!
! The laplacian looks like :
!   >    1       0                      ... 0
!   > -2(n-1)  (n-1)  -2(n-1)   0       ... 0
!   >    0    -2(n-1)  (n-1)  -2(n-1) 0 ... 0
!   >              ....
!   >              ....
!   >    0 ...              -2(n-1)  n-1  -2(n-1)
!   >    0 ...                        0     1
!
! The Right-hand-side member is defined by
!   $$RHS_i = - 4\pi^2 sin(2\pi x)$$
!
! The solution should be $$X_i = sin(2\pi x)$$
!
! The solution is stored in the file "result.plot" that can be ploted using
! > gnuplot> plot "./result.plot" u 2:3 t "reference", "./result.plot" u 2:4 t "Result, first degree of freedom"
!
!
! Usage:
!   > ./fmurge <size> <DofNbr>
!
! Authors:
!   Pascal JACQ    - jacq@labri.fr
!   Xavier LACOSTE - lacoste@labri.fr
!

PROGRAM main
  !$ use omp_lib, only : omp_get_num_threads
  IMPLICIT NONE
! will be replaced during compilation by replaceCOEF.sh

  INCLUDE "mpif.h"
  INCLUDE "murge.inc"
  ! MPI Data
  INTEGER(KIND=MURGE_INTS_KIND)          :: ierr
  INTEGER       :: Me, NTasks, required, provided
  ! MURGE Data


  INTEGER(KIND=MURGE_INTS_KIND)                         :: id, m, job, mode
  ! CSC Data
  INTEGER(KIND=MURGE_INTS_KIND)                         :: n, dof
  INTEGER(KIND=MURGE_INTL_KIND)                         :: nnzeros, edgenbr
  ! __COEF__ will be replaced during compilation by replaceCOEF.sh
  __COEF__                                                  :: val
  ! Local Data
  INTEGER(KIND=MURGE_INTS_KIND)                         :: localnodenbr,interior,taille
  INTEGER(KIND=MURGE_INTS_KIND), DIMENSION(:) , POINTER :: nodelist
  INTEGER(KIND=MURGE_INTS_KIND)                         :: root
  INTEGER(KIND=MURGE_INTS_KIND)                         :: base
  __COEF__, DIMENSION(:) , POINTER                          :: lrhs
  __COEF__, DIMENSION(:) , POINTER                          :: globrhs
  __COEF__, DIMENSION(:) , POINTER                          :: globrhs_recv
  __COEF__, DIMENSION(:) , POINTER                          :: globx
  __COEF__, DIMENSION(:) , POINTER                          :: globprod
  ! Other data
  __COEF__, DIMENSION(:) , POINTER                          :: expand
  __COEF__                                                  :: val2
  REAL(KIND=MURGE_REAL_KIND)                            :: prec
  REAL(KIND=8)                                          :: xmin, xmax
  INTEGER(KIND=MURGE_INTS_KIND)                         :: NArgs, i, j, k, myfirstrow, mylastrow, l
  CHARACTER(LEN=100)                                    :: args
  CHARACTER(LEN=20)                                     :: atester
  INTEGER(KIND=MURGE_INTS_KIND)                         :: solver
  INTEGER(KIND=MURGE_INTS_KIND)                         :: paramidx
  INTEGER(KIND=MURGE_INTS_KIND)                         :: zero=0
  INTEGER(KIND=MURGE_INTS_KIND)                         :: one=1
  INTEGER(KIND=MURGE_INTS_KIND)                         :: two=2
  INTEGER(KIND=MURGE_INTS_KIND)                         :: nb_threads

  NArgs = command_argument_count()
  root = -1
  base = 1

  required=MPI_THREAD_MULTIPLE
  call MPI_Init_thread(required,provided,ierr)

  CALL MPI_COMM_SIZE(MPI_Comm_world, NTasks,ierr )
  CALL MPI_COMM_RANK(MPI_Comm_world, Me, ierr)

  IF (NArgs >= 2) THEN
     CALL GETARG(1,args)
     READ(args,*) n
     ! degrees of freedom
     CALL GETARG(2,args)
     READ(args,*) dof
  ELSE
     IF (me == 0) THEN
	PRINT *, "Usage: ./fmurge <size> <DofNumber>"
     END IF
     CALL Abort()
  END IF

  xmin = 0.
  xmax = 1.

  ! Starting MURGE
  CALL MURGE_INITIALIZE(one, ierr)
  if (ierr /= MURGE_SUCCESS) call abort()
  id = 0

  ! Set Options
  prec = 1e-7
  !
  ! Call MURGE_Get_Solver(solver)
  !
  solver = MURGE_SOLVER_PASTIX




  IF (solver == MURGE_SOLVER_PASTIX) THEN
     CALL MURGE_SetDefaultOptions(id, zero, ierr)
     CALL MURGE_SetOptionINT(id, IPARM_VERBOSE, API_VERBOSE_NO, ierr)
     CALL MURGE_SetOptionINT(id, IPARM_MATRIX_VERIFICATION, API_YES, ierr)
     CALL MURGE_SetOptionINT(id, IPARM_FREE_CSCUSER, API_NO, ierr)
     nb_threads = 1
     !$OMP PARALLEL shared(nb_threads)
     !$ nb_threads = omp_get_num_threads()
     !$OMP END PARALLEL
     if (Me .eq. 0) then
	print *, "Running on", nb_threads,"threads and", NTasks, "MPI Tasks"
     end if
     CALL MURGE_SetOptionINT(id, IPARM_THREAD_NBR, nb_threads, ierr)
  ELSE IF (solver == MURGE_SOLVER_HIPS) THEN
#ifdef HIPS
     IF (method == 1) THEN
	CALL MURGE_SetDefaultOptions(id, HIPS_ITERATIVE, ierr)
     ELSE
	CALL MURGE_SetDefaultOptions(id, HIPS_HYBRID, ierr )
	CALL MURGE_SetOptionINT(id, HIPS_PARTITION_TYPE, zero, ierr)
	CALL MURGE_SetOptionINT(id, HIPS_DOMSIZE, domsize, ierr)
     END IF
     CALL MURGE_SetOptionINT(id, HIPS_SYMMETRIC, zero, ierr)
     CALL MURGE_SetOptionINT(id, HIPS_LOCALLY, zero, ierr)
     CALL MURGE_SetOptionINT(id, HIPS_ITMAX, itmax, ierr)
     CALL MURGE_SetOptionINT(id, HIPS_KRYLOV_RESTART, restart, ierr)
     CALL MURGE_SetOptionINT(id, HIPS_VERBOSE, verbose, ierr)
     CALL MURGE_SetOptionINT(id, HIPS_DOMNBR, NTasks, ierr)
     CALL MURGE_SetOptionINT(id, HIPS_CHECK_GRAPH, one, ierr)
     CALL MURGE_SetOptionINT(id, HIPS_CHECK_MATRIX, one, ierr)
#endif
  END IF
  CALL MURGE_SetOptionINT(id, MURGE_IPARAM_DOF, dof, ierr)
  CALL MURGE_SetOptionINT(id, MURGE_IPARAM_SYM, MURGE_BOOLEAN_FALSE, ierr)
  CALL MURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL, base, ierr)

  CALL MURGE_SetOptionREAL(id, MURGE_RPARAM_EPSILON_ERROR, PREC, ierr)
  ! Set the graph : all processor enter some edge of the
  ! graph that corresponds to non-zeros location in the matrix


  ! ****************************************
  ! ** Enter the matrix non-zero pattern  **
  ! ** you can use any distribution       **
  ! ****************************************

  ! this processor enters the A(myfirstrow:mylastrow, *)
  ! part of the matrix non-zero pattern
  IF (me == 0) THEN
     edgenbr = 3*n-4

     CALL MURGE_GRAPHBEGIN(id, n, edgenbr, ierr)

     ! Dirichlet boundary condition
     CALL MURGE_GRAPHEDGE(id, one, one, ierr)
     CALL MURGE_GRAPHEDGE(id, n, n, ierr)

     ! Interior
     DO i = 2, n-1
	DO j = -1,1
	   CALL MURGE_GRAPHEDGE(id, i, i+j, ierr)
	   !IF (j /= 0) THEN
	   !CALL MURGE_GRAPHEDGE(id, j+i, i, ierr)
	   !END IF
	END DO
     END DO
  ELSE
     edgenbr = 0
     CALL MURGE_GRAPHBEGIN(id, n, edgenbr,  ierr)
  END IF
  CALL MURGE_GRAPHEND(id, ierr)

  ! Get Local nodes
  CALL MURGE_GETLOCALNODENBR(id, localnodenbr, ierr)
  ALLOCATE(nodelist(localnodenbr))

  CALL MURGE_GETLOCALNODELIST(id, nodelist, ierr)

  mode = 0

  ! job = 0 : initialize the matrix with this coefficient
  ! job = 1 : these coefficients are added to the matrix
  job = 0

  ! compute the number of non-zeros;
  nnzeros = 0
  DO m = 1, localnodenbr
     i = nodelist(m)
     IF (i == 1 .OR. i == n) THEN
	! Boundaries
	nnzeros = nnzeros + 1
     ELSE
	! Interior
	DO k = -1, 1
	   nnzeros = nnzeros + 1
	END DO
     END IF
  END DO

  ! We are using dof so a non zero is in fact a block of size dof**2
  nnzeros = nnzeros * dof**2

  ! You can enter only coefficient (i,j) that are in A(nodelist, nodelist)
  ! on this processor

  ! We enter the lower and upper triangular part of the matrix so sym = 0

  ! expand is the identity matrix of size 'dof' stored by line
  ALLOCATE(expand(dof**2))

  k = 1
  expand = 0.
  DO i = 1,dof
     DO j = 1,dof
	IF (i == j) expand(k) = 1.
!        expand(k)= k
	k = k + 1
     END DO
  END DO

  CALL MURGE_ASSEMBLYBEGIN(id, n, nnzeros, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW, &
       MURGE_ASSEMBLY_FOOL, MURGE_BOOLEAN_FALSE, ierr)
  !$OMP PARALLEL default(none)&
  !$OMP private(m, i, val, ierr, k) &
  !$OMP shared(expand, localnodenbr, n, xmin, xmax, nodelist, id)
  !$OMP DO
  DO m = 1, localnodenbr
     i = nodelist(m)
     IF (i == 1 .OR. i == n) THEN
	! Boundaries
	CALL GetCoef(val,i,i,xmin,xmax,n)
	CALL MURGE_ASSEMBLYSETNODEVALUES(id, i, i, val*expand, ierr)
	!print *, i, i, val*expand
     ELSE
	DO k = -1,1
	   CALL GetCoef(val,i+k,i,xmin,xmax,n)
	   CALL MURGE_ASSEMBLYSETNODEVALUES(id, i, i+k, val*expand, ierr)
	   !print *, i, i+k, val*expand
	END DO
     END IF
  END DO
  !$OMP END DO
  !$OMP END PARALLEL
  CALL MURGE_ASSEMBLYEND(id, ierr)


  ! We expand the rhs
  ALLOCATE(lrhs(localnodenbr*dof))
  ALLOCATE(globrhs(n*dof))
  globrhs = 0.0
  k = 1
  DO m = 1,localnodenbr
     CALL GetRhs(val,nodelist(m),xmin,xmax,n)
     globrhs((nodelist(m)-1)*dof+1:(nodelist(m)-1)*dof+dof) = val
     lrhs(k:k+dof-1) = val
     DO l = 1, dof
	!print *, "rhs", (nodelist(m)-1)*dof+l, lrhs((m-1)*dof+l)
     END DO
     k = k + dof
  END DO
  ALLOCATE(globrhs_recv(n*dof))
  if (__IS_COMPLEX__) then ! __IS_COMPLEX__will be replaced during compilation by replaceCOEF.sh
     if (MURGE_COEF_KIND==4) then
	CAll MPI_Allreduce(globrhs, globrhs_recv, n*dof, MPI_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
     else
	CAll MPI_Allreduce(globrhs, globrhs_recv, n*dof, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD, ierr)
     end if
  else
     if (MURGE_COEF_KIND==4) then
	CAll MPI_Allreduce(globrhs, globrhs_recv, n*dof, MPI_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
     else
	CAll MPI_Allreduce(globrhs, globrhs_recv, n*dof, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
     end if
  end if
  DEALLOCATE(globrhs)
  CALL MURGE_SETLOCALRHS(id, lrhs, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW, ierr)


  ! Get the global solution
  ALLOCATE(globx(n*dof))
  CALL MURGE_GETGLOBALSOLUTION(id, globx, root, ierr)

  CALL MURGE_SETGLOBALRHS(id, globx, -one, MURGE_ASSEMBLY_OVW, ierr)
  ALLOCATE(globprod(n*dof))
  CALL MURGE_GETGLOBALPRODUCT(id, globprod, -one, ierr)

  print *, "||AX - B||/||AX||  :", SQRT(SUM((globprod(:) - globrhs_recv(:))**2)/SUM(globprod(:)**2))

  ! Store in a file
  IF (Me == 0) THEN
     CALL store(globx,xmin,xmax,dof)
  END IF

  ! I'm Free
  CALL MURGE_CLEAN(id, ierr)
  CALL MURGE_FINALIZE(ierr)
  CALL MPI_FINALIZE(ierr)

  DEALLOCATE(expand, nodelist, lrhs, globx, globprod, globrhs_recv)

CONTAINS

  !
  ! Subroutine: GetCoef
  !
  ! Computes the value for a given coefficient.
  !
  ! Parameters:
  !   val  - Value to set
  !   i    - Row of the coefficient.
  !   j    - Column of the coefficient.
  !   xmin - Minimum value of the interval.
  !   xmax - Maximum value of the interval.
  !   n    - Number of points in the interval.
  SUBROUTINE GetCoef(val,i,j,xmin,xmax,n)
    __COEF__        , INTENT(OUT) :: val
    INTEGER(KIND=MURGE_INTS_KIND)        , INTENT(IN)  :: i, j, n
    REAL(KIND=8), INTENT(IN)  :: xmin, xmax

    REAL(KIND=8) :: dx_1

    dx_1 = (n-1.)/ (xmax - xmin)

    val = 0.

    IF (i==j) THEN
       IF (i==1 .OR. i == n) THEN
	  ! Boundary Condition (Dirichlet)
	  val = 1
       ELSE
	  ! Interior diagonnal part
	  val = -2 * dx_1
       END IF
    ELSE
       val = dx_1
    END IF
  END SUBROUTINE GetCoef


  !
  ! Subroutine: GetRhs
  !
  ! computes the value of a coefficient of the Right-hand-side member.
  !
  ! Parameters:
  !   val  - Value to set.
  !   i    - Index of the value.
  !   xmin - Minimum value of the interval.
  !   xmax - Maximum value of the interval.
  !   n    - Number of points in the interval.
  SUBROUTINE GetRhs(val,i,xmin,xmax,n)
    __COEF__                , INTENT(OUT) :: val
    INTEGER(KIND=MURGE_INTS_KIND)                , INTENT(IN)  :: i, n
    REAL(KIND=8)           , INTENT(IN)  :: xmin, xmax

    REAL(KIND=8) :: dx , x,Pi
    INTEGER :: c

    dx = (xmax - xmin) / (n-1.)
    x = xmin + (i-1)*dx

    Pi =ACOS(-1.)

    val = -4*Pi**2*SIN(2*Pi*x)
    ! Boundary Condition (Dirichlet)
    IF (i == n .OR. i==1) THEN
       val = 0.
    ELSE
       val = dx*val
    END IF

  END SUBROUTINE GetRhs

  !
  ! Subroutine: store
  !
  ! Write the solution into a file result.plot.
  !
  ! The file contains :
  ! > k x sin(2 \pi x) sol(k:k+dofnbr-1)
  ! Where k goes from 1 to n*dofnbr, dofnbr by dofnbr and x goes from
  ! xmin to xmax, with a step of (xmax - xmin) / (n-1).
  !
  ! Parameters:
  !   sol  - The solution of the problem
  !   xmin - The minimum value of the interval.
  !   xmax - The maximum value of the interval.
  !   dof  - The Number of degree of freedom.
  !
  SUBROUTINE store(sol,xmin,xmax,dof)
    __COEF__,         DIMENSION(:), INTENT(IN) :: sol
    REAL(KIND=8),               INTENT(IN) :: xmin,xmax
    INTEGER(KIND=MURGE_INTS_KIND),                       INTENT(IN) :: dof
    REAL(KIND=8)                           :: x,dx,Pi2,s
    INTEGER(KIND=MURGE_INTL_KIND)                                   :: i,j,k, n
    CHARACTER(len=100)                     :: ecriture

    n = SIZE(sol) / dof
    x = xmin
    dx = (xmax - xmin) / (n-1.)
    Pi2 = 2.*ACOS(-1.)

    WRITE(ecriture,*) '(i9,X,',dof+2,'(E15.8,X))'

    OPEN(unit = 20, file="result.plot")
    k = 1
    DO i = 1, n
       WRITE(20,ecriture) k,x,SIN(Pi2*x),sol(k:k+dof-1)
       k = k + dof
       x = x + dx
    END DO
    CLOSE(20)
  END SUBROUTINE store
END PROGRAM main
