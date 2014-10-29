  ! File: utils.f90
  !
  ! Utilities for Fortran exemple
  ! Print usage (see <usage>)
  ! Get option from commant line (see <get_option>)
  ! Reads the matrix (see <read_matrix>)
  
#include "pastix_fortran.h"
#define RSA       1
#define LAPLACIAN 2
#define STR_SIZE 64
  
  Module utils
    Implicit None
  Contains 

    ! Subroutine: usage
    !
    ! Prints usage
    !
    Subroutine usage()
      Character(len=STR_SIZE) :: string
      Call getarg(0,string)
      Print *, "Usage : " , Trim(string) ," [option] "
      Print *, "   options : "
      Print *, "     -rsa    [filename], default driver RSA (use Fortran) "
      Print *, "     -lap    <integer>           generate a laplacian of size <integer>"
      Print *, "     -t      <integer>           define thread number"
      Print *, "     -v      <integer>           verbose level (1, 2 or 3)"
      Print *, "     -h                          print this help"
    End Subroutine usage

    ! Subroutine: get_option
    !
    ! Gets options from command line
    ! 
    ! Parameters:
    !    driver_num - integer giving driver to use to read the matrix
    !    filename   - Path to the file to read
    !    nbthread   - Number of threads in PaStiX
    !    verbose    - Verbose mode 1, 2 or 3
    !    matsize    - size of the matrix wanted
    ! 
    Subroutine get_option(driver_num, filename, nbthread, verbose, matsize)
      Integer,                Intent(out) :: driver_num
      Character(len=STR_SIZE),Intent(out) :: filename
      pastix_int_t,           Intent(out) :: nbthread
      pastix_int_t,           Intent(out) :: verbose
      pastix_int_t,           Intent(out) :: matsize
      Integer                             :: narg
      Integer                             :: i 
      Character(len=STR_SIZE)             :: string1
      Character(len=STR_SIZE)             :: string2

      verbose = API_VERBOSE_NO

      narg = command_argument_count()
      If (narg == 0) Then
         Call usage()
         Stop
      End If

      i = 1
      nbthread = 1;
      Do While(i < narg +1)
         Call get_command_argument(i,string1)
         Select Case(string1)
         Case ("-lap")
            driver_num = LAPLACIAN
            Call get_command_argument(i+1,string2)
            Read (string2,*) matsize
            i = i + 1
         Case ("-rsa")
            driver_num = RSA
            Call get_command_argument(i+1,filename)
            i = i + 1
            If (filename(1:1) == '-') Then
               Call usage()
               Stop
            End If
         Case ("-t")
            Call get_command_argument(i+1,string2)
            Read (string2,*) nbthread
            i = i + 1
         Case ("-v")
            Call get_command_argument(i+1,string2)
            Read (string2,*) verbose
            verbose = verbose - 1
            i = i + 1
         Case ("-h")
            Call usage()
            Stop
         Case default
            Call usage()
            Stop 
         End Select

         i = i + 1
      End Do
    End Subroutine get_option

    ! Subroutine: read_matrix
    !
    ! reads the matrix, only RSA supported yet.
    ! (see <read_rsa>)
    !
    ! Parameters:
    !    driver_type - integer corresponding to driver to use to read the matrix
    !    filename    - Path to the file to read
    !    n           - size of the matrix
    !    ia          - Index of first element of each column in *ja* and *vals*
    !    ja          - Column of each element of the matrix.
    !    val         - Value of each element of the matrix.
    !    rhs         - right-hand-side member 
    !    type        - type of the matrix
    !    rhstype     - type of the right-hand-side member 
    !    pastix_comm - MPI communicator
    !    ierr        - return value
    !
    Subroutine read_matrix(driver_type, filename, &
         n, ia, ja, val, rhs, type, rhstype, pastix_comm, ierr)

#ifndef FORCE_NOMPI
      include 'mpif.h'
#endif
      Integer,                                    Intent(in)    :: driver_type 
      Character(len=STR_SIZE),                    Intent(in)    :: filename
      pastix_int_t,                               Intent(inout) :: n
      pastix_int_t,    Dimension(:), Allocatable, Intent(out)   :: ia  ! Index of 1st element of each column in ja and avals
      pastix_int_t,    Dimension(:), Allocatable, Intent(out)   :: ja  ! Row of each element
      pastix_float_t,  Dimension(:), Allocatable, Intent(out)   :: val ! Value of each element
      pastix_float_t,  Dimension(:), Allocatable, Intent(out)   :: rhs ! right hand side
      Character(len=4),                           Intent(out)   :: type
      Character(len=4),                           Intent(out)   :: rhstype
      Integer,                                    Intent(in)    :: pastix_comm
      Integer,                                    Intent(out)   :: ierr
      Integer                                                   :: mpid
      pastix_int_t                                              :: nnzero
      Integer                                                   :: i
      Integer                                                   :: j

      ierr = 0
#ifndef FORCE_NOMPI
      call MPI_Comm_rank(pastix_comm,mpid, ierr)
#else
      mpid = 0
#endif
      If (ierr == 1) Then
         return 
      End If

      If (mpid == 0) Then

         Select Case(driver_type)
         Case(RSA)
            Call read_rsa(filename, n, nnzero,  ia, ja, val, type, rhstype, ierr)
         Case(LAPLACIAN)
            Call genlaplacian(n, nnzero,  ia, ja, val, rhs, type, rhstype, ierr)
         Case default
            ierr = 1
            Return
         End Select

         !         print *, n, "columns", nnzero, "non zeros"  


         If (driver_type /= LAPLACIAN) Then
            allocate(rhs(n))
            Do i = 1, n 
               Do j = ia(i), ia(i+1)-1
                  rhs(ja(j)) = val(j)
                  If (type(2:2) == 'S') Then
                     rhs(i) = val(j)
                  End If
               End Do
            End Do
         End If
      End if

#ifndef FORCE_NOMPI
      call MPI_Bcast(n,1,MPI_PASTIX_INT,0,pastix_comm,ierr)
      call MPI_Bcast(nnzero,1,MPI_PASTIX_INT,0,pastix_comm,ierr)
#endif

      If (mpid /= 0) Then
         allocate(ia(n+1))
         allocate(ja(nnzero))
         allocate(val(nnzero))
         allocate(rhs(n))
      End If
#ifndef FORCE_MPI
      !      print *, mpid, "ia"
      call MPI_Bcast(ia, n+1,    MPI_PASTIX_INT,   0, pastix_comm, ierr)
      !      print *, mpid, ia(1)
      !      print *, mpid, "ja"
      call MPI_Bcast(ja, nnzero, MPI_PASTIX_INT,   0, pastix_comm, ierr)
      !      print *, mpid, ja(1)
      !      print *, mpid, "val"
      call MPI_Bcast(val,nnzero, MPI_PASTIX_FLOAT, 0, pastix_comm, ierr)
      !      print *, mpid, val(1)
      !      print *, mpid, "rhs"
      call MPI_Bcast(rhs,n     , MPI_PASTIX_FLOAT, 0, pastix_comm, ierr)
      !      print *, mpid, "type"
      call MPI_Bcast(type,    3, MPI_CHARACTER,    0, pastix_comm, ierr)
#endif
    End Subroutine read_matrix

    ! Subroutine: read_rsa
    !
    ! Reads a matrix in RSA format
    ! 
    ! Parameters: 
    !    filename - Path to the file to read from.
    !    n        - Size of the matrix
    !    nnzero   - Number of non zeros
    !    ia       - Index of first element of each column in *ja* and *vals*
    !    ja       - Column of each element of the matrix.
    !    val      - Value of each element of the matrix.
    !    rhs      - right-hand-side member 
    !    type     - type of the matrix
    !    rhstype  - type of the right-hand-side member 
    !    ierr     - return value
    !
    Subroutine read_rsa(filename, n, nnzero, ia, ja, val, type, rhstype, ierr)
      Character(len=STR_SIZE),                    Intent(in)  :: filename
      pastix_int_t,                               Intent(out) :: n
      pastix_int_t,    Dimension(:), Allocatable, Intent(out) :: ia  ! Index of first element of each column in ja and avals
      pastix_int_t,    Dimension(:), Allocatable, Intent(out) :: ja  ! Row of each element
      pastix_float_t,  Dimension(:), Allocatable, Intent(out) :: val ! Value of each element
      Character(len=3),                           Intent(out) :: Type
      Character(len=3),                           Intent(out) :: rhstype
      pastix_int_t,                               Intent(out) :: nnzero
      Integer,                                    Intent(out) :: ierr
      Integer                                                 :: ncol
      Integer                                                 :: nrow
      Integer                                                 :: nnz
      Integer                                                 :: nrhs
      Character(len=8)                                        :: key
      Character(len=72)                                       :: title
      Integer                                                 :: tmp
      Integer         ,Dimension(:), Allocatable              :: tmpia
      Integer         ,Dimension(:), Allocatable              :: tmpja
      real*8          ,Dimension(:), Allocatable              :: tmpval
      real*8          ,Dimension(:), Allocatable              :: tmprhs
      Integer                                                 :: i 

      ierr = 0
      call read_rsa_header(filename, nrow, ncol, nnz, type, rhstype, ierr)
      nnzero = nnz

      If (nrow /= nrow .or. ierr == 1) Then 
         ierr = 1
         return
      end If
      n = ncol

      allocate(tmpia(ncol+1))
      allocate(tmpja(nnzero))
      allocate(tmpval(nnzero))

      tmp = 2; ! job without rhs
      call wreadmtc(nrow,nnz,tmp,filename,STR_SIZE,tmpval,tmpja,tmpia,tmprhs, &
           nrhs,rhstype,nrow,ncol,nnz,title,key,type,ierr)

      If (ierr == 1) then
         return
      end If

      allocate(ia(ncol+1))

      Do i = 1, ncol +1
         ia(i) = tmpia(i)
      End Do

      deallocate(tmpia)

      allocate(ja(nnzero))

      Do i = 1, nnzero
         ja(i) = tmpja(i)
      End Do

      deallocate(tmpja)

      allocate(val(nnzero))  
      Do i = 1, nnzero
         val(i) = tmpval(i)
      End Do

      deallocate(tmpval)



    End Subroutine read_rsa
    ! Subroutine: read_rsa_header
    !
    ! Reads a matrix header in RSA format
    ! 
    ! Parameters: 
    !    filename - Path to the file to read from.
    !    nrow     - Number of rows
    !    ncol     - Number of columns
    !    nnzero   - Number of non zeros
    !    type     - type of the matrix
    !    rhstype  - type of the right-hand-side member 
    !    ierr     - return value
    !
    Subroutine read_rsa_header(filename, nrow, ncol, nnzero, type, rhstype, ierr)
      Character(len=STR_SIZE),                    Intent(in)  :: filename
      Integer,                                    Intent(out) :: ncol
      Integer,                                    Intent(out) :: nrow
      Integer,                                    Intent(out) :: nnzero
      Character(len=3),                           Intent(out) :: Type
      Character(len=3),                           Intent(out) :: rhstype
      Integer,                                    Intent(out) :: ierr
      pastix_int_t,    Dimension(0)                           :: ia  ! Index of first element of each column in ja and avals
      pastix_int_t,    Dimension(0)                           :: ja  ! Row of each element
      pastix_float_t,  Dimension(0)                           :: val ! Value of each element
      pastix_float_t,  Dimension(0)                           :: rhs ! right hand side
      Integer                                                 :: tmp
      Integer                                                 :: nrhs
      Character(len=8)                                        :: key
      Character(len=72)                                       :: title
      ierr = 0
      tmp  = 0
      call wreadmtc(tmp,tmp,tmp,filename,STR_SIZE,val,ja,ia,rhs,nrhs, &
           rhstype,nrow,ncol, nnzero,title,key,type,ierr)

    End Subroutine read_rsa_header


!!$  Function: genlaplacien
!!$
!!$  Generate a laplacien of size *n*
!!$
!!$  Parameters:
!!$    n  	    - Size of the wanted matrix
!!$    nnzeros - Number of non zeros in the produced matrice
!!$    ia      - Index of first element of each column in *row* and *val*
!!$    ja 	    - Row of eah element				       
!!$    avals   - Value of each element				       
!!$    rhs     - Right-hand-side member
!!$    type    - Type of the matrix				       	 
!!$    rhstype - Type of the right hand side.			       	 

    subroutine genlaplacian(n, nnzeros, ia, ja, avals, rhs, type, rhstype,success)
      pastix_int_t,     intent(in)                              ::    n
      pastix_int_t,     intent(out)                             :: nnzeros
      pastix_int_t,     dimension(:) , allocatable, intent(out) :: ia,ja
      pastix_float_t,   dimension(:) , allocatable, intent(out) :: avals,rhs
      character(len=4)         , intent(out)                    :: type
      character(len=4)         , intent(out)                    :: rhstype
      integer                                                   :: success
                                                               
      pastix_int_t                                              :: i, j

      nnzeros = 3*n - 2
      ! Allocating
      allocate( ia     (n+1)    )
      allocate( ja     (nnzeros))
      allocate( avals  (nnzeros))
      allocate( rhs    (n)      )

      ! Building ia, ja and avals and rhs
      j=1
      do i = 1, n
         ia(i) = j
         ! /* ONLY triangular inferior matrix */
         ! /*       if (i != 0) */
         ! /* 	{ */
         ! /* 	  (*ja)[j]    = i; */
         ! /* 	  (*avals)[j] = -1; */
         ! /* 	  j++; */
         ! /* 	} */
         ja(j)    = i
         avals(j) = 2
         j=j+1
         if (i /= n) then
            ja(j)    = i+1
            avals(j) = -1.
            j=j + 1
         end if

         rhs(i) = 0
      end do
      ia(n+1) = j
      rhs(1)  = 1
      rhs(n)  = 1

      ! type and rhstype
      type         = "RSA"
      type(4:4)    = char(0)
      rhstype(1:1) = char(0)

      success = 0
    end subroutine genlaplacian



  End Module utils

