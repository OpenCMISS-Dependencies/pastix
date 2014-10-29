!******************************************************************************
!* File: murge_pastix.inc                                                     *
!*                                                                            *
!* Fortran function declaration for Murge function used only in PaStiX.       *
!*                                                                            *
!* Authors:                                                                   *
!*   Xavier Lacoste - xavier.lacoste@inria.fr                                 *
!*                                                                            *
!******************************************************************************


!******************************************************************************
!* Function: MURGE_ANALYZE                                                    *
!*                                                                            *
!* Perform matrix analyse:                                                    *
!*   - Compute a new ordering of the unknows                                  *
!*   - Compute the symbolic factorisation of the matrix                       *
!*   - Distribute column blocks and computation on processors                 *
!*                                                                            *
!* Parameters:                                                                *
!*   id - Solver instance identification number.                              *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned succesfuly.                     *
!*   MURGE_ERR_ORDER     - If function the graph is not built.                *
!*   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
!*                                                                            *
!******************************************************************************

   SUBROUTINE MURGE_ANALYZE(ID, IERROR)
     INTS, INTENT(IN)  :: ID
     INTS, INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ANALYZE

!******************************************************************************
!* Function: MURGE_FACTORIZE                                                  *
!*                                                                            *
!* Perform matrix factorisation.                                              *
!*                                                                            *
!* Parameters:                                                                *
!*   id - Solver instance identification number.                              *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned succesfuly.                     *
!*   MURGE_ERR_ORDER     - If function the graph is not built.                *
!*   MURGE_ERR_PARAMETER - If *murge_id* is not a valid ID.                   *
!*                                                                            *
!******************************************************************************

   SUBROUTINE MURGE_FACTORIZE(ID, IERROR)
     INTS, INTENT(IN)  :: ID
     INTS, INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_FACTORIZE


!******************************************************************************
!* Function: MURGE_SETORDERING                                                *
!*                                                                            *
!* Set a personal ordering to perform factorization.                          *
!*                                                                            *
!* Parameters:                                                                *
!*   id - Solver instance identification number.                              *
!*   permutation - Permutation to perform factorization.                      *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned succesfuly.                     *
!*   MURGE_ERR_ORDER     - If graph has not been built.                       *
!*   MURGE_ERR_PARAMETER - If permutation is not a valid memory area.         *
!*                                                                            *
!******************************************************************************

  SUBROUTINE MURGE_SETORDERING(ID, PERMUTATION, IERROR)
     INTS,               INTENT(IN)  :: ID
     INTS, DIMENSION(0), INTENT(IN)  :: PERMUTATION
     INTS,               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_SETORDERING

!******************************************************************************
!* Function: MURGE_GETLOCALPRODUCT                                            *
!*                                                                            *
!* Perform the product A * X.                                                 *
!*                                                                            *
!* The vector must have been given trough <MURGE_SetLocalRHS> or              *
!* <MURGE_SetGlobalRHS>.                                                      *
!*                                                                            *
!* Parameters:                                                                *
!*   id - Solver instance identification number.                              *
!*   x  - Array in which the local part of the product will be stored.        *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_ERR_ORDER  - If values have not been set.                          *
!*                                                                            *
!******************************************************************************

  SUBROUTINE MURGE_GETLOCALPRODUCT(ID, X, IERROR)
     INTS,               INTENT(IN)  :: ID
     COEF, DIMENSION(0), INTENT(OUT) :: X
     INTS,               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_GETLOCALPRODUCT

!******************************************************************************
!* Function: MURGE_GETGLOBALPRODUCT                                           *
!*                                                                            *
!* Perform the product A * X.                                                 *
!*                                                                            *
!* The vector must have been given trough <MURGE_SetLocalRHS> or              *
!* <MURGE_SetGlobalRHS>.                                                      *
!*                                                                            *
!* Parameters:                                                                *
!*   id   - Solver instance identification number.                            *
!*   x    - Array in which the product will be stored.                        *
!*   root - Rank of the process which will own the product at end of call,    *
!*          use -1 for all processes.                                         *
!* Returns:                                                                   *
!*   MURGE_ERR_ORDER  - If values have not been set.                          *
!*                                                                            *
!******************************************************************************
  SUBROUTINE MURGE_GETGLOBALPRODUCT(ID, X, root, IERROR)
     INTS,               INTENT(IN)  :: ID
     COEF, DIMENSION(0), INTENT(OUT) :: X
     INTS,               INTENT(IN)  :: root
     INTS,               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_GETGLOBALPRODUCT

!******************************************************************************
!* Function: MURGE_PRODUCTSETLOCALNODENBR                                     *
!*                                                                            *
!*                                                                            *
!*                                                                            *
!* The vector must have been given trough <MURGE_SetLocalRHS> or              *
!* <MURGE_SetGlobalRHS>.                                                      *
!*                                                                            *
!* Parameters:                                                                *
!*   id   - Solver instance identification number.                            *
!*   n    - Number of local nodes.                                            *
!* Returns:                                                                   *
!*   MURGE_ERR_ORDER  - If values have not been set.                          *
!*                                                                            *
!******************************************************************************
  SUBROUTINE MURGE_PRODUCTSETLOCALNODENBR(ID, n, IERROR)
     INTS,               INTENT(IN)  :: ID
     INTS,               INTENT(IN)  :: n
     INTS,               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_PRODUCTSETLOCALNODENBR

!******************************************************************************
!* Function: MURGE_PRODUCTSETGLOBALNODENBR                                    *
!*                                                                            *
!*                                                                            *
!*                                                                            *
!* The vector must have been given trough <MURGE_SetLocalRHS> or              *
!* <MURGE_SetGlobalRHS>.                                                      *
!*                                                                            *
!* Parameters:                                                                *
!*   id   - Solver instance identification number.                            *
!*   N    - Number of global nodes.                                           *
!* Returns:                                                                   *
!*   MURGE_ERR_ORDER  - If values have not been set.                          *
!*                                                                            *
!******************************************************************************
  SUBROUTINE MURGE_PRODUCTSETGLOBALNODENBR(ID, N, IERROR)
     INTS,               INTENT(IN)  :: ID
     INTS,               INTENT(IN)  :: N
     INTS,               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_PRODUCTSETGLOBALNODENBR

!******************************************************************************
!* Function: MURGE_PRODUCTSETLOCALNODELIST                                    *
!*                                                                            *
!*                                                                            *
!*                                                                            *
!* The vector must have been given trough <MURGE_SetLocalRHS> or              *
!* <MURGE_SetGlobalRHS>.                                                      *
!*                                                                            *
!* Parameters:                                                                *
!*   id   - Solver instance identification number.                            *
!*   l2g  - List of local nodes.                                              *
!* Returns:                                                                   *
!*   MURGE_ERR_ORDER  - If values have not been set.                          *
!*                                                                            *
!******************************************************************************

  SUBROUTINE MURGE_PRODUCTSETLOCALNODELIST(ID, l2g, IERROR)
     INTS,               INTENT(IN)  :: ID
     INTS, DIMENSION(0), INTENT(IN)  :: l2g
     INTS,               INTENT(OUT) :: IERROR
  END SUBROUTINE MURGE_PRODUCTSETLOCALNODELIST

!******************************************************************************
!* Function: MURGE_FORCENOFACTO                                               *
!*                                                                            *
!* Prevent Murge from running factorisation even if matrix has changed.       *
!*                                                                            *
!* Parameters:                                                                *
!*   id - Solver instance identification number.                              *
!* Returns:                                                                   *
!*   MURGE_SUCCESS                                                            *
!*                                                                            *
!******************************************************************************!
   SUBROUTINE MURGE_FORCENOFACTO(ID, IERROR)
     INTS,               INTENT(IN)  :: ID
     INTS,               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_FORCENOFACTO

   SUBROUTINE MURGE_SETLOCALNODELIST(ID, N, LIST, IERROR)
     INTS,               INTENT(IN)  :: ID, N
     INTS, DIMENSION(N), INTENT(IN)  :: LIST
     INTS,               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_SETLOCALNODELIST

!******************************************************************************
!* Function: MURGE_AssemblySetSequence                                        *
!*                                                                            *
!* Create a sequence of entries to build a matrix and store it for being      *
!* reused.                                                                    *
!*                                                                            *
!* Parameters:                                                                *
!*   id      - Solver instance identification number.                         *
!*   coefnbr - Number of entries.                                             *
!*   ROWs    - List of rows in the sequence.                                  *
!*   COLs    - List of columns in the sequence.                               *
!*   op      - Operation to perform for coefficient which appear              *
!*             several tim (see <MURGE_ASSEMBLY_OP>).                         *
!*   op2     - Operation to perform when a coefficient is set by              *
!*             two different processors (see <MURGE_ASSEMBLY_OP>).            *
!*   mode    - Indicates if user ensure he will respect solvers distribution  *
!*             (see <MURGE_ASSEMBLY_MODE>).                                   *
!*   nodes   - 0 entries are entered value by value,                          *
!*             1 entries are entries node by node.                            *
!*   id_seq  - Sequence ID.                                                   *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned successfully.                   *
!*   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
!*   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
!*   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
!*                         *op*, *mode*, *sym*, or *coefnbr* are not valid.   *
!******************************************************************************/
   SUBROUTINE MURGE_ASSEMBLYSETSEQUENCE (ID, NNZ, ROWs, COLs, OP, OP2, MODE, NODES, ID_SEQ, IERROR)
     INTS,                     INTENT(IN)  :: ID, OP, OP2, MODE, NODES
     INTL,                     INTENT(IN)  :: NNZ
     INTS, DIMENSION(NNZ),     INTENT(IN)  :: ROWs, COLs
     INTS,                     INTENT(OUT) :: ID_SEQ, IERROR
   END SUBROUTINE MURGE_ASSEMBLYSETSEQUENCE
          

!******************************************************************************
!* MURGE_AssemblySetSequence                                                  *
!*                                                                            *
!* Assembly the matrix using a stored sequence.                               *
!*                                                                            *
!* Parameters:                                                                *
!*   id      - Solver instance identification number.                         *
!*   id_seq  - Sequence ID.                                                   *
!*   values  - Values to insert in the CSC.                                   *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned successfully.                   *
!*   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
!*   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
!*   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
!*                         *id_seq* or *values* are not valid.                *
!******************************************************************************/
   SUBROUTINE MURGE_ASSEMBLYUSESEQUENCE(ID, ID_SEQ, VALUES, IERROR)
     INTS,               INTENT(IN)  :: ID, ID_SEQ 
     COEF, DIMENSION(0), INTENT(IN)  :: VALUES
     INTS,               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ASSEMBLYUSESEQUENCE
      
!******************************************************************************
!* Function: MURGE_AssemblyDeleteSequence                                     *
!*                                                                            *
!* Destroy an assembly sequence                                               *
!*                                                                            *
!*   id      - Solver instance identification number.                         *
!*   id_seq  - Sequence ID.                                                   *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned successfully.                   *
!*   MURGE_ERR_ORDER     - If graph hasn't been built before.                 *
!*   MURGE_ERR_ALLOCATE  - If Allocation didn't worked.                       *
!*   MURGE_ERR_PARAMETER - If *id* is not in solver arrays range, or          *
!*                         *id_seq* is not valid.                             *
!******************************************************************************/
   SUBROUTINE MURGE_ASSEMBLYDELETESEQUENCE(ID, ID_SEQ, IERROR)
     INTS,               INTENT(IN)  :: ID, ID_SEQ 
     INTS,               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_ASSEMBLYDELETESEQUENCE

!******************************************************************************
!* Function: MURGE_SetDropNodes                                               *
!*                                                                            *
!* Set a mask to ignore entries in the graph/matrix.                          *
!*                                                                            *
!* Ignore an entry if dropmask[col] and dropmask[lig]                         *
!*                                                                            *
!*   id       - Solver instance identification number.                        *
!*   nodenbr  - Number of nodes in the global matrix.                         *
!*   dropmask - Mask to indicate if we have to ignore entries.                *
!*                                                                            *
!* Returns:                                                                   *
!*   MURGE_SUCCESS       - If function runned successfully.                   *
!******************************************************************************/
   SUBROUTINE MURGE_SETDROPNODES(ID, NODENBR, DROPMASK, IERROR)
     INTS,               INTENT(IN) :: ID, NODENBR
     INTS, DIMENSION(0), INTENT(IN) :: DROPMASK
     INTS,               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_SETDROPNODES

   SUBROUTINE MURGE_COLGETNONZEROSNBR(ID, ROW, NNZNBR, IERROR)
     INTS,               INTENT(IN) :: ID, ROW
     INTS,               INTENT(INOUT) :: NNZNBR
     INTS,               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_COLGETNONZEROSNBR

   SUBROUTINE MURGE_COLGETNONZEROSIDX(ID, ROW, NNZIDX, IERROR)
     INTS,               INTENT(IN) :: ID, ROW
     INTS, DIMENSION(0), INTENT(INOUT) :: NNZIDX
     INTS,               INTENT(OUT) :: IERROR
   END SUBROUTINE MURGE_COLGETNONZEROSIDX

