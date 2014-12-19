/*
 * Program: murge-mm
 *
 * Example which reads a file "matrix.mm" and solves
 * the system using *Murge* interface.
 *
 * The right-hand-side member used is such as the solution will be 1.
 *
 * The given matrix must be in a IJV Matrix Market format, with 1 as baseval.
 *
 * This example is independant of the Murged solver used.
 *
 * Authors:
 *    X. LACOSTE
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#ifdef TYPE_COMPLEX
#include <complex.h>
#endif /* TYPE_COMPLEX */
#include "murge.h"
#include "mmio.h"

/* Section: Macros */
/*
 * Macro: MURGE_CALL
 *
 * Execute *call* fonction and test the return value.
 *
 * Parameters:
 *   call - The murge function to execute.
 */
#define MURGE_CALL(call)                                    \
  {                                                         \
    INTS ret;                                               \
    ret = call;                                             \
    if (ret != MURGE_SUCCESS)                               \
      {                                                     \
        fprintf(stderr, "%s:%d error in murge call : %d\n", \
                __FILE__, __LINE__, (int)ret);              \
        exit(1);                                            \
      }                                                     \
  }


/* Section: Functions */

/*
 * Function: main
 *
 * Main function of <murge-dist> example, solving laplacian.
 *
 * Initialize MURGE and default solver options,
 *
 * Enter the graph in ditributed mode,
 * edge by edge using <MURGE_GraphBegin>, <MURGE_GraphEdge>
 * and <MURGE_GraphEnd>.
 *
 * Get the new computed distribution using <MURGE_GetLocalUnknownNbr>
 * and <MURGE_GetLocalUnknownList>.
 *
 * Fill the ditributed matrix, value by value using <MURGE_AssemblyBegin>,
 * <MURGE_AssemblySetValue> and <MURGE_AssemblyEnd>.
 *
 * Set the right hand side member using <MURGE_SetLocalRHS>.
 *
 * Solve the problem and get the local solution using <MURGE_GetLocalSolution>.
 *
 */
int main(int argc, char *argv[])
{
  COEF   *rhs      = NULL;    /* Local right-hand-side member                */
  COEF   *xx       = NULL;    /* Local solution                              */
  INTS    globalN;            /* Size  of the problem                        */
  INTL    lnnz     = 0;       /* Local number of non zeros                   */
  INTS    nrow;               /* Local number rows                           */
  INTL    i;                  /* iterator                                    */
  int     proc_id;            /* MPI process ID                              */
  int     nproc;              /* Number of MPI processes                     */
  INTS    id;                 /* Solver instance id                          */
  int     baseval  = 1;       /* base 0 or 1, used to define the matrix      */
  INTS    sym;                /* indicate if the matrix pattern is symmetric */
  int     required;           /* MPI thread level required                   */
  int     provided;           /* MPI thread level provided                   */
  FILE   *file     = NULL;
  double  re, im;
  long    row, col;
  int     tmp1, tmp2, tmp3;
  MM_typecode matcode;

  /** Init MPI environment **/
  required = MPI_THREAD_MULTIPLE;
  provided = -1;
  MPI_Init_thread(&argc, &argv, required, &provided);
  if (provided != required)
    {
      switch (provided)
        {
        case MPI_THREAD_SINGLE:
          printf("MPI_Init_thread level = MPI_THREAD_SINGLE\n");
          break;
        case MPI_THREAD_FUNNELED:
          printf("MPI_Init_thread level = MPI_THREAD_FUNNELED\n");
          break;
        case MPI_THREAD_SERIALIZED:
          printf("MPI_Init_thread level = MPI_THREAD_SERIALIZED\n");
          break;
        case MPI_THREAD_MULTIPLE:
          printf("MPI_Init_thread level = MPI_THREAD_MULTIPLE\n");
          break;
        default:
          printf("MPI_Init_thread level = ???\n");
        }
    }
  MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);


  if (proc_id == 0)
    {
      file = fopen ("matrix.mm","r");
      if (file==NULL)
        {
          fprintf(stderr,"cannot load matrix.mm\n");
          exit(1);
        }

      if (mm_read_banner(file, &matcode) != 0)
        {
          fprintf(stderr,"Could not process Matrix Market banner.\n");
          exit(1);
        }

      sym = MURGE_BOOLEAN_FALSE;
      if (mm_is_symmetric(matcode))
        sym = MURGE_BOOLEAN_TRUE;

      /* find out size of sparse matrix .... */

      if (mm_read_mtx_crd_size(file, &tmp1, &tmp2, &tmp3) !=0)
        exit(1);
      nrow    = tmp1;
      globalN = tmp2;
      lnnz    = tmp3;

    }

  MPI_Bcast(&sym, sizeof(INTS), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&nrow, sizeof(INTS), MPI_BYTE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&globalN, sizeof(INTS), MPI_BYTE, 0, MPI_COMM_WORLD);
  /***************************************/
  /* Initialize MURGE for 1 problem      */
  /***************************************/
  MURGE_CALL(MURGE_Initialize(1));

  id = 0; /** id of the linear system **/

  /***************************************/
  /* Initialize Default solver options   */
  /***************************************/
  MURGE_CALL(MURGE_SetDefaultOptions(id, 0));

  /***************************************/
  /* Set the communicator                */
  /***************************************/
  MURGE_CALL(MURGE_SetCommunicator(id, MPI_COMM_WORLD));

  /**********************************/
  /* Read the matrix from file      */
  /**********************************/

  MURGE_CALL(MURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL, (INTS)baseval));

  MURGE_CALL(MURGE_SetOptionINT(id, MURGE_IPARAM_SYM, sym));

  /***************************************************/
  /* ENTER THE GRAPH : PARALLEL INTERFACE            */
  /* Every processor deals with a part of the graph  */
  /***************************************************/

  /***************************************************/
  /*                                                 */
  /* ENTER THE MATRIX GRAPH : Distributed INTERFACE  */
  /*                                                 */
  /***************************************************/
  if (proc_id == 0)
    {
      MURGE_CALL(MURGE_GraphBegin(id, globalN, lnnz));
      if (mm_is_complex(matcode))
        {
#ifndef TYPE_COMPLEX
          fprintf(stderr,
                  "WARNING : Matrix should not be complex."
                  " Only real part will be taken.\n");
#endif
          for (i=0; i<(lnnz); i++)
            {
              fscanf(file,"%ld %ld %lg %lg\n", &row, &col, &re, &im);
              MURGE_CALL(MURGE_GraphEdge(id, (INTS)row, (INTS)col));
            }
        }
      else
        {
          for (i=0; i<lnnz; i++)
            {
              fscanf(file,"%ld %ld %lg\n", &row, &col, &re);
              MURGE_CALL(MURGE_GraphEdge(id, (INTS)row, (INTS)col));
            }
        }
      MURGE_CALL(MURGE_GraphEnd(id));
      fclose(file);
    }
  else
    {
      MURGE_CALL(MURGE_GraphBegin(id, globalN, 0));
      MURGE_CALL(MURGE_GraphEnd(id));
    }

  /***************************************************/
  /*                                                 */
  /* ENTER THE MATRIX COEFS : Distributed INTERFACE  */
  /***************************************************/
  rhs = (COEF*) malloc(globalN*sizeof(COEF));
  memset(rhs, 0, globalN*sizeof(COEF));
  im = 0;

  if (proc_id == 0)
    {
      file = fopen ("matrix.mm","r");
      if (file==NULL)
        {
          fprintf(stderr,"cannot load matrix.mm\n");
          exit(1);
        }

      if (mm_read_banner(file, &matcode) != 0)
        {
          fprintf(stderr,"Could not process Matrix Market banner.\n");
          exit(1);
        }

      if (mm_read_mtx_crd_size(file, &tmp1, &tmp2, &tmp3) !=0)
        exit(1);

      MURGE_CALL(MURGE_AssemblyBegin(id, globalN, lnnz,
                                     MURGE_ASSEMBLY_OVW,
                                     MURGE_ASSEMBLY_OVW,
                                     MURGE_ASSEMBLY_FOOL, sym));
      if (mm_is_complex(matcode))
        {
          for (i=0; i<(lnnz); i++)
            {
              fscanf(file,"%ld %ld %lg %lg\n", &row, &col, &re, &im);

#ifdef TYPE_COMPLEX
              MURGE_CALL(MURGE_AssemblySetValue(id,
                                                (INTS)row, (INTS)col, re+I*im));
              rhs[row-baseval] += re+I*im;
#else
              MURGE_CALL(MURGE_AssemblySetValue(id, (INTS)row, (INTS)col, re));
              rhs[row-baseval] += re;
#endif
            }
        }
      else
        {
          for (i=0; i<lnnz; i++)
            {
              fscanf(file,"%ld %ld %lg\n", &row, &col, &re);
#ifdef TYPE_COMPLEX
              MURGE_CALL(MURGE_AssemblySetValue(id,
                                                (INTS)row, (INTS)col, re+I*im));
              rhs[row-baseval] += re+I*im;
#else
              MURGE_CALL(MURGE_AssemblySetValue(id, (INTS)row, (INTS)col, re));
              rhs[row-baseval] += re;
#endif
            }
        }
      fclose(file);
      MURGE_CALL(MURGE_AssemblyEnd(id));

    }
  else
    {
      MURGE_CALL(MURGE_AssemblyBegin(id, globalN, 0, 
                                     MURGE_ASSEMBLY_OVW,
                                     MURGE_ASSEMBLY_OVW,
                                     MURGE_ASSEMBLY_FOOL, sym));
      MURGE_CALL(MURGE_AssemblyEnd(id));
    }

  /****************************************************/
  /* Set the global rhs                               */
  /* rhs is only significant on the master processor  */
  /****************************************************/
  MURGE_CALL(MURGE_SetGlobalRHS(id, rhs, 0, MURGE_ASSEMBLY_OVW));

  /****************************************************/
  /* Get the global solution on processor 0           */
  /* Original ordering                                */
  /****************************************************/
  xx = (COEF *)malloc(sizeof(COEF)*globalN);
  MURGE_CALL(MURGE_GetLocalSolution(id, xx));

  /***************************************************/
  /* Free Solver internal structures for problem id  */
  /***************************************************/
  MURGE_CALL(MURGE_Clean(id));

  /************************************/
  /* Free Solver internal structures  */
  /************************************/
  MURGE_CALL(MURGE_Finalize());

  free(xx);
  free(rhs);

  /** End MPI **/
  MPI_Finalize();

  return 0;
}
