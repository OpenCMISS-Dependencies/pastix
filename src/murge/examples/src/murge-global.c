/*
  Program: murge-global

  Example which build a global laplacian and solves 
  the system using *Murge* interface.
  
  This example is independant of the Murged solver used.
 
  The solved problem is very simple :

  >  2 -1  0  0 ... 0   1   1
  > -1  2 -1  0 ... 0   0   1
  >   ....            x . = .
  >  0 ...    2 -1  0   .   .
  >  0 ...   -1  2 -1   0   1
  >  0 ...    0 -1  2   1   1
  
  Authors:
     M. FAVERGE, J. GAIDAMOUR, P. HENON , X. LACOSTE
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#ifdef TYPE_COMPLEX
#include <complex.h>
#endif /* TYPE_COMPLEX */
#include "murge.h"

/* Section: Macros */

/*
  Macro: MURGE_CALL
  
  Execute *call* fonction and test the return value.
  
  Parameters:
    call - The murge function to execute.
*/
#define MURGE_CALL(call)						\
  {									\
    INTS ret;								\
    ret = call;								\
    if (ret != MURGE_SUCCESS)						\
    fprintf(stderr, "%s:%d error in murge call : %d\n",			\
	    __FILE__, __LINE__, (int)ret);				\
  }									

/* Section: Defines */

/* Define: DIAG_VALUE       
   The value of the diagonal elements of the matrix
*/
#define DIAG_VALUE         2

/* Define: EXTRA_DIAG_VALUE 
   The value of the extra diagonal elements of the matrix.
*/
#define EXTRA_DIAG_VALUE  -1

/* Section: Functions */

/*
  Function: main

  Main function of <murge-dist> example, solving laplacian.

  Initialize MURGE and default solver options,

  Enter the graph in centralised mode,
  edge by edge using <MURGE_GraphGlobalIJV>

  Fill the ditributed matrix, in one <MURGE_MatrixGlobalIJV> call.

  Set the right hand side member using <MURGE_SetGlobalRHS>.
  
  Solve the problem and get the local solution using <MURGE_GetGlobalSolution>.
  
*/
int main(int argc, char *argv[])
{  
  INTS *rows    = NULL; /* Global list of rows                         */
  INTS *cols    = NULL;	/* Global list of columns                      */
  COEF *values  = NULL;	/* Global list of values                       */
  COEF *xx      = NULL;	/* Global solution                             */
  COEF *rhs     = NULL;	/* Global right-hand-side member               */
  INTL  nnz;		/* Gloabl number of non zeros                  */
  INTS  globalN = 1000;	/* Global number of unknowns                   */
  INTS  sym;	/* 1 if the pattern of the matrix is symmetric */
  int   nproc;		/* Number of MPI processes                     */
  int   proc_id;	/* MPI process ID                              */
  INTS  id;		/* Solver instance ID                          */
  INTS  index;		/* iterator                                    */
  INTS  i;		/* iterator                                    */
  int   baseval = 1;    /* Numerotation starting index                 */
  int   required;           /* MPI thread level required                   */
  int   provided;           /* MPI thread level provided                   */

  sym = MURGE_BOOLEAN_TRUE;

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
  /* Build the graph of the matrix  */
  /**********************************/
  if (proc_id == 0)
    {

      if (sym == 1)
	nnz = 2*globalN - 1;
      else
	nnz = 3*globalN - 2;

	rows = malloc(nnz*sizeof(INTS));
	cols = malloc(nnz*sizeof(INTS));

	index = 0;
	for (i = 0; i < globalN; i++)
	  {
	    if (sym == 0 && i != 0) 
	      { 
		rows[index] = i + baseval;
		cols[index] = i - 1 + baseval;
		index ++;
	      }
	    rows[index] = i + baseval;
	    cols[index] = i + baseval;
	    index ++;
	    if (i != globalN - 1)
	      {
		rows[index] = i + baseval;
		cols[index] = i + 1 + baseval;
		index ++;
	      }
	  }
    }

  MURGE_CALL(MURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL, (INTS)baseval));

  /***************************************************/
  /*                                                 */
  /* ENTER THE MATRIX GRAPH : CENTRALIZED INTERFACE  */
  /*                                                 */
  /***************************************************/
  MURGE_CALL(MURGE_GraphGlobalIJV(id, globalN, nnz, rows, cols, 0));

  /***************************************************/
  /*                                                 */
  /* ENTER THE MATRIX GOEFS : CENTRALIZED INTERFACE  */
  /*                                                 */
  /***************************************************/
  if (proc_id == 0)
    {
	values = malloc(nnz*sizeof(COEF));
	index = 0;
	for (i = 0; i < globalN; i++)
	  {
	    if (sym == 0 && i != 0) 
	      { 
		values[index] = EXTRA_DIAG_VALUE;
		index ++;
	      }
	    values[index] = DIAG_VALUE;
	    index ++;
	    if (i != globalN - 1)
	      {
		values[index] = EXTRA_DIAG_VALUE;
		index ++;
	      }
	  }
    }
  MURGE_CALL(MURGE_MatrixGlobalIJV(id, globalN, nnz, rows, cols, values, 
				   0, MURGE_ASSEMBLY_OVW, sym));

  if(proc_id == 0)
    {
      free(rows);
      free(cols);
      free(values);
    }
     
  /****************************************/
  /* Set the right hand side (from proc 0)*/
  /****************************************/
  if(proc_id == 0)
    {
      rhs = (COEF *)malloc(globalN*sizeof(COEF));
      rhs[0]         = DIAG_VALUE - EXTRA_DIAG_VALUE;
      for (i = 1; i < globalN -1; i++)
	rhs[i] = DIAG_VALUE - 2 * EXTRA_DIAG_VALUE;
      rhs[globalN-1] = DIAG_VALUE - EXTRA_DIAG_VALUE;
    }
  else
    rhs = NULL;
  
 
  /****************************************************/
  /* Set the global rhs                               */
  /* rhs is only significant on the master processor  */
  /****************************************************/ 
  MURGE_CALL(MURGE_SetGlobalRHS(id, rhs, 0, MURGE_ASSEMBLY_OVW));

  /****************************************************/
  /* Get the global solution on processor 0           */
  /* Original ordering                                */
  /****************************************************/ 
  if(proc_id == 0)
    xx = (COEF *)malloc(sizeof(COEF)*globalN);
  
  MURGE_CALL(MURGE_GetGlobalSolution(id, xx, 0));
 
  /***************************************************/
  /* Free Solver internal structures for problem id  */
  /***************************************************/
  MURGE_CALL(MURGE_Clean(id));

  /************************************/
  /* Free Solver internal structures  */
  /************************************/
  MURGE_CALL(MURGE_Finalize());

  if(proc_id == 0)
    {
      free(xx);
      free(rhs);
    }

  /** End MPI **/
  MPI_Finalize();
  
  return 0;
}
