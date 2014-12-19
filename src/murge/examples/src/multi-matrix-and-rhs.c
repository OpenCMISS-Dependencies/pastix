/*
  Program: murge-matrix-and-rhs

  Example which build several distributed laplacian and solves 
  the system using *Murge* interface.
  
  This example is independant of the Murged solver used.
 
  The solved problem are very simple :

  >  2 -1  0  0 ... 0   x_0   1
  > -1  2 -1  0 ... 0   x_1   0
  >   ....            x .   = .
  >  0 ...    2 -1  0   .     .
  >  0 ...   -1  2 -1   .     0
  >  0 ...    0 -1  2   x_n   1

  >  4 -1  0  0 ... 0   x_0   6
  > -1  4 -1  0 ... 0   x_1   4
  >   ....            x .   = .
  >  0 ...    4 -1  0   .     .
  >  0 ...   -1  4 -1   .     4
  >  0 ...    0 -1  4   x_n   6
  
  etc..

  Authors:
     M. FAVERGE, J. GAIDAMOUR, P. HENON , X. LACOSTE
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#ifdef TYPE_COMPLEX
#include <complex.h>
#endif /* TYPE_COMPLEX */
#include "murge.h"

/* Section: Macros */

/* 
   Macro: MIN
   
   Compute the minimum of two values.

   Parameters:
     x - The first value.
     y - The second value.
*/
#define MIN(x,y) (((x)<(y))?(x):(y))
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
      {									\
	fprintf(stderr, "%s:%d error in %s : %d\n",		\
		__FILE__, __LINE__,#call, (int)ret);			\
	return 1;							\
      }									\
  }									

/* Section: Defines */

/* Define: DIAG_VALUE       
   The value of the diagonal elements of the matrix
   Define: EXTRA_DIAG_VALUE 
   The value of the extra diagonal elements of the matrix.
*/
#define DIAG_VALUE         2
#define EXTRA_DIAG_VALUE  -1
#define MATR_NBR           3
#define RHS_NBR            3
/* Section: Functions */

/*
  Function: main

  Main function of <murge-dist> example, solving laplacian.

  Initialize MURGE and default solver options,

  Enter the graph in ditributed mode, 
  edge by edge using <MURGE_GraphBegin>, <MURGE_GraphEdge>
  and <MURGE_GraphEnd>.

  Get the new computed distribution using <MURGE_GetLocalUnknownNbr> 
  and <MURGE_GetLocalUnknownList>.

  Fill the ditributed matrix, value by value using <MURGE_AssemblyBegin>,
  <MURGE_AssemblySetValue> and <MURGE_AssemblyEnd>.

  Set the right hand side member using <MURGE_SetLocalRHS>.
  
  Solve the problem and get the local solution using <MURGE_GetLocalSolution>.
  
*/
int main(int argc, char *argv[])
{  
  COEF   *rhs      = NULL;    /* Local right-hand-side member                */
  COEF   *xx       = NULL;    /* Local solution                              */
  INTS   *loc2glob = NULL;    /* Local to global column numbers              */
  INTS    globalN  = 1000;    /* Size  of the problem                        */
  INTL    lnnz;               /* Local number of non zeros                   */
  INTS    localn;             /* Local number of unkowns                     */
  INTL    start;              /* First local column                          */
  INTL    end;                /* Last local column                           */
  INTL    i, matr, irhs;      /* iterator                                    */
  int     proc_id;            /* MPI process ID                              */
  int     nproc;              /* Number of MPI processes                     */
  INTS    id;                 /* Solver instance id                          */
  int     baseval  = 0;       /* base 0 or 1, used to define the matrix      */
  INTS    sym;                /* indicate if the matrix pattern is symmetric */
  int     required;           /* MPI thread level required                   */
  int     provided;           /* MPI thread level provided                   */
  double norme1, norme1_recv;
  double norme2, norme2_recv;
    
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
  /* Read the matrix from file      */
  /**********************************/
  start  = proc_id*((globalN + 1)/nproc);
  end    = MIN(((proc_id+1)*((globalN + 1)/nproc))-1, globalN-1);
  if (proc_id == nproc-1)
    end = globalN -1;
  localn = end - start + 1;

  MURGE_CALL(MURGE_SetOptionINT(id, MURGE_IPARAM_BASEVAL, (INTS)baseval));

  MURGE_CALL(MURGE_SetOptionINT(id, MURGE_IPARAM_SYM, sym));
  if (sym == MURGE_BOOLEAN_TRUE)
    lnnz = 2*localn;
  else
    lnnz = 3*localn;

  if (sym == MURGE_BOOLEAN_FALSE && start == 0)
    lnnz --;
  if (end == globalN-1)
    lnnz --;
  /***************************************************/
  /* ENTER THE GRAPH : PARALLEL INTERFACE            */
  /* Every processor deals with a part of the graph  */
  /***************************************************/
  fprintf(stdout, "Processor %ld enter edges [%ld %ld] of the graph \n", (long)proc_id, (long)start, (long)end);

  /***************************************************/
  /*                                                 */
  /* ENTER THE MATRIX GRAPH : Distributed INTERFACE  */
  /*                                                 */
  /***************************************************/
  MURGE_CALL(MURGE_GraphBegin(id, globalN, lnnz));

  for (i = start; i < end+1; i++)
    {
      if (sym == 0 && i != 0) 
	MURGE_CALL(MURGE_GraphEdge(id, i + baseval, i - 1 + baseval));

      MURGE_CALL(MURGE_GraphEdge(id, i + baseval, i + baseval));

      if (i != globalN - 1)
	MURGE_CALL(MURGE_GraphEdge(id, i + baseval, i + 1 + baseval));
    }
  MURGE_CALL(MURGE_GraphEnd(id));

  /***************************************************/
  /*                                                 */
  /*           RETRIEVE LOCAL UNKNOWN LIST           */
  /***************************************************/

  MURGE_CALL(MURGE_GetLocalUnknownNbr (id, &localn));

  if (NULL == (loc2glob = (INTS*)malloc(localn*sizeof(INTS))))
    {
      fprintf(stderr, "Error: couldn't allocate\n");
      return EXIT_FAILURE;
    }
  MURGE_CALL(MURGE_GetLocalUnknownList(id, loc2glob));

  if (sym == 1)
    lnnz = 2*localn;
  else
    lnnz = 3*localn;
  for (i = 0; i < localn; i++)
    if (( sym == 0 && loc2glob[i] == baseval) ||
	loc2glob[i] == globalN - 1 + baseval)
      lnnz --;

  for (matr = 1; matr < MATR_NBR+1; matr++)
    {
      double extra_diag_value = EXTRA_DIAG_VALUE;
      double diag_value = matr*DIAG_VALUE;
      /***************************************************/
      /*                                                 */
      /* ENTER THE MATRIX COEFS : Distributed INTERFACE  */
      /***************************************************/
      MURGE_CALL(MURGE_AssemblyBegin(id, globalN, lnnz, 
				     MURGE_ASSEMBLY_OVW, 
				     MURGE_ASSEMBLY_OVW, 
				     MURGE_ASSEMBLY_FOOL, sym)); 
      
      for (i = 0; i < localn; i++)
	{
	  if ((sym == 0) && (loc2glob[i] != baseval)) 
	    MURGE_CALL(MURGE_AssemblySetValue(id, loc2glob[i], 
					      loc2glob[i] - 1, 
					      extra_diag_value )); 
      
	  MURGE_CALL(MURGE_AssemblySetValue(id, loc2glob[i], 
					    loc2glob[i], diag_value)); 
	  
	  if (loc2glob[i] != globalN + baseval - 1)
	    MURGE_CALL(MURGE_AssemblySetValue(id, loc2glob[i], 
					      loc2glob[i] + 1, 
					      extra_diag_value)); 
      
	}
      MURGE_CALL(MURGE_AssemblyEnd(id));
     
      for (irhs = 1; irhs < RHS_NBR+1; irhs++)
	{
	  rhs = (COEF*) malloc(localn*sizeof(COEF));
	  for (i=0;i<localn;i++)
	    {
	      if (loc2glob[i] == baseval || loc2glob[i] == globalN - 1 + baseval) 
		{
		  rhs[i] = irhs*(diag_value + extra_diag_value);
		}
	      else
		{
		  rhs[i] = irhs*(diag_value + 2 * extra_diag_value);
		}
	    }
      
	  /****************************************************/
	  /* Set the global rhs                               */
	  /* rhs is only significant on the master processor  */
	  /****************************************************/ 
	  MURGE_CALL(MURGE_SetLocalRHS(id, rhs, MURGE_ASSEMBLY_OVW, MURGE_ASSEMBLY_OVW));
	  
	  /****************************************************/
	  /* Get the global solution on processor 0           */
	  /* Original ordering                                */
	  /****************************************************/ 
	  xx = (COEF *)malloc(sizeof(COEF)*localn);
	  MURGE_CALL(MURGE_GetLocalSolution(id, xx));
	  
	  norme1= norme2 = 0;
	  for (i= 0; i < localn; i++)
	    {
	      /*fprintf(stdout, "xx[%ld] = %lg (expected %lg)\n", (long)loc2glob[i], (double)xx[i], (double)irhs);*/
	      norme1 += (double)((xx[i] - irhs)*(xx[i] - irhs));
	      norme2 += (double)(irhs * irhs);
	    } 
	  MPI_Reduce(&norme1, &norme1_recv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  MPI_Reduce(&norme2, &norme2_recv, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	  if (proc_id == 0)
	    {
	      fprintf(stdout, "Precision : ||solution - expected||/||expected|| = %.20lg\n", sqrt(norme1_recv/norme2_recv)); 
	      
	    }
	}
      free(xx);
      free(rhs);
    }
  /***************************************************/
  /* Free Solver internal structures for problem id  */
  /***************************************************/
  MURGE_CALL(MURGE_Clean(id));
  
  /************************************/
  /* Free Solver internal structures  */
  /************************************/
  MURGE_CALL(MURGE_Finalize());

  free(loc2glob);
  loc2glob = NULL;
  
  
  /** End MPI **/
  MPI_Finalize();
  
  return 0;
}
