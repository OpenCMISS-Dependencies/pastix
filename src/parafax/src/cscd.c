/* Copyright 2008 BORDEAUX I UNIVERSITY & INRIA 
**
** This file is part of the PaStiX parallel sparse matrix package.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
#include <sys/time.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef FORCE_NOMPI
#include "nompi.h"
#else
#include <mpi.h>
#endif

#include "common_pastix.h"
#ifdef WITH_SCOTCH
#ifdef    DISTRIBUTED
#include "ptscotch.h"
#else  /* DISTRIBUTED */
#include "scotch.h"
#endif /* DISTRIBUTED */
#endif
#include "cscd.h"


/* 
   Concentre la CSCd sur le noeud maitre 0 
   Not used, not sure it works, it shouldn't 
   works if proc 0 doesn't own frist columns... (XL)
*/
int cscd_Fusion(CSCD * cscd_graph) 
{ 
  PASTIX_FLOAT     *buffloat = NULL;  
  PASTIX_INT       *bufint   = NULL; 
  PASTIX_INT       *vect, *vectrcv;
  MPI_Status status;
  int        myrank, nbproc;
  PASTIX_INT        i, j, edgeglb=0, vertglb=0;

  MPI_Comm_rank( (cscd_graph)->mpi_comm, &myrank );
  MPI_Comm_size( (cscd_graph)->mpi_comm, &nbproc );

  if (!(vect    = (PASTIX_INT *)memAlloc(sizeof(PASTIX_INT)*2*(nbproc+1))))
    MALLOC_ERROR("vect");
  if (!(vectrcv = (PASTIX_INT *)memAlloc(sizeof(PASTIX_INT)*3*(nbproc+1))))
    MALLOC_ERROR("vectrecv"); 

  vect    = memSet(vect,   0,sizeof(PASTIX_INT)*2*(nbproc+1));
  vectrcv = memSet(vectrcv,0,sizeof(PASTIX_INT)*3*(nbproc+1));

  vect[myrank]          = CSCD_VERTLOC(cscd_graph);
  vect[nbproc + myrank] = CSCD_EDGELOC(cscd_graph);
  vect[2*nbproc]        = CSCD_VERTLOC(cscd_graph);
  vect[2*nbproc+1]      = CSCD_EDGELOC(cscd_graph);

  MPI_Allreduce(vect, vectrcv, 2*(nbproc+1), COMM_INT, MPI_MAX, (cscd_graph)->mpi_comm);

  if (myrank == 0)
    {
      PASTIX_INT   *ia;
      PASTIX_INT   *loc2glb;
      PASTIX_INT   *ja;
      PASTIX_FLOAT *rhs;
      PASTIX_FLOAT *aval;

      for(i=0; i<nbproc; i++) 
	{
	  edgeglb += vectrcv[i+nbproc]; 
	  vertglb += vectrcv[i];
	  /* Decalage de 2 valeurs a cause du max sommet et du max arete */
	  vectrcv[2*nbproc+i+2] = ((i==0) ? 0 : vectrcv[2*nbproc+i+1] + vectrcv[i-1]);
	}

      ia      = (PASTIX_INT *)  memAlloc(sizeof(PASTIX_INT)  *vertglb);
      loc2glb = (PASTIX_INT *)  memAlloc(sizeof(PASTIX_INT)  *vertglb);
      ja      = (PASTIX_INT *)  memAlloc(sizeof(PASTIX_INT)  *edgeglb);
      rhs     = (PASTIX_FLOAT *)memAlloc(sizeof(PASTIX_FLOAT)*vertglb);
      aval    = (PASTIX_FLOAT *)memAlloc(sizeof(PASTIX_FLOAT)*edgeglb);
      
      bufint   = memAlloc(sizeof(PASTIX_INT)  *(vectrcv[2*nbproc+1]+1));
      buffloat = memAlloc(sizeof(PASTIX_FLOAT)* vectrcv[2*nbproc+1]);
      
      for(i=0; i<nbproc; i++) 
	{
	  if (i==0) 
	    {
	      for (j=0; j<vectrcv[0]; j++) 
		{
		  ia[j]  = CSCD_IA(cscd_graph)[j];
		  rhs[j] = CSCD_RHS(cscd_graph)[j];
		}
 	      ia[vectrcv[0]] = CSCD_IA(cscd_graph)[vectrcv[0]];
	      for (j=0; j<vectrcv[nbproc]; j++) 
		{
		  ja[j]   = CSCD_JA(cscd_graph)[j];
		  aval[j] = CSCD_AVAL(cscd_graph)[j];
		}
	      for (j=0; j<vertglb; j++)
		loc2glb[j] = j+1;
	    }
	  else 
	    {
	      MPI_Recv(bufint, vectrcv[i]+1, COMM_INT, i, 1, cscd_graph->mpi_comm, &status);
	      for (j=0; j<=vectrcv[i]; j++)
		ia[j+vectrcv[2*nbproc+i+2]] = bufint[j] + ia[vectrcv[2*nbproc+i+2]]-1;
	      
	      MPI_Recv(bufint, vectrcv[nbproc+i], COMM_INT, i, 1, cscd_graph->mpi_comm, &status);
	      for(j=0;j<vectrcv[nbproc+i];j++)
		ja[ia[vectrcv[2*nbproc+i+2]]+j-1] = bufint[j];
	      
	      MPI_Recv(buffloat, vectrcv[i], COMM_FLOAT, i, 1, cscd_graph->mpi_comm, &status);
	      for (j=0; j<vectrcv[i]; j++)
		rhs[j+vectrcv[2*nbproc+i+2]] = buffloat[j];
	      
	      MPI_Recv(buffloat, vectrcv[nbproc+i], COMM_FLOAT, i, 1, cscd_graph->mpi_comm, &status);
	      for(j=0;j<vectrcv[nbproc+i];j++)
		aval[ia[vectrcv[2*nbproc+i+2]]+j-1] = buffloat[j];
	    }
	}
      cscd_Init(cscd_graph, loc2glb, vertglb, edgeglb, ia, ja, aval, rhs, cscd_graph->mpi_comm);
    }
  else 
    {
      MPI_Send(CSCD_IA(cscd_graph),  CSCD_VERTLOC(cscd_graph)+1, COMM_INT,   0, 1, cscd_graph->mpi_comm);
      MPI_Send(CSCD_JA(cscd_graph),  CSCD_EDGELOC(cscd_graph),   COMM_INT,   0, 1, cscd_graph->mpi_comm);
      MPI_Send(CSCD_RHS(cscd_graph), CSCD_VERTLOC(cscd_graph),   COMM_FLOAT, 0, 1, cscd_graph->mpi_comm);
      MPI_Send(CSCD_AVAL(cscd_graph),CSCD_EDGELOC(cscd_graph),   COMM_FLOAT, 0, 1, cscd_graph->mpi_comm);
    }
  
  if (myrank != 0) 
    {
      CSCD_VERTLOC(cscd_graph) = 0;
      CSCD_EDGELOC(cscd_graph) = 0;
    }

  if (buffloat) memFree_null(buffloat);
  if (bufint)   memFree_null(bufint);
  if (vectrcv)  memFree_null(vectrcv);
  if (vect)     memFree_null(vect);
  return 0;
}

/* Chaque processus appelle cette fonction en passant la cscd a remplir avec le bon */
/* communicateur, 0 doit posseder la CSC complete centralisee */
/* Not used, works ? (XL) */
int cscd_Explode(CSCD * cscd_graph, MPI_Comm mpi_comm) 
{
  MPI_Status status;
  PASTIX_INT       *ia;
  PASTIX_INT       *ja;   
  PASTIX_FLOAT     *aval;
  PASTIX_FLOAT     *rhs;
  PASTIX_INT       *loc2glb;
  PASTIX_INT       *nbvertedge;
  int        myrank, nbproc;
  PASTIX_INT        i;
  PASTIX_INT        edgeloc, vertloc;

  MPI_Comm_rank( mpi_comm, &myrank);
  MPI_Comm_size( mpi_comm, &nbproc);

  if (!(nbvertedge = (PASTIX_INT *)memAlloc(3*nbproc*sizeof(PASTIX_INT))))
    MALLOC_ERROR("nbvertedge"); 

  if (myrank == 0) 
    {
      PASTIX_INT moyenne;
      PASTIX_INT sum = 0;

      if (nbproc > CSCD_VERTLOC(cscd_graph)) 
	{
	  errorPrint("Explode, %d sommets, %d procs", (int)CSCD_VERTLOC(cscd_graph), (int)nbproc);
	  return EXIT_FAILURE;
	}

      moyenne = (PASTIX_INT) CSCD_VERTLOC(cscd_graph)/nbproc;
      for (i=0; i<nbproc; i++) 
	{
	  if (i<CSCD_VERTLOC(cscd_graph)%nbproc)
	    nbvertedge[i] = moyenne + 1;
	  else 
	    nbvertedge[i] = moyenne;
	}
      
      for (i=0; i<nbproc; i++) 
	{
	  nbvertedge[nbproc+i]   = CSCD_IA(cscd_graph)[sum+nbvertedge[i]] - CSCD_IA(cscd_graph)[sum];
	  nbvertedge[2*nbproc+i] = sum;
	  sum += nbvertedge[i];
	}
    }
  
  MPI_Bcast( &nbvertedge[0], 3*nbproc, COMM_INT, 0, mpi_comm );

  vertloc = nbvertedge[myrank]; 
  edgeloc = nbvertedge[nbproc+myrank];

  ia      = (PASTIX_INT *)  memAlloc(sizeof(PASTIX_INT)  *(vertloc+1));
  ja      = (PASTIX_INT *)  memAlloc(sizeof(PASTIX_INT)  *edgeloc);
  aval    = (PASTIX_FLOAT *)memAlloc(sizeof(PASTIX_FLOAT)*edgeloc);
  rhs     = (PASTIX_FLOAT *)memAlloc(sizeof(PASTIX_FLOAT)*vertloc);
  loc2glb = (PASTIX_INT *)  memAlloc(sizeof(PASTIX_INT)  *vertloc);

  for(i=0; i<vertloc; i++)
    loc2glb[i] = nbvertedge[2*nbproc+myrank] + i+1;
  
  if (myrank == 0) 
    {
      for (i=1; i<nbproc; i++) 
	{
	  MPI_Send(&(CSCD_IA(cscd_graph)[nbvertedge[2*nbproc+i]]), 
		   nbvertedge[i],        COMM_INT,    i, 1, mpi_comm);
      
	  MPI_Send(&(CSCD_JA(cscd_graph)[CSCD_IA(cscd_graph)[nbvertedge[2*nbproc+i]]-1]), 
		   nbvertedge[nbproc+i], COMM_INT,    i, 1, mpi_comm);

	  MPI_Send(&(CSCD_RHS(cscd_graph)[nbvertedge[2*nbproc+i]]), 
		   nbvertedge[i],        COMM_FLOAT, i, 1, mpi_comm);
      
	  MPI_Send(&(CSCD_AVAL(cscd_graph)[CSCD_IA(cscd_graph)[nbvertedge[2*nbproc+i]]-1]), 
		   nbvertedge[nbproc+i], COMM_FLOAT, i, 1, mpi_comm);
	}
      
      for(i=0; i<=vertloc; i++)
	ia[i] = CSCD_IA(cscd_graph)[i];
      for(i=0; i<edgeloc; i++)
	ja[i] = CSCD_JA(cscd_graph)[i];
      for(i=0; i<vertloc; i++)
	rhs[i] = CSCD_RHS(cscd_graph)[i];
      for(i=0; i<edgeloc; i++)
	aval[i] = CSCD_AVAL(cscd_graph)[i];
    }
  else 
    {
      MPI_Recv(ia,   vertloc, COMM_INT,    0, 1, mpi_comm, &status);
      MPI_Recv(ja,   edgeloc, COMM_INT,    0, 1, mpi_comm, &status);
      MPI_Recv(rhs,  vertloc, COMM_FLOAT, 0, 1, mpi_comm, &status);
      MPI_Recv(aval, edgeloc, COMM_FLOAT, 0, 1, mpi_comm, &status);

      for(i=vertloc-1; i>=0; i--)
	ia[i] = ia[i] - ia[0] + 1;
    }

  for(i=0; i<vertloc; i++)
    loc2glb[i] = nbvertedge[2*nbproc+myrank] + i + 1;
  
  cscd_Init(cscd_graph, loc2glb, vertloc, edgeloc, ia, ja, aval, rhs, mpi_comm);
  
  ia[vertloc] = edgeloc+1;
  
  return 0;
}

/* 
   Function: csc_Free
   
   Free csc structure

   Parameters:
     csc_graph - CSC structure to free.
*/
int csc_Free( CSC * csc_graph ) 
{
  if ( csc_graph != NULL ) 
    {
      if ( csc_graph->ia   != NULL ) memFree_null(csc_graph->ia);
      if ( csc_graph->ja   != NULL ) memFree_null(csc_graph->ja);
      if ( csc_graph->aval != NULL ) memFree_null(csc_graph->aval);
      if ( csc_graph->rhs  != NULL ) memFree_null(csc_graph->rhs);
      return EXIT_SUCCESS;
    }
  return EXIT_FAILURE;
}

/* 
   Function: cscd_Free
   
   Free a cscd structure.

   Parameters:
     cscd_graph - CSCD structure to free.

*/
int
cscd_Free( CSCD * cscd_graph ) 
{
  if (cscd_graph != NULL) 
    {

      if (cscd_graph->csc != NULL) 
	if (csc_Free(cscd_graph->csc) == EXIT_SUCCESS ) 
	  {
	    memFree_null(cscd_graph->csc);
	    memFree_null(cscd_graph->loc2glb);
	  }

#ifdef DISTRIBUTED
      if (cscd_graph->dgrafptr != NULL)  /* liberation de cscd_graph->dgrafptr */
	memFree_null(cscd_graph->dgrafptr);
#endif
    }
  return 0;
}


/* 
   Function: csc_Init

   Fill-in the csc structure.

   Parameters:
     csc_graph     - CSC structure to fill.
     vertlocnumber - Number of vertices/unknowns/Matrix size.
     edgelocnumber - Number of edge/Nnzeros
     ia            - starting index of each column in ja and aval.
     ja            - row of each element of the matrix.
     aval          - value of each element.
     rhs           - Right hand side
*/ 
int 
csc_Init(  CSC     *  csc_graph,
	   PASTIX_INT        vertlocnbr,
	   PASTIX_INT        edgelocnbr,
	   PASTIX_INT     *  ia,
	   PASTIX_INT     *  ja,
	   PASTIX_FLOAT   *  aval,
	   PASTIX_FLOAT   *  rhs) 
{
  
  csc_graph->vertlocnbr = vertlocnbr;
  csc_graph->edgelocnbr = edgelocnbr;
  csc_graph->ia         = NULL;
  csc_graph->ja         = NULL;
  csc_graph->aval       = NULL;
  csc_graph->rhs        = NULL;

  if (vertlocnbr > 0)
    {
      if (NULL == (csc_graph->ia = (PASTIX_INT *)memAlloc((vertlocnbr+1)*sizeof(PASTIX_INT))))
	MALLOC_ERROR("ia");
      memcpy(csc_graph->ia,ia,(vertlocnbr+1)*sizeof(PASTIX_INT));
    }
  if (edgelocnbr > 0)
    {
      if (NULL == (csc_graph->ja = (PASTIX_INT *)memAlloc(edgelocnbr*sizeof(PASTIX_INT))))
	MALLOC_ERROR("ja");
      memcpy(csc_graph->ja,ja,edgelocnbr*sizeof(PASTIX_INT));
    }

  /* Possibilite de passer aval = NULL */
  if (aval != NULL && edgelocnbr > 0)
    {
      if (NULL == (csc_graph->aval = (PASTIX_FLOAT *)memAlloc(edgelocnbr*sizeof(PASTIX_FLOAT))))
	MALLOC_ERROR("aval");
      memcpy(csc_graph->aval,aval,edgelocnbr*sizeof(PASTIX_FLOAT));
    }

  /* Possibilite de passer  rhs = NULL */
  if (rhs != NULL && vertlocnbr > 0)
    {
      if (NULL == (csc_graph->rhs = (PASTIX_FLOAT *)memAlloc(vertlocnbr*sizeof(PASTIX_FLOAT))))
	MALLOC_ERROR("rhs");
      memcpy(csc_graph->rhs,rhs,vertlocnbr*sizeof(PASTIX_FLOAT));
    }  
  return 0;
}
/* 
   Function: cscd_Init

   Fill-in the cscd structure.

   Parameters:
     cscd_graph    - CSCD structure to fill.
     loc2glb       - global index of each local column.
     vertlocnumber - Number of vertices/unknowns/Matrix size.
     edgelocnumber - Number of edge/Nnzeros
     ia            - starting index of each column in ja and aval.
     ja            - row of each element of the matrix.
     aval          - value of each element.
     rhs           - Right hand side.
     mpi_comm      - MPI communicator.
*/ 
int cscd_Init( CSCD     *  cscd_graph,
	       PASTIX_INT      *  loc2glb,
	       PASTIX_INT         vertlocnbr,
	       PASTIX_INT         edgelocnbr,
	       PASTIX_INT      *  ia, 	 
	       PASTIX_INT      *  ja, 
	       PASTIX_FLOAT    *  aval,
	       PASTIX_FLOAT    *  rhs,
	       MPI_Comm    mpi_comm)
{


  PASTIX_INT ret;
  cscd_graph->loc2glb  = NULL;
#ifdef DISTRIBUTED
  cscd_graph->dgrafptr = NULL;
#endif
  cscd_graph->mpi_comm = mpi_comm;
  cscd_graph->j_cst    = 10;
  cscd_graph->j_size   = edgelocnbr;

  /*  if(((*cscd_graph)      = (CSCD *)  memAlloc(sizeof(CSCD)))   ==  NULL  ) return -1; */
  if(NULL == (cscd_graph->csc = (CSC *)   memAlloc(sizeof(CSC )))) 
    MALLOC_ERROR("csc");


  /* Initialise la structure de CSC interne au CSCD */
  ret = csc_Init(cscd_graph->csc, vertlocnbr,	edgelocnbr, ia, ja, aval, rhs);
  if (ret != EXIT_SUCCESS) errorPrint("csc_Init returned : %d \n", (int)CSCD_VERTLOC(cscd_graph) );

  if (vertlocnbr > 0)
    {
      if (NULL == (cscd_graph->loc2glb = memAlloc(sizeof(PASTIX_INT)*vertlocnbr)))
	MALLOC_ERROR("loc2glb");
      memcpy(cscd_graph->loc2glb,loc2glb,sizeof(PASTIX_INT)*vertlocnbr);
    }

#ifdef DISTRIBUTED
  if (NULL == (cscd_graph->dgrafptr = memAlloc(sizeof(SCOTCH_Dgraph))))
    MALLOC_ERROR("dgrafptr");
  /* Initialise le dgraph en lui assignant un communicateur MPI */
  SCOTCH_dgraphInit(cscd_graph->dgrafptr, mpi_comm);
#endif /* DISTRIBUTED */

  return 0;
}


/* 
   Function: cscd_getGloballocal
   
   returns the global column of a local column.
   
   Parameters: 
     cscd_graph - CSCD structure.
     i          - local column number.
   
*/
PASTIX_INT 
cscd_getGloballocal(CSCD * cscd_graph, PASTIX_INT i) 
{
  return CSCD_LOC2GLB(cscd_graph)[i];
}


/*
  Function: cscd_indIfLoc
  
  Compute the local index corresponding to a global column number/vertex of the graph.
  
  Parameters:
    cscd_graph - graph
    vert_ind   - global vertex/column number
 */
/* renvoie -1 si le sommet au bout de l'arete est distant, l'indice local sinon */
PASTIX_INT cscd_indIfLoc(CSCD * cscd_graph, PASTIX_INT vert_ind) 
{
  PASTIX_INT min, max, pivot;

  min = 1; 
  max = CSCD_VERTLOC(cscd_graph);
  pivot = min + (max - min +1)/2;
  if (vert_ind > CSCD_LOC2GLB(cscd_graph)[max-1] || vert_ind < CSCD_LOC2GLB(cscd_graph)[min -1])
    return -1;

  while (CSCD_LOC2GLB(cscd_graph)[pivot-1] != vert_ind && max != min)
    {
      if (CSCD_LOC2GLB(cscd_graph)[pivot-1] < vert_ind)
	{
	  min = pivot;
	}
      else
	{
	  if (max == pivot)
	    max = pivot -1;
	  else
	    max = pivot;
	  
	}
      pivot = min + (max - min +1)/2;
    }

  if (CSCD_LOC2GLB(cscd_graph)[pivot-1] != vert_ind)
    return -1;
  else
    return pivot;
}

/*
  Function: cscd_add_diag

  
 */ 
int cscd_add_diag(CSCD * cscd_graph) 
{
  PASTIX_INT i;
  for (i=0; i<CSCD_VERTLOC(cscd_graph); i++)
    cscd_Addedge(cscd_graph, i, CSCD_LOC2GLB(cscd_graph)[i], 0);
  return 0;
}


int cscd_sort(CSCD * cscd_graph) 
{
  PASTIX_INT i;
  print_debug(DBG_CSCD, "->cscd_sort\n");
  for (i=0; i<CSCD_VERTLOC(cscd_graph); i++)
    cscd_quicksort(CSCD_JA(cscd_graph), 
		   CSCD_AVAL(cscd_graph), 
		   CSCD_IA(cscd_graph)[i]-1, 
		   CSCD_IA(cscd_graph)[i+1]-2);
  print_debug(DBG_CSCD, "<-cscd_sort\n");
  return 0;
}


void cscd_quicksort(PASTIX_INT *inttab, PASTIX_FLOAT *valtab, PASTIX_INT start, PASTIX_INT end)
{
  if(start < end) 
    {
      PASTIX_INT boundary;
      boundary = cscd_partition(inttab, valtab, start, end, start);
      if(boundary == -1)
	return;
      cscd_quicksort(inttab, valtab, start, boundary-1);
      cscd_quicksort(inttab, valtab, boundary+1, end);
    }
}
   
   
int cscd_partition(PASTIX_INT * inttab, PASTIX_FLOAT * valtab, PASTIX_INT start, PASTIX_INT end, PASTIX_INT pivot)
{
  PASTIX_INT   up;
  PASTIX_INT   down;
  PASTIX_INT   tmp;
  PASTIX_FLOAT tmpval;
  PASTIX_INT   cursor;
  up   = start;
  down = end;
  
  while(up<down)
    {
      /** Move UP to  the first node > PIVOT **/
      cursor = up;
      while( cursor <= end && inttab[cursor] <= inttab[pivot])
	cursor++;
      
      if(cursor <=end)
	up = cursor;
      else
	{
	  if(inttab[end] == inttab[pivot])
	    up = end;
	}

      /** Move DOWN to  the first node <= PIVOT **/
      cursor = down;
      while(cursor >= start  && inttab[cursor] > inttab[pivot] )
	cursor--;
      
      if(cursor >= start)
	down = cursor;


      if(up<down)
	{
	  /* Exchange up and down */
	  tmp          = inttab[down]; 
     	  inttab[down] = inttab[up]; 
	  inttab[up]   = tmp; 
	  if (valtab) 
	    {
	      tmpval       = valtab[down];
	      valtab[down] = valtab[up];	
	      valtab[up]   = tmpval;
	    }
	}
    }
  
  /*exchange value in down and pivindex */
  tmp           = inttab[down];
  inttab[down]  = inttab[pivot];
  inttab[pivot] = tmp;
  if (valtab) 
    {
      tmpval        = valtab[down];
      valtab[down]  = valtab[pivot];
      valtab[pivot] = tmpval;
    }

  return down;
}


/* Function: cscd_Sym 

   symetrize the graph structure.

   Parameters:
     cscd_graph - graph to symetrize
     iparm6     - determine what to do with second occurence of a coefficient.

 **/
int cscd_Sym(CSCD * cscd_graph, PASTIX_INT iparm6) 
{
  PASTIX_INT   * vectsend      = NULL;  
  PASTIX_INT   * vectsize      = NULL;
  PASTIX_INT   * vectsizercv   = NULL;
  PASTIX_FLOAT * vectvalues    = NULL;
  PASTIX_FLOAT * valuesvect    = NULL;
  PASTIX_FLOAT * valuesvectrcv = NULL;
  PASTIX_INT   * edgevect      = NULL;
  PASTIX_INT   * edgevectrcv   = NULL;
  PASTIX_INT     i, j, jj;
  PASTIX_INT     pred;
  PASTIX_INT     write; 
  int     myrank;
  int     nbproc;
  PASTIX_INT     pos_write  = 0;
  PASTIX_INT     count      = 0;  /* count contient le nombre d'arete a envoyer */
  /* on alloue à la moitié des paires d'entier que l'on peut avoir au max à envoyer */
  PASTIX_INT     size       = (CSCD_EDGELOC(cscd_graph))/2;
  
  MPI_Comm_rank( cscd_graph->mpi_comm, &myrank);
  MPI_Comm_size( cscd_graph->mpi_comm, &nbproc);

  if (CSCD_AVAL(cscd_graph) != NULL)
    vectvalues = memAlloc(size*sizeof(PASTIX_FLOAT));

  vectsend = memAlloc(size*2*sizeof(PASTIX_INT));/* vectsend contient le vecteur d'arete a envoyer */

  /* PHASE DE CONSTRUCTION */
  /* determine la taille de Ti et recopie des aretes a distribuer */
  for (i = 0; i < CSCD_VERTLOC(cscd_graph); i++)
    for (j = CSCD_IA(cscd_graph)[i]-1; j < CSCD_IA(cscd_graph)[i+1]-1; j++) 
      { 
	if (cscd_indIfLoc(cscd_graph, CSCD_JA(cscd_graph)[j]) == -1) 
	  {
	    if (count >= size) 
	      {
		size      += (size+1)/2; /* +1 au cas ou size = 0 */
		vectsend   = memRealloc(vectsend,  size*2*sizeof(PASTIX_INT));
		if (CSCD_AVAL(cscd_graph) != NULL)
		  vectvalues = memRealloc(vectvalues,size*sizeof(PASTIX_FLOAT));
	    }
	    /* 	  fprintf(stdout, "[P%d] %d-%d  arete ajoutee...\n", myrank, CSCD_LOC2GLB(cscd_graph)[i], CSCD_JA(cscd_graph)[j]);  */
	    vectsend[2*count]   = CSCD_JA(cscd_graph)[j];
	    vectsend[2*count+1] = CSCD_LOC2GLB(cscd_graph)[i];
	    if (CSCD_AVAL(cscd_graph) != NULL)
	      vectvalues[count]   = CSCD_AVAL(cscd_graph)[j];
	    count++;
	  }
      }
  
/*   fprintf(stdout, "[P%d] nb d'arete a echanger : %d\n", (int)myrank, (int)count); */


  /* PHASE D'ECHANGE */
  vectsize    = memAlloc(sizeof(PASTIX_INT)*(nbproc+1));
  vectsizercv = memAlloc(sizeof(PASTIX_INT)*(nbproc+1));
  vectsize    = memSet(vectsize,   0,sizeof(PASTIX_INT)*nbproc+1);
  vectsizercv = memSet(vectsizercv,0,sizeof(PASTIX_INT)*nbproc+1);
  vectsize[myrank] = count; /* on va indiquer a son successeur que l'on occupera count place */
  vectsize[nbproc] = count; /* ce champ va sommer le nombre d'arete a echanger.              */
  
  MPI_Allreduce(vectsize, vectsizercv, nbproc+1, COMM_INT, MPI_SUM, cscd_graph->mpi_comm);
  
  edgevect    = memAlloc(sizeof(PASTIX_INT)*2*vectsizercv[nbproc]); 
  edgevectrcv = memAlloc(sizeof(PASTIX_INT)*2*vectsizercv[nbproc]); 
  edgevect    = memSet(edgevect,   0,sizeof(PASTIX_INT)*2*vectsizercv[nbproc]);
  edgevectrcv = memSet(edgevectrcv,0,sizeof(PASTIX_INT)*2*vectsizercv[nbproc]);
  for (i=0; i<myrank; i++)
    pos_write +=vectsizercv[i]; /* on connait l'endroit ou l'on doit ecrire. */

  memcpy(&edgevect[2*pos_write]  ,vectsend  ,2*count*sizeof(PASTIX_INT));

  /* le tableau edgevect ne contient que des 0, sauf aux endroits ou je devais mettre mes aretes */
  MPI_Allreduce(edgevect, edgevectrcv, 2*vectsizercv[nbproc],COMM_INT, MPI_SUM, cscd_graph->mpi_comm);

  if (CSCD_AVAL(cscd_graph) != NULL) {
    valuesvect     = memAlloc(sizeof(PASTIX_FLOAT)*vectsizercv[nbproc]);
    valuesvectrcv  = memAlloc(sizeof(PASTIX_FLOAT)*vectsizercv[nbproc]);
    valuesvect     = memSet(valuesvect   ,   0,sizeof(PASTIX_FLOAT)*vectsizercv[nbproc]);
    valuesvectrcv  = memSet(valuesvectrcv,   0,sizeof(PASTIX_FLOAT)*vectsizercv[nbproc]);
    memcpy(&valuesvect[pos_write],vectvalues,count*sizeof(PASTIX_FLOAT));
    MPI_Allreduce(valuesvect, valuesvectrcv, vectsizercv[nbproc],COMM_FLOAT, MPI_SUM, cscd_graph->mpi_comm);  
  }

  /* tout le monde possede l'ensemble des aretes de chacun des halos */

  
  /* PHASE D'ANALYSE */
  
  /* chaque proc va devoir lire le edgevect */
  /* sommet local qui va gagner une arete... */

  /* On compte de combien on doit faire grossir la cscd */
  count = 0; 
  for (i = 0; i < CSCD_VERTLOC(cscd_graph); i++)
    for (j = CSCD_IA(cscd_graph)[i]-1; j < CSCD_IA(cscd_graph)[i+1]-1; j++) 
      if (cscd_indIfLoc(cscd_graph, CSCD_JA(cscd_graph)[j]) != -1 &&
	  CSCD_JA(cscd_graph)[j] != CSCD_LOC2GLB(cscd_graph)[i])
	count ++;
  for (i=0; i<vectsizercv[nbproc]; i++) 
    {
      pred = cscd_indIfLoc    (cscd_graph, edgevectrcv[2*i]);
      if (pred>=0)  /* pred est local, un seul proc peut repondre vrai a ce test */
	{
	  write=0;
	  for(j = CSCD_IA(cscd_graph)[pred-1]-1; j < CSCD_IA(cscd_graph)[pred]-1; j++) 
	    {
	      /*fprintf(stderr, "[P%d] arete %d -- %d \n", myrank,edgevectrcv[2*i], CSCD_JA(cscd_graph)[j]); */
	      if ( CSCD_JA(cscd_graph)[j] == edgevectrcv[2*i+1] )
		{ 
		  write = -1; /* l'arete existe deja */
		  break;
		}
	    }
	  if (write==0 || iparm6) /* iparm 6 a 1 oblige la recopie des aretes */
	    { 
/* 	      fprintf(stderr, "[P%d] Ajout de %ld -- %ld\n", (int)myrank,(long)edgevectrcv[2*i], (long)edgevectrcv[2*i+1] ); */
	      count++;
	      /* cscd_Addedge(cscd_graph, pred, edgevectrcv[2*i+1]); */
	    }
	}
    }

  cscd_graph->j_size  = CSCD_EDGELOC(cscd_graph)+count;
  CSCD_JA(cscd_graph) = memRealloc(CSCD_JA(cscd_graph),(CSCD_EDGELOC(cscd_graph)+count)*sizeof(PASTIX_INT));
  if (CSCD_AVAL(cscd_graph) != NULL)
    CSCD_AVAL(cscd_graph) = memRealloc(CSCD_AVAL(cscd_graph),(CSCD_EDGELOC(cscd_graph)+count)*sizeof(PASTIX_FLOAT));
 
  for (i = 0; i < CSCD_VERTLOC(cscd_graph); i++)
    for (j = CSCD_IA(cscd_graph)[i]-1; j < CSCD_IA(cscd_graph)[i+1]-1; j++) 
      {
	pred = cscd_indIfLoc(cscd_graph, CSCD_JA(cscd_graph)[j]);
	if (pred != -1 && CSCD_JA(cscd_graph)[j] != CSCD_LOC2GLB(cscd_graph)[i])
	  {
	    write=0;
	    for(jj = CSCD_IA(cscd_graph)[pred-1]-1; jj < CSCD_IA(cscd_graph)[pred]-1; jj++) 
	      {
		if ( CSCD_LOC2GLB(cscd_graph)[i] == CSCD_JA(cscd_graph)[jj] )
		  {
		    write = -1; /* l'arete existe deja */
		    break;
		  }
	      }
	    if (write == 0)
	      {
/* 		fprintf(stdout, "[P%d] Ajout de %ld --- %ld\n",   */
/* 			(int)myrank,(long)pred, (long)CSCD_LOC2GLB(cscd_graph)[i]);  */
		if (CSCD_AVAL(cscd_graph) != NULL)
		  cscd_Addedge(cscd_graph,  pred, CSCD_LOC2GLB(cscd_graph)[i],valuesvectrcv[i]); 
		else
		  cscd_Addedge(cscd_graph,  pred, CSCD_LOC2GLB(cscd_graph)[i], 0); 
	      }
	  }
      }
   for (i=0; i<vectsizercv[nbproc]; i++) 
    {
      pred = cscd_indIfLoc    (cscd_graph, edgevectrcv[2*i]);
      if (pred>=0)  /* pred est local, un seul proc peut repondre vrai a ce test */
	{
	  write=0;
	  for(j = CSCD_IA(cscd_graph)[pred-1]-1; j < CSCD_IA(cscd_graph)[pred]-1; j++) 
	    {
	      /*fprintf(stderr, "[P%d] arete %d -- %d \n", myrank,edgevectrcv[2*i], CSCD_JA(cscd_graph)[j]); */
	      if ( CSCD_JA(cscd_graph)[j] == edgevectrcv[2*i+1] )
		{ 
		  write = -1; /* l'arete existe deja */
		  break;
		}
	    }

	  if (write==0 || iparm6) /* iparm 6 a 1 oblige la recopie des aretes */
	    { 
/* 	      if (CSCD_LOC2GLB(cscd_graph)[pred-1] == 12111 &&  edgevectrcv[2*i+1] == 5000) */
/* 		fprintf(stdout, "[P%d] Ajout de %ld -- %ld\n", */
/* 			(int)myrank,(long)pred, (long)edgevectrcv[2*i+1] ); */
	      if (CSCD_AVAL(cscd_graph) != NULL)
		cscd_Addedge(cscd_graph, pred, edgevectrcv[2*i+1],valuesvectrcv[i]); 
	      else
		cscd_Addedge(cscd_graph, pred, edgevectrcv[2*i+1],0); 
	    }
	}
    }

   if (NULL != CSCD_AVAL(cscd_graph))
     {
       memFree_null(vectvalues); 
       memFree_null(valuesvect);
       memFree_null(valuesvectrcv);
     }
   memFree_null(vectsend);
   memFree_null(vectsize);
   memFree_null(vectsizercv);
   memFree_null(edgevect);
   memFree_null(edgevectrcv);
   return 0;
}

/* 
   Function: cscd_Addedge

   Add an edge between i and j with value 0 if aval is allocated.
   After that, rows are not sorted anymore.

   Parameters: 
     cscd_graph - CSCD
     i          - vertex/column (local)
     j          - edge/row 
     value      - value if aval of cscd_graph is not NULL
 */
/** cscd_Addedge ajoute une arete a la structure entre le sommet de 
 ** numerotation local i et le sommet de numerotation global j
 ** On supposera que cette arete n'existe pas. Aucune verification 
 ** n'est faite en ce sens. La fonction appelante le fait.
 **/
int cscd_Addedge(CSCD * cscd_graph, PASTIX_INT i, PASTIX_INT j, PASTIX_FLOAT value) 
{
  PASTIX_INT k;

  if ( CSCD_EDGELOC(cscd_graph) >= cscd_graph->j_size) /* besoin d'une realloc :) */
    {
      cscd_graph->j_cst  *= 2;
      cscd_graph->j_size += cscd_graph->j_cst;
      
      CSCD_JA(cscd_graph) = (PASTIX_INT *) memRealloc(CSCD_JA(cscd_graph), sizeof(PASTIX_INT)*(cscd_graph->j_size));
      if (CSCD_AVAL(cscd_graph))
	CSCD_AVAL(cscd_graph) = (PASTIX_FLOAT *)memRealloc(CSCD_AVAL(cscd_graph), sizeof(PASTIX_FLOAT)*(cscd_graph->j_size));
    }
  
/*   fprintf(stdout, " ADDING %d %d %d %d\n", i, j,CSCD_VERTLOC(cscd_graph) , CSCD_EDGELOC(cscd_graph)); */
  /* la recopie se fait par decalage des valeurs aux frontieres */
  for(k=CSCD_VERTLOC(cscd_graph); k>i-1; k--)
    {
/*       fprintf(stdout, " ADDING %d %d : CSCD_JA[%d] = CSCD_JA[%d] = %d\n", i, j, */
/* 	      CSCD_IA(cscd_graph)[k]-1, */
/* 	      CSCD_IA(cscd_graph)[k-1]-1,  */
/* 	      CSCD_JA(cscd_graph)[CSCD_IA(cscd_graph)[k-1]-1]); */

      CSCD_JA(cscd_graph)[CSCD_IA(cscd_graph)[k]-1] = CSCD_JA(cscd_graph)[CSCD_IA(cscd_graph)[k-1]-1];
      if (CSCD_AVAL(cscd_graph)) 
	CSCD_AVAL(cscd_graph)[CSCD_IA(cscd_graph)[k]-1] = CSCD_AVAL(cscd_graph)[CSCD_IA(cscd_graph)[k-1]-1];
      CSCD_IA(cscd_graph)[k]++; 
    }

/*   fprintf(stdout, "CSCD_JA[%d] = %d\n",CSCD_IA(cscd_graph)[i-1]-1, j); */
  CSCD_JA(cscd_graph)[CSCD_IA(cscd_graph)[i-1]-1] = j;
  if (CSCD_AVAL(cscd_graph)) 
    CSCD_AVAL(cscd_graph)[CSCD_IA(cscd_graph)[i-1]-1] = value;
  CSCD_EDGELOC(cscd_graph)++;
/*   fprintf(stdout, " ADDING2 %d %d %d %d\n", i, j,CSCD_VERTLOC(cscd_graph) , CSCD_EDGELOC(cscd_graph)); */
  return 0;
}


/* Function: cscd_Convert_dgraph 
   
   Symetrize the cscd then convert it to dgraph.
   After that you can call dgraph_build.

   Symetrise la structure puis la convertie en
   dgraph. On apelle alors le dgraph_build sur le dgraph obtenu

   Parameters:
     cscd_graph - csc to convert to dgraph
     iparm6     - 
 */
int cscd_Convert_dgraph(CSCD * cscd_graph, PASTIX_INT iparm6) 
{
  int myrank;
  PASTIX_INT dgraph_ret, cscd_ret;
  PASTIX_INT vertglbmax, vertlocmax;

  vertlocmax = CSCD_VERTLOC(cscd_graph);

  MPI_Allreduce(&vertlocmax, &vertglbmax, 1, COMM_INT, MPI_MAX, cscd_graph->mpi_comm);

  MPI_Comm_rank( cscd_graph->mpi_comm, &myrank);

  /* Symetriser la structure CSCd */
  cscd_ret = cscd_Sym(cscd_graph, iparm6);

  /* cscd_ret = cscd_add_diag(cscd_graph); */

  cscd_ret = cscd_sort(cscd_graph);

  /*   PASTIX_INT * edgegsttab = memAlloc(sizeof(PASTIX_INT)*CSCD_EDGELOC(cscd_graph)); */
  /*   if (edgegsttab == NULL) fprintf(stderr, "Error : memAlloc\n"); */

  /* Remplir la structure Dgraph */
#ifdef DISTRIBUTED
  dgraph_ret = 
    SCOTCH_dgraphBuild (
			cscd_graph->dgrafptr,
			1,                        /* baseval */
			CSCD_VERTLOC(cscd_graph),
			vertglbmax,               /* Maximum number of local vertices     */
			CSCD_IA(cscd_graph),
			CSCD_IA(cscd_graph)+1,
			NULL,                     /* Local vertex load array (if any)     */
			/*cscd_getLoc2Glb(cscd_graph),*/
			NULL,                     /* Local vertex label array (if any)    */
			CSCD_EDGELOC(cscd_graph),
			CSCD_EDGELOC(cscd_graph),
			CSCD_JA(cscd_graph),      /* Local edge array                     */
			NULL,                     /* Ghost edge array (if any); not const */
			NULL );
  
  if (dgraph_ret != 0) 
    fprintf(stderr, "Error : Convert CSCd/Dgraph > dgraphBuild\n");


  /*memFree_null(edgegsttab);*/
  
  dgraph_ret = SCOTCH_dgraphGhst(cscd_graph->dgrafptr);
  if (dgraph_ret != 0) 
    fprintf(stderr, "Error : Dgraph ghost\n");
#endif  
  return 0;
}


/** cscd_Addvertices ajoute un ensemble de sommets a la csc locale, ce qui implique de connaitre
 ** la position de ces sommets par rapport aux autres, leurs numeros globaux, leurs nombres d'arete
 ** ainsi que leurs vecteurs d'arete
 ** On suppose que les sommets sont ordonnes suivant leur numero locaux.
 **/
int cscd_Addvertices(CSCD * cscd_graph, PASTIX_INT * loc2glb, PASTIX_INT * ia, PASTIX_INT * ja, PASTIX_FLOAT * aval, PASTIX_FLOAT * rhs, PASTIX_INT nbedge, PASTIX_INT nbvert) 
{
  MPI_Comm mpi_comm = cscd_graph->mpi_comm;
  PASTIX_FLOAT *rhsf;
  PASTIX_FLOAT *avalf;
  PASTIX_INT   *iaf;
  PASTIX_INT   *jaf;
  PASTIX_INT   *loc2glbf;
  PASTIX_INT    i,j;
  PASTIX_INT    vali = -1, valj = -1;
  PASTIX_INT    vert_tot   = nbvert + CSCD_VERTLOC(cscd_graph);
  PASTIX_INT    edge_tot   = nbedge + CSCD_EDGELOC(cscd_graph);

  iaf      = (PASTIX_INT   *)memAlloc(sizeof(PASTIX_INT)  *(vert_tot+1));
  jaf      = (PASTIX_INT   *)memAlloc(sizeof(PASTIX_INT)  *edge_tot);
  loc2glbf = (PASTIX_INT   *)memAlloc(sizeof(PASTIX_INT)  *vert_tot);
  rhsf     = (PASTIX_FLOAT *)memAlloc(sizeof(PASTIX_FLOAT)*vert_tot);
  avalf    = (PASTIX_FLOAT *)memAlloc(sizeof(PASTIX_FLOAT)*edge_tot);
  
  iaf[0]        = 1;
  iaf[vert_tot] = edge_tot+iaf[0];

  for(i=0,j=0;i+j<vert_tot;) 
    {
      if (i >= CSCD_VERTLOC(cscd_graph) && j < nbvert) 
	{
	  valj = loc2glb[j]; /* n'ayant plus de sommet dans i, on finira d'ecrire j */
	  vali = valj+1;     /* pour cela, vali sera toujours superieur a valj      */
	}
      
      if (i < CSCD_VERTLOC(cscd_graph) && j >= nbvert) 
	{
	  vali = CSCD_LOC2GLB(cscd_graph)[i]; /* n'ayant plus de sommet dans j, on finira d'ecrire i */
	  valj = vali+1;                      /* pour cela, valj sera toujours superieur a vali      */
	}
      
      if (i < CSCD_VERTLOC(cscd_graph) && j < nbvert) 
	{
	  vali = CSCD_LOC2GLB(cscd_graph)[i]; valj = loc2glb[j];
	}

      /* ajout du sommet de cscd */
      if (vali < valj) 
	{
	  loc2glbf[i+j] = CSCD_LOC2GLB(cscd_graph)[i];
	  iaf[i+j+1]    = iaf[i+j] + CSCD_IA(cscd_graph)[i+1] - CSCD_IA(cscd_graph)[i];
	  rhsf[i+j]     = CSCD_RHS(cscd_graph)[i];
	  
	  memCpy(&(jaf[iaf[i+j]-1]), &(CSCD_JA(cscd_graph)[CSCD_IA(cscd_graph)[i]-1]), 
		 (CSCD_IA(cscd_graph)[i+1]-CSCD_IA(cscd_graph)[i])*sizeof(PASTIX_INT) );
	  memCpy(&(avalf[iaf[i+j]-1]), &(CSCD_AVAL(cscd_graph)[CSCD_IA(cscd_graph)[i]-1]), 
		 (CSCD_IA(cscd_graph)[i+1]-CSCD_IA(cscd_graph)[i])*sizeof(PASTIX_FLOAT));
	  i++;
	}
      /* ajout du sommet depuis le pool de sommets */
      else 
	{
	  loc2glbf[i+j] = loc2glb[j];
	  iaf[i+j+1]    = iaf[i+j] + ia[j+1] - ia[j];
	  rhsf[i+j]     = rhs[j];
	  
	  memCpy(&(jaf[iaf[i+j]-1]),   &(ja[ia[j]-1]),   (ia[j+1] - ia[j])*sizeof(PASTIX_INT));
	  memCpy(&(avalf[iaf[i+j]-1]), &(aval[ia[j]-1]), (ia[j+1] - ia[j])*sizeof(PASTIX_FLOAT));
	  j++;
	}
    }

  cscd_Free(cscd_graph);
  cscd_Init(cscd_graph, loc2glbf, vert_tot, edge_tot, iaf, jaf, avalf, rhsf, mpi_comm);
  
  return 0;
}

/** rank : rang du noeud ou va aller la cscd en construction
 ** exchange_tab : donne les le rang auquel appartient le sommet local
 **/ 
int cscd_Remvertices(CSCD * cscd_graph, CSCD * cscd_target, PASTIX_INT * exchange_tab, PASTIX_INT rank) 
{
  MPI_Comm mpi_comm = cscd_graph->mpi_comm;
  PASTIX_INT     *iaf;
  PASTIX_INT     *jaf;
  PASTIX_INT     *loc2glbf;
  PASTIX_FLOAT   *rhsf     = NULL;
  PASTIX_FLOAT   *avalf    = NULL;
  PASTIX_INT     *iaex;
  PASTIX_INT     *jaex;
  PASTIX_INT     *loc2glbex;
  PASTIX_FLOAT   *rhsex    = NULL;
  PASTIX_FLOAT   *avalex   = NULL;
  PASTIX_INT      i, j, vert_tot, edge_tot;
  PASTIX_INT      vert_ex  = 0;
  PASTIX_INT      edge_ex  = 0;
  /* vert_tot, sommet restant locaux, vert_ex, sommets a echanger. idem pour edge... */

  for(i=0; i<CSCD_VERTLOC(cscd_graph); i++) 
    {
      if(exchange_tab[i] == rank) 
	{
	  vert_ex ++;
	  edge_ex += CSCD_IA(cscd_graph)[i+1] - CSCD_IA(cscd_graph)[i];
	}
    }

  vert_tot = CSCD_VERTLOC(cscd_graph) - vert_ex; 
  edge_tot = CSCD_EDGELOC(cscd_graph) - edge_ex;

  iaf       = (PASTIX_INT   *)memAlloc(sizeof(PASTIX_INT)   * (vert_tot+1));
  jaf       = (PASTIX_INT   *)memAlloc(sizeof(PASTIX_INT)   * edge_tot);
  loc2glbf  = (PASTIX_INT   *)memAlloc(sizeof(PASTIX_INT)   * vert_tot);
  if (CSCD_RHS(cscd_graph))  
    rhsf    = (PASTIX_FLOAT *)memAlloc(sizeof(PASTIX_FLOAT) * vert_tot); 
  if (CSCD_AVAL(cscd_graph)) 
    avalf   = (PASTIX_FLOAT *)memAlloc(sizeof(PASTIX_FLOAT) * edge_tot);
  
  iaex      = (PASTIX_INT   *)memAlloc(sizeof(PASTIX_INT)   * (vert_ex+1));
  jaex      = (PASTIX_INT   *)memAlloc(sizeof(PASTIX_INT)   * edge_ex);
  loc2glbex = (PASTIX_INT   *)memAlloc(sizeof(PASTIX_INT)   * vert_ex);
  if (CSCD_RHS(cscd_graph))
    rhsex   = (PASTIX_FLOAT *)memAlloc(sizeof(PASTIX_FLOAT) * vert_ex);
  if (CSCD_AVAL(cscd_graph)) 
    avalex  = (PASTIX_FLOAT *)memAlloc(sizeof(PASTIX_FLOAT) * edge_ex);
  
  iaf[0]  = 1; 
  iaex[0] = 1;
  iaf[vert_tot] = edge_tot; 
  iaex[vert_ex] = edge_ex;
  
  for(i=0, j=0; i+j<CSCD_VERTLOC(cscd_graph); ) 
    {
      /* i, indice du vecteur d'echange */
      if(exchange_tab[i+j] == rank) 
	{
	  iaex[i+1] = iaex[i] + CSCD_IA(cscd_graph)[i+j+1] - CSCD_IA(cscd_graph)[i+j];
	  loc2glbex[i] = CSCD_LOC2GLB(cscd_graph)[i+j];
	  if (CSCD_RHS(cscd_graph))
	    rhsex[i] = CSCD_RHS(cscd_graph)[i+j];

	  memCpy(&(jaex[iaex[i]-1]), &(CSCD_JA(cscd_graph)[CSCD_IA(cscd_graph)[i+j]-1]), 
		 (CSCD_IA(cscd_graph)[i+j+1]-CSCD_IA(cscd_graph)[i+j])*sizeof(PASTIX_INT));
	  if (CSCD_AVAL(cscd_graph))
	    memCpy(&(avalex[iaex[i]-1]), &(CSCD_AVAL(cscd_graph)[CSCD_IA(cscd_graph)[i+j]-1]), 
		   (CSCD_IA(cscd_graph)[i+j+1]-CSCD_IA(cscd_graph)[i+j])*sizeof(PASTIX_FLOAT));
	  i++;
	}
      /* j, indice du vecteur local */
      else 
	{
	  iaf[j+1]    = iaf[j] + CSCD_IA(cscd_graph)[i+j+1] - CSCD_IA(cscd_graph)[i+j];
	  loc2glbf[j] = CSCD_LOC2GLB(cscd_graph)[i+j];
	  if (CSCD_RHS(cscd_graph))
	    rhsf[j] = CSCD_RHS(cscd_graph)[i+j];
	  memCpy(&(jaf[iaf[j]-1]), &(CSCD_JA(cscd_graph)[CSCD_IA(cscd_graph)[i+j]-1]), 
		 (CSCD_IA(cscd_graph)[i+j+1]-CSCD_IA(cscd_graph)[i+j])*sizeof(PASTIX_INT));
	  if (CSCD_AVAL(cscd_graph))
	    memCpy(&(avalf[iaf[j]-1]), &(CSCD_AVAL(cscd_graph)[CSCD_IA(cscd_graph)[i+j]-1]), 
		   (CSCD_IA(cscd_graph)[i+j+1]-CSCD_IA(cscd_graph)[i+j])*sizeof(PASTIX_FLOAT));
	  j++;
	}
    }

  cscd_Init(cscd_target, loc2glbex, vert_ex,  edge_ex,  iaex, jaex, avalex, rhsex, mpi_comm);
  cscd_Free(cscd_graph);
  cscd_Init(cscd_graph,  loc2glbf,  vert_tot, edge_tot, iaf,  jaf,  avalf,  rhsf,  mpi_comm);
  
  return 0;
}




int cscd_Redistribute(CSCD * cscd_graph, PASTIX_INT * permtab) 
{

  /** On suppose que la perm tab donne la localisation de chaque sommet.
   ** Si ce n'etait pas le cas, avec le vecteur de renumerotation, 
   ** et le vecteur donnant le partitionnement global, on est en mesure de le determiner **/

  /* Calcul eventuel du vecteur de permutation */

  /* Necessite de savoir a qui on doit envoyer nos sommets. */
  /* Il serait plus aise de savoir qui va nous envoyer des sommets */
  
  /* Au cas ou. On construit un vecteur de taille nbproc^2 */
  /* On compte le nombre de sommet que l'on doit envoyer a chaque proc */
  /* On met ce nombre dans le vecteur en position i*nbproc+j,  */
  /* ou i correspond au proc local, et j le proc distant */
  /* un piti MPI_Allreduce(MPI_SUM/MPI_MAX)  */
  /* et hop, chacun peut savoir qui va lui parler et cb de sommet il va envoyer */
  
  /* Une fois que l'on a connaissance de qui va nous envoyer des messages,  */
  /* on initialise des receptions non bloquantes que l'on place dans une table de requete */
  
  /* On commence a decouper la CSCD locale en fonction de chaque processeur */
  
  /* On envoie chaque vecteur a son nouveau proprio.  */
  /* Pour cela une fonction gerant l'envoie de tous les vecteurs doit etre ecrite. */
  
  /* On se met en WaitAny sur la table de requete */
  
  /* Des qu'une requete est completee, on recupere tous les vecteurs en appelant une fonction */
  /* qui va se charger des comm et de l'appel a Addvertices(...) */
  
  
  /* On reboucle sur le WaitAny tant que l'ensemble des requetes n'est pas complete */
    
  return 0;
}



/*
PASTIX_INT 
cscd_Sendint(PASTIX_INT * vect, PASTIX_INT size, PASTIX_INT cible) {

  PASTIX_INT i;
  PASTIX_INT size_max = 1024*1024*1024/sizeof(PASTIX_INT); // 1Go par buffer, exemple...






  return 0;
}


PASTIX_INT
cscd_Senddouble(PASTIX_FLOAT * vect, PASTIX_INT size, PASTIX_INT cible) {

  PASTIX_INT i;
  PASTIX_INT size_max = 1024*1024*1024/sizeof(PASTIX_FLOAT); // 1Go par buffer, exemple...
  
  PASTIX_INT nbpacket = (PASTIX_INT) size/size_max +1;
  
  MPI_Send();

  



  return 0;
}
*/


/* Le cscd_graph ne contient que le squelette de la csc. Aucune valeur dedans
   cscd_val contient plus de sommet que la cscd_graph ainsi que toutes les valeurs associes.
   Le but est d'inserer les conditions aux limites (rhs+aval) dans la premier cscd.
   On mettra des zeros dans les valeurs qui ont ete cree par la symetrisation de la matrice.
*/
int cscd_Addvalues(CSCD * cscd_graph, CSCD * cscd_val) 
{
  PASTIX_INT i_graph, i_val;
  PASTIX_INT j_graph, j_val;

  if (!CSCD_AVAL(cscd_graph)) CSCD_AVAL(cscd_graph) = memAlloc(sizeof(PASTIX_FLOAT) * CSCD_EDGELOC(cscd_graph));
  if (!CSCD_RHS(cscd_graph))  CSCD_RHS(cscd_graph)  = memAlloc(sizeof(PASTIX_FLOAT) * CSCD_VERTLOC(cscd_graph));

  CSCD_AVAL(cscd_graph) = memSet(CSCD_AVAL(cscd_graph), 0, sizeof(PASTIX_FLOAT)*CSCD_EDGELOC(cscd_graph));
  CSCD_RHS(cscd_graph)  = memSet(CSCD_RHS(cscd_graph) , 0, sizeof(PASTIX_FLOAT)*CSCD_VERTLOC(cscd_graph));

  for(i_val=0, i_graph=0; i_graph < CSCD_VERTLOC(cscd_graph) && i_val < CSCD_VERTLOC(cscd_val); )
    {
      
      while (CSCD_LOC2GLB(cscd_graph)[i_graph] < CSCD_LOC2GLB(cscd_val)[i_val] 
	     && i_graph < CSCD_VERTLOC(cscd_graph))
	i_graph++;
      if (i_graph < CSCD_VERTLOC(cscd_graph))
	while (CSCD_LOC2GLB(cscd_graph)[i_graph] > CSCD_LOC2GLB(cscd_val)[i_val] 
	       && i_val < CSCD_VERTLOC(cscd_val))
	  i_val++;
      
      if (i_graph < CSCD_VERTLOC(cscd_graph) && i_val < CSCD_VERTLOC(cscd_val))
	if (CSCD_LOC2GLB(cscd_graph)[i_graph] == CSCD_LOC2GLB(cscd_val)[i_val]) 
	  {
	    /* copie du rhs ! */
	    /* fprintf(stderr, "%d \n", i_graph); */
	    CSCD_RHS(cscd_graph)[i_graph] = CSCD_RHS(cscd_val)[i_val];
	    
	    /* parcours du vecteur */
	    for(j_val = CSCD_IA(cscd_val)[i_val]-1, j_graph=CSCD_IA(cscd_graph)[i_graph]-1; 
		(j_graph < CSCD_IA(cscd_graph)[i_graph+1]-1) && (j_val < CSCD_IA(cscd_val)[i_val+1]-1); ) 
	      {
		while (CSCD_JA(cscd_graph)[j_graph] < CSCD_JA(cscd_val)[j_val] 
		       && j_graph < CSCD_IA(cscd_graph)[i_graph+1]-1)
		  j_graph++;

		if (j_graph < CSCD_IA(cscd_graph)[i_graph+1]-1)
		  while (CSCD_JA(cscd_val)[j_val] < CSCD_JA(cscd_graph)[j_graph] 
			 && j_val < CSCD_IA(cscd_val)[i_val+1]-1)
		    j_val++;
		
		if (j_graph < CSCD_IA(cscd_graph)[i_graph+1]-1 && j_val < CSCD_IA(cscd_val)[i_val+1]-1) 
		  {
		    CSCD_AVAL(cscd_graph)[j_graph] = CSCD_AVAL(cscd_val)[j_val];
		    j_graph++; j_val++;
		  }
	      }
	    
	    i_graph++; i_val++;      
	  }
    }
  return 0;
}


