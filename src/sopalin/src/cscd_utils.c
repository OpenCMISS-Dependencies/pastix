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
/*
 * File: cscd_utils.c
 *
 * Several operations on CSCD.
 *
 */
#ifdef FORCE_NOMPI
#  include "nompi.h"
#else
#  include <mpi.h>
#endif

#include "common_pastix.h"
#include "tools.h"
#include "cscd_utils.h"
#include "cscd_utils_intern.h"

#define TAG_SIZE   1
#define TAG_COL    2
#define TAG_ROW    3
#define TAG_VAL    4
#define TAG_L2G    5
#define TAG_RHS    6
#define TAG_SIZE2  7
#define TAG_COUPLE 8

/*
 *  Function: csc_dispatch
 *
 *  Distribute a CSC to a CSCD
 *
 *  Parameters:
 *     gN                - global number of columns
 *     gcolptr           - global starting index of each column in *grows* and
 *                         *gavals*.
 *     grows             - global rows of each element.
 *     gavals            - global values of each element.
 *     gperm             - global permutation tabular.
 *     ginvp             - global reverse permutation tabular.
 *     lN                - local number of columns (output).
 *     lcolptr           - starting index of each local column (output).
 *     lrowptr           - row number of each local element (output).
 *     lavals            - values of each local element (output).
 *     lrhs              - local part of the right hand side (output).
 *     lperm             - local part of the permutation tabular (output).
 *     loc2glob          - global numbers of local columns (before permutation).
 *     dispatch          - choose how to dispatch the csc
 *     pastix_comm       - PaStiX MPI communicator.
 */
#ifndef csc_dispatch
#  error "redefinition required"
#endif
void csc_dispatch(PASTIX_INT  gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow, PASTIX_FLOAT *  gavals,
                  PASTIX_FLOAT *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
                  PASTIX_INT *lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow, PASTIX_FLOAT ** lavals,
                  PASTIX_FLOAT ** lrhs, PASTIX_INT ** lperm,
                  PASTIX_INT **loc2glob, int dispatch, MPI_Comm pastix_comm)
{
  PASTIX_INT i;
  int rank, comm_size;
  PASTIX_INT index;
  PASTIX_INT rowindex;
  (void)pastix_comm;
  (void)ginvp;

  MPI_Comm_rank(pastix_comm,&rank);
  MPI_Comm_size(pastix_comm,&comm_size);

  /* Split proportionnaly the columns number on each proc */
  *lN = gN / comm_size;
  if (rank < (gN%comm_size))
    (*lN)++;

  /* Ici on a des malloc parce que le free est externe */
  MALLOC_EXTERN((*lcolptr), (*lN + 1), PASTIX_INT);
  MALLOC_EXTERN((*lrhs),    (*lN),     PASTIX_FLOAT);
  if (gperm != NULL)
    MALLOC_EXTERN((*lperm), (*lN), PASTIX_INT);

  MALLOC_EXTERN((*loc2glob), (*lN), PASTIX_INT);

  index    = 0;
  rowindex = 1;

  switch (dispatch)
    {
    case CSC_DISP_CYCLIC:
      {
        for (i = rank; i < gN; i+=comm_size)
          {
            (*loc2glob)[index] = i+1;
            (*lcolptr )[index] = rowindex;
            (*lrhs    )[index] = grhs[i];
            if (gperm != NULL)
              {
                (*lperm   )[index] = gperm[i];
                /*           (*linvp   )[index] = ginvp[i];  */
              }
            rowindex += gcolptr[i+1] - gcolptr[i];
            index++;
          }
      }
      (*lcolptr )[index] = rowindex;
      break;
    case CSC_DISP_SIMPLE:
    default:
      {
        PASTIX_INT ideb;
        PASTIX_INT ifin;
        ideb = rank * (gN / comm_size) + MIN(gN%comm_size, rank);
        ifin = ideb + (*lN);

        for (i = ideb; i < ifin; i++)
          {
            (*loc2glob)[index] = i+1;
            (*lcolptr )[index] = rowindex;
            (*lrhs    )[index] = grhs[i];
            if (gperm != NULL)
              {
                (*lperm   )[index] = gperm[i];
                /*           (*linvp   )[index] = ginvp[i];  */
              }
            rowindex += gcolptr[i+1] - gcolptr[i];
            index++;
          }
      }
    }

  ASSERT(index == (*lN), MOD_UNKNOWN);

  (*lcolptr )[*lN] = rowindex;
  if ((*lcolptr)[(*lN)]-1 > 0)
    {
      MALLOC_EXTERN(*lrow,   (*lcolptr)[(*lN)]-1, PASTIX_INT);
      MALLOC_EXTERN(*lavals, (*lcolptr)[(*lN)]-1, PASTIX_FLOAT);
    }
  else
    {
      *lrow   = NULL;
      *lavals = NULL;
    }

  for (i = 0; i < *lN; i++)
    {
      PASTIX_INT iG = (*loc2glob)[i];

      memcpy(&(*lrow)[(*lcolptr)[i]-1],
             &grow[gcolptr[iG-1]-1],
             (gcolptr[iG] - gcolptr[iG-1])* sizeof(PASTIX_INT));
      memcpy(&(*lavals)[(*lcolptr)[i]-1],
             &gavals[gcolptr[iG-1]-1],
             (gcolptr[iG] - gcolptr[iG-1])* sizeof(PASTIX_FLOAT));
    }
}

/*
 * Function: csc_cyclic_distribution
 *
 * Distribute the CSC cyclicaly.
 *
 * Parameters:
 *   column      - column number to distribute
 *   columnnbr   - Number of colmuns.
 *   pastix_comm - PaStiX MPI communicator
 *
 * Return:
 *   owner of the column (column%commSize)
 */
#ifndef csc_cyclic_distribution
#  error "redefinition required"
#endif
PASTIX_INT csc_cyclic_distribution(PASTIX_INT column, PASTIX_INT columnnbr, MPI_Comm pastix_comm)
{
  int commSize;
  (void)pastix_comm;
  (void)columnnbr;

  MPI_Comm_size(pastix_comm,&commSize);

  return column%commSize;
}

/*
 * Function: csc_simple_distribution
 *
 * Distribute the CSC.
 * First columns are for first proc and so on.
 *
 * Parameters:
 *   column      - column number to distribute
 *   columnnbr   - Number of colmuns.
 *   pastix_comm - PaStiX MPI communicator
 *
 * Return:
 *   owner of the column (column*commSize/columnnbr)
 */
PASTIX_INT csc_simple_distribution(PASTIX_INT column, PASTIX_INT columnnbr, MPI_Comm pastix_comm)
{
  int commSize;
  (void)pastix_comm;

  MPI_Comm_size(pastix_comm,&commSize);

  return (PASTIX_INT)floor((double)column*(double)commSize/(double)columnnbr);
}



/*
 * Function: cscd_build_g2l
 *
 * Construct global to local tabular containing local number of global columns
 * if one column is local, and -owner if column is not local.
 *
 * For i in 0, gN
 *    g2l[i] = i local number if i is local
 *    g2l[i] = -p if p is the owner of the column i
 *
 * Parameters:
 *   n        - Number of local columns
 *   colptr   - Starting index of each columns in *ja*
 *   rows     - Row of each element.
 *   values   - Value of each element.
 *   l2g      - global number of each local column.
 *   correct  - Flag indicating if we can correct the symmetry.
 *   dof      - Number of degrees of freedom.
 *   comm     - MPI communicator
 */
#ifndef cscd_build_g2l
#  error "redefinition required"
#endif
int cscd_build_g2l(PASTIX_INT       ncol,
                   PASTIX_INT      *loc2glob,
                   MPI_Comm  comm,
                   PASTIX_INT      *gN,
                   PASTIX_INT     **g2l)
{
  int commRank;
  int commSize;
  int i,j;
  PASTIX_INT ig;
  PASTIX_INT nrecv;
  PASTIX_INT *l2grecv = NULL;
  (void)comm;

  MPI_Comm_rank(comm, &commRank);
  MPI_Comm_size(comm, &commSize);

  if (*gN == -1) /* Needed for MURGE_ProductSetGlobalNodeNbr */
    MPI_Allreduce(&ncol, gN, 1, COMM_INT, MPI_SUM, comm);

  MALLOC_INTERN(*g2l, *gN, PASTIX_INT);
  for(ig=0; ig < *gN; ig++)
    (*g2l)[ig] = -commSize;

  for (i = 0; i < commSize; i++)
    {
      if (commRank == i)
        {
          nrecv   = ncol;
          l2grecv = loc2glob;

          MPI_Bcast(&nrecv,      1, COMM_INT, i, comm);
          if (nrecv > 0)
            {
              MPI_Bcast(l2grecv, nrecv, COMM_INT, i, comm);

              for(j = 0; j < nrecv; j++)
                {
                  (*g2l)[l2grecv[j]-1] = j+1;
                }
            }
        }
      else
        {
          MPI_Bcast(&nrecv,      1, COMM_INT, i, comm);
          if (nrecv > 0)
            {
              MALLOC_INTERN(l2grecv, nrecv, PASTIX_INT);
              MPI_Bcast(l2grecv, nrecv, COMM_INT, i, comm);

              for(j = 0; j < nrecv; j++)
                {
                  /* With Murge Product, one column can be on 2 processors
                   * If so it's always marked in g2l as local on those processors.
                   * On other processors we don't care who it belongs too.
                   */
                  if ((*g2l)[l2grecv[j]-1] == -commSize)
                    (*g2l)[l2grecv[j]-1] = -i;
                }
              memFree_null(l2grecv);
            }
        }
    }
  for(ig=0; ig < *gN; ig++) {
    if ((*g2l)[ig] == -commSize)
      fprintf(stdout, "((*g2l)[%ld] = %ld\n", (long)ig, (long)((*g2l)[ig]));
    ASSERT((*g2l)[ig] != -commSize, MOD_SOPALIN);
  }

  return EXIT_SUCCESS;
}

/*
 * Function: cscd_checksym
 *
 * Check if the CSCD graph is symetric.
 *
 * For all local column C,
 *
 * For all row R in the column C,
 *
 * If the row number R correspond to a local column,
 * We look in column R if we have the row number C.
 *
 * Else, we verify the column R belong to one proc
 *
 * We had (ROW, l2g(COL)) (Fortran numbering) to
 * list of the couples to send to the owner.
 *
 * If we are not allowed to correct the matrix,
 * we check that we have a local symmetry and return EXIT_FAILURE
 * if it was not good on one processor.
 *
 * If Local symmetry is good, or we are allowed to correct the matrix,
 * we exchange couples (J,I) that processors should have and
 * we check that this couples are present.
 *
 * Parameters:
 *   n        - Number of local columns
 *   colptr   - Starting index of each columns in *ja*
 *   rows     - Row of each element.
 *   values   - Value of each element.
 *   l2g      - global number of each local column.
 *   correct  - Flag indicating if we can correct the symmetry.
 *   alloc    - indicate if allocation on CSC uses internal malloc.
 *   dof      - Number of degrees of freedom.
 *   comm     - MPI communicator
 */
#ifndef cscd_checksym
#  error "redefinition required"
#endif
int cscd_checksym(PASTIX_INT      n,
                  PASTIX_INT     *colptr,
                  PASTIX_INT    **rows,
                  PASTIX_FLOAT  **values,
                  PASTIX_INT     *l2g,
                  int      correct,
                  int      alloc,
                  int      dof,
                  MPI_Comm comm)
{
  int            commSize;
  int            proc;
  int            commRank;
  PASTIX_INT            gN         = -1;
  PASTIX_INT            i,j,k,l;
  int            found;
  int            sym;
  int            all_sym;
  PASTIX_INT         *  g2l        = NULL;
  PASTIX_INT         *  nbtosend   = NULL;
  PASTIX_INT            toaddsize[2];
  PASTIX_INT         *  toadd      = NULL;
  PASTIX_INT         *  tmpcolptr  = NULL;
  PASTIX_INT         *  tmprows    = NULL;
  PASTIX_INT         *  newcolptr  = NULL;
  PASTIX_INT         *  newrows    = NULL;
  PASTIX_FLOAT       *  newvalues  = NULL;
  PASTIX_INT         *  tosendsize = NULL;
  PASTIX_INT         ** tosend     = NULL;
  PASTIX_INT         *  torecvsize = NULL;
  PASTIX_INT         ** torecv     = NULL;
  MPI_Request *  requests   = NULL;

#ifndef FORCE_NOMPI
  MPI_Status     status;
#endif
  PASTIX_INT            column;

  MPI_Comm_size(comm,&commSize);
  MPI_Comm_rank(comm,&commRank);

  /* Build the g2l vector */
  cscd_build_g2l(n, l2g, comm, &gN, &g2l);

  MALLOC_INTERN(nbtosend,   commSize, PASTIX_INT);
  MALLOC_INTERN(tosendsize, commSize, PASTIX_INT);
  MALLOC_INTERN(tosend,     commSize, PASTIX_INT*);
  for (i = 0; i < commSize; i++)
    {
      nbtosend[i]   = 0;
      tosendsize[i] = n / (2*commSize);
      if (i != commRank)
        MALLOC_INTERN(tosend[i], 2*tosendsize[i], PASTIX_INT);
    }

  toaddsize[0] = 0;
  toaddsize[1] = 0;

  /* For all local column C,
     For all row R in the column C,

     If the row number R correspond to a local column,
     We look in column R if we have the row number C.

     Else, we verify the column R belong to one proc

     We had (ROW, l2g(COL)) (Fortran numbering) to
     list of the couples to send to the owner.
  */
  sym = 1;
  for (i = 0; i < n && sym == 1; i++)
    for (j = colptr[i]-1; (j < colptr[i+1]-1) && (sym == 1); j++)
      {
        /* look if *rows[j] column is local */
        k = g2l[(*rows)[j]-1];
        if (k > 0)
          {
            found = 0;
            for (l = (colptr)[k-1]-1; l < (colptr)[k-1+1]-1; l++)
              {
                if (l2g[i] == (*rows)[l])
                  {
                    found = 1;
                    break;
                  }
              }
            if (found == 0)
              {
                if (correct == API_NO)
                  sym = 0;
                else
                  {
                    if (toaddsize[1] == 0)
                      {
                        toaddsize[1] = n/(2*commSize);
                        MALLOC_INTERN(toadd, 2*toaddsize[1], PASTIX_INT);
                      }
                    if (toaddsize[0] >= toaddsize[1])
                      {
                        toaddsize[1] += toaddsize[1]/2 + 1;
                        if (NULL ==
                            (toadd =
                             (PASTIX_INT*)memRealloc(toadd,
                                              2*toaddsize[1]*sizeof(PASTIX_INT))))
                          MALLOC_ERROR("toadd");

                      }
                    toadd[2*toaddsize[0]]     = (*rows)[j];
                    toadd[2*toaddsize[0] + 1] = l2g[i];
                    /*  fprintf(stdout, "Adding %ld %ld %ld, %ld\n",
                        (long)l2g[i], (long) i, (long) n, (long)(*rows)[j]); */
                    toaddsize[0]++;
                  }
              }
          }
        /* k is not local */
        else
          {
            proc = -k;
            if (proc >= commSize)
              errorPrint("Missing column in CSCD (%ld)", (long)((*rows)[j]));

            /* owner is processor proc */
            if (nbtosend[proc] >= tosendsize[proc])
              {
                tosendsize[proc] += tosendsize[proc]/2 +1;
                if (NULL ==
                    (tosend[proc] =
                     (PASTIX_INT*)memRealloc(tosend[proc],
                                      2*tosendsize[proc]*sizeof(PASTIX_INT))))
                  MALLOC_ERROR("tosend[proc]");
                if (nbtosend[proc] >= tosendsize[proc])
                  errorPrint("cscd_checksym : should'nt happen");
              }
            tosend[proc][2*nbtosend[proc]]   = (*rows)[j];
            tosend[proc][2*nbtosend[proc]+1] = l2g[i];
            nbtosend[proc]++;
          }
      }

  memFree_null(tosendsize);
  /*
   * We check that we have a local symmetry.
   */
  all_sym = 0;
  MPI_Allreduce(&sym,&all_sym,1,MPI_INT,MPI_SUM,comm);
  if (all_sym != commSize) {
    memFree_null(nbtosend);
    for (i = 0; i < commSize; i++) {
      if (i != commRank) {
        memFree_null(tosend[i]);
      }
    }
    memFree_null(tosend);
    return EXIT_FAILURE;
  }
  /*
   * If Local symmetry id good
   */
  MALLOC_INTERN(torecvsize, commSize, PASTIX_INT);
  MALLOC_INTERN(requests, commSize, MPI_Request);

  /* We Send and receive list of couples each processor should have */
  for (i = 0; i < commSize; i++) {
    if (i != commRank) {
      MPI_Isend(&nbtosend[i], 1, COMM_INT, i, i+100*commRank,
                comm, &requests[i]);
    }
  }
  for (i = 0; i < commSize; i++) {
    if (i != commRank) {
      MPI_Recv(&torecvsize[i], 1, COMM_INT, i, commRank + i * 100,
               comm, &status);
    }
  }
  for (i = 0; i < commSize; i++) {
    if (i != commRank) {
      MPI_Wait(&requests[i], &status);
    }
  }

  MALLOC_INTERN(torecv, commSize, PASTIX_INT*);

  for (i = 0; i < commSize; i++) {
    torecv[i] = NULL;
    if (i != commRank) {
      MALLOC_INTERN(torecv[i], torecvsize[i]*2, PASTIX_INT);
    }
  }
  for (i = 0; i < commSize; i++) {
    if (i != commRank) {
      MPI_Isend(tosend[i], 2*nbtosend[i], COMM_INT, i , i+100*commRank,
                comm, &requests[i]);
    }
  }
  for (i = 0; i < commSize; i++) {
    if (i != commRank) {
      MPI_Recv(torecv[i], 2*torecvsize[i], COMM_INT, i, commRank + i * 100,
               comm, &status);
    }
  }
  for (i = 0; i < commSize; i++) {
    if (i != commRank) {
      MPI_Wait(&requests[i], &status);
    }
  }

  /* Free sent datas and requests */
  memFree_null(requests);
  memFree_null(nbtosend);
  for (i = 0; i < commSize; i++) {
    if (i != commRank) {
      memFree_null(tosend[i]);
    }
  }
  memFree_null(tosend);

  /* Check that received couples are present */
  for (i = 0; i < commSize && sym == 1; i++) {
    if (i != commRank && sym == 1)
      {
        for (j = 0; j < torecvsize[i]; j++)
          {
            /* look for column torecv[i][2*j] */
            column = g2l[torecv[i][2*j]-1];
            if (column < 0)
              {
                errorPrint("Shoulnd't happen");
                sym = 0;
                break;
              }

            /* look for row torecv[i][2*j+1] in column k */
            found = 0;
            k = column -1;
            for (l = (colptr)[k]-1; l < (colptr)[k+1]-1; l++)
              {
                if ((*rows)[l] == torecv[i][2*j+1])
                  {
                    found =1;
                    break;
                  }
              }

            if (found == 0)
              {
                if (correct == API_NO)
                  {
                    sym = 0;
                    break;
                  }
                else
                  {
                    if (toaddsize[1] == 0)
                      {
                        toaddsize[1] = n/(2*commSize);
                        MALLOC_INTERN(toadd, 2*toaddsize[1], PASTIX_INT);
                      }
                    if (toaddsize[0] >= toaddsize[1])
                      {
                        toaddsize[1] += toaddsize[1]/2 + 1;
                        if (NULL ==
                            (toadd =
                             (PASTIX_INT*)memRealloc(toadd,
                                              2*toaddsize[1]*sizeof(PASTIX_INT))))
                          MALLOC_ERROR("toadd");

                      }
                    toadd[2*toaddsize[0]]     = torecv[i][2*j];
                    toadd[2*toaddsize[0] + 1] = torecv[i][2*j+1];
                    toaddsize[0]++;
                  }

              }
          }
      }
  }
  memFree_null(torecvsize);
  for (i = 0; i < commSize; i++)
    if (i != commRank && torecv[i] != NULL )
      memFree_null(torecv[i]);
  memFree_null(torecv);


  all_sym = 0;
  MPI_Allreduce(&sym,&all_sym,1,MPI_INT,MPI_SUM,comm);
  if (all_sym != commSize)
    return EXIT_FAILURE;

  if (correct == API_YES && toaddsize[0] > 0)
    {

      MALLOC_INTERN(tmpcolptr, n + 1, PASTIX_INT);

      /* Build tmpcolptr

         tmpcolptr[i+1] will contain the number of element of
         the column i
      */
      for (i = 0; i <  n+1; i++)
        tmpcolptr[i] = 0;
      for (i = 0; i <  toaddsize[0]; i++)
        {
          tmpcolptr[g2l[toadd[2*i]-1]]++;
        }

      /* tmpcolptr[i] will contain the index of
         first element of column i
      */
      tmpcolptr[0] = 1;
      for (i = 0; i <  n; i++)
        {
          tmpcolptr[i+1] = tmpcolptr[i] + tmpcolptr[i+1];
        }

      if (tmpcolptr[n] - 1 != toaddsize[0])
        {
          errorPrint("Wrong CSCd construction. %ld %ld",
                     (long) tmpcolptr[n], (long)toaddsize[0]);
          return EXIT_FAILURE;
        }

      MALLOC_INTERN(tmprows, tmpcolptr[n] - 1, PASTIX_INT);
      for (i = 0; i <  toaddsize[0]; i++)
        {
          tmprows[tmpcolptr[g2l[toadd[2*i]-1]-1]-1] = toadd[2*i+1];
          tmpcolptr[g2l[toadd[2*i]-1]-1]++;
        }

      memFree_null(toadd);
      /* restore tmpcolptr */
      for (i = 0; i <  n; i++)
        {
          tmpcolptr[n - i] = tmpcolptr[n - i - 1];
        }
      tmpcolptr[0] = 1;

      if (tmpcolptr[n] - 1 != toaddsize[0])
        {
          errorPrint("Wrong CSCd construction %ld %ld",
                     (long) tmpcolptr[n], (long)toaddsize[0]);
          return EXIT_FAILURE;
        }

      /* On ajoute les 2 cscd, Allocation externe */
      if (values != NULL)
        {
          cscd_addlocal_int(n,   colptr,   *rows   , *values   ,  l2g,
                            n,   tmpcolptr,  tmprows,  NULL,      l2g,
                            &n, &newcolptr, &newrows, &newvalues,
                            &add_two_floats,
                            dof, alloc);
        }
      else
        {
          cscd_addlocal_int(n,   colptr,   *rows   ,   NULL,    l2g,
                            n,   tmpcolptr,  tmprows,  NULL,    l2g,
                            &n, &newcolptr, &newrows,  NULL,
                            &add_two_floats, dof, alloc);
        }
      memFree_null(tmpcolptr);
      memFree_null(tmprows);
      memcpy(colptr, newcolptr, (n+1)*sizeof(PASTIX_INT));
      if (alloc == API_NO)
        {
          free(newcolptr);
          free(*rows);
        }
      else
        {
          memFree_null(newcolptr);
          memFree_null(*rows);
        }
      if (values != NULL)
        {
          if (alloc == API_NO)
            {
              free(*values);
            }
          else
            {
              memFree_null(*values);
            }
        }
      *rows   = newrows;
      if (values != NULL)
        *values = newvalues;

    }
  memFree_null(g2l);
  return EXIT_SUCCESS;
}


/*
 * Function: cscd_symgraph
 *
 * Symetrize the graph of the given CSCD.
 *
 * External function.
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - Starting index of each columns in *ja* and *a*
 *   ja          - Row of each element.
 *   a           - Values of each element.
 *   newn        - New number of local columns
 *   newia       - Starting index of each columns in *newja* and *newa*
 *   newja       - Row of each element.
 *   newa        - Values of each element.
 *   l2g         - global number of each local column.
 *   comm        - MPI communicator.
 */
#ifndef cscd_symgraph
#  error "redefinition required"
#endif
int cscd_symgraph(PASTIX_INT      n, PASTIX_INT *      ia, PASTIX_INT *        ja, PASTIX_FLOAT *     a,
                  PASTIX_INT * newn, PASTIX_INT **  newia, PASTIX_INT **    newja, PASTIX_FLOAT ** newa,
                  PASTIX_INT *  l2g, MPI_Comm comm)
{
  return cscd_symgraph_int(n,    ia,    ja,    a,
                           newn, newia, newja, newa,
                           l2g,  comm,  API_NO);
}

/*
 * Function: cscd_symgraph_int
 *
 * Symetrize the graph of the given CSCD.
 *
 * Internal function.
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - Starting index of each columns in *ja* and *a*
 *   ja          - Row of each element.
 *   a           - Values of each element.
 *   newn        - New number of local columns
 *   newia       - Starting index of each columns in *newja* and *newa*
 *   newja       - Row of each element.
 *   newa        - Values of each element.
 *   l2g         - global number of each local column.
 *   comm        - MPI communicator.
 *   malloc_flag - flag to indicate if function call is intern to pastix
 *                or extern.
 */

#ifndef cscd_symgraph_int
#  error "redefinition required"
#endif
int cscd_symgraph_int(PASTIX_INT      n, PASTIX_INT *      ia, PASTIX_INT *        ja, PASTIX_FLOAT *     a,
                      PASTIX_INT * newn, PASTIX_INT **  newia, PASTIX_INT **    newja, PASTIX_FLOAT ** newa,
                      PASTIX_INT *  l2g, MPI_Comm comm, int malloc_flag)
{
  int            commSize;
  int            proc;
  int            commRank;
  PASTIX_INT            gN     = -1;
  PASTIX_INT            i,j,k,l;
  int            found;
  PASTIX_INT         *  nbtosend   = NULL;
  PASTIX_INT         *  tosendsize = NULL;
  PASTIX_INT         ** tosend     = NULL;
  PASTIX_INT         *  torecvsize = NULL;
  PASTIX_INT         ** torecv     = NULL;
  MPI_Request *  requests   = NULL;
#ifndef FORCE_NOMPI
  MPI_Status     status;
#endif
  PASTIX_INT         *  toadd      = NULL;
  PASTIX_INT         *  added      = NULL;
  PASTIX_INT            tmpn;
  PASTIX_INT         *  tmpia      = NULL;
  PASTIX_INT         *  tmpja      = NULL;
  PASTIX_INT         *  tmpl2g     = NULL;
  PASTIX_INT         *  g2l        = NULL;
  PASTIX_INT            colnum;
  PASTIX_INT            rownum;

  MPI_Comm_size(comm,&commSize);
  MPI_Comm_rank(comm,&commRank);

  print_debug(DBG_CSCD,"->cscd_symgraph\n");

  /* Build the g2l vector */
  cscd_build_g2l(n, l2g, comm, &gN, &g2l);

  MALLOC_INTERN(nbtosend,   commSize, PASTIX_INT);
  MALLOC_INTERN(tosendsize, commSize, PASTIX_INT);
  MALLOC_INTERN(tosend,     commSize, PASTIX_INT*);
  for (i = 0; i < commSize; i++)
    {
      nbtosend[i]   = 0;
      tosendsize[i] = n/(2*commSize);
      if (i != commRank)
        MALLOC_INTERN(tosend[i], 2*tosendsize[i], PASTIX_INT);
    }

  MALLOC_INTERN(toadd, n, PASTIX_INT);
  for (i = 0; i < n; i++)
    toadd[i] = 0;

  /* For all local column C,
     For all row R in the column C,

     If the row number R correspond to a local column,
     We look in column R if we have the row number C.

     Else, we verify the column R belong to one proc

     We had (ROW, l2g(COL)) (Fortran numbering) to
     list of the couples to send to the owner.
  */
  /* Pour chaque colonne i, pour chaque ligne ja[j]*/
  for (i = 0; i < n; i++)
    for (j = ia[i]-1; j < ia[i+1]-1; j++)
      {
        /* look if ja[j] column is local */
        k = g2l[ja[j]-1];
        if (k > 0) /* local */
          {
            /* si ja[j] correspond à une colonne locale
               on parcours les lignes de la colonne qui correspond à ja[j]
               si on ne trouve pas la colonne global[i], il faut l'ajouter.
            */
            k     = k-1;
            found = 0;

            for (l = ia[k]-1; l < ia[k+1]-1; l++)
              {
                if (l2g[i] == ja[l])
                  {
                    found = 1;
                    break; /* Back to the j loop */
                  }
              }
            if (found == 0)
              {
                toadd[k] ++;
              }
          }
        else
          {
            /*
              si ja[j] n'est pas locale, on envoie
              ja[j], l2g[i] au processeur qui la possède.
            */
            proc = -k;
            if (proc >= commSize)
              errorPrint("Missing column in CSCD (%ld)", (long)(ja[j]));

            if (nbtosend[proc] >= tosendsize[proc])
              {
                tosendsize[proc] += tosendsize[proc]/2 +1;
                if (NULL ==
                    (tosend[proc] =
                     (PASTIX_INT*)memRealloc(tosend[proc],
                                      2*tosendsize[proc]*sizeof(PASTIX_INT))))
                  MALLOC_ERROR("tosend[proc]");
                if (nbtosend[proc] >= tosendsize[proc])
                  errorPrint("cscd_checksym : should'nt happen");
              }
            tosend[proc][2*nbtosend[proc]]   = ja[j];
            tosend[proc][2*nbtosend[proc]+1] = l2g[i];
            nbtosend[proc]++;
          }
      }

  memFree_null(tosendsize);

  MALLOC_INTERN(torecvsize, commSize, PASTIX_INT);
  MALLOC_INTERN(requests,   commSize, MPI_Request);

  for (i = 0; i < commSize; i++) {
    if (i != commRank) {
      MPI_Isend(&nbtosend[i],         /* buffer            */
                1, COMM_INT, i,       /* count, type, dest */
                TAG_SIZE2,            /* tag               */
                comm, &requests[i]);  /* comm, request     */
    }
  }

  for (i = 0; i < commSize; i++) {
    if (i != commRank) {
      MPI_Recv(&torecvsize[i],        /* buffer            */
               1, COMM_INT, i,                /* count, type, dest */
               TAG_SIZE2,             /* tag               */
               comm, &status);        /* comm, request     */
    }
  }

  for (i = 0; i < commSize; i++) {
    if (i != commRank) {
      MPI_Wait(&requests[i], &status);
    }
  }

  MALLOC_INTERN(torecv, commSize, PASTIX_INT*);

  for (i = 0; i < commSize; i++) {
    if (i != commRank && torecvsize[i] > 0) {
      MALLOC_INTERN(torecv[i], torecvsize[i]*2, PASTIX_INT);
    }
  }

  for (i = 0; i < commSize; i++) {
    if (i != commRank && nbtosend[i] > 0) {
      MPI_Isend(tosend[i],            /* buffer            */
                2*nbtosend[i],        /* count             */
                COMM_INT, i ,         /* type, dest        */
                TAG_COUPLE,           /* tag               */
                comm, &requests[i]);  /* comm, request     */
    }
  }

  for (i = 0; i < commSize; i++) {
    if (i != commRank && torecvsize[i] > 0) {
      MPI_Recv(torecv[i],             /* buffer            */
               2*torecvsize[i],       /* count             */
               COMM_INT, i,           /* type, dest        */
               TAG_COUPLE,            /* tag               */
               comm, &status);        /* comm, request     */
    }
  }

  for (i = 0; i < commSize; i++) {
    if (i != commRank)
      {
        if (nbtosend[i] > 0)
          {
            MPI_Wait(&requests[i], &status);
          }
        memFree_null(tosend[i]);
      }
  }
  memFree_null(requests);
  /* on ajoute les elements recus dans les elements a ajouter */
  for (i = 0; i < commSize ; i++){
    if (i != commRank)
      for (j = 0; j < torecvsize[i]; j++)
        {
          colnum = g2l[torecv[i][2*j]-1];
          toadd[colnum-1] ++;
        }
  }

  memFree_null(nbtosend);
  memFree_null(tosend);

  tmpn  = n;
  MALLOC_INTERN(tmpia, tmpn+1, PASTIX_INT);
  tmpia[0] = 1;
  for (i = 0; i < n ; i++) {
    tmpia[i+1] = tmpia[i] + toadd[i];
  }
  memFree_null(toadd);

  MALLOC_INTERN(tmpja, tmpia[tmpn]-1, PASTIX_INT);

  MALLOC_INTERN(added, n, PASTIX_INT);
  for (i = 0; i < n; i++)
    added[i] = 0;

  /* On ajoute les elements recus */
  tmpl2g = l2g;
  for (i = 0; i < commSize ; i++){
    if (i != commRank) {
      for (j = 0; j < torecvsize[i]; j++)
        {
          colnum = g2l[torecv[i][2*j]-1];
          rownum =        torecv[i][2*j+1];
          tmpja[tmpia[colnum-1]-1+added[colnum -1]] = rownum;
          added[colnum-1] ++;
        }
    }
  }

  for (i = 0; i < commSize; i++) {
    if (i != commRank && torecvsize[i] > 0) {
      memFree_null(torecv[i]);
    }
  }
  memFree_null(torecv);
  memFree_null(torecvsize);

  /* on ajoute les symetriques locaux */
  for (i = 0; i < n; i++) {
    for (j = ia[i]-1; j < ia[i+1]-1; j++)
      {

        /* look if ja[j] column is local */
        k = g2l[ja[j]-1];

        if (k > 0) /* local */
          {
            found = 0;
            for (l = ia[k-1]-1; l < ia[k]-1; l++)
              {
                if (l2g[i] == ja[l])
                  {
                    found = 1;
                    break;
                  }
              }
            if (found == 0)
              {
                colnum = k-1;
                rownum = l2g[i];
                tmpja[tmpia[colnum]-1+added[colnum]] = rownum;
                added[colnum] ++;
              }
          }
      }
  }
  memFree_null(g2l);
  memFree_null(added);
  /* On ajoute les 2 cscd */
  cscd_addlocal_int(n   , ia   , ja   , a   , l2g,
                    tmpn, tmpia, tmpja, NULL, tmpl2g,
                    newn, newia, newja, newa, &add_two_floats, 1, malloc_flag);
  memFree_null(tmpia);
  memFree_null(tmpja);

  print_debug(DBG_CSCD,"<-cscd_symgraph\n");

  return EXIT_SUCCESS;
}

/*
 * Function: cscd_addlocal
 *
 * Interface to <cscd_addlocal_int>, for external usage.
 * Add second cscd to first cscd into third cscd (unallocated)
 * Only adds columns from the second CSCD which belongs to the first one.
 *
 * Parameters:
 *   n           - First cscd size
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   l2g         - local 2 global column numbers for first cscd
 *   addn        - CSCD to add size
 *   addia       - CSCD to add starting index of each column in
 *                 *addja* and *adda*
 *   addja       - Row of each element in second CSCD
 *   adda        - value of each cscd in second CSCD (can be NULL -> add 0)
 *   addl2g      - local 2 global column numbers for second cscd
 *   newn        - new cscd size (same as first)
 *   newia       - CSCD to add starting index of each column in *newja*
 *                 and *newwa*
 *   newja       - Row of each element in third CSCD
 *   newa        - value of each cscd in third CSCD
 *   OP          - Operation to manage common CSCD coefficients.
 *   dof         - Number of degrees of freedom.
 */
#ifndef cscd_addlocal
#  error "redefinition required"
#endif
int cscd_addlocal(PASTIX_INT   n   , PASTIX_INT *  ia   , PASTIX_INT *  ja   , PASTIX_FLOAT *  a   ,
                  PASTIX_INT * l2g,
                  PASTIX_INT   addn, PASTIX_INT *  addia, PASTIX_INT *  addja, PASTIX_FLOAT *  adda,
                  PASTIX_INT * addl2g,
                  PASTIX_INT * newn, PASTIX_INT ** newia, PASTIX_INT ** newja, PASTIX_FLOAT ** newa,
                  CSCD_OPERATIONS_t OP, int dof)
{

  switch(OP)
    {
    case CSCD_ADD:
      return cscd_addlocal_int(n,    ia,    ja,    a,    l2g,
                               addn, addia, addja, adda, addl2g,
                               newn, newia, newja, newa,
                               &add_two_floats, dof, API_NO);
    case CSCD_KEEP:
      return cscd_addlocal_int(n,    ia,    ja,    a,    l2g,
                               addn, addia, addja, adda, addl2g,
                               newn, newia, newja, newa, &keep_first,
                               dof, API_NO);
#ifndef TYPE_COMPLEX
    case CSCD_MAX:
      return cscd_addlocal_int(n,    ia,    ja,    a,    l2g,
                               addn, addia, addja, adda, addl2g,
                               newn, newia, newja, newa, &get_max, dof, API_NO);
    case CSCD_MIN:
      return cscd_addlocal_int(n,    ia,    ja,    a,    l2g,
                               addn, addia, addja, adda, addl2g,
                               newn, newia, newja, newa, &get_min, dof, API_NO);
#endif
    case CSCD_OVW:
      return cscd_addlocal_int(n,    ia,    ja,    a,    l2g,
                               addn, addia, addja, adda, addl2g,
                               newn, newia, newja, newa, &keep_last,
                               dof, API_NO);

    default:
      return EXIT_FAILURE;
    }
}



#define searchInList(n, list, list_size, index){  \
    PASTIX_INT min = 0;                                  \
    PASTIX_INT max = list_size-1;                        \
    if  (list_size == 0) {                        \
      index = -1;                                 \
    }                                             \
    else {                                        \
      while (1) {                                 \
        PASTIX_INT cursor = (max +min)/2;                \
        if ((list)[cursor] == n)                  \
          {                                       \
            index = cursor;                       \
            break;                                \
          }                                       \
        else {                                    \
          if (max <= min) {                       \
            index = -1;                           \
            break;                                \
          }                                       \
          if ((list)[cursor] > n)                 \
            {                                     \
              max = cursor-1;                     \
            }                                     \
          else {                                  \
            if ((list)[cursor] < n)               \
              {                                   \
                min = cursor+1;                   \
              }                                   \
          }                                       \
        }                                         \
      }                                           \
    }                                             \
  }



/*
 * Function: cscd_addlocal_int
 *
 * Add second cscd to first cscd into third cscd (unallocated)
 * Only adds columns from the second CSCD which belongs to the first one.
 * All this work is localy done.
 *
 * Parameters:
 *   n           - First cscd size
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   l2g         - local 2 global column numbers for first cscd
 *   addn        - CSCD to add size
 *   addia       - CSCD to add starting index of each column in *addja*
 *                 and *adda*
 *   addja       - Row of each element in second CSCD
 *   adda        - value of each cscd in second CSCD (can be NULL -> add 0)
 *   addl2g      - local 2 global column numbers for second cscd
 *   newn        - new cscd size (same as first)
 *   newia       - CSCD to add starting index of each column in *newja*
 *                 and *newwa*
 *   newja       - Row of each element in third CSCD
 *   newa        - value of each cscd in third CSCD
 *   add_fct     - Fonction indicating what to do if a coefficient is in
 *                 the two cscd.
 *   dof         - Number of degrees of freedom.
 *   malloc_flag - flag to indicate if function call is intern to pastix
 *                 or extern.
 */

#ifndef cscd_addlocal_int
#  error "redefinition required"
#endif
int cscd_addlocal_int(PASTIX_INT   n   , PASTIX_INT *  ia   , PASTIX_INT *  ja   , PASTIX_FLOAT *  a   ,
                      PASTIX_INT * l2g,
                      PASTIX_INT   addn, PASTIX_INT *  addia, PASTIX_INT *  addja, PASTIX_FLOAT *  adda,
                      PASTIX_INT * addl2g,
                      PASTIX_INT * newn, PASTIX_INT ** newia, PASTIX_INT ** newja, PASTIX_FLOAT ** newa,
                      PASTIX_FLOAT (*add_fct)(PASTIX_FLOAT , PASTIX_FLOAT), int dof,
                      int malloc_flag)
{
  PASTIX_INT   i,j,k,l;
  PASTIX_INT   i2;
  int   baseval = ia[0];
  int   iterdof;
  int   val_flag;
  (*newia) = NULL;
  (*newja) = NULL;
  if (newa != NULL) (*newa = NULL);
  /* A can be null, Adda too but newa has to be filled with zeros...
   * => Test on newa
   */
  if (newa != NULL)
    val_flag = 1;
  else
    val_flag = 0;

  *newn = n;
  MALLOC_INTOREXTERN((*newia), (n+1), PASTIX_INT, malloc_flag);

  i2 = 0;
  (*newia)[0]=baseval;
  /* Construction of newia */
  for (i = 0; i < n; i++)
    {
      (*newia)[i+1] = (*newia)[i] + ia[i+1] - ia[i];
      /* On cherche si la colonne l2g[i] existe aussi dans addl2g,
         si oui on la comptabilise */
      searchInList(l2g[i], addl2g, addn, i2);
      if (i2>=0)
        {
          /* Add coefficient from addia wich are not in ia */
          for (j = addia[i2] - baseval; j < addia[i2+1] - baseval; j++)
            {
              PASTIX_INT index;
              PASTIX_INT *searchTab = &(ja[ia[i] - baseval]);
              PASTIX_INT searchTabSize = ia[i+1] - ia[i];
              searchInList(addja[j], searchTab, searchTabSize, index);
              if (index < 0)
                {
                  (*newia)[i+1]++;
                }
            }

        }
    }

  /* Allocation of newja */
  print_debug(DBG_CSCD, "allocation %ld\n", (long)(*newia)[n]-baseval);
  if (((*newia)[n]-baseval) > 0)
    MALLOC_INTOREXTERN((*newja), ((*newia)[n]-baseval), PASTIX_INT, malloc_flag);
  if (val_flag == 1) {
    MALLOC_INTOREXTERN((*newa), ((*newia)[n]-baseval)*dof*dof, PASTIX_FLOAT,
                       malloc_flag);

  }
  if (val_flag == 1) {
      memset(*newa,  0, ((*newia)[n]-baseval)*dof*dof*sizeof(PASTIX_FLOAT));
  }
  if (n > 0) {
    memset(*newja, 0, ((*newia)[n]-baseval)*sizeof(PASTIX_INT));
  }

  PASTIX_INT maxj,maxk,maxl;
  for (i = 0; i < addn; i++) {
    if(addia[i+1] - addia[i])
      intSort1asc1(&addja[addia[i]-baseval], addia[i+1] - addia[i]);
  }

  for (i = 0; i < n; i++)
    {
      searchInList(l2g[i], addl2g, addn, i2);
      if (i2>=0)
        {
          j = ia[i] - baseval;
          maxj = ia[i+1] - baseval;
          k = addia[i2] - baseval;
          maxk = addia[i2+1] - baseval;
          l = (*newia)[i] - baseval;
          maxl = (*newia)[i+1] - baseval;
          /* loop of size maxl <= maxk + makj */
          while (l < maxl)
            {
              /* while the new column is not filled */
              if ( j >= maxj ||
                   ( k < maxk &&
                     ja[j] > addja[k]))
                {
                  /* If nothing in ja
                     or addja is not empty and addja[k] must be before ja[j]
                  */

                  (*newja)[l] = addja[k];
                  if (val_flag == 1)
                    {
                      if (adda != NULL)
                        for (iterdof = 0; iterdof < dof*dof; iterdof ++)
                          {
                            (*newa)[l*dof*dof+iterdof] =
                              adda[k*dof*dof+iterdof];
                          }
                      else
                        for (iterdof = 0; iterdof < dof*dof; iterdof ++)
                          {
                            (*newa)[l*dof*dof+iterdof] = 0.0;
                          }
                    }
                  l++;
                  k++;
                }
              else
                {

                  /* else, There are element in ja
                     and addja is empty or ja[j] must
                     be before addja[j]
                  */
                  (*newja)[l] = ja[j];
                  if (val_flag == 1)
                    for (iterdof = 0; iterdof < dof*dof; iterdof ++)
                      {
                        if (a != NULL)
                          (*newa)[l*dof*dof+iterdof] = a[j*dof*dof+iterdof];
                        else
                          (*newa)[l*dof*dof+iterdof] = 0.0;
                      }
                  /* coef in both CSCs */
                  if (addja != NULL && k < maxk && ja[j] == addja[k])
                    {
                      if (val_flag == 1)
                        for (iterdof = 0; iterdof < dof*dof; iterdof ++)
                          {
                            (*newa)[l*dof*dof+iterdof] =
                              add_fct((*newa)[l*dof*dof+iterdof],
                                      adda[k*dof*dof+iterdof]);
                          }
                      k++;
                    }
                  l++;
                  j++;
                }
            }
        }
      else
        for (j = 0; j < ia[i+1] - ia[i]; j++)
          {
            (*newja)[(*newia)[i]-1+j] = ja[ia[i]-1+j];
            if (val_flag == 1)
              for (iterdof = 0; iterdof < dof*dof; iterdof ++)
                {
                  if (a != NULL)
                    (*newa)[((*newia)[i]-1+j)*dof*dof+iterdof] =
                      a[(ia[i]-1+j)*dof*dof+iterdof];
                  else
                    (*newa)[((*newia)[i]-1+j)*dof*dof+iterdof] = 0.0;
                }
          }
    }

  return EXIT_SUCCESS;
}


/*
 * Function: cscd_redispatch_scotch
 *
 * Redistribute the columns to have first columns on first proc for Scotch
 *
 * - Checks if the matrix is not allready well ditributed.
 * - Build all new loc2globs on all processors.
 * - Create new CSC for each processor and send it.
 * - Merge all new CSC to the new local CSC with <cscd_addlocal_int>.
 *
 *
 * UNUSED
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   l2g         - local 2 global column numbers for first cscd
 *   dn          - Number of local columns
 *   dia         - First cscd starting index of each column in *ja* and *a*
 *   dja         - Row of each element in first CSCD
 *   da          - value of each cscd in first CSCD (can be NULL)
 *   l2g         - local 2 global column numbers for first cscd
 *   comm        - MPI communicator
 *
 * Returns:
 *   EXIT_SUCCESS if already well distributed, 2 if redistributed
 *
 */

#ifndef cscd_redispatch_scotch
#  error "redefinition required"
#endif
int cscd_redispatch_scotch(PASTIX_INT   n, PASTIX_INT *   ia, PASTIX_INT *   ja, PASTIX_FLOAT *   a,
                           PASTIX_INT *   l2g,
                           PASTIX_INT *dn, PASTIX_INT ** dia, PASTIX_INT ** dja, PASTIX_FLOAT ** da,
                           PASTIX_INT ** dl2g,
                           MPI_Comm comm)
{
#ifdef FORCE_NOMPI
  (void)n; (void)ia; (void)ja; (void)a; (void)l2g;
  (void)dn; (void)dia; (void)dja; (void)da; (void)dl2g;
  (void)comm;
  return EXIT_SUCCESS;
#else /* FORCE_NOMPI */
  PASTIX_INT i,j;
  int           OK = 0;
  int           OKRecv = 0;
  PASTIX_INT           gn;
  PASTIX_INT           tmpn;
  PASTIX_INT         * tmpia  = NULL;
  PASTIX_INT         * tmpja  = NULL;
  PASTIX_FLOAT       * tmpa   = NULL;
  PASTIX_INT         * tmpia2 = NULL;
  PASTIX_INT         * tmpja2 = NULL;
  PASTIX_FLOAT       * tmpa2  = NULL;
  PASTIX_INT           start;
  PASTIX_INT           end;
  PASTIX_INT           size;
  PASTIX_INT         * newloc2globssize = NULL;
  PASTIX_INT        ** newloc2globs     = NULL;
  PASTIX_INT        ** colptr2send      = NULL;
  PASTIX_INT        ** rows2send        = NULL;
  PASTIX_FLOAT      ** values2send      = NULL;
  MPI_Request * requests_size    = NULL;
  MPI_Request * requests_col     = NULL;
  MPI_Request * requests_row     = NULL;
  MPI_Request * requests_val     = NULL;
  MPI_Status    status;
  MPI_Status    status2;
  PASTIX_INT           toreceive;
  int           commSize;
  int           rank;

  MPI_Comm_size(comm, &commSize);
  MPI_Comm_rank(comm, &rank);
  if (commSize == 1)
    goto simply_copy;

  /* Check if matrix is not allready correct for scotch */
  /* PT-Scotch needs consecutives column blocs */
  for (i = 0; i < n-1; i++)
    if (l2g[i] != l2g[i+1] -1)
      OK = 1;

  MPI_Allreduce(&OK, &OKRecv, 1, MPI_INT, MPI_SUM, comm);
  /* If it is correct, simply copy it */
  if (OKRecv == 0)
    {
    simply_copy:
      *dn = n;
      MALLOC_INTERN(*dia, n+1, PASTIX_INT);
      memcpy(*dia,ia,(n+1)*sizeof(PASTIX_INT));
      MALLOC_INTERN(*dl2g, n, PASTIX_INT);
      memcpy(*dl2g,l2g,n*sizeof(PASTIX_INT));
      MALLOC_INTERN(*dja, ia[n]-1, PASTIX_INT);
      memcpy(*dja,ja,(ia[n]-1)*sizeof(PASTIX_INT));
      if (NULL != a)
        {
          MALLOC_INTERN(*da, ia[n]-1, PASTIX_FLOAT);
          memcpy(*da,a,(ia[n]-1)*sizeof(PASTIX_FLOAT));
        }

      return EXIT_SUCCESS;
    }

  /* Build all newloc2globs on all processors */
  MALLOC_INTERN(newloc2globssize, commSize, PASTIX_INT );
  MALLOC_INTERN(newloc2globs,     commSize, PASTIX_INT*);
  gn = 0;
  MPI_Allreduce(&n, &gn, 1, COMM_INT, MPI_SUM, comm);
  for (i = 0; i < commSize; i++)
    {
      start = (PASTIX_INT)floor((double)i*(double)gn/(double)commSize)+1;
      end   = (PASTIX_INT)floor((double)(i+1)*(double)gn/(double)commSize);
      newloc2globssize[i] = end - start +1;

      MALLOC_INTERN(newloc2globs[i], newloc2globssize[i], PASTIX_INT);

      for (j = 0; j < newloc2globssize[i]; j++)
        newloc2globs[i][j] = j+start;
    }

  /* Create one CSCD for each proc */
  MALLOC_INTERN(colptr2send, commSize, PASTIX_INT*);
  MALLOC_INTERN(rows2send,   commSize, PASTIX_INT*);
  MALLOC_INTERN(values2send, commSize, PASTIX_FLOAT*);

  for (i = 0; i < commSize; i++)
    {
      colptr2send[i] = NULL;
      rows2send[i]   = NULL;
      values2send[i] = NULL;
    }

  MALLOC_INTERN(requests_size, commSize, MPI_Request);
  MALLOC_INTERN(requests_col,  commSize, MPI_Request);
  MALLOC_INTERN(requests_row,  commSize, MPI_Request);
  MALLOC_INTERN(requests_val,  commSize, MPI_Request);


  /* Sending informations to all proc */
  for (i = 0; i < commSize; i++)
    {
      MALLOC_INTERN(colptr2send[i], newloc2globssize[i]+1, PASTIX_INT);
      for (j = 0; j < newloc2globssize[i]+1; j++)
        colptr2send[i][j] = 1;

      /* Adding old CSCD to the CSCD to send */
      cscd_addlocal_int(newloc2globssize[i], colptr2send[i], rows2send[i],
                        values2send[i], newloc2globs[i],
                        n                  , ia            , ja          ,
                        a             , l2g,
                        &tmpn              , &tmpia        , &tmpja      ,
                        &tmpa         , &add_two_floats, 1, API_YES);

      if (newloc2globssize[i] != tmpn)
        errorPrint("newloc2globssize[i] != tmpn\n");
      memFree_null(colptr2send[i]);
      if (rows2send[i]   != NULL)
        memFree_null(rows2send[i]);
      if (values2send[i] != NULL)
        memFree_null(values2send[i]);

      colptr2send[i] = tmpia;
      rows2send[i]   = tmpja;
      values2send[i] = tmpa;

      tmpia = NULL;
      tmpja = NULL;
      tmpa  = NULL;

      /* Sending the CSCD built */
      if (i != rank)
        {
          size = colptr2send[i][newloc2globssize[i]]-1;
          MPI_Isend(&size, 1, COMM_INT, i, TAG_SIZE, comm, &requests_size[i]);
          if (size > 0)
            {
              MPI_Isend(colptr2send[i], newloc2globssize[i]+1, COMM_INT,    i,
                        TAG_COL, comm, &requests_col[i]);
              MPI_Isend(rows2send[i],   size,                  COMM_INT,    i,
                        TAG_ROW, comm, &requests_row[i]);
              if (a != NULL)
                MPI_Isend(values2send[i], size,                  COMM_FLOAT, i,
                          TAG_VAL, comm, &requests_val[i]);
            }
        }
    }

  toreceive = commSize -1;

  while(toreceive > 0)
    {
      MPI_Recv(&size, 1, COMM_INT, MPI_ANY_SOURCE, TAG_SIZE, comm, &status);
      if (size > 0)
        {
          /* Adding contribution CSCD to local CSCD */
          MALLOC_INTERN(tmpia, newloc2globssize[rank]+1, PASTIX_INT);
          MPI_Recv( tmpia, (int)(newloc2globssize[rank]+1), COMM_INT,
                    status.MPI_SOURCE, TAG_COL, comm, &status2);

          MALLOC_INTERN(tmpja, size, PASTIX_INT);
          MPI_Recv( tmpja, size, COMM_INT,  status.MPI_SOURCE, TAG_ROW, comm,
                    &status2);
          if (a != NULL)
            {
              MALLOC_INTERN(tmpa, size, PASTIX_FLOAT);
              MPI_Recv( tmpa, size, COMM_FLOAT, status.MPI_SOURCE, TAG_VAL,
                        comm, &status2);
            }
          cscd_addlocal_int(newloc2globssize[rank], colptr2send[rank],
                            rows2send[rank], values2send[rank],
                            newloc2globs[rank],
                            newloc2globssize[rank], tmpia            ,
                            tmpja          , tmpa             ,
                            newloc2globs[rank],
                            &newloc2globssize[rank],&tmpia2          ,
                            &tmpja2        , &tmpa2           ,
                            &add_two_floats, 1, API_YES);

          memFree_null(colptr2send[rank]);
          memFree_null(tmpia);
          if (rows2send[rank]   != NULL)
            memFree_null(rows2send[rank]);
          memFree_null(tmpja);
          if (values2send[rank] != NULL)
            memFree_null(values2send[rank]);
          if (NULL != tmpa)
            memFree_null(tmpa);

          colptr2send[rank] = tmpia2;
          rows2send[rank]   = tmpja2;
          values2send[rank] = tmpa2;
        }
      toreceive--;
    }
  for (i = 0; i < commSize; i++)
    {
      if (i != rank)
        {
          size = colptr2send[i][newloc2globssize[i]]-1;
          memFree_null(newloc2globs[i]);
          MPI_Wait(&requests_size[i], &status);
          if (size > 0)
            {
              MPI_Wait(&requests_col[i], &status);
              MPI_Wait(&requests_row[i], &status);
              memFree_null(rows2send[i]);
              if (a != NULL)
                {
                  MPI_Wait(&requests_val[i], &status);
                  memFree_null(values2send[i]);
                }
            }
          memFree_null(colptr2send[i]);
        }
    }
  *dn   = newloc2globssize[rank];
  *dl2g = newloc2globs[rank];
  *dia  = colptr2send[rank];
  *dja  = rows2send[rank];
  if ( a != NULL)
    *da   = values2send[rank];

  memFree_null(newloc2globssize);
  memFree_null(newloc2globs);
  memFree_null(colptr2send);
  memFree_null(rows2send);
  memFree_null(values2send);
  memFree_null(requests_size);
  memFree_null(requests_col);
  memFree_null(requests_row);
  memFree_null(requests_val);

  return 2;
#endif
}


/*
 * Function: cscd_redispatch
 *
 * Redistribute the first cscd into a new one using *dl2g*.
 *
 * - gather all new loc2globs on all processors.
 * - allocate *dia*, *dja* and *da*.
 * - Create new CSC for each processor and send it.
 * - Merge all new CSC to the new local CSC with <cscd_addlocal_int>.
 *
 * If communicator size is one, check that n = dn and
 * l2g = dl2g and simply create a copy of the first cscd.
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   rhs         - right-hand-side member corresponding to the first CSCD
 *                (can be NULL)
 *   nrhs        - number of right-hand-side.
 *   l2g         - local 2 global column numbers for first cscd
 *   dn          - Number of local columns
 *   dia         - New cscd starting index of each column in *ja* and *a*
 *   dja         - Row of each element in new CSCD
 *   da          - value of each cscd in new CSCD
 *   rhs         - right-hand-side member corresponding to the new CSCD
 *   dl2g        - local 2 global column numbers for new cscd
 *   comm        - MPI communicator
 *
 * Returns:
 *   NO_ERR           - If all goes well
 *   BADPARAMETER_ERR - If commsize = 1 and *n* != *dn* or *l2g* != *dl2g*.
 */
#ifndef cscd_redispatch
#  error "redefinition required"
#endif
int cscd_redispatch(PASTIX_INT   n, PASTIX_INT *   ia, PASTIX_INT *   ja, PASTIX_FLOAT *   a,
                    PASTIX_FLOAT *  rhs,  PASTIX_INT nrhs, PASTIX_INT *   l2g,
                    PASTIX_INT  dn, PASTIX_INT ** dia, PASTIX_INT ** dja, PASTIX_FLOAT ** da,
                    PASTIX_FLOAT ** drhs, PASTIX_INT *  dl2g,
                    MPI_Comm comm, PASTIX_INT dof)
{
  return cscd_redispatch_int(n,  ia,  ja,  a,  rhs,  nrhs, l2g,
                             dn, dia, dja, da, drhs, dl2g,
                             API_NO, comm, dof);
}
/*
 * Function: cscd_redispatch_int
 *
 * Redistribute the first cscd into a new one using *dl2g*.
 *
 * - gather all new loc2globs on all processors.
 * - allocate *dia*, *dja* and *da*.
 * - Create new CSC for each processor and send it.
 * - Merge all new CSC to the new local CSC with <cscd_addlocal_int>.
 *
 * If communicator size is one, check that n = dn and
 * l2g = dl2g and simply create a copy of the first cscd.
 *
 * Suppose that values and rhs are NOT NULL, or NULL  on all processors.
 * Thus, it doesn't works if all matrix is on one pocessor and we want to
 * distribute it with values...
 *
 * Parameters:
 *   n           - Number of local columns
 *   ia          - First cscd starting index of each column in *ja* and *a*
 *   ja          - Row of each element in first CSCD
 *   a           - value of each cscd in first CSCD (can be NULL)
 *   rhs         - right-hand-side member corresponding to the first CSCD
 *                 (can be NULL)
 *   nrhs        - number of right-hand-side.
 *   l2g         - local 2 global column numbers for first cscd
 *   dn          - Number of local columns
 *   dia         - New cscd starting index of each column in *ja* and *a*
 *   dja         - Row of each element in new CSCD
 *   da          - value of each cscd in new CSCD
 *   rhs         - right-hand-side member corresponding to the new CSCD
 *   dl2g        - local 2 global column numbers for new cscd
 *   malloc_flag - Internal (API_YES) or external (API_NO) malloc use.
 *   comm        - MPI communicator
 *
 * Returns:
 *   EXIT_SUCCESS - If all goes well
 *   EXIT_FAILURE - If commsize = 1 and *n* != *dn* or *l2g* != *dl2g*.
 */
#ifndef cscd_redispatch_int
#  error "redefinition required"
#endif
int cscd_redispatch_int(PASTIX_INT   n, PASTIX_INT *   ia, PASTIX_INT *   ja, PASTIX_FLOAT *   a,
                        PASTIX_FLOAT *  rhs,  PASTIX_INT nrhs, PASTIX_INT *   l2g,
                        PASTIX_INT  dn, PASTIX_INT ** dia, PASTIX_INT ** dja, PASTIX_FLOAT ** da,
                        PASTIX_FLOAT ** drhs, PASTIX_INT *  dl2g,
                        int  malloc_flag, MPI_Comm comm, PASTIX_INT dof)
{
#ifdef FORCE_NOMPI
  (void)n; (void)ia; (void)ja; (void)a; (void)rhs; (void)nrhs; (void)l2g;
  (void)dn; (void)dia; (void)dja; (void)da; (void)drhs; (void)dl2g;
  (void)malloc_flag; (void)comm; (void)dof;
  return NO_ERR;
#else /* FORCE_NOMPI */
  PASTIX_INT i,j;
  PASTIX_INT           globidx;
  PASTIX_INT           locidx;
  PASTIX_INT           tmpn;
  PASTIX_INT         * tmpia  = NULL;
  PASTIX_INT         * tmpja  = NULL;
  PASTIX_FLOAT       * tmpa   = NULL;
  PASTIX_FLOAT       * tmprhs = NULL;
  PASTIX_INT         * tmpia2 = NULL;
  PASTIX_INT         * tmpja2   = NULL;
  PASTIX_FLOAT       * tmpa2  = NULL;
  PASTIX_INT           size;
  PASTIX_INT         * sizep                 = NULL;
  PASTIX_INT         * newloc2globssize      = NULL;
  PASTIX_INT         * newloc2globssize_recv = NULL;
  PASTIX_INT        ** newloc2globs          = NULL;
  PASTIX_INT        ** colptr2send           = NULL;
  PASTIX_INT        ** rows2send             = NULL;
  PASTIX_FLOAT      ** values2send           = NULL;
  PASTIX_FLOAT      ** rhs2send              = NULL;
  MPI_Request * requests_size         = NULL;
  MPI_Request * requests_col          = NULL;
  MPI_Request * requests_row          = NULL;
  MPI_Request * requests_val          = NULL;
  MPI_Request * requests_rhs          = NULL;
  MPI_Status    status;
  MPI_Status    status2;
  PASTIX_INT           toreceive;
  int           commSize;
  int           rank;
  int           flags[2];
  int           flags_recv[2];

  if (NULL != a)
    flags[0] = 1;
  else
    flags[0] = 0;

  if (NULL != rhs)
    flags[1] = 1;
  else
    flags[1] = 0;

  MPI_Allreduce(flags, flags_recv, 2, MPI_INT, MPI_LOR, comm);
  MPI_Comm_size(comm, &commSize);
  MPI_Comm_rank(comm, &rank);
  if (commSize == 1)
    {
      if (n != dn)
        {
          return BADPARAMETER_ERR;
        }
      for (i = 0; i < n; i++)
        if (dl2g[i] != l2g[i])
          {
            return BADPARAMETER_ERR;
          }

      MALLOC_INTOREXTERN(*dia, n+1, PASTIX_INT, malloc_flag);
      memcpy(*dia,ia,(n+1)*sizeof(PASTIX_INT));
      MALLOC_INTOREXTERN(*dja, ia[n]-1, PASTIX_INT, malloc_flag);
      memcpy(*dja,ja,(ia[n]-1)*sizeof(PASTIX_INT));
      if (flags_recv[0] == 1)
        {
          MALLOC_INTOREXTERN(*da, (ia[n]-1)*dof*dof, PASTIX_FLOAT, malloc_flag);
          memcpy(*da,a,(ia[n]-1)*sizeof(PASTIX_FLOAT)*dof*dof);
        }
      if (flags_recv[1] == 1)
        {
          MALLOC_INTOREXTERN(*drhs, n*nrhs*dof, PASTIX_FLOAT, malloc_flag);
          memcpy(*drhs,rhs,n*nrhs*sizeof(PASTIX_FLOAT)*dof);
        }
      return NO_ERR;
    }

  /* Build all newloc2globs on all processors */
  MALLOC_INTERN(newloc2globssize, commSize + 1, PASTIX_INT );

  for (i = 0; i < rank+1; i++)
    newloc2globssize[i] = 0;
  for (i = rank+1; i < commSize+1; i++)
    newloc2globssize[i] = dn;

  /* Construct newloc2globssize
   * It will contain start of each dl2g in newloc2globs
   */
  MALLOC_INTERN(newloc2globssize_recv, commSize + 1, PASTIX_INT );
  MPI_Allreduce(newloc2globssize, newloc2globssize_recv, commSize+1,
                COMM_INT, MPI_SUM, comm);
  memFree_null(newloc2globssize);
  newloc2globssize = newloc2globssize_recv;

  /* Construct newloc2globs
   * It will contains all dl2g in one tabular
   */
  MALLOC_INTERN(newloc2globs, commSize, PASTIX_INT*);
  for (i = 0; i < commSize; i++)
    {
      if (rank != i)
        {
          newloc2globs[i] = NULL;
          MALLOC_INTERN(newloc2globs[i],
                        newloc2globssize[i+1] - newloc2globssize[i], PASTIX_INT);
          for (j = 0; j < newloc2globssize[i+1] - newloc2globssize[i]; j++)
            newloc2globs[i][j] = 0;
        }
      else
        {
          newloc2globs[i] = dl2g;
        }
      MPI_Bcast(newloc2globs[i], newloc2globssize[i+1] - newloc2globssize[i],
                COMM_INT, i, comm);
    }
  /* Create one CSCD for each proc */
  MALLOC_INTERN(colptr2send, commSize, PASTIX_INT*);
  MALLOC_INTERN(rows2send,   commSize, PASTIX_INT*);
  if (flags_recv[0] == 1)
    MALLOC_INTERN(values2send, commSize, PASTIX_FLOAT*);
  if (flags_recv[1] == 1)
    MALLOC_INTERN(rhs2send,    commSize, PASTIX_FLOAT*);

  for (i = 0; i < commSize; i++)
    {
      colptr2send[i] = NULL;
      rows2send[i]   = NULL;
      if (flags_recv[0] == 1)
        values2send[i] = NULL;
      if (flags_recv[1] == 1)
        rhs2send[i]    = NULL;
    }

  MALLOC_INTERN(sizep,         commSize, PASTIX_INT);
  MALLOC_INTERN(requests_size, commSize, MPI_Request);
  MALLOC_INTERN(requests_col,  commSize, MPI_Request);
  MALLOC_INTERN(requests_row,  commSize, MPI_Request);
  if (flags_recv[0] == 1)
    MALLOC_INTERN(requests_val,  commSize, MPI_Request);
  if (flags_recv[1] == 1)
    MALLOC_INTERN(requests_rhs,  commSize, MPI_Request);

  /* Sending informations to all proc */
  for (i = 0; i < commSize; i++)
    {
      MALLOC_INTERN(colptr2send[i],
                    newloc2globssize[i+1] - newloc2globssize[i] + 1, PASTIX_INT);
      for (j = 0; j < newloc2globssize[i+1] - newloc2globssize[i] + 1; j++)
        colptr2send[i][j] = 1;

      /* Adding old CSCD to the CSCD to send */
      if (i == rank)
        {
          cscd_addlocal_int(newloc2globssize[i+1] - newloc2globssize[i],
                            colptr2send[i], rows2send[i],
                            ((flags_recv[0]==0)?(NULL):(values2send[i])),
                            newloc2globs[i],
                            n, ia, ja, ((flags_recv[0]==0)?(NULL):(a))   , l2g,
                            &tmpn, &tmpia, &tmpja,
                            ((flags_recv[0]==0)?(NULL):(&tmpa)),
                            &add_two_floats,
                            dof, malloc_flag);
        }
      else
        {
          cscd_addlocal_int(newloc2globssize[i+1] - newloc2globssize[i],
                            colptr2send[i], rows2send[i],
                            ((flags_recv[0]==0)?(NULL):(values2send[i])),
                            newloc2globs[i],
                            n, ia, ja, ((flags_recv[0]==0)?(NULL):(a)), l2g,
                            &tmpn, &tmpia, &tmpja,
                            ((flags_recv[0]==0)?(NULL):(&tmpa)),
                            &add_two_floats,
                            dof, API_YES);
        }
      if (newloc2globssize[i+1] - newloc2globssize[i] != tmpn)
        errorPrint("newloc2globssize[i] != tmpn\n");
      memFree_null(colptr2send[i]);

      if (flags_recv[1] == 1 && i != rank)
        {
          PASTIX_INT newrhs_size = (newloc2globssize[i+1] - newloc2globssize[i]);
          MALLOC_INTERN(rhs2send[i],
                        dof*nrhs*newrhs_size,
                        PASTIX_FLOAT);
          memset(rhs2send[i], 0,
                 dof*nrhs*newrhs_size*sizeof(PASTIX_FLOAT));
          for (j = 0; j < newloc2globssize[i+1] - newloc2globssize[i]; j++)
            {
              int iter_rhs;
              for (iter_rhs = 0; iter_rhs < nrhs; iter_rhs++)
                {
                  globidx = newloc2globs[i][j];
                  searchInList(globidx, l2g, n, locidx);
                  if (locidx >= 0)
                    {
                      PASTIX_INT k;
                      for (k = 0; k< dof; k++)
                        rhs2send[i][dof*(j+iter_rhs*newrhs_size)+k] =
                          rhs[dof*(locidx+iter_rhs*n)+k];
                    }
                }
            }
        }

      colptr2send[i] = tmpia;
      rows2send[i]   = tmpja;
      if (flags_recv[0] == 1)
        values2send[i] = tmpa;

      tmpia = NULL;
      tmpja = NULL;
      if (flags_recv[0] == 1)
        tmpa  = NULL;

      /* Sending the CSCD built */
      if (i != rank)
        {
          size     = newloc2globssize[i+1] - newloc2globssize[i];
          sizep[i] = colptr2send[i][size]-1;
          MPI_Isend(&sizep[i], 1, COMM_INT, i, TAG_SIZE, comm,
                    &requests_size[i]);
          if (sizep[i] > 0)
            {
              MPI_Isend(colptr2send[i], size+1,
                        COMM_INT, i, TAG_COL, comm, &requests_col[i]);
              MPI_Isend(rows2send[i],   sizep[i],
                        COMM_INT, i, TAG_ROW, comm, &requests_row[i]);
              if (flags_recv[0] == 1)
                MPI_Isend(values2send[i], sizep[i]*dof*dof,
                          COMM_FLOAT, i, TAG_VAL, comm, &requests_val[i]);
            }
          if (flags_recv[1] == 1)
            {
              PASTIX_INT newrhs_size = (newloc2globssize[i+1] - newloc2globssize[i]);
              MPI_Isend(rhs2send[i], nrhs*newrhs_size*dof,
                        COMM_FLOAT, i, TAG_RHS, comm, &requests_rhs[i]);
            }

        }
    }

  toreceive = commSize -1;
  /* Receive and add to local CSCD */
  if (flags_recv[1] == 1)
    {
      PASTIX_INT iter_rhs;
      MALLOC_INTOREXTERN(*drhs, dn*nrhs*dof, PASTIX_FLOAT, malloc_flag);
      MALLOC_INTERN(tmprhs, dn*nrhs*dof, PASTIX_FLOAT);
      for (i = 0; i < dn*nrhs*dof; i++)
        (*drhs)[i] = 0;

      /* Initialize local RHS */
      for (iter_rhs = 0; iter_rhs < nrhs; iter_rhs++)
        {
          for (j = 0; j < dn; j++)
            {
              globidx = newloc2globs[rank][j];
              searchInList(globidx, l2g, n, locidx);
              if (locidx >= 0)
                {
                  PASTIX_INT k;
                  for (k = 0; k< dof; k++)
                    (*drhs)[dof*(j + iter_rhs*dn)+k] = rhs[dof*(locidx + iter_rhs*n)+k];
                }
            }
        }
    }
  while(toreceive > 0)
    {
      size = 0;
      MPI_Recv( &size, 1, COMM_INT, MPI_ANY_SOURCE, TAG_SIZE, comm, &status);
      if (size > 0)
        {
          /* Adding contribution CSCD to local CSCD */
          MALLOC_INTERN(tmpia,
                        newloc2globssize[rank+1] - newloc2globssize[rank]+1,
                        PASTIX_INT);
          MPI_Recv( tmpia, (newloc2globssize[rank+1] -
                            newloc2globssize[rank] + 1),
                    COMM_INT,
                    status.MPI_SOURCE, TAG_COL, comm, &status2);
          MALLOC_INTERN(tmpja, size, PASTIX_INT);
          MPI_Recv( tmpja, size, COMM_INT,  status.MPI_SOURCE,
                    TAG_ROW, comm, &status2);
          if (flags_recv[0] == 1)
            {
              MALLOC_INTERN(tmpa, dof*dof*size, PASTIX_FLOAT);
              MPI_Recv( tmpa, dof*dof*size, COMM_FLOAT, status.MPI_SOURCE,
                        TAG_VAL, comm, &status2);
            }
          else
            {
              tmpa = NULL;
            }

          cscd_addlocal_int(newloc2globssize[rank+1] - newloc2globssize[rank],
                            colptr2send[rank], rows2send[rank],
                            ((flags_recv[0]==0)?(NULL):(values2send[rank])),
                            newloc2globs[rank],
                            newloc2globssize[rank+1] - newloc2globssize[rank],
                            tmpia            , tmpja          ,
                            tmpa             , newloc2globs[rank],
                            &tmpn            , &tmpia2          ,
                            &tmpja2        ,
                            (flags_recv[0]==1)?(&tmpa2):NULL, &add_two_floats , dof,
                            malloc_flag);


          memFree_null(tmpia);
          memFree_null(tmpja);
          if (NULL != tmpa)
            memFree_null(tmpa);
          FREE_NULL_INTOREXT(colptr2send[rank], malloc_flag);
          if (rows2send[rank]   != NULL)
            FREE_NULL_INTOREXT(rows2send[rank], malloc_flag);
          if (flags_recv[0]==1)
            if (values2send[rank] != NULL)
              FREE_NULL_INTOREXT(values2send[rank], malloc_flag);
          colptr2send[rank] = tmpia2;
          rows2send[rank]   = tmpja2;
          if (flags_recv[0]==1)
            values2send[rank] = tmpa2;
        }

      if (flags_recv[1] == 1)
        {
          PASTIX_INT iter_rhs;
          MPI_Recv(tmprhs, dof*dn*nrhs, COMM_FLOAT, MPI_ANY_SOURCE, TAG_RHS,
                   comm, &status);
          for (iter_rhs = 0; iter_rhs < nrhs; iter_rhs++)
            {
              for (j = 0; j < dn; j++)
                {
                  PASTIX_INT k;
                  for (k = 0; k< dof; k++)
                    (*drhs)[dof*(j + iter_rhs*dn)+k] += tmprhs[dof*(j + iter_rhs*dn)+k];
                }
            }
        }
      toreceive--;
    }

  if (flags_recv[1] == 1)
    memFree_null(tmprhs);

  for (i = 0; i < commSize; i++)
    {
      if (i != rank)
        {
          size = colptr2send[i][newloc2globssize[i+1] - newloc2globssize[i]]-1;
          if (newloc2globs[i])
            memFree_null(newloc2globs[i]);
          MPI_Wait(&requests_size[i], &status);
          if (size > 0)
            {
              MPI_Wait(&requests_col[i], &status);
              MPI_Wait(&requests_row[i], &status);
              memFree_null(rows2send[i]);
              if (flags_recv[0] == 1)
                {
                  MPI_Wait(&requests_val[i], &status);
                  memFree_null(values2send[i]);
                }
            }
          memFree_null(colptr2send[i]);
          if (rows2send[i]   != NULL)
            memFree_null(rows2send[i]);
          if (flags_recv[0] == 1)
            if (values2send[i] != NULL)
              memFree_null(values2send[i]);

          if (flags_recv[1] == 1)
            {
              MPI_Wait(&requests_rhs[i], &status);
              memFree_null(rhs2send[i]);
            }

        }
    }

  MPI_Barrier(comm);

  *dia  = colptr2send[rank];
  *dja  = rows2send[rank];
  if ( flags_recv[0] == 1)
    *da   = values2send[rank];

  memFree_null(newloc2globssize);
  memFree_null(newloc2globs);
  memFree_null(colptr2send);
  memFree_null(rows2send);
  if (flags_recv[0] == 1)
    memFree_null(values2send);
  if (flags_recv[1] == 1)
    memFree_null(rhs2send);
  memFree_null(sizep);
  memFree_null(requests_size);
  memFree_null(requests_col);
  memFree_null(requests_row);
  if (flags_recv[0] == 1)
    memFree_null(requests_val);
  if (flags_recv[1] == 1)
    memFree_null(requests_rhs);

  return NO_ERR;
#endif
}


/*
   Function: cscd2csc

   Transform a cscd to a csc.
   colptr2, row2, avals2, rhs2, perm2, invp2 are allocated here.

   External function, allocation are not of the internal type.

   Parameters:
   lN          - number of local column.
   lcolptr     - starting index of each local column in row and avals.
   lrow        _ row number of each local element.
   lavals      - values of each local element.
   lrhs        - local part of the right hand side.
   lperm       - local part of the permutation tabular.
   linvp       - Means nothing, to suppress.
   gN          - global number of columns (output).
   gcolptr     - starting index of each column in row2 and avals2 (output).
   grow        - row number of each element (output).
   gavals      - values of each element (output).
   grhs        - global right hand side (output).
   gperm       - global permutation tabular (output).
   ginvp       - global reverse permutation tabular (output).
   loc2glob    - global number of each local column.
   pastix_comm - PaStiX MPI communicator.

*/
#ifndef cscd2csc
#  error "redefinition required"
#endif
void  cscd2csc(PASTIX_INT  lN, PASTIX_INT *  lcolptr, PASTIX_INT * lrow, PASTIX_FLOAT * lavals,
               PASTIX_FLOAT * lrhs, PASTIX_INT * lperm, PASTIX_INT * linvp,
               PASTIX_INT *gN, PASTIX_INT ** gcolptr, PASTIX_INT **grow, PASTIX_FLOAT **gavals,
               PASTIX_FLOAT **grhs, PASTIX_INT **gperm, PASTIX_INT **ginvp,
               PASTIX_INT *loc2glob, MPI_Comm pastix_comm, PASTIX_INT ndof)
{
  cscd2csc_int(lN, lcolptr, lrow, lavals,
               lrhs, lperm, linvp,
               gN, gcolptr, grow, gavals,
               grhs, gperm, ginvp,
               loc2glob, pastix_comm, ndof, API_NO);
}
/*
 *   Function: cscd2csc_int
 *
 *   Transform a cscd to a csc.
 *   colptr2, row2, avals2, rhs2, perm2, invp2 are allocated here.
 *
 *   Parameters:
 *     lN          - number of local column.
 *     lcolptr     - starting index of each local column in row and avals.
 *     lrow        _ row number of each local element.
 *     lavals      - values of each local element.
 *     lrhs        - local part of the right hand side.
 *     lperm       - local part of the permutation tabular.
 *     linvp       - Means nothing, to suppress.
 *     gN          - global number of columns (output).
 *     gcolptr     - starting index of each column in row2 and avals2 (output).
 *     grow        - row number of each element (output).
 *     gavals      - values of each element (output).
 *     grhs        - global right hand side (output).
 *     gperm       - global permutation tabular (output).
 *     ginvp       - global reverse permutation tabular (output).
 *     loc2glob    - global number of each local column.
 *     pastix_comm - PaStiX MPI communicator.
 *     intern_flag - Decide if malloc will use internal or external macros.
 *
 */
#ifndef cscd2csc_int
#  error "redefinition required"
#endif
void  cscd2csc_int(PASTIX_INT  lN, PASTIX_INT *  lcolptr, PASTIX_INT * lrow, PASTIX_FLOAT * lavals,
                   PASTIX_FLOAT * lrhs, PASTIX_INT * lperm, PASTIX_INT * linvp,
                   PASTIX_INT *gN, PASTIX_INT ** gcolptr, PASTIX_INT **grow, PASTIX_FLOAT **gavals,
                   PASTIX_FLOAT **grhs, PASTIX_INT **gperm, PASTIX_INT **ginvp,
                   PASTIX_INT *loc2glob, MPI_Comm pastix_comm, PASTIX_INT ndof, int intern_flag)
{
  PASTIX_INT      i,j;
  int      rank;
  int      commSize;
  PASTIX_INT      index;
  PASTIX_INT    * AllLocalN   = NULL;
  PASTIX_INT   ** AllColptr   = NULL;
  PASTIX_INT   ** AllRow      = NULL;
  PASTIX_INT   ** AllLoc2glob = NULL;
  PASTIX_FLOAT ** AllAvals    = NULL;
  PASTIX_FLOAT ** AllRhs      = NULL;
  PASTIX_INT   ** AllPerm     = NULL;
  /*   PASTIX_INT   ** AllInvp; */
  int      proc;
  (void)pastix_comm; (void)linvp;

  MPI_Comm_rank(pastix_comm,&rank);
  MPI_Comm_size(pastix_comm,&commSize);

  MALLOC_INTERN(AllLocalN,   commSize, PASTIX_INT );
  MALLOC_INTERN(AllColptr,   commSize, PASTIX_INT*);
  MALLOC_INTERN(AllRow,      commSize, PASTIX_INT*);
  MALLOC_INTERN(AllLoc2glob, commSize, PASTIX_INT*);

  if (gavals != NULL)
    MALLOC_INTERN(AllAvals, commSize, PASTIX_FLOAT*);

  if (grhs != NULL)
    MALLOC_INTERN(AllRhs, commSize, PASTIX_FLOAT*);

  if (gperm != NULL)
    MALLOC_INTERN(AllPerm, commSize, PASTIX_INT*);

  for (i = 0; i < commSize; i++)
    {
      if (i == rank)
        {
          AllLocalN[i]   = lN;
          AllColptr[i]   = lcolptr;
          AllRow[i]      = lrow;
          if (gavals != NULL)
            AllAvals[i]    = lavals;
          if (grhs != NULL)
            AllRhs[i]      = lrhs;
          if (gperm != NULL)
            AllPerm[i]     = lperm;
          AllLoc2glob[i] = loc2glob;
        }

      MPI_Bcast(&AllLocalN[i] , 1, COMM_INT, i, pastix_comm);
      if (rank != i)
        {
          MALLOC_INTERN(AllColptr[i],   AllLocalN[i]+1, PASTIX_INT);
          MALLOC_INTERN(AllLoc2glob[i], AllLocalN[i],   PASTIX_INT);
          if (grhs != NULL)
            MALLOC_INTERN(AllRhs[i], ndof*AllLocalN[i], PASTIX_FLOAT);

          if (gperm != NULL)
            MALLOC_INTERN(AllPerm[i], AllLocalN[i], PASTIX_INT);
        }

      MPI_Bcast(AllColptr[i], AllLocalN[i]+1, COMM_INT  , i, pastix_comm);
      if (rank != i)
        {
          MALLOC_INTERN(AllRow[i], AllColptr[i][AllLocalN[i]]-1, PASTIX_INT);
          if (gavals != NULL)
            MALLOC_INTERN(AllAvals[i], ndof*ndof*(AllColptr[i][AllLocalN[i]]-1), PASTIX_FLOAT);
        }

      MPI_Bcast(AllRow[i], AllColptr[i][AllLocalN[i]]-1,
                COMM_INT, i, pastix_comm);
      MPI_Bcast(AllLoc2glob[i], AllLocalN[i], COMM_INT, i, pastix_comm);
      if (gperm != NULL) {
        MPI_Bcast(AllPerm[i], AllLocalN[i], COMM_INT, i, pastix_comm);
      }

      if (grhs != NULL) {
        MPI_Bcast(AllRhs[i], ndof*AllLocalN[i], COMM_FLOAT, i, pastix_comm);
      }
      if (gavals != NULL) {
        MPI_Bcast(AllAvals[i], ndof*ndof*(AllColptr[i][AllLocalN[i]]-1),
                  COMM_FLOAT, i, pastix_comm);
      }
    }


  *gN = 0;
  for (i = 0; i < commSize; i++)
    {
      *gN += AllLocalN[i];
    }

  MALLOC_INTOREXTERN(*gcolptr, ((*gN)+1), PASTIX_INT, intern_flag);
  /* Fill in gcolptr */
  (*gcolptr)[0] = 1;
  for (i = 0; i < (*gN); i++)
    {
      for (proc = 0; proc < commSize; proc++)
        {
          searchInList(i+1, AllLoc2glob[proc], AllLocalN[proc], index);
          if ( index >= 0 )
            {
              (*gcolptr)[i+1] = (*gcolptr)[i] +
                AllColptr[proc][index+1] -
                AllColptr[proc][index];
              break;
            }

        }
    }
  MALLOC_INTOREXTERN(*grow, (*gcolptr)[(*gN)]-1, PASTIX_INT, intern_flag);
  if (gavals != NULL)
    MALLOC_INTOREXTERN(*gavals, ndof*ndof*((*gcolptr)[(*gN)]-1), PASTIX_FLOAT, intern_flag);
  if (grhs != NULL)
    MALLOC_INTOREXTERN(*grhs, *gN*ndof, PASTIX_FLOAT, intern_flag);
  if (gperm != NULL)
    {
      MALLOC_INTOREXTERN(*gperm, *gN, PASTIX_INT, intern_flag);
      MALLOC_INTOREXTERN(*ginvp, *gN, PASTIX_INT, intern_flag);
    }

  /* Fill-in grow, gavals, grhs, gperm and ginvp*/
  for (proc = 0; proc < commSize; proc++)
    {
      for (i = 0; i < AllLocalN[proc]; i++)
        {
          memcpy(&(*grow)[(*gcolptr)[AllLoc2glob[proc][i]-1]-1],
                 &AllRow[proc][AllColptr[proc][i]-1],
                 (AllColptr[proc][i+1] - AllColptr[proc][i])*sizeof(PASTIX_INT));
          if (gavals != NULL)
            memcpy(&(*gavals)[(*gcolptr)[AllLoc2glob[proc][i]-1]-1],
                   &AllAvals[proc][AllColptr[proc][i]-1],
                   ndof*ndof*(AllColptr[proc][i+1] - AllColptr[proc][i])*sizeof(PASTIX_FLOAT));
          if (grhs != NULL)
            for (j = 0; j < ndof; j++)
              (*grhs)[ndof*(AllLoc2glob[proc][i]-1)+j]  = AllRhs[proc][ndof*i+j];
          if (gperm != NULL)
            {
              (*gperm)[AllLoc2glob[proc][i]-1] = AllPerm[proc][i];
              /* (*ginvp)[AllLoc2glob[proc][i]-1] = AllInvp[proc][i]; */

            }
        }
    }
  if (gperm != NULL)
    {
      for (i = 0; i < *gN; i++)
        (*ginvp)[(*gperm)[i]] = i;
    }
  for (i = 0; i < commSize; i++)
    {
      if (rank != i)
        {
          memFree_null(AllColptr[i]);
          memFree_null(AllRow[i]);
          if (gavals != NULL)
            memFree_null(AllAvals[i]);
          if (grhs != NULL)
            memFree_null(AllRhs[i]);
          memFree_null(AllLoc2glob[i]);
          if (gperm != NULL)
            {
              memFree_null(AllPerm[i]);
              /*            memFree_null(AllInvp[i]); */
            }
        }
    }
  memFree_null(AllLocalN);
  memFree_null(AllColptr);
  memFree_null(AllRow);
  if (gavals != NULL)
    memFree_null(AllAvals);
  if (grhs != NULL)
    memFree_null(AllRhs);
  memFree_null(AllLoc2glob);
  if (gperm != NULL)
    {
      memFree_null(AllPerm);
      /*       memFree_null(AllInvp); */
    }

}

/*
 *   Function: csc2cscd
 *
 *   Transform a csc to a cscd.
 *   Allocate the CSCD.
 *   If grhs == NULL forget right hand side part.
 *   If gperm == NULL forget permutation and reverse permutation part.
 *
 *
 *   TODO: USE_EXTERNAL MALLOCS
 *
 *   Parameters:
 *     gN       - global number of columns
 *     gcolptr  - global starting index of each column in grows ans gavals.
 *     grows    - global rows of each element.
 *     gavals   - global values of each element.
 *     gperm    - global permutation tabular.
 *     ginvp    - global reverse permutation tabular.
 *     lN       - local number of columns.
 *     lcolptr  - starting index of each local column.
 *     lrowptr  - row number of each local element.
 *     lavals   - values of each local element.
 *     lrhs     - local part of the right hand side (output).
 *     lperm    - local part of the permutation tabular (output).
 *     linvp    - local part of the reverse permutation tabular (output).
 *     loc2glob - global numbers of local columns (before permutation).
 */
#ifndef csc2cscd
#  error "redefinition required"
#endif
void  csc2cscd(PASTIX_INT gN, PASTIX_INT *  gcolptr, PASTIX_INT *  grow,
               PASTIX_FLOAT *  gavals, PASTIX_FLOAT *  grhs, PASTIX_INT *  gperm, PASTIX_INT *  ginvp,
               PASTIX_INT lN, PASTIX_INT ** lcolptr, PASTIX_INT ** lrow,
               PASTIX_FLOAT ** lavals, PASTIX_FLOAT ** lrhs, PASTIX_INT ** lperm, PASTIX_INT ** linvp,
               PASTIX_INT *loc2glob)
{
  PASTIX_INT   i;
  (void)gN; (void)ginvp; (void)linvp;

  if (grhs != NULL)
    {
      MALLOC_INTERN(*lrhs, lN, PASTIX_FLOAT);
      for (i = 0; i < lN; i++)
        (*lrhs )[i] = grhs [loc2glob[i]];
    }
  if (gperm != NULL)
    {
      MALLOC_INTERN(*lperm, lN, PASTIX_INT);
      for (i = 0; i < lN; i++)
        {
          (*lperm)[i] = gperm[loc2glob[i]];
        }
    }

  MALLOC_INTERN(*lcolptr, lN+1, PASTIX_INT);
  (*lcolptr)[0] = 1;
  for (i = 0; i < lN; i++)
    {
      (*lcolptr)[i+1] = (*lcolptr)[i] +
        gcolptr[loc2glob[i]] - gcolptr[loc2glob[i]-1];
    }

  MALLOC_INTERN(*lrow, (*lcolptr)[lN]-1, PASTIX_INT);
  if (gavals != NULL)
    MALLOC_INTERN(*lavals, (*lcolptr)[lN]-1, PASTIX_FLOAT);

  for (i = 0; i < lN; i++)
    {
      memcpy(&((*lrow)  [(*lcolptr)[i]-1]),
             &(  grow  [gcolptr[loc2glob[i]-1]-1]),
             (gcolptr[loc2glob[i]] - gcolptr[loc2glob[i]-1])*sizeof (PASTIX_INT));
      if (gavals != NULL)
        memcpy(&((*lavals)[(*lcolptr)[i]-1]),
               &(  gavals[gcolptr[loc2glob[i]-1]]),
               (gcolptr[loc2glob[i]] - gcolptr[loc2glob[i]-1])*sizeof(PASTIX_FLOAT));
    }


}

/*
 *  Function: cscd_noDiag
 *
 *  Removes diagonal elements from a CSCD.
 *  *ja* and *a* can be reallocated to
 *  ia[n]-1 elements after this call.
 *
 *
 *  Parameters:
 *    n           - Number of local columns
 *    ia          - First cscd starting index of each column in *ja* and *a*
 *    ja          - Row of each element in first CSCD
 *    a           - value of each cscd in first CSCD (can be NULL)
 *    l2g         - local 2 global column numbers for first cscd
 *
 */
#ifndef cscd_noDiag
#  error "redefinition required"
#endif
int cscd_noDiag(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT * a, PASTIX_INT * l2g)
{
  PASTIX_INT     i;
  PASTIX_INT     j;
  PASTIX_INT     previous_index = 0;
  PASTIX_INT     index          = 0;

  print_debug(DBG_CSCD, "->cscd_noDiag\n");
  /* for all local column i */
  for(i = 0; i < n; i++)
    {
      /* for all line in the colupmn */
      for (j = ia[i]; j < ia[i+1]; j++)
        {
          /* if the element is not on diag */
          if (ja[j-1] != l2g[i] )
            {
              /* we add it to ja */
              ja[index] = ja[j-1];
              if (NULL != a)
                a[index] = a[j-1];
              index++;
            }
        }
      /* we change ia[i] to it's new value */
      ia[i] = previous_index+1;
      previous_index = index;
    }
  ia[n] = index+1;
  print_debug(DBG_CSCD, "<-cscd_noDiag\n");
  return EXIT_SUCCESS;
}


/*
 *  Function: cscd_save
 *
 *  save a distributed csc to disk.
 *  files are called $(filename) and $(filename)$(RANK)
 *  if filename is NULL then filename = cscd_matrix.
 *
 *  file filename contains the number of processors/files
 *  on first line. Then each line contain the name of each file
 *  (here $(filename)$(RANK)).
 *
 *
 *
 *  Parameters:
 *    n           - Number of local columns
 *    ia          - First cscd starting index of each column in *ja* and *a*
 *    ja          - Row of each element in first CSCD
 *    a           - value of each cscd in first CSCD (can be NULL)
 *    rhs         - Right hand side.
 *    l2g         - local 2 global column numbers for first cscd
 *    dof         - Number of degrees of freedom
 *    filename    - name of the files.
 *    comm        - MPI communicator
 */
#ifndef cscd_save
#  error "redefinition required"
#endif
int cscd_save(PASTIX_INT          n,
              PASTIX_INT        * ia,
              PASTIX_INT        * ja,
              PASTIX_FLOAT      * a,
              PASTIX_FLOAT      * rhs,
              PASTIX_INT        * l2g,
              int          dof,
              const char * filename,
              MPI_Comm     comm)
{
  FILE * outfile;
  char   file[64];
  char   file2[64];
  int    myrank, nbproc;
  PASTIX_INT    i, j, k, l, vertglb = 0;
  (void)comm;

  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &nbproc);

  MPI_Allreduce(&n, &vertglb, 1, COMM_INT, MPI_SUM, comm);

  if (filename == NULL)
    sprintf(file, "cscd_matrix");
  else
    sprintf(file, "%s", filename);

  if (myrank == 0)
    {
      outfile = fopen(file, "w+");
      fprintf(outfile, "%d %d %d\n", (int)nbproc, (a != NULL)?1:0,
              (rhs != NULL)?1:0);
      for(i=0; i<nbproc; i++)
        fprintf(outfile, "%s%d\n", file, (int)i);
      fclose(outfile);
    }

  sprintf(file2, "%s%d", file, myrank);

  outfile = fopen(file2, "w+");
  fprintf(outfile, "%ld %ld\n", (long)vertglb, (long)n);
  fprintf(outfile, "%ld\n", (long)ia[n]-1);

  /* Copie de Loc2Glob */
  if (l2g != NULL)
    {
      for (i=0; i<n; i++)
        {
          fprintf(outfile, "%ld ", (long)l2g[i]);
          if (i%4 == 3) fprintf(outfile, "\n");
        }
      if ((i-1)%4 !=3) fprintf(outfile, "\n");
    }

  /* Copie de IA */
  for (i=0; i<n+1; i++)
    {
      fprintf(outfile, "%ld ", (long)ia[i]);
      if (i%4 == 3) fprintf(outfile, "\n");
    }
  if ((i-1)%4 !=3) fprintf(outfile, "\n");

  /* Copie de JA */
  for (i=0; i<ia[n]-1; i++)
    {
      fprintf(outfile, "%ld ", (long)ja[i]);
      if (i%4 == 3) fprintf(outfile, "\n");
    }
  if ((i-1)%4 !=3) fprintf(outfile, "\n");

  /* Copie de Avals */

  if (a != NULL)
    {
      for (i=0; i<(ia[n]-1)*dof*dof; i++)
        {
#ifdef CPLX
          fprintf(outfile, "%lg %lg ", creal(a[i]), cimag(a[i]));
          if (i%2 == 1) fprintf(outfile, "\n");
#else
          fprintf(outfile, "%lg ", a[i]);
          if (i%4 == 3) fprintf(outfile, "\n");
#endif
        }
#ifdef CPLX
      if ((i-1)%2 !=3) fprintf(outfile, "\n");
#else
      if ((i-1)%4 !=3) fprintf(outfile, "\n");
#endif
    }
  /* Copie de RHS */
  if (rhs != NULL)
    {
      for (i=0; i<n*dof; i++)
        {
#ifdef CPLX
          fprintf(outfile, "%lg %lg ", creal(rhs[i]), cimag(rhs[i]));
          if (i%2 == 1) fprintf(outfile, "\n");
#else
          fprintf(outfile, "%lg ", rhs[i]);
          if (i%4 == 3) fprintf(outfile, "\n");
#endif
        }
#ifdef CPLX
      if ((i-1)%2 !=3) fprintf(outfile, "\n");
#else
      if ((i-1)%4 !=3) fprintf(outfile, "\n");
#endif
    }
  else
    {
      fprintf(outfile, "rhs vaut NULL\n");
    }

  fclose(outfile);

  sprintf(file2, "%s_ijv%d", file, myrank);
  outfile = fopen(file2, "w+");

  for (i = 0; i < n; i++)
    for (j = ia[i]-1; j < ia[i+1]-1; j++)
      {
        for (k = 0; k < dof; k++)/* column dof */
          {
            for(l = 0; l < dof; l++) /* line dof */
              {
                if (l2g != NULL)
                  {
                    if (a != NULL)
                      {
#ifdef CPLX
#ifdef DUMP_BY_NODE
                        fprintf(outfile, "%ld (%ld)\t%ld (%ld)\t%.12lg\t%.12lg\n",
                                (long)l2g[i], (long)(k+1),
                                (long)ja[j],  (long)(l+1),
                                (double)creal(a[j*dof*dof+k*dof+l]),
                                (double)cimag(a[j*dof*dof+k*dof+l]));
#else
                        fprintf(outfile, "%ld\t%ld\t%.12lg\t%.12lg\n",
                                (long)(l2g[i]-1)*dof + (k+1),
                                (long)(ja[j]-1)*dof + (l+1),
                                (double)creal(a[j*dof*dof+k*dof+l]),
                                (double)cimag(a[j*dof*dof+k*dof+l]));
#endif
#else

#ifdef DUMP_BY_NODE
                        fprintf(outfile, "%ld (%ld)\t%ld (%ld)\t%.12lg\n",
                                (long)l2g[i], (long)(k+1),
                                (long)ja[j],  (long)(l+1),
                                (double)a[j*dof*dof+k*dof+l]);
#else
                        fprintf(outfile, "%ld\t%ld\t%.12lg\n",
                                (long)(l2g[i]-1)*dof + (k+1),
                                (long)(ja[j]-1)*dof + (l+1),
                                (double)a[j*dof*dof+k*dof+l]);
#endif
#endif
                      }
                    else
                      {
                        fprintf(outfile, "%ld\t%ld\n",
                                (long)((l2g[i]-1)*dof+k+1),
                                (long)((ja[j]-1)*dof+l+1));
                      }
                  }
                else
                  {
                    if (a != NULL)
                      {
#ifdef CPLX
                        fprintf(outfile, "%ld\t%ld\t%.12lg\t%.12lg\n",
                                (long)((i)*dof+k+1),
                                (long)((ja[j]-1)*dof+l+1),
                                (double)creal(a[j*dof*dof+k*dof+l]),
                                (double)cimag(a[j*dof*dof+k*dof+l]));
#else
                        fprintf(outfile, "%ld\t%ld\t%.12lg\n",
                                (long)((i)*dof+k+1),
                                (long)((ja[j]-1)*dof+l+1),
                                (double)a[j*dof*dof+k*dof+l]);
#endif
                      }
                    else
                      {
                        fprintf(outfile, "%ld\t%ld\n",
                                (long)((i)*dof+k+1),
                                (long)((ja[j]-1)*dof+l+1));
                      }
                  }
              }
          }

      }
  fclose(outfile);
  if (rhs != NULL)
    {
      sprintf(file2, "%s_rhs%d", file, myrank);
      outfile = fopen(file2, "w");
      for (i=0; i<n; i++)
        {

          for (j = 0; j < dof; j++)
            {
              if (l2g != NULL)
                {
#ifdef CPLX
                  fprintf(outfile, "%ld %lg %lg\n", (long)((l2g[i]-1)*dof+j+1),
                          (double)creal(rhs[i*dof+j]),
                          (double)cimag(rhs[i*dof+j]));
#else
                  fprintf(outfile, "%ld %lg\n", (long)((l2g[i]-1)*dof+j+1),
                          (double)rhs[i*dof+j]);
#endif
                }
              else
                {
#ifdef CPLX
                  fprintf(outfile, "%ld %lg %lg\n", (long)(i*dof+j+1),
                          (double)creal(rhs[i*dof+j]),
                          (double)cimag(rhs[i*dof+j]));
#else
                  fprintf(outfile, "%ld %lg\n", (long)(i*dof+j+1),
                          (double)rhs[i*dof+j]);
#endif
                }
            }
        }

      fclose(outfile);
    }

  return 0;
}

/*
 *  Function: cscd_load
 *
 *  Loads a distributed csc from disk.
 *  if filename is NULL then filename = cscd_matrix.
 *
 *  Parameters:
 *    n           - Number of local columns
 *    ia          - First cscd starting index of each column in *ja* and *a*
 *    ja          - Row of each element in first CSCD
 *    a           - value of each cscd in first CSCD (can be NULL)
 *    rhs         - Right hand side.
 *    l2g         - local 2 global column numbers for first cscd
 *    filename    - name of the files.
 *    comm        - MPI communicator
 */
#ifndef cscd_load
#  error "redefinition required"
#endif
int cscd_load(PASTIX_INT *n, PASTIX_INT ** ia, PASTIX_INT ** ja, PASTIX_FLOAT ** a, PASTIX_FLOAT ** rhs,
              PASTIX_INT ** l2g, const char * filename, MPI_Comm mpi_comm)
{
  FILE * infile;
  char   file[64];
  int    myrank, nbproc;
  int    nbprocexpected;
  PASTIX_INT    i;
  int    valbool = 0;
  int    rhsbool = 0;
  long   tmp1,tmp2,tmp3,tmp4;
  double   tmpflt1,tmpflt2,tmpflt3,tmpflt4;
#ifdef TYPE_COMPLEX
  double   tmpflt5,tmpflt6,tmpflt7,tmpflt8;
#endif
  long   nnz;
  (void)mpi_comm;

  MPI_Comm_rank(mpi_comm, &myrank);
  MPI_Comm_size(mpi_comm, &nbproc);

  if (filename == NULL)
    sprintf(file, "cscd_matrix");
  else
    sprintf(file, "%s", filename);

  if (myrank == 0)
    {
      infile = fopen(file, "r");
      if (3 != fscanf(infile, "%ld %ld %ld\n", &tmp1, &tmp2, &tmp3)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }
      nbprocexpected = tmp1;
      valbool        = tmp2;
      rhsbool        = tmp3;
      if (nbprocexpected != nbproc)
        {
          errorPrint("Number of processors should be %d\n",nbprocexpected);
          return EXIT_FAILURE;
        }
      for(i=0; i<myrank; i++)
        if (1 != fscanf(infile, "%s\n", file)){
          errorPrint("CSCD badly formated");
          return EXIT_FAILURE;
        }

      fclose(infile);
    }
  fprintf(stdout,"%s", file);
  infile = fopen(file, "r");
  if (2 != fscanf(infile, "%ld %ld\n", &tmp1, &tmp2)){
    errorPrint("CSCD badly formated");
    return EXIT_FAILURE;
  }

  *n      = tmp2;
  if (1 != fscanf(infile, "%ld\n", &nnz)){
    errorPrint("CSCD badly formated");
    return EXIT_FAILURE;
  }


  /* Copie de Loc2Glob */
  *l2g = NULL;
  MALLOC_INTERN(*l2g, *n, PASTIX_INT);
  for (i=0; i<*n-4+1; i+=4)
    {
      if (4 != fscanf(infile, "%ld %ld %ld %ld", &tmp1, &tmp2, &tmp3, &tmp4)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*l2g)[i  ] = tmp1;
      (*l2g)[i+1] = tmp2;
      (*l2g)[i+2] = tmp3;
      (*l2g)[i+3] = tmp4;
    }
  switch (*n - i )
    {
    case 3:
      if (3 != fscanf(infile, "%ld %ld %ld", &tmp1, &tmp2, &tmp3)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*l2g)[i  ] = tmp1;
      (*l2g)[i+1] = tmp2;
      (*l2g)[i+2] = tmp3;
      break;
    case 2:
      if (2 != fscanf(infile, "%ld %ld", &tmp1, &tmp2)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*l2g)[i  ] = tmp1;
      (*l2g)[i+1] = tmp2;
      break;
    case 1:
      if (1 != fscanf(infile, "%ld", &tmp1)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*l2g)[i  ] = tmp1;
      break;
    }

  /* Copie de IA */
  *ia = NULL;
  MALLOC_INTERN(*ia, *n+1, PASTIX_INT);
  for (i=0; i<*n+1+1-4; i+=4)
    {
      if (4 != fscanf(infile, "%ld %ld %ld %ld", &tmp1, &tmp2, &tmp3, &tmp4)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*ia)[i  ] = tmp1;
      (*ia)[i+1] = tmp2;
      (*ia)[i+2] = tmp3;
      (*ia)[i+3] = tmp4;
    }
  switch (*n +1 - i)
    {
    case 3:
      if (3 != fscanf(infile, "%ld %ld %ld", &tmp1, &tmp2, &tmp3)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*ia)[i  ] = tmp1;
      (*ia)[i+1] = tmp2;
      (*ia)[i+2] = tmp3;
      break;
    case 2:
      if (2 != fscanf(infile, "%ld %ld", &tmp1, &tmp2)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*ia)[i  ] = tmp1;
      (*ia)[i+1] = tmp2;
      break;
    case 1:
      if (1 != fscanf(infile, "%ld", &tmp1)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*ia)[i  ] = tmp1;
      break;
    }


  /* Copie de JA */
  (*ja) = NULL;
  MALLOC_INTERN(*ja, nnz, PASTIX_INT);
  for (i=0; i<nnz+1-4; i+=4)
    {
      if (4 != fscanf(infile, "%ld %ld %ld %ld", &tmp1, &tmp2, &tmp3, &tmp4)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*ja)[i  ] = tmp1;
      (*ja)[i+1] = tmp2;
      (*ja)[i+2] = tmp3;
      (*ja)[i+3] = tmp4;
    }

  switch (nnz - i)
    {
    case 3:
      if (3 != fscanf(infile, "%ld %ld %ld", &tmp1, &tmp2, &tmp3)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*ja)[i  ] = tmp1;
      (*ja)[i+1] = tmp2;
      (*ja)[i+2] = tmp3;
      break;
    case 2:
      if (2 != fscanf(infile, "%ld %ld", &tmp1, &tmp2)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*ja)[i  ] = tmp1;
      (*ja)[i+1] = tmp2;
      break;
    case 1:
      if (1 != fscanf(infile, "%ld", &tmp1)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*ja)[i  ] = tmp1;
      break;
    }

  /* Copie de Avals */
  if (valbool)
    {
      (*a) = NULL;
      MALLOC_INTERN(*a, nnz, PASTIX_FLOAT);

      for (i=0; i<nnz+1-4; i+=4)
        {
#ifdef TYPE_COMPLEX
          if (8 != fscanf(infile, "%lf %lf %lf %lf %lf %lf %lf %lf",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4,
                          &tmpflt5, &tmpflt6, &tmpflt7, &tmpflt8)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*a)[i  ] = (PASTIX_FLOAT)(tmpflt1 + I * tmpflt2);
          (*a)[i+1] = (PASTIX_FLOAT)(tmpflt3 + I * tmpflt4);
          (*a)[i+2] = (PASTIX_FLOAT)(tmpflt5 + I * tmpflt6);
          (*a)[i+3] = (PASTIX_FLOAT)(tmpflt7 + I * tmpflt8);
#else
          if (4 != fscanf(infile, "%lf %lf %lf %lf",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*a)[i  ] = (PASTIX_FLOAT)tmpflt1;
          (*a)[i+1] = (PASTIX_FLOAT)tmpflt2;
          (*a)[i+2] = (PASTIX_FLOAT)tmpflt3;
          (*a)[i+3] = (PASTIX_FLOAT)tmpflt4;
#endif
        }
      switch (nnz - i )
        {
        case 3:
#ifdef TYPE_COMPLEX
          if (6 != fscanf(infile, "%lf %lf %lf %lf %lf %lf",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4,
                          &tmpflt5, &tmpflt6)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*a)[i  ] = (PASTIX_FLOAT)(tmpflt1 + I * tmpflt2);
          (*a)[i+1] = (PASTIX_FLOAT)(tmpflt3 + I * tmpflt4);
          (*a)[i+2] = (PASTIX_FLOAT)(tmpflt5 + I * tmpflt6);
#else
          if (3 != fscanf(infile, "%lf %lf %lf",
                          &tmpflt1, &tmpflt2, &tmpflt3)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*a)[i  ] = (PASTIX_FLOAT)tmpflt1;
          (*a)[i+1] = (PASTIX_FLOAT)tmpflt2;
          (*a)[i+2] = (PASTIX_FLOAT)tmpflt3;
#endif
          break;
        case 2:
#ifdef TYPE_COMPLEX
          if (4 != fscanf(infile, "%lf %lf %lf %lf",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*a)[i  ] = (PASTIX_FLOAT)(tmpflt1 + I * tmpflt2);
          (*a)[i+1] = (PASTIX_FLOAT)(tmpflt3 + I * tmpflt4);
#else
          if (2 != fscanf(infile, "%lf %lf",
                          &tmpflt1, &tmpflt2)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*a)[i  ] = (PASTIX_FLOAT)tmpflt1;
          (*a)[i+1] = (PASTIX_FLOAT)tmpflt2;
#endif
          break;
        case 1:
#ifdef TYPE_COMPLEX
          if (2 != fscanf(infile, "%lf %lf",
                          &tmpflt1, &tmpflt2)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*a)[i  ] = (PASTIX_FLOAT)(tmpflt1 + I * tmpflt2);
#else
          if (1 != fscanf(infile, "%lf",
                          &tmpflt1)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*a)[i  ] = (PASTIX_FLOAT)tmpflt1;
#endif
          break;
        }
    }
  /* Copie de RHS */
  if (rhsbool)
    {
      (*rhs) = NULL;
      MALLOC_INTERN(*rhs, *n, PASTIX_FLOAT);
      for (i=0; i<*n+1-4; i+=4)
        {
          if (4 != fscanf(infile, "%lf %lf %lf %lf",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }

          (*rhs)[i  ] = (PASTIX_FLOAT)tmpflt1;
          (*rhs)[i+1] = (PASTIX_FLOAT)tmpflt2;
          (*rhs)[i+2] = (PASTIX_FLOAT)tmpflt3;
          (*rhs)[i+3] = (PASTIX_FLOAT)tmpflt4;
        }

      switch (*n - i)
        {
        case 3:
          if (3 != fscanf(infile, "%lf %lf %lf",
                          &tmpflt1, &tmpflt2, &tmpflt3)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }

          (*rhs)[i  ] = (PASTIX_FLOAT)tmpflt1;
          (*rhs)[i+1] = (PASTIX_FLOAT)tmpflt2;
          (*rhs)[i+2] = (PASTIX_FLOAT)tmpflt3;
          break;
        case 2:
          if (2 != fscanf(infile, "%lf %lf", &tmpflt1, &tmpflt2)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }

          (*rhs)[i  ] = (PASTIX_FLOAT)tmpflt1;
          (*rhs)[i+1] = (PASTIX_FLOAT)tmpflt2;
          break;
        case 1:
          if (1 != fscanf(infile, "%lf", &tmpflt1)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }

          (*rhs)[i  ] = (PASTIX_FLOAT)tmpflt1;
          break;
        }

    }
  fclose(infile);

  return 0;
}
