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
  File: csc_utils.c

  Contains functions to manage CSC

 */
#include "common_pastix.h"
#include "csc_utils.h"

#include "sparRow.h"
#include "sort_row.h"

/*
   Function: cmp_colrow

   Used for qsort to sort arrays of PASTIX_INT following their first element.

   Returns the difference between the first element of *p1* and
   the first element of *p2*

   Parameters:
     p1 - the first array to compare
     p2 - the second array to compare
*/
int
cmp_colrow(const void *p1, const void *p2)
{
  return ((* (PASTIX_INT * const *) p1) - (*(PASTIX_INT * const *) p2));
}
/*
  Function: csc_symgraph


  Modify the CSC to a symetric graph one.
  Don't use it on a lower symetric CSC
  it would give you all the CSC upper + lower.

  External function.

  Parameters:
    n     - Number of columns/vertices
    ia    - Starting index of each column in *ja* and *a*
    ja    - Row index of each element
    a     - Value of each element,can be NULL
    newn  - New number of column
    newia - Starting index of each column in *ja* and *a*
    newja - Row index of each element
    newa  - Value of each element,can be NULL

 */
int csc_symgraph(PASTIX_INT n,     PASTIX_INT * ia,    PASTIX_INT * ja,    PASTIX_FLOAT * a,
                 PASTIX_INT *newn, PASTIX_INT **newia, PASTIX_INT **newja, PASTIX_FLOAT **newa)
{
  return csc_symgraph_int(n, ia, ja, a, newn, newia, newja, newa, API_NO);
}

/*
  Function: csc_symgraph_int


  Modify the CSC to a symetric graph one.
  Don't use it on a lower symetric CSC
  it would give you all the CSC upper + lower.

  Parameters:
    n           - Number of columns/vertices
    ia          - Starting index of each column in *ja* and *a*
    ja          - Row index of each element
    a           - Value of each element,can be NULL
    newn        - New number of column
    newia       - Starting index of each column in *ja* and *a*
    newja       - Row index of each element
    newa        - Value of each element,can be NULL
    malloc_flag - flag to indicate if function call is intern to pastix or extern.
 */
int csc_symgraph_int (PASTIX_INT n,     PASTIX_INT * ia,    PASTIX_INT * ja,    PASTIX_FLOAT * a,
                      PASTIX_INT *newn, PASTIX_INT **newia, PASTIX_INT **newja, PASTIX_FLOAT **newa,
                      int malloc_flag)
{
  PASTIX_INT * nbrEltCol = NULL; /* nbrEltCol[i] = Number of elt to add in column i */
  PASTIX_INT * cia       = NULL; /* ia of diff between good CSC and bad CSC */
  PASTIX_INT * cja       = NULL; /* ja of diff between good CSC and bad CSC */
  PASTIX_INT   nbr2add;          /* Number of elt to add */
  PASTIX_INT   itercol, iterrow, iterrow2; /* iterators */
  PASTIX_INT   l = ia[n] -1;
  PASTIX_INT   newl;

  /* Ncol=Nrow don't need change */
  *newn = n;

  MALLOC_INTERN(nbrEltCol, n, PASTIX_INT);
  /* !! Need check for malloc */

  /* Init nbrEltCol */
  for (itercol=0; itercol<n; itercol++)
    {
      nbrEltCol[itercol]=0;
    }

  /* Compute number of element by col to add for correcting the CSC */
  for (itercol=0; itercol<n; itercol++)
    {
      for (iterrow=ia[itercol]-1; iterrow<ia[itercol+1]-1; iterrow++)
        {
          if (ja[iterrow] != (itercol+1))
            {
              /* Not diagonal elt */
              /* So we have a (i,j) and we are looking for a (j,i) elt */
              /* i = itercol+1, j=ja[iterrow] */
              int rowidx=ja[iterrow]-1;
              int flag=0;

              for (iterrow2=ia[rowidx]-1; iterrow2<ia[rowidx+1]-1; iterrow2++)
                {
                  if (ja[iterrow2] == itercol+1)
                    {
                      /* Ok we found (j,i) so stop this madness */
                      flag = 1;
                      break;
                    }
                }

              if (flag==0)
                {
                  /* We never find (j,i) so increase nbrEltCol[j] */
                  (nbrEltCol[ja[iterrow]-1])++;
                }
            }
        }
    }

  /* Compute number of element to add */
  /* And cia=ia part of csc of element to add */
  /* kind of a diff between the corrected one and the original CSC */
  MALLOC_INTERN(cia, n+1, PASTIX_INT);
  /* !! Need checking good alloc) */
  nbr2add=0;
  for (itercol=0;itercol<n;itercol++)
    {
      cia[itercol]=nbr2add;
      nbr2add += nbrEltCol[itercol];
    }
  cia[n]=nbr2add;
  /*fprintf(stderr, "nbr of elt to add %ld\n", nbr2add);*/

  if (nbr2add != 0)
    {
      /* Build cja */
      /* like cia, cja is ja part of diff CSC */
      MALLOC_INTERN(cja, nbr2add, PASTIX_INT);
      /* !! again we need check of memAlloc */

      /* We walkthrough again the csc */
      for (itercol=0;itercol<n;itercol++)
        {
          for (iterrow=ia[itercol]-1;iterrow<ia[itercol+1]-1;iterrow++)
            {
              if (ja[iterrow] != itercol+1)
                {
                  /* we find (i,j) need to find (j,i) */
                  int rowidx=ja[iterrow]-1;
                  int flag=0;

                  for (iterrow2=ia[rowidx]-1;iterrow2<ia[rowidx+1]-1;iterrow2++)
                    {
                      if (ja[iterrow2] == itercol+1)
                        {
                          /* find (j,i) */
                          flag=1;
                          break;
                        }
                    }

                  if (flag==0)
                    {
                      /* We don't find (j,i) so put in diff CSC (cia,cja,0) */
                      PASTIX_INT index=ja[iterrow]-1;
                      /* cia[index] = index to put in cja the elt */
                      cja[cia[index]] = itercol+1;
                      (cia[index])++;
                    }
                }
            }
        }

      /* Restore cia */
      cia[0]=0;
      for (itercol=0;itercol<n;itercol++)
        {
          cia[itercol+1]=cia[itercol]+nbrEltCol[itercol];
        }

      memFree_null(nbrEltCol);

      /* Build corrected csc */
      newl = l+nbr2add;

      if (malloc_flag == API_NO)
        {
          /* Ici on a des malloc car le free est externe */
          MALLOC_EXTERN(*newia, n+1, PASTIX_INT);
          MALLOC_EXTERN(*newja, newl, PASTIX_INT);
          if (a != NULL)
            MALLOC_EXTERN(*newa, newl, PASTIX_FLOAT);
        }
      else
        {
          /* Ici on a des memAlloc car le free est interne */
          MALLOC_INTERN(*newia, n+1, PASTIX_INT);
          MALLOC_INTERN(*newja, newl, PASTIX_INT);
          if (a != NULL)
            MALLOC_INTERN(*newa, newl, PASTIX_FLOAT);
        }
      iterrow2 = 0; /* iterator of the CSC diff */
      for (itercol=0; itercol<n; itercol++)
        {
          (*newia)[itercol] = ia[itercol]+iterrow2;
          for (iterrow=ia[itercol]-1;iterrow<ia[itercol+1]-1;iterrow++)
            {
              /* we add the new elt with respect of order in row */
              while ((iterrow2<cia[itercol+1]) &&
                     (ja[iterrow] > cja[iterrow2]))
                {
                  /* we have elt(s) to add with a row lower than ja[iterrow] */
                  (*newja)[iterrow+iterrow2]=cja[iterrow2];
                  if (a != NULL)
                    (*newa)[iterrow+iterrow2]=0.;
                  iterrow2++;
                }

              /* Put the elt from the origin CSC */
              (*newja)[iterrow+iterrow2] = ja[iterrow];
              if (a != NULL)
                (*newa)[iterrow+iterrow2] = a[iterrow];
            }

          /* Since we put elt with a row lower than elt in the origin CSC */
          /* We could have some elt to add after the last elt in the column */
          while(iterrow2<cia[itercol+1])
            {
              (*newja)[iterrow+iterrow2]=cja[iterrow2];
              if (a != NULL)
                (*newa)[iterrow+iterrow2]=0.;
              iterrow2++;
            }
        }

      (*newia)[n]=ia[n]+iterrow2;
      memFree_null(cja);

    }
  else
    {
      /* No correction to do */
      memFree_null(nbrEltCol);
      newl = l;
      if (malloc_flag == API_NO)
        {
          /* ici on a des mallocs car le free est externe */
          MALLOC_EXTERN(*newia, n+1, PASTIX_INT);
          MALLOC_EXTERN(*newja, l, PASTIX_INT);
          if (a != NULL)
            MALLOC_EXTERN(*newa, l, PASTIX_FLOAT);
        }
      else
        {
          /* ici on a des memAllocs car le free est interne */
          MALLOC_INTERN(*newia, n+1, PASTIX_INT);
          MALLOC_INTERN(*newja, l,   PASTIX_INT);
          if (a != NULL)
            MALLOC_INTERN(*newa, l, PASTIX_FLOAT);
        }
      memcpy((*newia), ia, (n+1)*sizeof(PASTIX_INT));
      memcpy((*newja), ja, l * sizeof(PASTIX_INT));
      if (a != NULL)
        memcpy((*newa) , a , l * sizeof(PASTIX_FLOAT));
    }
  memFree_null(cia);

  return EXIT_SUCCESS;
}


/**
    Function: csc_noDiag

    Supress diagonal term.
    After this call, *ja* can be reallocated to *ia[n] -1*.

    Parameters:
      n  - size of the matrix.
      ia - Index in *ja* and *a* of the first element of each column
      ja - row of each element
      a  - value of each element, can be set to NULL

    Returns:
      ia and ja tabulars modified.
*/
void csc_noDiag(PASTIX_INT baseval, PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a)
{
  PASTIX_INT i, j;
  PASTIX_INT indj;
  PASTIX_INT *old_ia = NULL;

  MALLOC_INTERN(old_ia, n+1, PASTIX_INT)
  memcpy(old_ia, ia, sizeof(PASTIX_INT)*(n+1));

  ASSERT(ia[0]==baseval,MOD_SOPALIN);

  indj = 0;
  /*fprintf(stdout, "NNZ with diag = %ld \n", ia[n]);*/

  for(i=0;i<n;i++)
    {
      /* ia[i] = number of column already counted */
      ia[i] = indj+baseval;
      /* for each row number in each column i */
      for(j=old_ia[i];j<old_ia[i+1];j++)
        /* if element is not diagonal
           we add it in ja and we count it */
        if(ja[j-baseval] != i+baseval)
          {
            ja[indj] = ja[j-baseval];
            if (a != NULL)
              a[indj] = a [j -baseval];
            indj++;
          }
    }
  ia[n] = indj+baseval;

  /*fprintf(stdout, "NNZ without diag = %ld \n", ia[n]);*/
  memFree_null(old_ia);

}

/*
  Function: csc_check_doubles

  Check if the csc contains doubles and if correct if asked

  Assumes that the CSC is sorted.

  Assumes that the CSC is Fortran numeroted (base 1)

  Parameters:
    n         - Size of the matrix.
    colptr    - Index in *rows* and *values* of the first element of each column
    rows      - row of each element
    values    - value of each element (Can be NULL)
    dof       - Number of degrees of freedom
    flag      - Indicate if user wants correction (<API_BOOLEAN>)
    flagalloc - indicate if allocation on CSC uses internal malloc.

  Returns:
    API_YES - If the matrix contained no double or was successfully corrected.
    API_NO  - Otherwise.
*/
int csc_check_doubles(PASTIX_INT      n,
                      PASTIX_INT   *  colptr,
                      PASTIX_INT   ** rows,
                      PASTIX_FLOAT ** values,
                      int      dof,
                      int      flag,
                      int      flagalloc)
{
  PASTIX_INT     i,j,k,d;
  int     doubles = 0;
  PASTIX_INT   * tmprows = NULL;
  PASTIX_FLOAT * tmpvals = NULL;
  PASTIX_INT     index = 0;
  PASTIX_INT     lastindex = 0;

  ASSERT(values == NULL || dof > 0, MOD_SOPALIN);
  ASSERT(flag == API_NO || flag == API_YES, MOD_SOPALIN);
  ASSERT(colptr[0] == 1, MOD_SOPALIN);
  ASSERT(n >= 0, MOD_SOPALIN);

  for (i = 0; i < n; i++)
    {
      for (j = colptr[i]-1; j < colptr[i+1]-1; j = k)
        {
          (*rows)[index]   = (*rows)[j];
          if (values != NULL)
            for (d = 0; d < dof*dof; d++)
              (*values)[index*dof*dof+d] = (*values)[j*dof*dof+d];

          k = j+1;
          while (k < colptr[i+1]-1 && (*rows)[j] == (*rows)[k])
            {
              if (flag == API_NO)
                return API_NO;
              if (values != NULL)
                for (d = 0; d < dof*dof; d++)
                  (*values)[index*dof*dof+d] += (*values)[k*dof*dof+d];
              doubles++;
              k++;
            }
          index++;
        }

      colptr[i] = lastindex+1;
      lastindex = index;
    }
  if (flag == API_NO)
    return API_YES;
  ASSERT(index == colptr[n]-1-doubles, MOD_SOPALIN);
  if (doubles > 0)
    {
      colptr[n] = lastindex+1;
      if (flagalloc == API_NO)
        {
          MALLOC_EXTERN(tmprows, lastindex, PASTIX_INT);
          if (values != NULL)
            MALLOC_EXTERN(tmpvals, lastindex*dof*dof, PASTIX_FLOAT);
        }
      else
        {
          MALLOC_INTERN(tmprows, lastindex, PASTIX_INT);
          if (values != NULL)
            MALLOC_INTERN(tmpvals, lastindex*dof*dof, PASTIX_FLOAT);
        }

      memcpy(tmprows, *rows,   lastindex*sizeof(PASTIX_INT));
      if (values != NULL)
        memcpy(tmpvals, *values, lastindex*dof*dof*sizeof(PASTIX_FLOAT));
      if (flagalloc == API_NO)
        {
          free(*rows);
          if (values != NULL)
            free(*values);
        }
      else
        {
          memFree_null(*rows);
          if (values != NULL)
            memFree_null(*values);
        }
      *rows   = tmprows;
      if (values != NULL)
        *values = tmpvals;
    }
  return API_YES;

}

/*
  Function: csc_checksym

    Check if the CSC graph is symetric.

    For all local column C,

    For all row R in the column C,

    We look in column R if we have the row number C.

    If we can correct we had missing non zeros.

    Assumes that the CSC is Fortran numbered (1 based).

    Assumes that the matrix is sorted.

  Parameters:
    n        - Number of local columns
    colptr   - Starting index of each columns in *ja*
    rows     - Row of each element.
    values   - Value of each element.
    correct  - Flag indicating if we can correct the symmetry.
    alloc    - indicate if allocation on CSC uses internal malloc.
    dof      - Number of degrees of freedom.
*/
int csc_checksym(PASTIX_INT      n,
                 PASTIX_INT     *colptr,
                 PASTIX_INT    **rows,
                 PASTIX_FLOAT  **values,
                 int      correct,
                 int      alloc,
                 int      dof)
{
  PASTIX_INT            i,j,k,l,d;
  PASTIX_INT            index1;
  PASTIX_INT            index2;
  int            found;
  PASTIX_INT            toaddsize;
  PASTIX_INT            toaddcnt;
  PASTIX_INT         *  toadd      = NULL;
  PASTIX_INT         *  tmpcolptr  = NULL;
  PASTIX_INT         *  tmprows    = NULL;
  PASTIX_FLOAT       *  tmpvals    = NULL;

  /* For all local column C,
     For all row R in the column C,

     If the row number R correspond to a local column,
     We look in column R if we have the row number C.

     Else,
   */
  toaddcnt  = 0;
  toaddsize = 0;
  for (i = 0; i < n; i++)
    {
      for (j = (colptr)[i]-1; j < (colptr)[i+1]-1; j++)
        {
          if ((*rows)[j] != i+1)
            {
              /* not in diagonal */
              k = (*rows)[j];
              found = 0;
              for (l = (colptr)[k-1]-1; l < (colptr)[k-1+1]-1; l++)
                {
                  if (i+1 == (*rows)[l])
                    {
                      found = 1;
                      break;
                    }
                  if (i+1 < (*rows)[l])
                    {
                      /* The CSC is sorted */
                      found = 0;
                      break;
                    }
                }
              if (found == 0)
                {
                  if (correct == API_NO)
                    return EXIT_FAILURE;
                  else
                    {
                      if (toaddsize == 0)
                        {
                          toaddsize = n/2;
                          MALLOC_INTERN(toadd, 2*toaddsize, PASTIX_INT);
                        }
                      if (toaddcnt >= toaddsize)
                        {
                          toaddsize += toaddsize/2 + 1;
                          if (NULL ==
                              (toadd =
                               (PASTIX_INT*)memRealloc(toadd,
                                                2*toaddsize*sizeof(PASTIX_INT))))
                            MALLOC_ERROR("toadd");
                        }
                      toadd[2*toaddcnt]     = (*rows)[j];
                      toadd[2*toaddcnt + 1] = i+1;
                      /* fprintf(stdout, "Adding %ld, %ld\n", (long)(i+1), (long)(*rows)[j]); */
                      toaddcnt++;
                    }
                }
            }
        }
    }

  if (toaddcnt > 0)
    {

      intSort2asc1(toadd, toaddcnt);
      /* Correct here is API_YES, otherwise we would have return EXIT_FAILURE
         Or toaddcnt == 0*/
      MALLOC_INTERN(tmpcolptr, n + 1, PASTIX_INT);
      if (alloc == API_NO)
        {
          MALLOC_EXTERN(tmprows, colptr[n]-1 + toaddcnt, PASTIX_INT);
          if (values != NULL)
            {
              MALLOC_EXTERN(tmpvals, colptr[n]-1 + toaddcnt, PASTIX_FLOAT);
            }
        }
      else
        {
          MALLOC_INTERN(tmprows, colptr[n]-1 + toaddcnt, PASTIX_INT);
          if (values != NULL)
            {
              MALLOC_INTERN(tmpvals, colptr[n]-1 + toaddcnt, PASTIX_FLOAT);
            }
        }
      /* Build tmpcolptr

         tmpcolptr[i+1] will contain the number of element of
         the column i
       */
      index1 = 0;
      index2 = 0;
      for (i = 0; i <  n; i++)
        {
          tmpcolptr[i] = index2+1;
          for (j = colptr[i]-1; j < colptr[i+1]-1; j++)
            {
              if (index1 < toaddcnt &&
                  (toadd[2*index1] == i+1) &&
                  (toadd[2*index1+1] < (*rows)[j]))
                {
                  tmprows[index2] = toadd[2*index1+1];
                  if (values != NULL)
                    {
                      for (d = 0; d < dof*dof ; d++)
                        tmpvals[index2*dof*dof+d] = 0.0;
                    }
                  index1++;
                  j--; /* hack do not increment j this step of the loop */
                }
              else
                {
                  tmprows[index2] = (*rows)[j];
                  if (values != NULL)
                    {
                      for (d = 0; d < dof*dof ; d++)
                        tmpvals[index2*dof*dof+d] = (*values)[j*dof*dof+d];
                    }
                }
              index2++;
            }

          while(index1 < toaddcnt && toadd[2*index1] == i+1)
            {
              tmprows[index2] = toadd[2*index1+1];
              if (values != NULL)
                {
                  for (d = 0; d < dof*dof ; d++)
                    tmpvals[index2*dof*dof+d] = 0.0;
                }
              index1++;
              index2++;
            }
        }
      tmpcolptr[n] = index2+1;
      ASSERT((tmpcolptr[n] - 1) == (colptr[n] - 1 + toaddcnt), MOD_SOPALIN);

      memcpy(colptr, tmpcolptr, (n+1)*sizeof(PASTIX_INT));
      memFree_null(tmpcolptr);
      memFree_null(toadd);
      if (alloc == API_NO)
        {
          free(*rows);
          if (values != NULL)
            {
              free(*values);
            }
        }
      else
        {
          memFree_null(*rows);
          if (values != NULL)
            {
              memFree_null(*values);
            }
        }
      *rows   = tmprows;
      if (values != NULL)
        {
          *values = tmpvals;
        }
    }
  return EXIT_SUCCESS;
}


/*
  Function: CSC_colPerm

  Performs column permutation on a CSC

  Parameters:
    n     - Size of the matrix.
    ia    - Index of first element of each column in *ia* and *a*
    ja    - Rows of non zeros of the matrix.
    a     - Values of non zeros of the matrix.
    cperm - Permutation to perform
*/
void CSC_colPerm(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_INT *cperm)
{
  PASTIX_INT i, k;
  PASTIX_INT   *newja = NULL;
  PASTIX_INT   *newia = NULL;
  PASTIX_FLOAT *newa  = NULL;
  int numflag, numflag2;

  numflag = ia[0];
  numflag2 = 1;
  for(i=0;i<n;i++)
    if(cperm[i] == 0)
      {
        numflag2 = 0;
        break;
      }

  if(numflag2 != numflag)
    {
      errorPrint("CSC_colPerm: rperm not in same numbering than the CSC.");
      exit(-1);
    }


  if(numflag == 1)
    {
      CSC_Fnum2Cnum(ja, ia, n);
      for(i=0;i<n;i++)
        cperm[i]--;
    }

  MALLOC_INTERN(newia, n+1,   PASTIX_INT);
  MALLOC_INTERN(newja, ia[n], PASTIX_INT);
  MALLOC_INTERN(newa,  ia[n], PASTIX_FLOAT);


  newia[0] = 0;
  for(i=0;i<n;i++)
    {
#ifdef DEBUG_KASS
      ASSERT(cperm[i]>=0 && cperm[i] < n, MOD_KASS);
#endif
      newia[cperm[i]+1] = ia[i+1]-ia[i];
    }

#ifdef DEBUG_KASS
  for(i=1;i<=n;i++)
    ASSERT(newia[i] >0, MOD_KASS);
#endif

  for(i=1;i<=n;i++)
    newia[i] += newia[i-1];

#ifdef DEBUG_KASS
  ASSERT(newia[n] == ia[n], MOD_KASS);
#endif


  for(i=0;i<n;i++)
    {
      k = cperm[i];
#ifdef DEBUG_KASS
      ASSERT(newia[k+1]-newia[k] == ia[i+1]-ia[i], MOD_KASS);
#endif
      memcpy(newja + newia[k], ja + ia[i], sizeof(PASTIX_INT)*(ia[i+1]-ia[i]));
      memcpy(newa + newia[k], a + ia[i], sizeof(PASTIX_FLOAT)*(ia[i+1]-ia[i]));

    }

  memCpy(ia, newia, sizeof(PASTIX_INT)*(n+1));
  memCpy(ja, newja, sizeof(PASTIX_INT)*ia[n]);
  memCpy(a, newa, sizeof(PASTIX_FLOAT)*ia[n]);

  memFree(newia);
  memFree(newja);
  memFree(newa);
  if(numflag == 1)
    {
      CSC_Cnum2Fnum(ja, ia, n);
      for(i=0;i<n;i++)
        cperm[i]++;
    }
}


/*
  Function: CSC_colScale

  Moved from kass, only used in MC64
*/
void CSC_colScale(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_FLOAT *dcol)
{
  PASTIX_INT i, j;
  int numflag;
  PASTIX_FLOAT d;
  numflag = ia[0];

  if(numflag == 1)
    CSC_Fnum2Cnum(ja, ia, n);

  for(i=0;i<n;i++)
    {
      d = dcol[i];
      for(j=ia[i];j<ia[i+1];j++)
        {
          /***@@@ OIMBE DSCAL **/
          a[j] *= d;
        }
    }

  if(numflag == 1)
    CSC_Cnum2Fnum(ja, ia, n);
}

/*
  Function: CSC_rowScale

  Moved from kass, only used in MC64
*/
void CSC_rowScale(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_FLOAT *drow)
{
  PASTIX_INT i, j;
  int numflag;
  numflag = ia[0];

  if(numflag == 1)
    CSC_Fnum2Cnum(ja, ia, n);

  for(i=0;i<n;i++)
    {
      for(j=ia[i];j<ia[i+1];j++)
        {
#ifdef DEBUG_KASS
          ASSERT(ja[j]>0 && ja[j] <n, MOD_KASS);
#endif
          a[j] *= drow[ja[j]];
        }
    }

  if(numflag == 1)
    CSC_Cnum2Fnum(ja, ia, n);
}

/*
 * CSC_sort:
 *
 * Sort CSC columns
 *
 * Parameters:
 *   n  - Number of columns
 *   ia - Index of first element of each column in *ia*.
 *   ja - Rows of each non zeros.
 *   a  - Values of each non zeros.
*/
#ifndef CSC_sort
#error "This function must be renamed via preprocessor."
#endif
#include <assert.h>
void CSC_sort(PASTIX_INT n, PASTIX_INT *ia, PASTIX_INT *ja, PASTIX_FLOAT *a, PASTIX_INT ndof)
{
  PASTIX_INT i;
  int numflag;
  PASTIX_INT ndof2;
  void * sortptr[3];
  numflag = ia[0];
  if(numflag == 1)
    CSC_Fnum2Cnum(ja, ia, n);
  if (a != NULL)
    {

      for(i=0;i<n;i++)
        {
          sortptr[0] = &ja[ia[i]];
          sortptr[1] = &a[ia[i]*ndof2];
          ndof2 = ndof*ndof;
          sortptr[2] = &ndof2;
          qsortIntFloatAsc(sortptr, ia[i+1] - ia[i]);
        }

    }
  else
    {
      for(i=0;i<n;i++)
        intSort1asc1(&ja[ia[i]], ia[i+1] - ia[i]);

    }
  if(numflag == 1)
    CSC_Cnum2Fnum(ja, ia, n);
}

/*
  Function: CSC_Fnum2Cnum

  Convert CSC numbering from fortran numbering to C numbering.

  Parameters:
    ja - Rows of each element.
    ia - First index of each column in *ja*
    n  - Number of columns
*/
void CSC_Fnum2Cnum(PASTIX_INT *ja, PASTIX_INT *ia, PASTIX_INT n)
{
  PASTIX_INT i, j;
  for(i=0;i<=n;i++)
    ia[i]--;

  for(i=0;i<n;i++)
    for(j=ia[i];j<ia[i+1];j++)
      ja[j]--;

}

/*
  Function: CSC_Cnum2Fnum

  Convert CSC numbering from C numbering to Fortran numbering.

  Parameters:
    ja - Rows of each element.
    ia - First index of each column in *ja*
    n  - Number of columns
*/
void CSC_Cnum2Fnum(PASTIX_INT *ja, PASTIX_INT *ia, PASTIX_INT n)
{
  PASTIX_INT i, j;

  for(i=0;i<n;i++)
    for(j=ia[i];j<ia[i+1];j++)
      ja[j]++;

  for(i=0;i<=n;i++)
    ia[i]++;
}
/*
  Function: CSC_buildZerosAndNonZerosGraphs

  Separate a graph in two graphs, following
  wether the diagonal term of a column is null or not.

  Parameters:
    n, colptr, rows, values  - The initial CSC
    n_nz, colptr_nz, rows_nz - The graph of the non-null diagonal part.
    n_z, colptr_z, rows_z    - The graph of the null diagonal part.
    perm                     - Permutation to go from the first graph to
                               the one composed of the two graph concatenated.
    revperm                  - Reverse permutation tabular.
    criteria                 - Value beside which a number is said null.
*/
int CSC_buildZerosAndNonZerosGraphs(PASTIX_INT     n,
                                    PASTIX_INT    *colptr,
                                    PASTIX_INT    *rows,
                                    PASTIX_FLOAT  *values,
                                    PASTIX_INT    *n_nz,
                                    PASTIX_INT   **colptr_nz,
                                    PASTIX_INT   **rows_nz,
                                    PASTIX_INT    *n_z,
                                    PASTIX_INT   **colptr_z,
                                    PASTIX_INT   **rows_z,
                                    PASTIX_INT    *perm,
                                    PASTIX_INT    *revperm,
                                    double  criteria)
{
  PASTIX_INT  itercol;
  PASTIX_INT  iterrow;

  PASTIX_INT  ncoefszeros  = 0;
  PASTIX_INT  ncoefsnzeros = 0;
  PASTIX_INT  itercol_nz   = 0;
  PASTIX_INT  itercol_z    = 0;
  int  seen;
  PASTIX_INT  cntrows;

  for (itercol = 0; itercol <n; itercol++)
    {
      seen = 0;
      for (iterrow = colptr[itercol]-1; iterrow < colptr[itercol+1]-1; iterrow++)
        {
          if (itercol == rows[iterrow] -1 )
            {
              if (ABS_FLOAT(values[iterrow]) < criteria)
                {
                  (*n_z) ++;
                  ncoefszeros += colptr[itercol+1] - colptr[itercol];
                  seen = 1;
                }
              else
                {
                  (*n_nz) ++;
                  ncoefsnzeros += colptr[itercol+1] - colptr[itercol];
                  seen = 1;
                }
              break;
            }
        }
      if (colptr[itercol+1] == colptr[itercol]) /* empty column */
        {
          (*n_z)++;
          seen = 1;
        }
      if (seen == 0)/*  column without diag*/
        {
          (*n_z)++;
          ncoefszeros += colptr[itercol+1] - colptr[itercol];
        }
    }
  fprintf(stdout, "n_z %ld\n",  (long)*n_z);
  fprintf(stdout, "n_nz %ld\n", (long)*n_nz);
  ASSERT(*n_z+*n_nz == n, MOD_SOPALIN);
  if (*n_z == 0 || *n_nz == 0)
    return NO_ERR;

  for (itercol = 0; itercol <n; itercol++)
    {
      seen = 0;
      for (iterrow = colptr[itercol]-1; iterrow < colptr[itercol+1]-1; iterrow++)
        {
          if (itercol == rows[iterrow] -1 )
            {
              if (ABS_FLOAT(values[iterrow]) < criteria)
                {
                  perm[itercol] = (*n_nz) + itercol_z + 1;
                  itercol_z++;
                  seen = 1;
                }
              else
                {
                  perm[itercol] = itercol_nz + 1;
                  itercol_nz++;
                  seen = 1;
                }
            }
        }
      if (colptr[itercol] == colptr[itercol+1])
        { /* empty column */
          perm[itercol] = (*n_nz) + itercol_z + 1;
          itercol_z++;
          seen =1;
        }
      if (seen == 0)/*  column without diag*/
        {
          perm[itercol] = (*n_nz) + itercol_z + 1;
          itercol_z++;
        }
    }

  ASSERT(itercol_nz == *n_nz, MOD_SOPALIN);
  ASSERT(itercol_z  == *n_z, MOD_SOPALIN);
  for(itercol = 0; itercol < n; itercol++)
    revperm[perm[itercol]-1] = itercol + 1;

  MALLOC_INTERN(*colptr_nz, *n_nz + 1, PASTIX_INT);
  MALLOC_INTERN(*rows_nz, ncoefsnzeros, PASTIX_INT);
  cntrows = 0;
  for (itercol = 0; itercol <*n_nz; itercol++)
    {
      (*colptr_nz)[itercol] = cntrows + 1;
      for (iterrow = colptr[revperm[itercol]-1]-1; iterrow < colptr[revperm[itercol]-1+1]-1; iterrow++)
        {
          if (perm[rows[iterrow]-1] - 1 < *n_nz )
            {
              (*rows_nz)[cntrows] = perm[rows[iterrow]-1];
              cntrows++;
            }

        }

    }
  (*colptr_nz)[*n_nz] = cntrows+1;

  MALLOC_INTERN(*colptr_z, *n_z + 1, PASTIX_INT);
  MALLOC_INTERN(*rows_z, ncoefszeros, PASTIX_INT);
  cntrows = 0;
  for (itercol = 0; itercol <*n_z; itercol++)
    {
      (*colptr_z)[itercol] = cntrows + 1;
      for (iterrow = colptr[revperm[itercol+*n_nz]-1]-1; iterrow < colptr[revperm[itercol+*n_nz]-1+1]-1; iterrow++)
        {
          if (perm[rows[iterrow]-1]  > *n_nz )
            {
              (*rows_z)[cntrows] = perm[rows[iterrow]-1] - *n_nz;
              cntrows++;
            }

        }

    }
  (*colptr_z)[*n_z] = cntrows+1;

  return NO_ERR;
}

/*
  Function: CSC_isolate

  Isolate a list of unknowns at the end of the CSC.

  Parameters:
    n            - Number of columns.
    colptr       - Index of first element of each column in *ia*.
    rows         - Rows of each non zeros.
    n_isolate    - Number of unknow to isolate.
    isolate_list - List of unknown to isolate.
    perm         - permutation tabular.
    revperm      - reverse permutation tabular.
*/
int CSC_isolate(PASTIX_INT     n,
                PASTIX_INT    *colptr,
                PASTIX_INT    *rows,
                PASTIX_INT     n_isolate,
                PASTIX_INT    *isolate_list,
                PASTIX_INT    *perm,
                PASTIX_INT    *revperm)
{
  PASTIX_INT  itercol;
  PASTIX_INT  iterrow;

  PASTIX_INT  iter_isolate  = 0;
  PASTIX_INT  iter_non_isolate  = 0;
  PASTIX_INT *tmpcolptr = NULL;
  PASTIX_INT *tmprows   = NULL;

  if (n_isolate == 0)
    {
      errorPrintW("No schur complement\n");
      return NO_ERR;
    }
  intSort1asc1(isolate_list, n_isolate);

  for (itercol = 0; itercol <n; itercol++)
    {
      if (iter_isolate < n_isolate &&
          itercol == isolate_list[iter_isolate]-1)
        {
          revperm[n-n_isolate+iter_isolate] = itercol;
          iter_isolate++;
        }
      else
        {
          revperm[iter_non_isolate] = itercol;
          iter_non_isolate++;
        }
    }
  ASSERT(iter_non_isolate == n - n_isolate, MOD_SOPALIN);
  ASSERT(iter_isolate == n_isolate, MOD_SOPALIN);

  MALLOC_INTERN(tmpcolptr, n - n_isolate + 1,       PASTIX_INT);
  memset(tmpcolptr, 0, (n - n_isolate + 1)*sizeof(PASTIX_INT));

  for(itercol = 0; itercol < n; itercol++)
    perm[revperm[itercol]] = itercol;

  for(itercol = 0; itercol < n; itercol++)
    {
      ASSERT(perm[itercol] < n, MOD_SOPALIN);
      ASSERT(perm[itercol] > -1, MOD_SOPALIN);
    }

  tmpcolptr[0] = 1;
  for (itercol = 0; itercol <n; itercol++)
    {
      if (perm[itercol] < n - n_isolate)
        {
          for (iterrow = colptr[itercol]-1; iterrow < colptr[itercol+1]-1; iterrow ++)
            {
              /* Count edges in each column of the new graph */
              if (perm[rows[iterrow]-1] < n-n_isolate)
                {
                  tmpcolptr[perm[itercol]+1]++;
                }
            }
        }
    }

  for (itercol = 0; itercol <n - n_isolate; itercol++)
    tmpcolptr[itercol+1] += tmpcolptr[itercol];

  MALLOC_INTERN(tmprows,   tmpcolptr[n- n_isolate]-1, PASTIX_INT);
  for (itercol = 0; itercol <n; itercol++)
    {
      if (perm[itercol] < n - n_isolate)
        {
          for (iterrow = colptr[itercol]-1; iterrow < colptr[itercol+1]-1; iterrow ++)
            {
              /* Count edges in each column of the new graph */
              if (perm[rows[iterrow]-1] < n-n_isolate)
                {
                  tmprows[tmpcolptr[perm[itercol]]-1] = perm[rows[iterrow]-1]+1;
                  tmpcolptr[perm[itercol]]++;
                }
            }
        }
    }


  /* restore tmpcolptr */

  for (itercol = 1; itercol <n - n_isolate ; itercol++)
    {
      tmpcolptr[n - n_isolate - itercol] = tmpcolptr[n - n_isolate - itercol - 1];
    }
  tmpcolptr[0] = 1;



  ASSERT(colptr[n] >= tmpcolptr[n-n_isolate], MOD_SOPALIN);
  memcpy(colptr, tmpcolptr, (n-n_isolate + 1)*sizeof(PASTIX_INT));
  memcpy(rows,   tmprows,   (colptr[n-n_isolate]-1)*sizeof(PASTIX_INT));

  memFree_null(tmpcolptr);
  memFree_null(tmprows);

  return NO_ERR;
}

/*
  Function: csc_save

  Save a csc on disk.

  Parameters:
    n       - number of columns
    colptr  - First cscd starting index of each column in *ja* and *a*
    rows    - Row of each element in first CSCD
    values  - value of each cscd in first CSCD (can be NULL)
    dof     - Number of degrees of freedom
    outfile - Output stream.

  Return:
    NO_ERR

*/
int csc_save(PASTIX_INT      n,
             PASTIX_INT    * colptr,
             PASTIX_INT    * rows,
             PASTIX_FLOAT  * values,
             int      dof,
             FILE   * outfile)
{

  PASTIX_INT i;
  fprintf(outfile, "%ld %ld %d\n", (long)n, (long)dof, ((values == NULL)?0:1));
  /* Copie du colptr */
  for (i=0; i<n+1; i++)
    {
      fprintf(outfile, "%ld ", (long)colptr[i]);
      if (i%4 == 3) fprintf(outfile, "\n");
    }
  if ((i-1)%4 !=3) fprintf(outfile, "\n");

  /* Copie de JA */
  for (i=0; i<colptr[n]-1; i++)
    {
      fprintf(outfile, "%ld ", (long)rows[i]);
      if (i%4 == 3) fprintf(outfile, "\n");
    }
  if ((i-1)%4 !=3) fprintf(outfile, "\n");

  /* Copie de Avals */

  if (values != NULL)
    {
      for (i=0; i<(colptr[n]-1)*dof*dof; i++)
        {
#ifdef CPLX
          fprintf(outfile, "%lg %lg ", (double)(creal(values[i])), (double)(cimag(values[i])));
#else
          fprintf(outfile, "%lg ", (double)(values[i]));
#endif
          if (i%4 == 3) fprintf(outfile, "\n");
        }
      if ((i-1)%4 !=3) fprintf(outfile, "\n");
    }
  return NO_ERR;

}

/*
  Function: csc_load

  Load a csc from disk.

  Fill *n*, *colptr*, *rows*, *values* and *dof* from *infile*.

  Parameters:
    n       - number of columns
    colptr  - First cscd starting index of each column in *ja* and *a*
    rows    - Row of each element in first CSCD
    values  - value of each cscd in first CSCD (can be NULL)
    dof     - Number of degrees of freedom
    outfile - Output stream.

  Return:
    NO_ERR

*/
int csc_load(PASTIX_INT    *  n,
             PASTIX_INT    ** colptr,
             PASTIX_INT    ** rows,
             PASTIX_FLOAT  ** values,
             int    *  dof,
             FILE   *  infile)
{
  int  hasval;
  long tmp1, tmp2, tmp3, tmp4;
  double tmpflt1, tmpflt2, tmpflt3, tmpflt4;
#ifdef TYPE_COMPLEX
  double tmpflt5, tmpflt6, tmpflt7, tmpflt8;
#endif
  PASTIX_INT i;
  if (3 != fscanf(infile, "%ld %ld %d\n", &tmp1, &tmp2, &hasval)){
    errorPrint("CSCD badly formated");
    return EXIT_FAILURE;
  }
  *n   = (PASTIX_INT)tmp1;
  *dof = (int)tmp2;
  /* Copie de IA */
  *colptr = NULL;
  MALLOC_INTERN(*colptr, *n+1, PASTIX_INT);
  for (i=0; i<*n+1+1-4; i+=4)
    {
      if (4 != fscanf(infile, "%ld %ld %ld %ld", &tmp1, &tmp2, &tmp3, &tmp4)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*colptr)[i  ] = tmp1;
      (*colptr)[i+1] = tmp2;
      (*colptr)[i+2] = tmp3;
      (*colptr)[i+3] = tmp4;
    }
  switch (*n +1 - i)
    {
    case 3:
      if (3 != fscanf(infile, "%ld %ld %ld", &tmp1, &tmp2, &tmp3)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*colptr)[i  ] = tmp1;
      (*colptr)[i+1] = tmp2;
      (*colptr)[i+2] = tmp3;
      break;
    case 2:
      if (2 != fscanf(infile, "%ld %ld", &tmp1, &tmp2)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*colptr)[i  ] = tmp1;
      (*colptr)[i+1] = tmp2;
      break;
    case 1:
      if (1 != fscanf(infile, "%ld", &tmp1)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*colptr)[i  ] = tmp1;
      break;
    }


  /* Copie de JA */
  (*rows) = NULL;
  MALLOC_INTERN(*rows, (*colptr)[*n]-(*colptr)[0], PASTIX_INT);
  for (i=0; i< (*colptr)[*n]-(*colptr)[0]+1-4; i+=4)
    {
      if (4 != fscanf(infile, "%ld %ld %ld %ld", &tmp1, &tmp2, &tmp3, &tmp4)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*rows)[i  ] = tmp1;
      (*rows)[i+1] = tmp2;
      (*rows)[i+2] = tmp3;
      (*rows)[i+3] = tmp4;
    }

  switch ( (*colptr)[*n]-(*colptr)[0] - i)
    {
    case 3:
      if (3 != fscanf(infile, "%ld %ld %ld", &tmp1, &tmp2, &tmp3)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*rows)[i  ] = tmp1;
      (*rows)[i+1] = tmp2;
      (*rows)[i+2] = tmp3;
      break;
    case 2:
      if (2 != fscanf(infile, "%ld %ld", &tmp1, &tmp2)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*rows)[i  ] = tmp1;
      (*rows)[i+1] = tmp2;
      break;
    case 1:
      if (1 != fscanf(infile, "%ld", &tmp1)){
        errorPrint("CSCD badly formated");
        return EXIT_FAILURE;
      }

      (*rows)[i  ] = tmp1;
      break;
    }

  /* Copie de Avals */
  if (hasval)
    {
      (*values) = NULL;

      MALLOC_INTERN(*values,  (*colptr)[*n]-(*colptr)[0], PASTIX_FLOAT);

      for (i=0; i< (*colptr)[*n]-(*colptr)[0]+1-4; i+=4)
        {
#ifdef TYPE_COMPLEX
          if (8 != fscanf(infile, "%lg %lg %lg %lg %lg %lg %lg %lg",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4,
                          &tmpflt5, &tmpflt6, &tmpflt7, &tmpflt8)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (PASTIX_FLOAT)(tmpflt1 + I * tmpflt2);
          (*values)[i+1] = (PASTIX_FLOAT)(tmpflt3 + I * tmpflt4);
          (*values)[i+2] = (PASTIX_FLOAT)(tmpflt5 + I * tmpflt6);
          (*values)[i+3] = (PASTIX_FLOAT)(tmpflt7 + I * tmpflt8);
#else
          if (4 != fscanf(infile, "%lg %lg %lg %lg",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (PASTIX_FLOAT)tmpflt1;
          (*values)[i+1] = (PASTIX_FLOAT)tmpflt2;
          (*values)[i+2] = (PASTIX_FLOAT)tmpflt3;
          (*values)[i+3] = (PASTIX_FLOAT)tmpflt4;
#endif
        }
      switch ( (*colptr)[*n]-(*colptr)[0] - i )
        {
        case 3:
#ifdef TYPE_COMPLEX
          if (6 != fscanf(infile, "%lg %lg %lg %lg %lg %lg",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4,
                          &tmpflt5, &tmpflt6)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (PASTIX_FLOAT)(tmpflt1 + I * tmpflt2);
          (*values)[i+1] = (PASTIX_FLOAT)(tmpflt3 + I * tmpflt4);
          (*values)[i+2] = (PASTIX_FLOAT)(tmpflt5 + I * tmpflt6);
#else
          if (3 != fscanf(infile, "%lg %lg %lg",
                          &tmpflt1, &tmpflt2, &tmpflt3)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (PASTIX_FLOAT)tmpflt1;
          (*values)[i+1] = (PASTIX_FLOAT)tmpflt2;
          (*values)[i+2] = (PASTIX_FLOAT)tmpflt3;
#endif
          break;
        case 2:
#ifdef TYPE_COMPLEX
          if (4 != fscanf(infile, "%lg %lg %lg %lg",
                          &tmpflt1, &tmpflt2, &tmpflt3, &tmpflt4)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (PASTIX_FLOAT)(tmpflt1 + I * tmpflt2);
          (*values)[i+1] = (PASTIX_FLOAT)(tmpflt3 + I * tmpflt4);
#else
          if (2 != fscanf(infile, "%lg %lg",
                          &tmpflt1, &tmpflt2)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (PASTIX_FLOAT)tmpflt1;
          (*values)[i+1] = (PASTIX_FLOAT)tmpflt2;
#endif
          break;
        case 1:
#ifdef TYPE_COMPLEX
          if (2 != fscanf(infile, "%lg %lg",
                          &tmpflt1, &tmpflt2)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (PASTIX_FLOAT)(tmpflt1 + I * tmpflt2);
#else
          if (1 != fscanf(infile, "%lg",
                          &tmpflt1)){
            errorPrint("CSCD badly formated");
            return EXIT_FAILURE;
          }
          (*values)[i  ] = (PASTIX_FLOAT)tmpflt1;
#endif
          break;
        }
    }
  return NO_ERR;
}
