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
  File: compute_trsm.c

  Computation functions.

  Pierre Ramet    : fev 2003
  Mathieu Faverge
  Xavier Lacoste

*/

/*
 * Definition of the 3 kernels for the 3 different factorizations
 *
 *    m   - The number of rows of L (and U)
 *    n   - The number of columns of L (and U)
 *    dL  - Diagonal block of matrix L
 *    dU  - Diagonal block of matrix U
 *    ldd - Leading dimension of matrix dL (and dU)
 *    L   - extra diagonal block of L
 *    U   - extra diagonal block of U
 *    ldl - Leading dimension of matrix L (and U)
 */
#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU

static inline
void API_CALL(kernel_trsm)(PASTIX_INT m, PASTIX_INT n,
               PASTIX_FLOAT *dL, PASTIX_FLOAT *dU, PASTIX_INT ldd,
               PASTIX_FLOAT *L,  PASTIX_FLOAT *U,  PASTIX_INT ldl )
{
  SOPALIN_TRSM("R", "U", "N", "N", m, n,
               fun, dL, ldd, L, ldl);

  SOPALIN_TRSM("R", "U", "N", "U", m, n,
               fun, dU, ldd, U, ldl);
}

#  else

static inline
void API_CALL(kernel_trsm)(PASTIX_INT m, PASTIX_INT n,
               PASTIX_FLOAT *dL, PASTIX_INT ldd,
               PASTIX_FLOAT *L,  PASTIX_INT ldl )
{
  SOPALIN_TRSM("R", "L",
               "T",
               "N", m, n,
         fun, dL, ldd, L, ldl);
}
#  endif /* SOPALIN_LU */

#else

static inline
void API_CALL(kernel_trsm)(PASTIX_INT m, PASTIX_INT n,
               PASTIX_FLOAT *dL, PASTIX_INT ldd,
               PASTIX_FLOAT *L,  PASTIX_FLOAT *L2, PASTIX_INT ldl )
{
  PASTIX_INT k;
#ifdef HERMITIAN
  SOPALIN_TRSM("R", "L",
               "C",
               "U", m, n,
               fun, dL, ldd, L, ldl);
#else
  SOPALIN_TRSM("R", "L",
               "T",
               "U", m, n,
               fun, dL, ldd, L, ldl);
#endif

  for (k=0; k<n; k++)
  {
    PASTIX_FLOAT alpha;
# ifdef COMPUTE
    ASSERTDBG(dL[k+k*ldd] != 0., MOD_SOPALIN);
    alpha = fun / dL[k+k*ldd];
# endif
    SOPALIN_COPY(m, &(L[ k*ldl]), iun,
            &(L2[k*m  ]), iun);
    SOPALIN_SCAL(m, alpha, &(L[k*ldl]), iun);
  }
}

#endif /* CHOL_SOPALIN */

/*********************************************************************
 *
 * Trsm executed in 1D version, or in 1D+esp:
 *   all the blocks are grouped together to do only one TRSM
 *
 *   sopalin_data - pointer on data associated to the factorization
 *   me           - Thread id
 *   c            - Block column id
 *
 */
void API_CALL(factor_trsm1d)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT c)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  PASTIX_FLOAT *dL, *L;
#ifdef SOPALIN_LU
  PASTIX_FLOAT *dU, *U;
#endif
#ifndef CHOL_SOPALIN
  Thread_Data_t *thread_data = sopalin_data->thread_data[me];
  PASTIX_FLOAT *L2;
#endif
  PASTIX_INT dima, dimb, stride, offsetD, offsetED;
  (void)me;

  offsetD  = SOLV_COEFIND(SYMB_BLOKNUM(c));
  offsetED = SOLV_COEFIND(SYMB_BLOKNUM(c)+1);

  /* diagonal column block address */
  dL = &(SOLV_COEFTAB(c)[offsetD]);

  /* first extra-diagonal bloc in column block address */
  L  = &(SOLV_COEFTAB(c)[offsetED]);

  stride = SOLV_STRIDE(c);

  /* horizontal dimension */
  dima = SYMB_LCOLNUM(c) - SYMB_FCOLNUM(c) + 1;
  /* vertical dimension */
  dimb = stride - dima;

#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU
  dU = &(SOLV_UCOEFTAB(c)[offsetD ]);
  U  = &(SOLV_UCOEFTAB(c)[offsetED]);
  API_CALL(kernel_trsm)(dimb, dima, dL, dU, stride, L,  U, stride);
#  else
  API_CALL(kernel_trsm)(dimb, dima, dL,     stride, L,     stride);
#  endif /* SOPALIN_LU */
#else
  L2 = thread_data->maxbloktab1;
  ASSERTDBG(SOLV_COEFMAX >= dimb*dima, MOD_SOPALIN);
  API_CALL(kernel_trsm)(dimb, dima, dL,     stride, L, L2, stride);
#endif /* CHOL_SOPALIN */
}


/*********************************************************************
 *
 * Trsm executed in 2D version:
 *   all the blocks are grouped together to do only one TRSM
 *
 *   sopalin_data - pointer on data associated to the factorization
 *   me           - Thread id
 *   task         - Task id
 *
 */
void API_CALL(factor_trsm2d)(Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT task)
{
  SolverMatrix  *datacode    = sopalin_data->datacode;
  PASTIX_FLOAT *dL, *L, *L2;
#ifdef SOPALIN_LU
  PASTIX_FLOAT *dU, *U, *U2;
#endif
  PASTIX_INT c, b, dima, dimb, stride, offsetED, size;
  (void)me;

  c = TASK_CBLKNUM(task);
  b = TASK_BLOKNUM(task);

  offsetED = SOLV_COEFIND(b);

  /* Address of diagonal block */
  /* Even if the data is local dL = RTASK_COEFTAB */
  dL = (PASTIX_FLOAT *)RTASK_COEFTAB(task);

  /* Adress of extra-diagonal block */
  L = &(SOLV_COEFTAB(c)[offsetED]);

  /* horizontal dimension / also leading dimension of dL since dL is square */
  dima = SYMB_LCOLNUM(c) - SYMB_FCOLNUM(c) + 1;
  /* vertical dimension */
  dimb = SYMB_LROWNUM(b) - SYMB_FROWNUM(b) + 1;

  /* Leading dimension of b */
  stride = SOLV_STRIDE(c);

  /* Size of L2/U2 */
#ifdef SOPALIN_LU
  size = dima * dimb * 2;
#else
  size = dima * dimb;
#endif
  MALLOC_INTERN(L2, size, PASTIX_FLOAT);

  print_debug(DBG_SOPALIN_ALLOC, "alloc block coeff %x\n",
        (unsigned int)(intptr_t)L2);

  STATS_ADD(size);

  /*
   * Resolution of L_kk * A_jk^t (ga, gb, stride, dimb, dima);
   */
#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU
  dU = dL + dima*dima;
  U2 = L2 + dima*dimb;
  U  = &(SOLV_UCOEFTAB(c)[offsetED]);

  API_CALL(kernel_trsm)(dimb, dima, dL, dU, dima, L,  U, stride);

  SOPALIN_LACPY(dimb, dima, L, stride, L2, dimb);
  SOPALIN_LACPY(dimb, dima, U, stride, U2, dimb);
#  else
  API_CALL(kernel_trsm)(dimb, dima, dL,     dima, L,     stride);
  SOPALIN_LACPY(dimb, dima, L, stride, L2, dimb);
#  endif /* SOPALIN_LU */
#else
  API_CALL(kernel_trsm)(dimb, dima, dL,     dima, L, L2, stride);
#endif /* CHOL_SOPALIN */

  /* Save pointer on the copy for future update and send */
  MUTEX_LOCK(&(sopalin_data->mutex_task[task]));
  STASK_COEFTAB(task) = L2;
  MUTEX_UNLOCK(&(sopalin_data->mutex_task[task]));

  /* Free E2 tasks waiting on previous pointer */
  pthread_cond_broadcast(&(sopalin_data->cond_task[task]));
}
