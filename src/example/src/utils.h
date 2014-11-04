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
#ifndef UTILS_H
#define UTILS_H
#ifdef FORCE_NOMPI
#  define EXCHANGE_AX {}
#  define EXCHANGE_NORME {}
#  define IF_RANK_0 if (1)
#else
#  define EXCHANGE_AX(ax)                                               \
  {                                                                     \
    pastix_float_t * EAX_ax_rcv;                                        \
    EAX_ax_rcv = malloc(globn*sizeof(pastix_float_t));                  \
    MPI_Allreduce(ax, EAX_ax_rcv, globn, MPI_PASTIX_FLOAT, MPI_SUM,     \
                  MPI_COMM_WORLD);                                      \
    free(ax);                                                           \
    ax = EAX_ax_rcv;                                                    \
  }

#  define EXCHANGE_NORME(norme1, norme2)                        \
  {                                                             \
    pastix_float_t EN_norme1_rcv, EN_norme2_rcv;                \
    MPI_Allreduce(&norme1, &EN_norme1_rcv, 1, MPI_PASTIX_FLOAT, \
                  MPI_SUM, MPI_COMM_WORLD);                     \
    norme1 = EN_norme1_rcv;                                     \
    MPI_Allreduce(&norme2, &EN_norme2_rcv, 1, MPI_PASTIX_FLOAT, \
                  MPI_SUM, MPI_COMM_WORLD);                     \
    norme2 = EN_norme2_rcv;                                     \
  }
#  define IF_RANK_0  if (mpid == 0)
#endif

#define PRINT_RHS_REAL(st, rh, nn, rk, verbose) \
  {                                             \
    if (verbose >= 5) {                         \
      int PRHS_ii;                              \
      fprintf(stdout,"%s (Proc %d) : ",st, rk); \
      for (PRHS_ii= 0; PRHS_ii< nn; PRHS_ii++)  \
        fprintf(stdout,"%.3g ",rh[PRHS_ii]);    \
      fprintf(stdout,"\n");                     \
    }                                           \
  }
#define PRINT_RHS_CPLX(st, rh, nn, rk, verbose) \
  {                                             \
    if (verbose >= 5) {                         \
      int PRHS_ii;                              \
      fprintf(stdout,"%s (Proc %d) : ",st, rk); \
      for (PRHS_ii= 0; PRHS_ii< nn; PRHS_ii++)  \
        fprintf(stdout,"(%.3g %.3g) ",          \
                creal(rh[PRHS_ii]),             \
                cimag(rh[PRHS_ii]));            \
      fprintf(stdout,"\n");                     \
    }                                           \
  }

#ifdef TYPE_COMPLEX
#  define PRINT_RHS PRINT_RHS_CPLX
#else
#  define PRINT_RHS PRINT_RHS_REAL
#endif

#define CONJ_REAL(x)  (x)
#define CONJ_CPLX(x)  conjf(x)
#define CONJ_DCPLX(x) conj(x)
#ifdef TYPE_COMPLEX
#  ifdef PREC_DOUBLE
#    define CONJ CONJ_DCPLX
#  else
#    define CONJ CONJ_CPLX
#  endif
#else
#  define CONJ CONJ_REAL
#endif

#define CHECK_SOL(sol, rhs, nn, rk)                                     \
  {                                                                     \
    int CS_ii,CS_jj;                                                    \
    ax = malloc(nn*sizeof(pastix_float_t));                             \
    memset(ax, 0, nn*sizeof(pastix_float_t));                           \
    if (iparm[IPARM_TRANSPOSE_SOLVE] == API_NO)                         \
      {                                                                 \
        for (CS_ii= 0; CS_ii < nn; CS_ii++)                             \
          {                                                             \
            for (CS_jj = colptr[CS_ii]-1;                               \
                 CS_jj < colptr[CS_ii+1] - 1;                           \
                 CS_jj++)                                               \
              {                                                         \
                ax[rows[CS_jj]-1] += values[CS_jj]*sol[CS_ii];          \
                if ((MTX_ISSYM(type) == 1) &&                           \
                    (CS_ii != (rows[CS_jj]-1)))                         \
                  {                                                     \
                    ax[CS_ii] += values[CS_jj]*sol[rows[CS_jj]-1];      \
                  }                                                     \
                if ((MTX_ISHER(type) == 1) &&                           \
                    (CS_ii != (rows[CS_jj]-1)))                         \
                  {                                                     \
                    ax[CS_ii] = ax[CS_ii] +                             \
                      CONJ(values[CS_jj])*sol[rows[CS_jj]-1];           \
                  }                                                     \
              }                                                         \
          }                                                             \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (CS_ii= 0; CS_ii < nn; CS_ii++)                             \
          {                                                             \
            for (CS_jj = colptr[CS_ii]-1;                               \
                 CS_jj < colptr[CS_ii+1] - 1;                           \
                 CS_jj++)                                               \
              {                                                         \
                ax[CS_ii] += values[CS_jj]*sol[rows[CS_jj]-1];          \
                if ((MTX_ISSYM(type) == 1) &&                           \
                    (CS_ii != (rows[CS_jj]-1)))                         \
                  {                                                     \
                    ax[rows[CS_jj]-1] = ax[rows[CS_jj]-1] +             \
                      values[CS_jj]*sol[CS_ii];                         \
                  }                                                     \
                if ((MTX_ISHER(type) == 1) &&                           \
                    (CS_ii != (rows[CS_jj]-1)))                         \
                  {                                                     \
                    ax[rows[CS_jj]-1] = ax[rows[CS_jj]-1] +             \
                      CONJ(values[CS_jj])*sol[rows[CS_jj]-1];           \
                  }                                                     \
              }                                                         \
          }                                                             \
      }                                                                 \
      norme1= norme2 = 0;                                               \
      for (CS_ii= 0; CS_ii < nn; CS_ii++)                               \
      {                                                                 \
        norme1 += (double)((ax[CS_ii] -                                 \
                            rhs[CS_ii])*CONJ(ax[CS_ii] - rhs[CS_ii]));  \
        norme2 += (double)(rhs[CS_ii] * CONJ(rhs[CS_ii]));              \
      }                                                                 \
    if (rk == 0)                                                        \
      fprintf(stdout, "Precision : ||ax -b||/||b|| = %.20lg\n",         \
              sqrt(norme1/norme2));                                     \
    free(ax);                                                           \
  }

#define CHECK_DIST_SOL(colptr2, rows2, values2, rhs2, ncol2,            \
                       loc2glob2, globn, rhssaved_g)                    \
  {                                                                     \
                                                                        \
    pastix_int_t   * CDS_glob2loc;                                      \
    pastix_float_t * CDS_ax;                                            \
    pastix_float_t   CDS_norme1, CDS_norme2;                            \
    pastix_int_t     CDS_j, CDS_iter;                                   \
    pastix_float_t * CDS_sol_g, *CDS_sol_g_recv;                        \
    CDS_glob2loc = malloc(globn*sizeof(pastix_int_t));                  \
    for (CDS_iter = 0; CDS_iter < globn; CDS_iter++)                    \
      CDS_glob2loc[CDS_iter] = -1;                                      \
    for (CDS_iter = 0; CDS_iter < ncol2; CDS_iter++)                    \
      CDS_glob2loc[loc2glob2[CDS_iter]-1] = CDS_iter+1;                 \
                                                                        \
    CDS_sol_g = malloc(globn*sizeof(pastix_float_t));                   \
    memset(CDS_sol_g, 0, globn*sizeof(pastix_float_t));                 \
    CDS_sol_g_recv = malloc(globn*sizeof(pastix_float_t));              \
    for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)                     \
    {                                                                   \
      CDS_sol_g[loc2glob2[CDS_iter]-1] = rhs2[CDS_iter];                \
    }                                                                   \
    MPI_Allreduce(CDS_sol_g, CDS_sol_g_recv, globn,                     \
                  MPI_PASTIX_FLOAT, MPI_SUM,                            \
                  MPI_COMM_WORLD);                                      \
    free(CDS_sol_g);                                                    \
    CDS_sol_g =CDS_sol_g_recv;                                          \
                                                                        \
    CDS_ax = malloc(globn*sizeof(pastix_float_t));                      \
    memset(CDS_ax, 0, globn*sizeof(pastix_float_t));                    \
    if (iparm[IPARM_TRANSPOSE_SOLVE] == API_NO)                         \
      {                                                                 \
        for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)                 \
          {                                                             \
            for (CDS_j = colptr2[CDS_iter]-1;                           \
                 CDS_j < colptr2[CDS_iter+1] - 1; CDS_j++)              \
              {                                                         \
                CDS_ax[rows2[CDS_j]-1] +=                               \
                  values2[CDS_j]*CDS_sol_g[loc2glob2[CDS_iter]-1];      \
                if ((MTX_ISSYM(type) == 1) &&                           \
                    (loc2glob2[CDS_iter]-1 != (rows2[CDS_j]-1)))        \
                  {                                                     \
                    CDS_ax[loc2glob2[CDS_iter]-1] += values2[CDS_j]*    \
                      CDS_sol_g[rows2[CDS_j]-1];                        \
                  }                                                     \
                if ((MTX_ISHER(type) == 1) &&                           \
                    (loc2glob2[CDS_iter]-1 != (rows2[CDS_j]-1)))        \
                  {                                                     \
                    CDS_ax[loc2glob2[CDS_iter]-1] =                     \
                      CDS_ax[loc2glob2[CDS_iter]-1] +                   \
                      CONJ(values2[CDS_j])*                             \
                      CDS_sol_g[rows2[CDS_j]-1];                        \
                  }                                                     \
              }                                                         \
          }                                                             \
      }                                                                 \
    else                                                                \
      {                                                                 \
        for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)                 \
          {                                                             \
            for (CDS_j = colptr2[CDS_iter]-1;                           \
                 CDS_j < colptr2[CDS_iter+1] - 1; CDS_j++)              \
              {                                                         \
                CDS_ax[loc2glob2[CDS_iter]-1] +=                        \
                  values2[CDS_j]*CDS_sol_g[rows2[CDS_j]-1];             \
                if ((MTX_ISSYM(type) == 1) &&                           \
                    (loc2glob2[CDS_iter]-1 != (rows2[CDS_j]-1)))        \
                  {                                                     \
                    CDS_ax[rows2[CDS_j]-1] += values2[CDS_j]*           \
                      CDS_sol_g[loc2glob2[CDS_iter]-1];                 \
                  }                                                     \
                if ((MTX_ISHER(type) == 1) &&                           \
                    (loc2glob2[CDS_iter]-1 != (rows2[CDS_j]-1)))        \
                  {                                                     \
                    CDS_ax[rows2[CDS_j]-1] =                            \
                      CDS_ax[rows2[CDS_j]-1] +                          \
                      CONJ(values2[CDS_j])*                             \
                      CDS_sol_g[loc2glob2[CDS_iter]-1];                 \
                  }                                                     \
              }                                                         \
          }                                                             \
      }                                                                 \
    free(CDS_sol_g);                                                    \
    EXCHANGE_AX(CDS_ax);                                                \
    CDS_norme1= CDS_norme2 = 0;                                         \
    for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)                     \
    {                                                                   \
      CDS_norme1 +=                                                     \
        (double)( ( CDS_ax[loc2glob2[CDS_iter]-1] -                     \
                    rhssaved_g[loc2glob2[CDS_iter]-1] )*                \
                  CONJ(CDS_ax[loc2glob2[CDS_iter]-1] -                  \
                       rhssaved_g[loc2glob2[CDS_iter]-1]));             \
      CDS_norme2 +=                                                     \
        (double)((rhssaved_g[loc2glob2[CDS_iter]-1])*                   \
                 CONJ(rhssaved_g[loc2glob2[CDS_iter]-1]));              \
    }                                                                   \
    EXCHANGE_NORME(CDS_norme1, CDS_norme2);                             \
    IF_RANK_0 {                                                         \
      fprintf(stdout, "Precision : ||ax -b||/||b|| = %.20lg\n",         \
              sqrt(CDS_norme1/CDS_norme2));                             \
    }                                                                   \
    free(CDS_ax);                                                       \
    free(CDS_glob2loc);                                                 \
  }
#endif /* not UTILS_H */
