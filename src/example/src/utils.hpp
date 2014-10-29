#ifndef UTILS_H
#define UTILS_H
#ifdef FORCE_NOMPI
#  define EXCHANGE_AX {}
#  define EXCHANGE_NORME {}
#  define IF_RANK_0 if (1)
#else
#  define EXCHANGE_AX(ax)						\
  {									\
    PaStiX::pastix_float_t * EAX_ax_rcv;				\
    EAX_ax_rcv = new PaStiX::pastix_float_t [globn];			\
    MPI_Allreduce(ax, EAX_ax_rcv, globn, MPI_PASTIX_FLOAT,      \
                  MPI_SUM,                                              \
                  MPI_COMM_WORLD);					\
    delete [] ax;							\
    ax = EAX_ax_rcv;							\
  }

#  define EXCHANGE_NORME(norme1, norme2)                        \
  {                                                             \
    PaStiX::pastix_float_t EN_norme1_rcv, EN_norme2_rcv;        \
    MPI_Allreduce(&norme1, &EN_norme1_rcv, 1,                   \
                  MPI_PASTIX_FLOAT,                     \
                  MPI_SUM, MPI_COMM_WORLD);                     \
    norme1 = EN_norme1_rcv;                                     \
    MPI_Allreduce(&norme2, &EN_norme2_rcv, 1,                   \
                  MPI_PASTIX_FLOAT,                     \
                  MPI_SUM, MPI_COMM_WORLD);                     \
    norme2 = EN_norme2_rcv;                                     \
  }
#  define IF_RANK_0  if (mpid == 0)
#endif

#define PRINT_RHS_REAL(st, rh, nn, rk, verbose)		\
  {							\
    if (verbose >= 5) {					\
      int PRHS_ii;					\
      std::cout << st << " ( Proc "<< rk<< " ) : ";	\
      for (PRHS_ii= 0; PRHS_ii< nn; PRHS_ii++)		\
	std::cout << rh[PRHS_ii] << " ";		\
      std::cout << std::endl;				\
    }							\
  }
#define PRINT_RHS_CPLX(st, rh, nn, rk, verbose)	\
  PRINT_RHS_REAL(st, rh, nn, rk, verbose)

#ifdef TYPE_COMPLEX
#  define PRINT_RHS PRINT_RHS_CPLX
#else
#  define PRINT_RHS PRINT_RHS_REAL
#endif


#define CONJ_CPLX(x)  std::conj(x)
#define CONJ_DCPLX(x) std::conj(x)
#define CONJ_REAL(x)  (x)
#ifdef TYPE_COMPLEX
#  define CONJ CONJ_CPLX
#else
#  define CONJ CONJ_REAL
#endif

#define CHECK_SOL(sol, rhs, nn, rk)                                     \
  {                                                                     \
    int CS_ii,CS_jj;                                                    \
    ax = new PaStiX::pastix_float_t[nn];				\
    memset(ax, 0, nn*sizeof(PaStiX::pastix_float_t));			\
    if (iparm[PaStiX::IPARM_TRANSPOSE_SOLVE] == PaStiX::API_NO)		\
      {                                                                 \
        for (CS_ii= 0; CS_ii < nn; CS_ii++)                             \
          {                                                             \
            for (CS_jj = colptr[CS_ii]-1; CS_jj < colptr[CS_ii+1] - 1; CS_jj++) \
              {                                                         \
                ax[rows[CS_jj]-1] += values[CS_jj]*sol[CS_ii];          \
                if ((MTX_ISSYM(type) == 1) && (CS_ii != (rows[CS_jj]-1))) \
                  {                                                     \
                    ax[CS_ii] += values[CS_jj]*sol[rows[CS_jj]-1];      \
                  }                                                     \
                if ((MTX_ISHER(type) == 1) && (CS_ii != (rows[CS_jj]-1))) \
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
            for (CS_jj = colptr[CS_ii]-1; CS_jj < colptr[CS_ii+1] - 1; CS_jj++) \
              {                                                         \
                ax[CS_ii] += values[CS_jj]*sol[rows[CS_jj]-1];          \
                if ((MTX_ISSYM(type) == 1) && (CS_ii != (rows[CS_jj]-1))) \
                  {                                                     \
                    ax[rows[CS_jj]-1] = ax[rows[CS_jj]-1] +             \
                      values[CS_jj]*sol[CS_ii];                         \
                  }                                                     \
                if ((MTX_ISHER(type) == 1) && (CS_ii != (rows[CS_jj]-1))) \
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
	std::complex<double> cplx;					\
	cplx = (ax[CS_ii] - rhs[CS_ii]) * CONJ(ax[CS_ii] - rhs[CS_ii]);	\
        norme1 += (double)real(cplx);					\
	cplx = rhs[CS_ii] * CONJ(rhs[CS_ii]);				\
        norme2 += (double)real(cplx);					\
      }                                                                 \
    if (rk == 0)                                                        \
      {									\
	std::cout << "Precision : ||ax -b||/||b|| = ";			\
	std::cout << sqrt(norme1/norme2) << std::endl;			\
      }									\
    delete [] ax;							\
  }

#define CHECK_DIST_SOL(colptr2, rows2, values2, rhs2, ncol2,            \
                       loc2glob2, globn, rhssaved_g)                    \
  {                                                                     \
                                                                        \
    PaStiX::pastix_int_t   * CDS_glob2loc;				\
    PaStiX::pastix_float_t * CDS_ax;					\
    PaStiX::pastix_float_t   CDS_norme1, CDS_norme2;			\
    PaStiX::pastix_int_t     CDS_j, CDS_iter;				\
    PaStiX::pastix_float_t * CDS_sol_g, *CDS_sol_g_recv;		\
    CDS_glob2loc = new PaStiX::pastix_int_t[globn];			\
    for (CDS_iter = 0; CDS_iter < globn; CDS_iter++)                    \
      CDS_glob2loc[CDS_iter] = -1;                                      \
    for (CDS_iter = 0; CDS_iter < ncol2; CDS_iter++)                    \
      CDS_glob2loc[loc2glob2[CDS_iter]-1] = CDS_iter+1;                 \
                                                                        \
    CDS_sol_g = new PaStiX::pastix_float_t[globn];			\
    memset(CDS_sol_g, 0, globn*sizeof(PaStiX::pastix_float_t));		\
    CDS_sol_g_recv = new PaStiX::pastix_float_t[globn];			\
    for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)                     \
    {                                                                   \
      CDS_sol_g[loc2glob2[CDS_iter]-1] = rhs2[CDS_iter];                \
    }                                                                   \
    MPI_Allreduce(CDS_sol_g, CDS_sol_g_recv, globn,                     \
                  MPI_PASTIX_FLOAT, MPI_SUM,                    \
                  MPI_COMM_WORLD);                                      \
    delete [] CDS_sol_g;						\
    CDS_sol_g =CDS_sol_g_recv;                                          \
                                                                        \
    CDS_ax = new PaStiX::pastix_float_t[globn];				\
    memset(CDS_ax, 0, globn*sizeof(PaStiX::pastix_float_t));		\
    for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)                     \
    {                                                                   \
      for (CDS_j = colptr2[CDS_iter]-1; CDS_j < colptr2[CDS_iter+1] - 1; CDS_j++) \
      {                                                                 \
        CDS_ax[rows2[CDS_j]-1] +=                                       \
          values2[CDS_j]*CDS_sol_g[loc2glob2[CDS_iter]-1];              \
        if ((MTX_ISSYM(type) == 1) &&                                   \
            (loc2glob2[CDS_iter]-1 != (rows2[CDS_j]-1)))                \
        {                                                               \
          CDS_ax[loc2glob2[CDS_iter]-1] += values2[CDS_j]*              \
            CDS_sol_g[rows2[CDS_j]-1];                                  \
        }                                                               \
        if ((MTX_ISHER(type) == 1) &&                                   \
            (loc2glob2[CDS_iter]-1 != (rows2[CDS_j]-1)))                \
        {                                                               \
          CDS_ax[loc2glob2[CDS_iter]-1] =                               \
            CDS_ax[loc2glob2[CDS_iter]-1] +                             \
            CONJ(values2[CDS_j])*                                       \
            CDS_sol_g[rows2[CDS_j]-1];                                  \
        }                                                               \
      }                                                                 \
    }                                                                   \
    delete [] CDS_sol_g;						\
    EXCHANGE_AX(CDS_ax);                                                \
    CDS_norme1= CDS_norme2 = 0;                                         \
    for (CDS_iter= 0; CDS_iter < ncol2; CDS_iter++)                     \
    {                                                                   \
      std::complex<double> cplx;					\
      cplx = ( CDS_ax[loc2glob2[CDS_iter]-1] -				\
	       rhssaved_g[loc2glob2[CDS_iter]-1] )*			\
	CONJ(CDS_ax[loc2glob2[CDS_iter]-1] -				\
	     rhssaved_g[loc2glob2[CDS_iter]-1]);			\
      CDS_norme1 += (double)( real(cplx) );				\
      cplx = (rhssaved_g[loc2glob2[CDS_iter]-1])*			\
	CONJ(rhssaved_g[loc2glob2[CDS_iter]-1]);			\
      CDS_norme2 += (double)( real(cplx) );				\
    }                                                                   \
    EXCHANGE_NORME(CDS_norme1, CDS_norme2);                             \
    IF_RANK_0 {                                                         \
      std::cout << "Precision : ||ax -b||/||b|| = ";			\
      std::cout << sqrt(CDS_norme1/CDS_norme2) << std::endl;		\
    }                                                                   \
    delete [] CDS_ax;							\
    delete [] CDS_glob2loc;						\
}
#endif /* not UTILS_H */
