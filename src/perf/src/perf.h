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
#ifndef PERF_H
#define PERF_H

#define PERF_MODEL "AMD 6180  MKL"

#define OPS_SCAL(i)     (i)
#define OPS_COPY(i)     (i)
/*
#define OPS_AXPY(i)     (i)
#define OPS_SWAP(i)     (i)
#define OPS_GEMV(i,j)   ((i)*(j))
#define OPS_TRSV(i,j)   (((i)*(i-1)*(j))/2.)
#define OPS_GER(i,j)    ((i)*(j))
*/
#define OPS_GEMM(i,j,k) (2.*(i)*(j)*(k))
#define OPS_TRSM(i,j)   ((i)*(i-1)*(j))
#define OPS_PPF(i)      ((2.*(i)*(i)*(i)+3.*(i)*(i)-5.*(i))/6.)
#define OPS_GEAM(i,j)   ((i)*(j))


/**GEMM**/
#define GEMM_A  2.429169e-10
#define GEMM_B  2.724804e-10
#define GEMM_C  1.328900e-09
#define GEMM_D  1.148989e-07
#define GEMM_E  -2.704179e-10
#define GEMM_F  1.216278e-06
#define PERF_GEMM(i,j,k) (GEMM_A*(double)(i)*(double)(j)*(double)(k)+GEMM_B*(double)(i)*(double)(j)+GEMM_C*(double)(j)*(double)(k)+GEMM_D*(double)(i)+GEMM_E*(double)(j)+GEMM_F)


/**GEAM**/
#define GEAM_A  1.358111e-09
#define GEAM_B -4.416379e-09
#define GEAM_C 2.270780e-08
#define GEAM_D  -3.335563e-07
#define PERF_GEAM(i,j)   (GEAM_A*(double)(i)*(double)(j)+GEAM_B*(double)(i)+GEAM_C*(double)(j)+GEAM_D)

/**TRSM**/
#define TRSM_A 2.626177e-10
#define TRSM_B 3.976198e-08
#define TRSM_C 3.255168e-06
#define PERF_TRSM(i,j)   (TRSM_A*(double)(i)*(double)(i)*(double)(j)+TRSM_B*(double)(i)+TRSM_C)

/**POF**/
#define POF_A 2.439599e-11
#define POF_B 1.707006e-08
#define POF_C -1.469893e-07
#define POF_D 4.071507e-07
#define PERF_POF(i) (POF_A*(double)(i)*(double)(i)*(double)(i)+POF_B*(double)(i)*(double)(i)+POF_C*(double)(i)+POF_D)

/**PPF**/
#define PPF_A 2.439599e-11
#define PPF_B 1.707006e-08
#define PPF_C -1.469893e-07
#define PPF_D 4.071507e-07
#define PERF_PPF(i) (PPF_A*(double)(i)*(double)(i)*(double)(i)+PPF_B*(double)(i)*(double)(i)+PPF_C*(double)(i)+PPF_D)

/**SCAL**/
#define SCAL_A 4.371793e-10
#define SCAL_B 2.052399e-07
#define PERF_SCAL(i) (SCAL_A*(double)(i)+SCAL_B)

/**COPY**/
#define COPY_A 9.177969e-10
#define COPY_B 2.266129e-07
#define PERF_COPY(i) (COPY_A*(double)(i)+COPY_B)

/**AXPY**/
#define AXPY_A 4.620143e-10
#define AXPY_B 2.101008e-07
#define PERF_AXPY(i) (AXPY_A*(double)(i)+AXPY_B)

/**GEMV**/
#define GEMV_A  6.192657e-10
#define GEMV_B -2.884799e-09
#define GEMV_C 7.594831e-10
#define GEMV_D  3.575035e-07
#define PERF_GEMV(i,j)   (GEMV_A*(double)(i)*(double)(j)+GEMV_B*(double)(i)+GEMV_C*(double)(j)+GEMV_D)

/**TRSV**/
#define TRSV_A 3.224536e-10
#define TRSV_B 1.709178e-08
#define TRSV_C 1.947268e-07
#define PERF_TRSV(i) (TRSV_A*(double)(i)*(double)(i)+TRSV_B*(double)(i)+TRSV_C)

/* en octets ...
   TIME : entre threads */

/* en octets ...
   CLUSTER : entre noeuds */

/* en octets ...
   SHARED : entre MPI shared */

/* old version compatibility
#define TIME_BANDWIDTH    1.5e-9
#define TIME_STARTUP      5.2e-6
#define CLUSTER_BANDWIDTH 5.9e-10
#define CLUSTER_STARTUP   3.9e-6
   end old                  */

#define TIME_BANDWIDTH_1    0.0
#define TIME_STARTUP_1      1e-8
#define SHARED_BANDWIDTH_1  1.0e-10
#define SHARED_STARTUP_1    0.2e-6
#define CLUSTER_BANDWIDTH_1 3.0e-10
#define CLUSTER_STARTUP_1   3.0e-6

#define TIME_BANDWIDTH_2    0.0
#define TIME_STARTUP_2      1e-8
#define SHARED_BANDWIDTH_2  3.0e-10
#define SHARED_STARTUP_2    0.4e-6
#define CLUSTER_BANDWIDTH_2 6.0e-10
#define CLUSTER_STARTUP_2   6.0e-6

#define TIME_BANDWIDTH_4    0.0
#define TIME_STARTUP_4      1e-8
#define SHARED_BANDWIDTH_4  6.0e-10
#define SHARED_STARTUP_4    0.8e-6
#define CLUSTER_BANDWIDTH_4 9.0e-10
#define CLUSTER_STARTUP_4   9.0e-6

#define TIME_BANDWIDTH_8    0.0
#define TIME_STARTUP_8      1e-8
#define SHARED_BANDWIDTH_8  6.0e-10
#define SHARED_STARTUP_8    0.8e-6
#define CLUSTER_BANDWIDTH_8 9.0e-10
#define CLUSTER_STARTUP_8   0.0e-6

#define PENALTY_STARTUP     0.0
#define PENALTY_BANDWIDTH   0.0

#endif /* PERF_H */
