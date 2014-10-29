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
#include "common_pastix.h"
#include "redefine_functions.h"
#include "sparse_gemm.h"
#include "sopalin_compute.h"

int sparse_gemm_cpu( char * transa, char * transb,
                     int m, int n, int k,
                     PASTIX_FLOAT alpha,
                     const PASTIX_FLOAT * a, int lda,
                     const PASTIX_FLOAT * b, int ldb,
                     PASTIX_FLOAT beta,
                     PASTIX_FLOAT       * c, unsigned int ldc,
                     int blocknbr,  const int * blocktab,
                     int fblocknbr, const int * fblocktab,
                     PASTIX_FLOAT *work, int worksize)
{
  int col;
  int   C_index = 0;
  int   W_index = 0;
  const int * fblock;
  const int * lfblock = &(fblocktab[2*(fblocknbr-1)]);
  const int * block;
  const int * lblock  = &(blocktab[2*(blocknbr-1)]);
  (void)worksize;

  SOPALIN_GEMM( transa, transb,
                m, n, k,
                alpha,
                a, lda,
                b, ldb,
                0.,
                work, m);

  for (block = blocktab, fblock = fblocktab;
       block <= lblock;
       W_index += block[1]-block[0]+1, block+=2)
    {
      int size;
      while ((!(block[0] >= fblock[0] && block[1] <= fblock[1])) &&
             lfblock >= fblock)
        {
          C_index += fblock[1]-fblock[0]+1;
          fblock+=2;
        }
      if (lfblock < fblock)
        {
          errorPrint("block [%d, %d] not found in facing column.",
                     block[0], block[1]);
          return BADPARAMETER_ERR;
        }

      size = block[1]-block[0]+1;
      /* fprintf(stdout, "%d %d %d %g %d\n", C_index, block[0], fblock[0], work[W_index], size); */
      /* fprintf(stdout, "%g\n", c[C_index + block[0]-fblock[0]]); */
      for (col = 0; col < n; col++)
        SOPALIN_SCAL(size, beta, &(c[C_index + block[0]-fblock[0] + col*ldc]), 1);
      SOPALIN_GEAM("N", "N", size, n, 1.0,
                   &(work[W_index]), m,
                   &(c[C_index + block[0]-fblock[0]]),  ldc);
      /* fprintf(stdout, "%g\n", c[C_index + block[0]-fblock[0]]); */
    }
  return NO_ERR;
}
