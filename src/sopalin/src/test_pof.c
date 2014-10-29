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
#include "commun.h"

#define PASTIX_potrf_block API_CALL(PASTIX_potrf_block)

int main(int argc, char ** argv)
{
  PASTIX_INT i,j,m,n,k,f;
  int init,end,step;
  PASTIX_INT nbpivot;
  double timemoy=0;
  double critere=1e-12;
  FILE * res;
  char filename[100];
  
  init=SIZEINIT;
  end=SIZEEND;
  step=SIZEPAS;

  sprintf(filename, "result_pof_%d_%d_%d.txt",init,end,step);
  res = fopen(filename,"w"); 
  n=-1;
  k=-1;

  for (m=init;m<end+1;m+=step)
    {
      
      for (i=0;i<m;i++)
        for (j=0;j<m;j++)
          matppf1[i*SIZEMAX+j]=1.0;
      for (i=0;i<m;i++)
        matppf1[i*SIZEMAX+i]=2.0*(PASTIX_FLOAT)m;
      
      memcpy(matbuf1,matppf1,SIZEMAX*SIZEMAX*sizeof(PASTIX_FLOAT));
      
      for (i=0;i<ITER;i++)
	{
	  TEST_CLOCK_INIT
          PASTIX_potrf_block(matbuf1, m, SIZEMAX, &nbpivot, critere);
	  TEST_CLOCK_STOP

	    if (i!=0)
	      {
		timemoy+=TEST_CLOCK_GET;
	      }
	}
      timemoy=timemoy/(ITER-1);
      fprintf(res,"[%d,%d,%d] : %.9lf\n", m, n, k, timemoy); 
    }
  fclose(res);
  exit(0);
}



