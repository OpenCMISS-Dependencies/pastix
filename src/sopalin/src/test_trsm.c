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

int main(int argc, char ** argv)
{
  PASTIX_INT i,j,m,n,k,f;
  int init,end,step;
  double timemoy=0;
  FILE * res;
  char filename[100];

  init=SIZEINIT;
  end=SIZEEND;
  step=SIZEPAS;

  sprintf(filename, "result_trsm_%d_%d_%d.txt",init,end,step);
  res = fopen(filename,"w");

  for (i=0;i<SIZEMAX*SIZEMAX;i++)
    matppf1[i]=(PASTIX_FLOAT)(rand()%10-5);

  k=-1;
 
  for (m=init;m<end+1;m+=step)
    {
    for (n=init;n<end+1;n+=step)
      {
	
	memcpy(matbuf1,matppf1,SIZEMAX*SIZEMAX*sizeof(PASTIX_FLOAT));
	memcpy(matbuf2,matppf1,SIZEMAX*SIZEMAX*sizeof(PASTIX_FLOAT));
	

	for (i=0;i<ITER;i++)
	  {
	    TEST_CLOCK_INIT
	    SOPALIN_TRSM("Left","Lower triangular","No transpose","Unit triangular",m,n,1.0,matbuf1,SIZEMAX,matbuf2,SIZEMAX);
	    TEST_CLOCK_STOP
	      if (i!=0)
		{
		  timemoy+=TEST_CLOCK_GET;
		} 
	    
	  }
	timemoy=timemoy/(ITER-1);
	fprintf(res,"[%d,%d,%d] : %.9lf\n", m, n, k, timemoy);

      }
    }
  fclose(res);
  exit(0);
}

