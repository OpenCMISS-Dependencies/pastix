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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "nr.h"

//#define REAL float
#define SIZE 100000
#define ma 2



long int tabcoord[SIZE];
REAL tabx[SIZE];
REAL taby[SIZE];  
REAL sig[SIZE];
REAL afunc[ma+1];
int ia[ma+1];
REAL a[ma+1];


void funcs( REAL i, REAL afunc[ma+1], int ma2)
{
  
  afunc[1]=1;
  afunc[2]=tabcoord[(int)i];
  //printf("%f %f %f \n",afunc[0],afunc[1],afunc[2]);
}


int main(int argc, char * argv[])
{
  REAL total=0.0;
  REAL ecart=0.0;
/*   long double total=0.0; */
/*   long double ecart=0.0; */
  char *filename = argv[1];
  int k,i;
  FILE * res;
  FILE *out;
  FILE *perf;
  int ndat=atoi(argv[2]);
  REAL ** covar;
  REAL *chisq = (REAL *)malloc(sizeof(REAL));
  res = fopen(filename,"r");
  out = fopen("toto.txt","w");
  perf = fopen("perf.h","a");

 covar = (REAL**) malloc((ma+1) *sizeof(REAL*));
  for (i=0;i<ma+1;i++)
    covar[i]=(REAL*)malloc((ma+1) *sizeof(REAL));

  for (k=1;k<ndat+1;k++)
    {
      long int i;
      REAL tmpfloat;
      fscanf(res,"[%li,-1,-1] : %lf\n", &i, &tmpfloat);
      tabcoord[k]=i; 
      taby[k]=tmpfloat;
      sig[k]=1;
      tabx[k]=k;
      //fprintf(out,"%f %f %f\n",tabcoord[k].x1 ,tabcoord[k].x2 ,tabcoord[k].x3);
      //fprintf(out,"%f %f\n",tabx[k] , taby[k]);
    }
  for (k=1;k<ma+1;k++)
    ia[k]=1;
  

  

  printf("toto\n");
  lfit(tabx, taby, sig, ndat, a, ia, ma, covar, chisq, &funcs);
  
  printf("done\n");
  printf("coefficients : \n");
  for (k=1;k<ma+1;k++)  
    {    
      printf("%.12lf\n", a[k]);
      //total+=a[k];
    }



  //calcul de l'ecart type
  for (k=1;k<ndat+1;k++)
    {
      double abs=0.0;
      abs += a[1];
      abs += tabcoord[k]*a[2];
      fprintf(out,"k=%i ; calcul : %lf ; reel : %lf ; ", k, abs, taby[k]);
      abs = abs - taby[k];
      if (abs < 0)
	abs = - abs;
      fprintf(out,"ecart : %lf\n ", abs);

      total += abs;
      //printf("%f %f %f\n",tabcoord[k].x1 ,tabcoord[k].x2 ,tabcoord[k].x3);
    }
  

  fprintf(perf,"//**AXPY**//\n");
  fprintf(perf,"#define AXPY_A %e\n#define AXPY_B %e\n", a[2], a[1]);
  fprintf(perf,"#define PERF_AXPY(i) (AXPY_A*(double)(i)+AXPY_B)\n");
  fprintf(perf,"\n");
  


  printf("total %lf\n", total);
  ecart = total / ndat;  
  printf("ecart moyen %lf\n", ecart);

  fclose(perf); 
  fclose(out);
  fclose(res);
  return 0;
}
