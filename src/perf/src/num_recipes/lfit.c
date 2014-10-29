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
#define NRANSI
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"

//#define REAL float

void lfit(REAL x[], REAL y[], REAL sig[], int ndat, REAL a[], int ia[],
	int ma, REAL **covar, REAL *chisq, void (*funcs)(REAL, REAL [], int))
{
	void covsrt(REAL **covar, int ma, int ia[], int mfit);
	void gaussj(REAL **a, int n, REAL **b, int m);
	int i,j,k,l,m,mfit=0;
	REAL ym,wt,sum,sig2i,**beta,*afunc;

	

	beta=matrix(1,ma,1,1);
	afunc=vector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	if (mfit == 0) nrerror("lfit: no parameters to be fitted");
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=mfit;k++) covar[j][k]=0.0;
		beta[j][1]=0.0;
	}

	for (i=1;i<=ndat;i++) {
		(*funcs)(x[i],afunc,ma);
		ym=y[i];
		if (mfit < ma) {
			for (j=1;j<=ma;j++)
				if (!ia[j]) ym -= a[j]*afunc[j];
		}
		sig2i=1.0/SQR(sig[i]);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=afunc[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) covar[j][++k] += wt*afunc[m];
				beta[j][1] += ym*wt;
			}
		}
	}
	
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++)
			covar[k][j]=covar[j][k];
	printf("lfit : gaussj\n");	
	gaussj(covar,mfit,beta,1);
	printf("lfit1\n");
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) a[l]=beta[++j][1];
	printf("lfit2\n");
	*chisq=0.0;

	for (i=1;i<=ndat;i++) {
		(*funcs)(x[i],afunc,ma);
		for (sum=0.0,j=1;j<=ma;j++) sum += a[j]*afunc[j];
		*chisq += SQR((y[i]-sum)/sig[i]);
	}
	covsrt(covar,ma,ia,mfit);
	free_vector(afunc,1,ma);
	free_matrix(beta,1,ma,1,1);
}
#undef NRANSI
