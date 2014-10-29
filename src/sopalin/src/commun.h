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
#include <time.h>
#include <inttypes.h>
#include <assert.h>

#include "common_pastix.h"
#include "sopalin_compute.h"

#define SIZEEND  128
#define SIZEPAS  4
#define SIZEINIT 4

#define ITER 100

#define SIZEMAX 1000

PASTIX_FLOAT vecbuf1[SIZEMAX];
PASTIX_FLOAT vecbuf2[SIZEMAX];
PASTIX_FLOAT matbuf1[SIZEMAX][SIZEMAX];
PASTIX_FLOAT matbuf2[SIZEMAX][SIZEMAX];
PASTIX_FLOAT matbuf3[SIZEMAX][SIZEMAX];
PASTIX_FLOAT matppf1[SIZEMAX*SIZEMAX];
PASTIX_FLOAT matppf2[SIZEMAX*SIZEMAX];

void init_commun();
void init_vecteur(PASTIX_INT n,PASTIX_FLOAT *v,PASTIX_FLOAT new[SIZEMAX],PASTIX_INT inc);
void end_vecteur(PASTIX_INT n,PASTIX_FLOAT *v,PASTIX_FLOAT new[SIZEMAX],PASTIX_INT inc);
void init_matrice(char *t,PASTIX_INT n,PASTIX_INT m,PASTIX_FLOAT *a,PASTIX_FLOAT new[SIZEMAX][SIZEMAX],PASTIX_INT lda);
void end_matrice(char *t,PASTIX_INT n,PASTIX_INT m,PASTIX_FLOAT *a,PASTIX_FLOAT new[SIZEMAX][SIZEMAX],PASTIX_INT lda);

void init_vecteur(PASTIX_INT n,PASTIX_FLOAT *v,PASTIX_FLOAT new[SIZEMAX],PASTIX_INT inc)
{
  PASTIX_INT i;
  for (i=0;i<n;i++)
    new[i]=v[i*inc];
}

void end_vecteur(PASTIX_INT n,PASTIX_FLOAT *v,PASTIX_FLOAT new[SIZEMAX],PASTIX_INT inc)
{
  PASTIX_INT i;
  for (i=0;i<n;i++)
    v[i*inc]=new[i];
}

void init_matrice(char *t,PASTIX_INT n,PASTIX_INT m,PASTIX_FLOAT *a,PASTIX_FLOAT new[SIZEMAX][SIZEMAX],PASTIX_INT lda)
{
  PASTIX_INT i,j;
  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      if (*t=='N')
	new[j][i]=a[i+j*lda];
      else
	new[i][j]=a[i+j*lda];
}

void end_matrice(char *t,PASTIX_INT n,PASTIX_INT m,PASTIX_FLOAT *a,PASTIX_FLOAT new[SIZEMAX][SIZEMAX],PASTIX_INT lda)
{
  PASTIX_INT i,j;
  for (i=0;i<n;i++)
    for (j=0;j<m;j++)
      if (*t=='N')
	a[i+j*lda]=new[j][i];
      else
	a[i+j*lda]=new[i][j];
}

static Clock test_clk;

#define TEST_CLOCK_INIT {clockInit(&test_clk);clockStart(&test_clk);}
#define TEST_CLOCK_STOP {clockStop(&test_clk);}
#define TEST_CLOCK_GET  clockVal(&test_clk)
