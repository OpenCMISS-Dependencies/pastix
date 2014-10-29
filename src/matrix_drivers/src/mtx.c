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
#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <sys/types.h>

#include "mmio.h"
#include "iohb.h"

#ifdef FORCE_NOMPI
#include "pastix_nompi.h"
#else
#include <mpi.h>
#endif

#ifdef X_ARCHsun
#include <inttypes.h>
#endif

#include "pastix.h"
#ifdef TYPE_COMPLEX
#if (defined X_ARCHalpha_compaq_osf1)

#ifndef USE_CXX

#ifndef   _RWSTD_HEADER_REQUIRES_HPP
#include <complex>
#else  /* _RWSTD_HEADER_REQUIRES_HPP */
#include <complex.hpp>
#endif /* _RWSTD_HEADER_REQUIRES_HPP */

#define ABS_FLOAT(x) abs(x)
#endif /* USE_CXX */

#else /*X_ARCHalpha_compaq_osf1*/
#include <complex.h>
#ifdef    PREC_DOUBLE
#define ABS_FLOAT(x) cabs(x)
#else  /* PREC_DOUBLE */
#define ABS_FLOAT(x) cabsf(x)
#endif /* PREC_DOUBLE */

#endif /* X_ARCHalpha_compaq_osf1 */
#else /* TYPE_COMPLEX */


#ifdef PREC_DOUBLE
#define ABS_FLOAT(x) fabs(x)
#else /* PREC_DOUBLE */
#define ABS_FLOAT(x) fabsf(x)
#endif /* FORCEDOUBLE */

#endif /* TYPE_COMPLEX */

#include "mtx.h"


#define MOD_SI 0

#define memFree_null(x) {if (x ==NULL) {fprintf(stdout,"%s:%d freeing NULL\n",__FILE__,__LINE__);} free(x); x=NULL;}

#define ASSERT(expr,module){\
    if((expr) == 0){\
      fprintf(stderr,"error in assert (line=%d,file=%s)\n",__LINE__,__FILE__);\
      exit(EXIT_FAILURE);}}

#define EXIT(module,error) {\
    exit(EXIT_FAILURE);}

#define MALLOC_ERROR(x) {\
    fprintf(stderr,"\nERROR %s allocation (line=%d,file=%s)\n\n",(x),__LINE__,__FILE__); \
    EXIT(MOD_UNKNOWN,ALLOC_ERR);}


/* set to extract only symmetric part in CCC format */
#define SYMPART

#ifndef USE_NOFORTRAN
#if (defined X_ARCHpower_ibm_aix)
#define FORTRAN_CALL(nom) nom
#else
#define FORTRAN_CALL(nom) nom ## _
#endif
#else
#define FORTRAN_CALL(nom) 
#endif

void  FORTRAN_CALL(wreadmtc)
     (int * tmp1,int * tmp2,int *tmp3,const char *filename,int *len, 
      double * val, int * row, int * col, double *crhs, int * nrhs,
      char * RhsType, int * tmpNrow, int * tmpNcol, int * tmpNnzero,
      char * title,char * key, char * Type, int * ierr);

/* #ifdef __INTEL_COMPILER */
/* void FORTRAN_CALL(readheaderfdupros)(char const *filename, int * tmpint1, int * tmpint2); */
/* void FORTRAN_CALL(readfdupros)      (char const *filename, int * Nrow,    int * Nnzero,  */
/* 				     int64_t * col, int64_t * row, double * val, double * rhs); */
/* #endif */

/* couple row/column */
typedef struct couple_
{
  pastix_int_t i,j;
} couple;



void ijvReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type)
{
  char line[BUFSIZ];
  long temp1,temp2,temp3;

  Type[0] = 'R';
  Type[1] = 'U';
  Type[2] = 'A';
  Type[3] = '\0';

  /* ncol nrow nnzero */
  fgets(line,BUFSIZ,infile);

  sscanf(line, "%ld %ld %ld", &temp1,&temp2,&temp3);
  if (temp1!=temp2)
    {
      temp2=temp1;
      fgets(line,BUFSIZ,infile);
      sscanf(line, "%ld", &temp3);
    }

  *Nrow = (pastix_int_t)temp1;
  *Ncol = (pastix_int_t)temp2;
  *Nnzero = (pastix_int_t)temp3;
}

void ijvRead(char *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType, pastix_int_t flagsort)
{
  pastix_int_t iter=0;
  pastix_int_t iter2=0;
  pastix_int_t *tempcol;
  FILE *infile;
  pastix_int_t baseval;

#ifdef TYPE_COMPLEX
  fprintf(stderr,"\nWARNING: Non complex driver\n\n");
#endif

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(1*sizeof(char));
  (*RhsType)[0] = '\0';
  infile=fopen(filename,"r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      /* EXIT(MOD_SI,FILE_ERR); */
      exit(EXIT_FAILURE);
    }

  /* Lire Entete */
  ijvReadHeader(infile, Nrow, Ncol, Nnzero, *Type);
  printf("IJV: Nrow=%ld Ncol=%ld Nnzero=%ld\n",(long)*Nrow,(long)*Ncol,(long)*Nnzero);

  /* Allocation memoire */
  tempcol = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));

  if ((tempcol==NULL) || ((*row) == NULL) || ((*val) == NULL))
    {
      fprintf(stderr, "mtxRead : Not enough memory\n");
      /* EXIT(MOD_SI,OUTOFMEMORY_ERR); */
      exit(EXIT_FAILURE);
    }
  
  /* Remplissage */
  for (iter=0; iter<(*Nnzero); iter++)
    {
      long temp1,temp2;
      double tmpval;
      fscanf(infile, "%ld %ld %lg\n", &temp1, &temp2, &tmpval);
      (*val)[iter] = tmpval;
      tempcol[iter]=(pastix_int_t)temp1;
      (*row)[iter]=(pastix_int_t)temp2;
    }

  /* Trie sur col et row */
  if (flagsort)
    for (iter=0; iter<(*Nnzero); iter++)
      for (iter2=0; iter2<(*Nnzero)-iter-1; iter2++)
	if ( (tempcol[iter2+1] < tempcol[iter2]) )
	  {
	    /* swap iter2 et iter2+1 */
	    pastix_int_t swapc = tempcol[iter2+1];
	    pastix_int_t swapr = (*row)[iter2+1];
	    pastix_float_t swapv = (*val)[iter2+1];
	    
	    tempcol[iter2+1] = tempcol[iter2];
	    (*row)[iter2+1] = (*row)[iter2];
	    (*val)[iter2+1] = (*val)[iter2];
	    tempcol[iter2] = swapc;
	    (*row)[iter2] = swapr;
	    (*val)[iter2] = swapv;
	  }

  (*col) = (pastix_int_t *) malloc((*Nrow+1)*sizeof(pastix_int_t));
  baseval=1; /* Attention on base a 1 */
  iter2=baseval; 
  for (iter=0; iter<(*Nrow); iter++)
    {
      (*col)[iter] = iter2;
      while (tempcol[iter2-baseval] == iter+1)
	{
	  iter2++;
	}
    }
  (*col)[(*Nrow)] = iter2;

  memFree_null(tempcol);

  fclose(infile);
}

void rsaReadHeader(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type, char *RhsType)
{
  int tmp;
  int *col=NULL;
  int *row=NULL;
  char title[72+1];
  char key[8+1];
  int nrhs,len,ierr;
  double *val=NULL;
  double *crhs=NULL;
  int tmpNrow,tmpNcol,tmpNnzero;

  len=strlen(filename);
  tmp=0;
  FORTRAN_CALL(wreadmtc)
    (&tmp,&tmp,&tmp,filename,&len,val,row,col,crhs,&nrhs,
     RhsType,&tmpNrow,&tmpNcol,&tmpNnzero,title,key,Type,&ierr);
  if(ierr != 0) {
    fprintf(stderr, "cannot read matrix (job=0)\n");
  }
  *Nrow=(pastix_int_t)tmpNrow;
  *Ncol=(pastix_int_t)tmpNcol;
  *Nnzero=(pastix_int_t)tmpNnzero;
  Type[3]='\0';

  /*ASSERT(*Nrow==*Ncol,MOD_SI);*/
  if ((*Nrow==*Ncol) == 0)
    {
      fprintf(stderr,"ERROR : (*Nrow!=*Ncol)\n");
      exit(EXIT_FAILURE);
    }
}


#include<ctype.h>

/******************************************************************************/
/* void myupcase(char *S)                                                     */
/*                                                                            */
/* Met en majuscule une chaine caractere                                      */
/*                                                                            */
/* S : chaine de caractere                                                    */
/******************************************************************************/
void myupcase(char *S)
{
  pastix_int_t iter=0;

  while (S[iter] != '\0')
    {
      S[iter] = (char)toupper(S[iter]);
      iter++;
    }
}

void rsaRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType)
{
  int i,tmp;
  char title[72+1];
  char key[8+1];
  int nrhs,len,ierr;
  double *crhs=NULL;
  int tmpNrow,tmpNcol,tmpNnzero;
  int *tmpcol,*tmprow;
#ifdef TYPE_COMPLEX
  double * tmpval;
#endif /* TYPE_COMPLEX */
  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));

  rsaReadHeader(filename, Nrow, Ncol, Nnzero, *Type, *RhsType);
  printf("RSA: Nrow=%ld Ncol=%ld Nnzero=%ld\n",(long)*Nrow,(long)*Ncol,(long)*Nnzero);
  tmpNrow  =(int)*Nrow;
  tmpNcol  =(int)*Ncol;
  tmpNnzero=(int)*Nnzero;

  len = strlen(filename);

  *col   =(pastix_int_t*)  malloc(((*Nrow)+1)*sizeof(pastix_int_t));  ASSERT(*col  !=NULL,MOD_SI);
  tmpcol =(int*)           malloc((tmpNrow+1)*sizeof(int));           ASSERT(tmpcol!=NULL,MOD_SI);
  *row   =(pastix_int_t*)  malloc((*Nnzero)  *sizeof(pastix_int_t));  ASSERT(*row  !=NULL,MOD_SI);
  tmprow =(int*)           malloc(tmpNnzero  *sizeof(int));           ASSERT(tmprow!=NULL,MOD_SI);
  *val   =(pastix_float_t*)malloc((*Nnzero)  *sizeof(pastix_float_t));ASSERT(*val  !=NULL,MOD_SI);

  tmp  = 2;
  nrhs = 0;

#ifdef TYPE_COMPLEX
  tmpval=(double*)malloc((*Nnzero)*sizeof(double));
  fprintf(stderr,"\nWARNING: Non complex driver\n\n");

  FORTRAN_CALL(wreadmtc)
    (&tmpNrow,&tmpNnzero,&tmp,filename,&len,tmpval,tmprow,tmpcol,crhs,
     &nrhs,*RhsType,&tmpNrow,&tmpNcol,&tmpNnzero,title,key,*Type,&ierr);
  {
    int ii;
    for (ii = 0; ii < *Nnzero; ii++) 
      {
	(*val)[ii] = tmpval[ii] + I*0;
      }

  }

#else /* TYPE_COMPLEX */
  FORTRAN_CALL(wreadmtc)
    (&tmpNrow,&tmpNnzero,&tmp,filename,&len,*val,tmprow,tmpcol,crhs,
     &nrhs,*RhsType,&tmpNrow,&tmpNcol,&tmpNnzero,title,key,*Type,&ierr);
#endif /*TYPE_COMPLEX */

  (*RhsType)[0]='\0';
  myupcase(*Type);
  if(ierr != 0) {
    fprintf(stderr, "cannot read matrix (job=2)\n");
  }
  for (i=0; i<tmpNrow+1; i++) (*col)[i]=(pastix_int_t)(tmpcol[i]);
  for (i=0; i<tmpNnzero; i++) (*row)[i]=(pastix_int_t)(tmprow[i]);
  memFree_null(tmpcol);
  memFree_null(tmprow);
  *Nrow  =(pastix_int_t)tmpNrow;
  *Ncol  =(pastix_int_t)tmpNcol;
  *Nnzero=(pastix_int_t)tmpNnzero;
}


void HBRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, 
	    pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType)
{
  int i;
  int nrhs;
  int tmpNrow,tmpNcol,tmpNnzero;
  int *tmpcol,*tmprow;
  int Nrow2, Ncol2, Nnzero2;
  double * tmpval;
  int ierr;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));

  readHB_info(filename, &Nrow2, &Ncol2, &Nnzero2, Type, &nrhs);

  *Nrow = Nrow2;
  *Ncol = Ncol2;
  *Nnzero = Nnzero2;

/*   fprintf(stderr,"Matrix in file %s is %ld x %ld, with %ld nonzeros with type %s;\n", */
/* 	  filename, (long)*Nrow, (long)*Ncol, (long)*Nnzero, *Type); */
/*   fprintf(stderr,"%d right-hand-side(s) available.\n",nrhs); */

/*   printf("RSA: Nrow=%ld Ncol=%ld Nnzero=%ld\n",(long)*Nrow,(long)*Ncol,(long)*Nnzero); */
#ifdef TYPE_COMPLEX
  fprintf(stderr,"Warning: HBRead is a real matrix driver\n");
  exit(EXIT_FAILURE);
#endif

  tmpNrow=(int)*Nrow;
  tmpNcol=(int)*Ncol;
  tmpNnzero=(int)*Nnzero;


  *col=(pastix_int_t*)malloc((*Nrow+1)*sizeof(pastix_int_t));
  ASSERT(*col!=NULL,MOD_SI);
  tmpcol=(int*)malloc((tmpNrow+1)*sizeof(int));
  ASSERT(tmpcol!=NULL,MOD_SI);
  *row=(pastix_int_t*)malloc(*Nnzero*sizeof(pastix_int_t));
  ASSERT(*row!=NULL,MOD_SI);
  tmprow=(int*)malloc(tmpNnzero*sizeof(int));
  ASSERT(tmprow!=NULL,MOD_SI);
  *val=(pastix_float_t*)malloc(*Nnzero*sizeof(pastix_float_t));
  ASSERT(*val!=NULL,MOD_SI);

  nrhs=0;
#if (defined PREC_DOUBLE && !defined TYPE_COMPLEX)
  tmpval = *val;
#else
  tmpval = (double*)malloc(*Nnzero*sizeof(double));
#endif

  ierr = readHB_mat_double(filename, tmpcol, tmprow, tmpval);
  if(ierr == 0) {
    fprintf(stderr, "cannot read matrix (job=2)\n");
  }

#if (!defined PREC_DOUBLE || defined TYPE_COMPLEX)
  for (i = 0; i < *Nnzero; i++)
    (*val)[i] = (pastix_float_t)tmpval[i];
#endif
  (*RhsType)[0]='\0';
  myupcase(*Type);
  for (i=0;i<tmpNrow+1;i++) (*col)[i]=(pastix_int_t)(tmpcol[i]);
  for (i=0;i<tmpNnzero;i++) (*row)[i]=(pastix_int_t)(tmprow[i]);
  memFree_null(tmpcol);
  memFree_null(tmprow);
  *Nrow=(pastix_int_t)tmpNrow;
  *Ncol=(pastix_int_t)tmpNcol;
  *Nnzero=(pastix_int_t)tmpNnzero;
}

void diag_dominance(pastix_int_t baseval, pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a)
{
  pastix_int_t i, j;
  
  ASSERT(baseval == ia[0],MOD_SI);

  for(i=0;i<n;i++)
    for(j=ia[i];j<ia[i+1];j++)
      {
	if(ja[j-baseval] == i+baseval)
	  a[j-baseval] = 2.0*(ia[i+1]-ia[i])-1.0;
	else
	  a[j-baseval] = -0.8;
      }
}

void diag_unite(pastix_int_t baseval, pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja, pastix_float_t *a)
{
  pastix_int_t i, j;
  
  ASSERT(baseval == ia[0],MOD_SI);

  for(i=0;i<n;i++)
    for(j=ia[i];j<ia[i+1];j++)
      {
	if(ja[j-baseval] == i+baseval)
	  a[j-baseval] = 1.0;
	else
	  a[j-baseval] = 0.0;
      }
}

void no_diag(pastix_int_t baseval, pastix_int_t n, pastix_int_t *ia, pastix_int_t *ja)
{
  /***************************************/
  /* Supress diagonal term               */
  /* On entry:                           */
  /*    n, ia, ja : CSR matrix           */
  /* On return:                          */
  /*    ia, ja : CSR matrix without diag */
  /***************************************/
  pastix_int_t i, j;
  pastix_int_t indj;
  pastix_int_t *old_ia;

  old_ia = (pastix_int_t*)malloc(sizeof(pastix_int_t)*(n+1));
  memcpy(old_ia, ia, sizeof(pastix_int_t)*(n+1));
  
  ASSERT(ia[0]==baseval,MOD_SI);

  indj = 0;

  for(i=0;i<n;i++)
    {
      ia[i] = indj+baseval;
      for(j=old_ia[i];j<old_ia[i+1];j++)
	if(ja[j-baseval] != i+baseval)
	  {
	    ja[indj] = ja[j-baseval];
	    indj++;
	  }
    }
  ia[n] = indj+baseval;

  memFree_null(old_ia);
}

void dimsym(pastix_int_t n, pastix_int_t **ia, pastix_int_t **ja)
{
  pastix_int_t *add;
  pastix_int_t iter,iter2,iter3;
  pastix_int_t cumul=0;
  pastix_int_t *newia,*newja;

  add = (pastix_int_t *) malloc((n+1)*sizeof(pastix_int_t));
  
  for (iter=0; iter<(n+1); iter++)
    {
      add[iter] = 0;
    }

  for (iter=0; iter<n; iter++)
    {
      for (iter2=(*ia)[iter]; iter2<(*ia)[iter+1]; iter2++)
	{
	  if (((*ja)[iter2-1]-1) != iter)
	    {
	      add[((*ja)[iter2-1])]++;
	    }
	}
    }

  newia = (pastix_int_t *) malloc((n+1)*sizeof(pastix_int_t));

  for (iter=0; iter<(n+1); iter++)
    {
      cumul += add[iter];
      newia[iter] = (*ia)[iter] + cumul;
      add[iter] = 0;
    }

  newja = (pastix_int_t *) malloc(((newia)[n]-1)*sizeof(pastix_int_t));
  
  for (iter=0; iter<n; iter++)
    {
      iter3 = newia[iter]+add[iter]-1;
      for (iter2=(*ia)[iter]-1; iter2<(*ia)[iter+1]-1; iter2++)
	{
	  newja[iter3]=(*ja)[iter2];
	  if ( ( (*ja)[iter2]-1) != iter) 
	    {
	      newja[newia[(*ja)[iter2]-1]+add[(*ja)[iter2]-1]-1] = iter+1;
	      add[(*ja)[iter2]-1]++;
	    }
	  iter3++;
	}
    }

  memFree_null(add);
  memFree_null(*ia);
  memFree_null(*ja);

  *ia = newia;
  *ja = newja;  
}

void cccReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type)
{
  long temp1,temp2;
  fscanf(infile, "%ld %ld\n", &temp1, &temp2);
  *Nrow=(pastix_int_t)temp1;
  *Nnzero=(pastix_int_t)temp2;
  *Ncol = *Nrow;
  Type[0] = 'C';
  Type[1] = 'U';
  Type[2] = 'A';
  Type[3] = '\0';
#ifdef SYMPART
  Type[1] = 'S';
  *Nnzero=(*Nnzero-*Ncol)/2+*Ncol;
#endif
}

#ifndef TYPE_COMPLEX
#define I 0
#endif

void cccRead(char const * filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType)
{
  FILE *infile,*infile1,*infile2;
  pastix_int_t iter,size,i=0,ii=0;
  double temp1,temp2;

  temp1=0;
  temp2=0;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));
  (*RhsType)[0] = '\0';

  infile = fopen(filename, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "hfile");
      EXIT(MOD_SI,FILE_ERR);
    }
  cccReadHeader(infile, Nrow, Ncol, Nnzero, *Type);
  fclose(infile);

  printf("Nrow %ld Ncol %ld Nnzero %ld\n", (long)*Nrow, (long)*Ncol, (long)*Nnzero);
  
  (*col) = (pastix_int_t *) malloc((*Ncol+1)*sizeof(pastix_int_t));
  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
 
  if (((*col) == NULL) || ((*row) == NULL) || ((*val) == NULL))
    fprintf(stderr, "cccRead : Not enough memory for \n");


  infile = fopen("ifile", "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "ifile");
      EXIT(MOD_SI,FILE_ERR);
    }
  for (iter=0; iter<(*Ncol+1); iter++)
    {
      long temp;
      fscanf(infile, "%ld", &temp);
      (*col)[iter]=(pastix_int_t)temp;
    }
  fclose(infile);

#ifdef SYMPART
  size=2*(*Nnzero-*Ncol)+*Ncol;
#else
  size=*Nnzero;
#endif
  
  infile1 = fopen("jfile", "r");
  if (infile1==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "jfile");
      EXIT(MOD_SI,FILE_ERR);
    }
  infile2 = fopen("afile", "r");
  if (infile2==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "afile");
      EXIT(MOD_SI,FILE_ERR);
    }
  for (iter=0; iter<size; iter++)
    {
      long x;
      fscanf(infile1, "%ld", &x);
#ifdef SYMPART
      if (iter+1>=(*col)[ii+1])
	{
	  (*col)[ii+1]=i+1;
	  ii++;
	}
      if ((pastix_int_t)x>=ii+1) 
	{
	  (*row)[i] = (pastix_int_t)x;
	}
#else
      (*row)[iter] = (pastix_int_t)x;
#endif

      fscanf(infile2, "%lf %lf", &temp1, &temp2);
#ifdef SYMPART
      if ((pastix_int_t)x>=ii+1)
	{
#if (defined X_ARCHalpha_compaq_osf1)
#ifdef CPLX
	  (*val)[i] = PASTIX_FLOAT(temp1,temp2);
#else
	  (*val)[i] = temp1;
#endif
#else
	  (*val)[i] = temp1;
	  (*val)[i] += ((pastix_float_t)temp2)*I;
#endif
	  i++;
	}
#else
      (*val)[iter] = temp1;
      (*val)[iter] += ((pastix_float_t)temp2)*I;
#endif
    }
  fclose(infile1);
  fclose(infile2);

#ifdef SYMPART
  (*col)[*Ncol]=i+1;
  ASSERT(i==*Nnzero,MOD_SI);
#endif
}

void olafReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type)
{
  long temp1;
  fscanf(infile, "%ld\n", &temp1);
  *Nrow = (pastix_int_t)temp1;
  fscanf(infile, "%ld\n", &temp1);
  *Nnzero = (pastix_int_t)temp1;
  *Ncol = *Nrow;
  Type[0] = 'R';
  Type[1] = 'S';
  Type[2] = 'A';
  Type[3] = '\0';
}

void olafRead(char const * filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType, pastix_float_t **rhs)
{
  FILE *infile;
  pastix_int_t iter,size;
  long temp1;
  double temp2;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));
  (*RhsType)[0] = 'A';
  (*RhsType)[1] = 'A';
  (*RhsType)[2] = 'A';

  infile = fopen(filename, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "olafcsr");
      EXIT(MOD_SI,FILE_ERR);
    }
  olafReadHeader(infile, Nrow, Ncol, Nnzero, *Type);

  printf("Nrow %ld Ncol %ld Nnzero %ld\n", (long)*Nrow, (long)*Ncol, (long)*Nnzero);
  
  (*col) = (pastix_int_t *) malloc((*Ncol+1)*sizeof(pastix_int_t));
  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  (*rhs) = (pastix_float_t *) malloc((*Ncol)*sizeof(pastix_float_t));
 
  if (((*col) == NULL) || ((*row) == NULL) || ((*val) == NULL) || ((*rhs) == NULL))
    fprintf(stderr, "olafRead : Not enough memory for \n");

  for (iter=0; iter<(*Ncol+1); iter++)
    {
      fscanf(infile, "%ld", &temp1);
      (*col)[iter] = (pastix_int_t)temp1;
    }

  size=*Nnzero;
  
  for (iter=0; iter<size; iter++)
    {
      fscanf(infile, "%ld", &temp1);
      (*row)[iter] = (pastix_int_t)temp1;
    }

  for (iter=0; iter<size; iter++)
    {
      fscanf(infile, "%lf", &temp2);
      (*val)[iter] = temp2;
    }

  fclose(infile);

  infile = fopen("olafrhs", "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", "olafrhs");
      EXIT(MOD_SI,FILE_ERR);
    }

  for (iter=0; iter<(*Ncol); iter++)
    {
      fscanf(infile, "%lf", &temp2);
      (*rhs)[iter] = temp2;
    }

  fclose(infile);
}

/******************************************************************************/
/* void mysubstr(char *s, char *S, pastix_int_t pos, pastix_int_t len)                          */
/*                                                                            */
/* remplit s avec le contenu de S a partir de la position pos avec len element*/
/*                                                                            */
/******************************************************************************/
void mysubstr(char *s, const char *S, const pastix_int_t pos, const pastix_int_t len)
{
  pastix_int_t iter;
  for (iter=0; iter<len;iter++)
    {
      s[iter] = S[pos+iter];
    }
  s[len] = '\0';
}

/******************************************************************************/
/* void mysubstr2(char *fmt, char a, char b, pastix_int_t *val)                        */
/*                                                                            */
/* val recoit le nombre contenu dans fmt compris entre les caracteres a et b  */
/******************************************************************************/
void mysubstr2(const char *fmt, const char a, const char b, pastix_int_t *val)
{
  pastix_int_t len = strchr(fmt,b) - strchr(fmt,a) -1;
  char *tmp = (char *) malloc(len+1);
  if (tmp == NULL)
    {
      fprintf(stderr, "mysubstr2 : Not enough memory for tmp\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  mysubstr(tmp, fmt, strchr(fmt, a) - fmt +1, len);
  *val = atoi(tmp);
  memFree_null(tmp);
}

/******************************************************************************/
/* void hbParseRfmt(char *fmt, pastix_int_t *perline, pastix_int_t *width, pastix_int_t *prec, pastix_int_t *flag)*/
/*                                                                            */
/* Parse un format fortran pour les reels flottants                           */
/*                                                                            */
/* fmt : format fortran                                                       */
/* perline : nombre de valeurs par ligne                                      */
/* width : nombre de chiffre                                                  */
/* prec : precision                                                           */
/* flag : nombre de caractere inutile en debut de ligne                       */
/******************************************************************************/
void chbParseRfmt(char *fmt, pastix_int_t *perline, pastix_int_t *width, pastix_int_t *prec, char *flag)
{
  myupcase(fmt);
  
  if (strchr(fmt,'E') != NULL)
    {
      *flag = 'E';
    }
  else if (strchr(fmt,'D') != NULL)
    { 
      *flag = 'D';
    }
  else if (strchr(fmt,'F') != NULL)
    { 
      *flag = 'F';
    }

  if (strchr(fmt,'P') != NULL)
    {
      if (strchr(fmt,'P')[1] == ',')
	{
	  if (strchr(fmt,'P')[2] == *flag)
	    {
	      /* (3(1P,E24.16)) */
	      mysubstr2(fmt, '(', *flag, perline);
	    }
	  else
	    {
	      /* (1P,3E25.16E3) */
	      mysubstr2(fmt, ',', *flag, perline);
	    }
	}

      else
	mysubstr2(fmt, 'P', *flag, perline);
    }
  else
    mysubstr2(fmt, '(', *flag, perline);
  
  if ( strchr(fmt,'.') )
    {
      mysubstr2(fmt, '.', ')', prec);
      mysubstr2(fmt, *flag, '.', width);
    }
  else
    {
      mysubstr2(fmt, *flag, ')', width);
    }
}

/******************************************************************************/
/* void hbParseIfmt(char *fmt, pastix_int_t *perline, pastix_int_t *width, pastix_int_t *flag)           */
/*                                                                            */
/* Parse un format fortran pour les entiers                                   */
/*                                                                            */
/* fmt : format fortran                                                       */
/* perline : nombre de valeurs par ligne                                      */
/* width : nombre de chiffre                                                  */
/* flag : nombre de caractere inutile en debut de ligne                       */
/******************************************************************************/
void chbParseIfmt(char *fmt, pastix_int_t *perline, pastix_int_t *width, pastix_int_t *flag)
{
  myupcase(fmt);

  *flag=0;

  if ( strchr(fmt, ',') )
    {
      mysubstr2(fmt, ',', 'I', perline);
      mysubstr2(fmt, 'I', ')', width);
      if (strchr(fmt, 'X'))
	{
	  *flag=1;
	}
    }
  else
    {
      mysubstr2(fmt, '(', 'I', perline);
      mysubstr2(fmt, 'I', ')', width);
      if (strchr(fmt, 'X'))
	{
	  *flag=1;
	}
    }
}

void chbReadHeader(FILE *infile, char *Type, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t *Nrhs, char *Ptrfmt, char *Indfmt, char *Valfmt, char *Rhsfmt, pastix_int_t *Ptrcrd, pastix_int_t *Indcrd, pastix_int_t *Valcrd, pastix_int_t *Rhscrd,  char *RhsType)
{
  char line[BUFSIZ];
  long totcrd;
  long neltvl;
  long temp1,temp2,temp3,temp4;

  /* first line */
  /* 1-72 title 73-80 Key */
  fgets(line, BUFSIZ, infile);

  /* Seconde line */
  /* 1-14 totcrd 15-28 ptrcrd 29-42 indcrd 43-56 valcrd 57-70 rhscrd */
  fgets(line, BUFSIZ ,infile);
  sscanf(line, "%14ld%14ld%14ld%14ld%14ld", &totcrd, &temp1,&temp2,&temp3,&temp4);
  *Ptrcrd = temp1;
  *Indcrd = temp2;
  *Valcrd = temp3;
  *Rhscrd = temp4;

  /* Third line*/
  /* 1-3 mxtype 15-28 nrow 29-42 ncol 43-56 nnzero 57-70 neltvl */
  fgets(line, BUFSIZ, infile);
  sscanf(line, "%3c%ld%ld%ld%ld", Type,&temp1,&temp2,&temp3,&neltvl);
  *Nrow = temp1;
  *Ncol = temp2;
  *Nnzero = temp3;
  Type[3] = '\0';
  myupcase(Type);

  /* fourth line */
  /* 1-16 ptrfmt 17-32 indfmt 33-52 valfmt 53-72 rhsfmt */
  fgets(line, BUFSIZ, infile);
  sscanf(line, "%16c%16c%20c%20c", Ptrfmt, Indfmt, Valfmt, Rhsfmt);
  Ptrfmt[16] = '\0';
  Indfmt[16] = '\0';
  Valfmt[20] = '\0';
  Rhsfmt[20] = '\0';

  /* fifth line optional */
  /* 1  2 rhstyp 3  15-28 nrhs 29-42 nrhsix */
  if (*Rhscrd != 0)
    {
      fgets(line, BUFSIZ, infile);
      sscanf(line,"%3c%ld", RhsType, &temp1);
      *Nrhs = temp1;
      RhsType[3] = '\0';
    }
  else
    {
      RhsType[0] = '\0';
    }
}

void chbRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row,
	     pastix_float_t **val, char **Type, char **RhsType, pastix_float_t **rhs)
{
  FILE *infile;
  pastix_int_t Nrhs;
  char Ptrfmt[17];
  char Indfmt[17];
  char Valfmt[21];
  char Rhsfmt[21];
  pastix_int_t Ptrcrd, Indcrd, Valcrd, Rhscrd=0;
  pastix_int_t Valperline, Valwidth, Valprec;
  char Valflag;
  pastix_int_t Rhsperline, Rhswidth, Rhsprec;
  char Rhsflag;
  pastix_int_t Indperline, Indwidth, Indflag;
  pastix_int_t Ptrperline, Ptrwidth, Ptrflag;
  char *element;
  pastix_int_t count;
  pastix_int_t iter,item,colcur;
  char line[BUFSIZ];
  pastix_float_t * rhs2 = NULL;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));
  
  infile = fopen(filename, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }
  
  chbReadHeader(infile, *Type, Nrow, Ncol, Nnzero, &Nrhs, Ptrfmt, Indfmt, Valfmt, Rhsfmt, &Ptrcrd, &Indcrd, &Valcrd, &Rhscrd, *RhsType);
  printf("CHB: Nrow=%ld Ncol=%ld Nnzero=%ld\n",(long)*Nrow,(long)*Ncol,(long)*Nnzero);

  chbParseRfmt(Valfmt, &Valperline, &Valwidth, &Valprec, &Valflag);
  chbParseIfmt(Indfmt, &Indperline, &Indwidth, &Indflag);
  chbParseIfmt(Ptrfmt, &Ptrperline, &Ptrwidth, &Ptrflag);

  if (Rhscrd != 0)
    {
      chbParseRfmt(Rhsfmt, &Rhsperline, &Rhswidth, &Rhsprec, &Rhsflag);
    }

  element = (char *) malloc(Ptrwidth+1);
  if (element == NULL)
    {
      fprintf(stderr, "chbRead : Not enough memory for element\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  element[Ptrwidth] = '\0';

  (*col) = (pastix_int_t *) malloc((*Ncol+1)*sizeof(pastix_int_t));
  if ((*col) == NULL)
    {
      fprintf(stderr, "chbRead : Not enough memory for *col\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }

  count = 0;
  for (iter=0; iter<Ptrcrd; iter++)
    {
      fgets(line, BUFSIZ, infile);

      colcur = Ptrflag;

      for (item=0; item<Ptrperline; item++)
	{
	  if (count > (*Ncol))
	    break;

	  strncpy(element, line+colcur, Ptrwidth);
	  (*col)[count] = atoi(element);

	  /*
	    if ((iter==0) || (iter==1))
	    {
	    fprintf(stderr, "count %ld element %s col %ld\n",
	    count,element, (*col)[count]);
	    }
	  */
	  count++;
	  colcur += Ptrwidth;

	}
    }
  memFree_null(element);

  element = (char *) malloc(Indwidth+1);
  if (element == NULL)
    {
      fprintf(stderr, "chbRead : Not enough memory for element\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  element[Indwidth] = '\0';

  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  if ((*row) == NULL)
    {
      fprintf(stderr, "chbRead : Not enough memory for *row\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }

  count = 0;
  for (iter=0; iter<Indcrd; iter++)
    {
      fgets(line, BUFSIZ, infile);

      colcur = Indflag;
      for (item=0; item<Indperline; item++)
	{
	  if (count == (*Nnzero))
	    break;

	  strncpy(element, line+colcur, Indwidth);
	  (*row)[count] = atoi(element);

	  /*
	    if ((iter==0) || (iter==1))
	    {
	    fprintf(stderr, "count %ld element %s row %ld\n",
	    count, element, (*row)[count]);
	    }
	  */
	  count++;
	  colcur += Indwidth;
	}
    }
  memFree_null(element);


  element = (char *) malloc(Valwidth+1);
  if (element == NULL)
    {
      fprintf(stderr, "ChbRead : Not enough memory for element\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  element[Valwidth] = '\0';
  
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  if ((*val) == NULL)
    {
      fprintf(stderr, "chbRead : Not enough memory for *val\n");
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }

  count = 0;
  for (iter=0; iter<Valcrd; iter++)
    {
      fgets(line, BUFSIZ, infile);

      colcur = 0;
      for (item=0; item<Valperline; item++)
	{
#ifdef CPLX
	  if (MTX_ISCOM(*Type))
	    {
	      if (count == 2*(*Nnzero))
		break;
	    }
	  else
	    {
#endif /* CPLX */
	      if (count == (*Nnzero))
		break;
#ifdef CPLX
	    }
#endif /* CPLX */

	  strncpy(element, line+colcur, Valwidth);
	  if (Valflag == 'D')
	    if (strchr(element, 'D') != NULL)
	      *(strchr(element, 'D')) = 'E';

#ifdef CPLX
	  if (MTX_ISCOM(*Type))
	    {
	      if (count %2)
		{
#if (defined X_ARCHalpha_compaq_osf1)
		  (*val)[(count-1)/2] += pastix_float_t(0.0, atof(element));
#else
		  (*val)[(count-1)/2] += atof(element)*I;
#endif
		}
	      else
		{
		  (*val)[count/2] = atof(element);
		}
	    }
	  else
	    {
#endif /* CPLX */
	      (*val)[count] = atof(element);
#ifdef CPLX
	    }
#endif /* CPLX */
	  count++;
	  colcur += Valwidth;
	}
    }
  memFree_null(element);
  fprintf(stderr, "Bufsiz %d\n", BUFSIZ);

  if (Rhscrd != 0)
    {
      pastix_int_t nbline = Rhscrd;

      if (MTX_ISRHX(*RhsType))
	nbline = Rhscrd/2;

      element = (char *) malloc(Rhswidth+1);
      if (element == NULL)
	{
	  fprintf(stderr, "ChbRead : Not enough memory for element\n");
	  EXIT(MOD_SI,OUTOFMEMORY_ERR);
	}
      element[Rhswidth] = '\0';

      if (MTX_ISRHX(*RhsType))
	(*rhs) = (pastix_float_t *) malloc(2*(*Ncol)*sizeof(pastix_float_t));
      else
	(*rhs) = (pastix_float_t *) malloc((*Ncol)*sizeof(pastix_float_t));
      if ((*rhs) == NULL)
	{
	  fprintf(stderr, "chbRead : Not enough memory for *rhs\n");
          EXIT(MOD_SI,OUTOFMEMORY_ERR);
	}

      count = 0;
      for (iter=0; iter<nbline; iter++)
	{
	  fgets(line, BUFSIZ, infile);

	  colcur=0;
	  for (item=0; item<Rhsperline; item++)
	    {
#ifdef CPLX
	      if (MTX_ISCOM(*Type))
		{
		  if (count == 2*(*Ncol))
		    break;
		}
	      else
		{
#endif /* CPLX */
	      if (count == (*Ncol))
		break;
#ifdef CPLX
		}
#endif /* CPLX */

	      strncpy(element, line+colcur, Rhswidth);

#ifdef CPLX
	      if (MTX_ISCOM(*Type))
		{
		  if (count % 2)
		    {
#if (defined X_ARCHalpha_compaq_osf1)
		      (*rhs)[(count-1)/2] += PASTIX_FLOAT(0.0, atof(element));
#else
		      (*rhs)[(count-1)/2] += atof(element)*I;
#endif
		    }
		  else
		    {
		      (*rhs)[count/2] = atof(element);
		    }
		}
	      else
		{
#endif /* CPLX */
		  (*rhs)[count] = atof(element);
#ifdef CPLX
		}
#endif /* CPLX */
	      count++;
	      colcur += Rhswidth;
	    }
	}

      if (MTX_ISRHX(*RhsType))
	{
	  rhs2 = &(*rhs)[*Ncol];
	  if (rhs2 == NULL)
	    {
	      fprintf(stderr, "chbRead : Not enough memory for *rhs2\n");
	      EXIT(MOD_SI,OUTOFMEMORY_ERR);
	    }

	  count = 0;
	  for (iter=0; iter<nbline; iter++)
	    {
	      fgets(line, BUFSIZ, infile);

	      colcur = 0;
	      for (item=0; item<Rhsperline; item++)
		{
#ifdef CPLX
		  if (MTX_ISCOM(*Type))
		    {
		      if (count == 2*(*Ncol))
			break;
		    }
		  else
		    {
#endif /* CPLX */
		  if (count == (*Ncol))
		    break;
#ifdef CPLX
		    }
#endif /* CPLX */
		  strncpy(element, line+colcur, Rhswidth);
#ifdef CPLX
		  if (MTX_ISCOM(*Type))
		    {
		      if (count % 2)
			{
#if (defined X_ARCHalpha_compaq_osf1)
			  rhs2[(count-1)/2] += PASTIX_FLOAT(0.0, atof(element));
#else
			  rhs2[(count-1)/2] += atof(element)*I;
#endif
			}
		      else
			{
			  rhs2[count/2] = atof(element);
			}
		    }
		  else
		    {
#endif /* CPLX */
		  rhs2[count] = atof(element);
#ifdef CPLX
		    }
#endif /* CPLX */
		  count++;
		  colcur += Rhswidth;
		}
	    }
	}
      else
	{
	  rhs2 = NULL;
	}
      memFree_null(element);
    }
  else
    {
      (*rhs) = NULL;
      rhs2 = NULL;
    }
  
  fclose(infile);
}

void peerRead(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType, pastix_float_t **rhs)
{
  FILE *infile;
  pastix_int_t iterfile;
  char line[BUFSIZ];
  pastix_int_t rowlocal,rowglobal;
  pastix_int_t nzlocal,nzglobal;
  long filenamenumber;
  char **filenametab;
  const pastix_int_t nbreltperline=4; /* nbr of elt per line */
  long tempint1,tempint2,tempint3,tempint4;
  long double tempfloat1, tempfloat2, tempfloat3, tempfloat4;

  *Nnzero=0;
  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));

  (*Type)[0] = 'R';
  (*Type)[1] = 'U';
  (*Type)[2] = 'A';
  (*Type)[3] = '\0';
  (*RhsType)[0] = 'A';
  (*RhsType)[2] = 'A';
  (*RhsType)[3] = '\0';

  /* Read rsaname */
  infile = fopen(filename, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      exit(EXIT_FAILURE);
    }
  fgets(line, BUFSIZ, infile);
  sscanf(line, "%ld", &filenamenumber); /* Read number of filename */
  filenametab = (char **) malloc(filenamenumber*sizeof(char *));
  for (iterfile=0; iterfile<filenamenumber; iterfile++)
    {
      filenametab[iterfile] = (char *) malloc(64*sizeof(char));
      fgets(line, BUFSIZ, infile);
      sscanf(line, "%s", filenametab[iterfile]);
    }
  fclose(infile);
  
  /* Calcul nnz global */
  for (iterfile=0; iterfile<filenamenumber; iterfile++)
    {
      infile = fopen(filenametab[iterfile], "r");
      if (infile==NULL)
	{
	  fprintf(stderr,"cannot load %s\n", filenametab[iterfile]);
	  exit(EXIT_FAILURE);
	}
      fgets(line, BUFSIZ, infile);
      sscanf(line, "%ld%ld", &tempint1,&tempint2);
      *Nrow = tempint1;
      rowlocal = tempint2;
      fprintf(stderr, "Nrow %ld rowlocal %ld\n", (long) *Nrow, (long) rowlocal);
      fgets(line, BUFSIZ, infile);
      sscanf(line, "%ld", &tempint1);
      nzlocal = tempint1;
      fprintf(stderr, "nzlocal %ld\n", (long) nzlocal);
      fclose(infile);
    
      *Nnzero += nzlocal;
    }
  *Ncol = *Nrow;
  fprintf(stderr, "Nnzero global %ld\n", (long int) *Nnzero);

  /* memory alloc */
  (*col) = (pastix_int_t *) malloc((*Nrow+1)*sizeof(pastix_int_t));
  if ((*col) == NULL)
    fprintf(stderr, "peerRead : Not enough memory for *col\n");
  (*row) = (pastix_int_t *) malloc(*Nnzero*sizeof(pastix_int_t));
  if ((*row) == NULL)
    fprintf(stderr, "peerRead : Not enough memory for *row\n");
  (*val) = (pastix_float_t *) malloc(*Nnzero*sizeof(pastix_float_t));
  if ((*val) == NULL)
    fprintf(stderr, "peerRead : Not enough memory for *val\n");
  (*rhs) = (pastix_float_t *) malloc(*Nrow*sizeof(pastix_float_t));
  if ((*rhs) == NULL)
    fprintf(stderr, "peerRead : Not enough memory for *rhs\n");

  rowglobal=0;
  nzglobal=0;
  for (iterfile=0; iterfile<filenamenumber; iterfile++)
    {
      pastix_int_t iterelt;
      
      infile = fopen(filenametab[iterfile], "r");
      if (infile==NULL)
	{
	  fprintf(stderr,"cannot load %s\n", filenametab[iterfile]);
	  exit(EXIT_FAILURE);
	}
      fgets(line,BUFSIZ,infile);
      sscanf(line, "%ld%ld", &tempint1, &tempint2);
      *Nrow = tempint1;
      rowlocal = tempint2;
      fgets(line,BUFSIZ,infile);
      sscanf(line, "%ld", &tempint1);
      nzlocal = tempint1;

      /* read col */
      for (iterelt=0; iterelt<rowlocal+1+1-nbreltperline;iterelt++)
	{
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld %ld", &tempint1, &tempint2, &tempint3, &tempint4);
	  (*col)[iterelt+rowglobal]   = (pastix_int_t)tempint1+nzglobal;
	  (*col)[iterelt+rowglobal+1] = (pastix_int_t)tempint2+nzglobal;
	  (*col)[iterelt+rowglobal+2] = (pastix_int_t)tempint3+nzglobal;
	  (*col)[iterelt+rowglobal+3] = (pastix_int_t)tempint4+nzglobal;
	  iterelt+=3;
	}
      
      switch (rowlocal-iterelt+1)
	{
	case 1:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld",&tempint1);
	  (*col)[iterelt+rowglobal] += (pastix_int_t)tempint1+nzglobal;
	  iterelt++;
	  break;
	case 2:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld", &tempint1, &tempint2);
	  (*col)[iterelt+rowglobal]   = (pastix_int_t)tempint1+nzglobal;
	  (*col)[iterelt+rowglobal+1] = (pastix_int_t)tempint2+nzglobal;
	  iterelt+=2;
	  break;
	case 3:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld", &tempint1, &tempint2, &tempint3);
	  (*col)[iterelt+rowglobal]   = (pastix_int_t)tempint1+nzglobal;
	  (*col)[iterelt+rowglobal+1] = (pastix_int_t)tempint2+nzglobal;
	  (*col)[iterelt+rowglobal+2] = (pastix_int_t)tempint3+nzglobal;
	  iterelt+=3;
	  break;
	}


      /* read row */
      for (iterelt=0; iterelt<nzlocal+1-nbreltperline; iterelt++)
	{
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld %ld", &tempint1,&tempint2,&tempint3,&tempint4);
	  (*row)[iterelt+nzglobal]   = (pastix_int_t)tempint1;
	  (*row)[iterelt+nzglobal+1] = (pastix_int_t)tempint2;
	  (*row)[iterelt+nzglobal+2] = (pastix_int_t)tempint3;
	  (*row)[iterelt+nzglobal+3] = (pastix_int_t)tempint4;
	  iterelt+=3;
	}
      switch (nzlocal-iterelt)
	{
	case 1:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld",&tempint1);
	  (*row)[iterelt+nzglobal] = (pastix_int_t)tempint1;
	  iterelt++;
	  break;
	case 2:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld", &tempint1, &tempint2);
	  (*row)[iterelt+nzglobal]   = (pastix_int_t)tempint1;
	  (*row)[iterelt+nzglobal+1] = (pastix_int_t)tempint2;
	  iterelt+=2;
	  break;
	case 3:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld", &tempint1, &tempint2, &tempint3);
	  (*row)[iterelt+nzglobal]   = (pastix_int_t)tempint1;
	  (*row)[iterelt+nzglobal+1] = (pastix_int_t)tempint2;
	  (*row)[iterelt+nzglobal+2] = (pastix_int_t)tempint3;
	  iterelt+=3;
	  break;
	}

      /* read val */
      for (iterelt=0; iterelt<nzlocal+1-nbreltperline; iterelt++)
	{
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3,&tempfloat4);
	  (*val)[iterelt+nzglobal]   = (pastix_float_t)tempfloat1;
	  (*val)[iterelt+nzglobal+1] = (pastix_float_t)tempfloat2;
	  (*val)[iterelt+nzglobal+2] = (pastix_float_t)tempfloat3;
	  (*val)[iterelt+nzglobal+3] = (pastix_float_t)tempfloat4;
	  iterelt+=3;
	}
      switch (nzlocal-iterelt)
	{
	case 1:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf",&tempfloat1);
	  (*val)[iterelt+nzglobal] = (pastix_float_t)tempfloat1;
	  iterelt++;
	  break;
	case 2:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf",&tempfloat1,&tempfloat2);
	  (*val)[iterelt+nzglobal]   = (pastix_float_t)tempfloat1;
	  (*val)[iterelt+nzglobal+1] = (pastix_float_t)tempfloat2;
	  iterelt+=2;
	  break;
	case 3:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3);
	  (*val)[iterelt+nzglobal]   = (pastix_float_t)tempfloat1;
	  (*val)[iterelt+nzglobal+1] = (pastix_float_t)tempfloat2;
	  (*val)[iterelt+nzglobal+2] = (pastix_float_t)tempfloat3;
	  iterelt+=3;
	  break;
	}
      nzglobal += nzlocal;

      /* read rhs */
      for (iterelt=0; iterelt<rowlocal+1-nbreltperline; iterelt++)
	{
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3,&tempfloat4);
	  (*rhs)[iterelt+rowglobal]   = (pastix_float_t)tempfloat1;
	  (*rhs)[iterelt+rowglobal+1] = (pastix_float_t)tempfloat2;
	  (*rhs)[iterelt+rowglobal+2] = (pastix_float_t)tempfloat3;
	  (*rhs)[iterelt+rowglobal+3] = (pastix_float_t)tempfloat4;
	  iterelt+=3;
	}

      switch (rowlocal-iterelt)
	{
	case 1:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf",&tempfloat1);
	  (*rhs)[iterelt+rowglobal] = (pastix_float_t)tempfloat1;
	  iterelt++;
	  break;
	case 2:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf",&tempfloat1,&tempfloat2);
	  (*rhs)[iterelt+rowglobal]   = (pastix_float_t)tempfloat1;
	  (*rhs)[iterelt+rowglobal+1] = (pastix_float_t)tempfloat2;
	  iterelt++;
	  break;
	case 3:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3);
	  (*rhs)[iterelt+rowglobal]   = (pastix_float_t)tempfloat1;
	  (*rhs)[iterelt+rowglobal+1] = (pastix_float_t)tempfloat2;
	  (*rhs)[iterelt+rowglobal+2] = (pastix_float_t)tempfloat3;
	  iterelt++;
	  break;
	}
      rowglobal += rowlocal;

      fclose(infile);
    }
}

void peerRead2(char const *filename, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, pastix_int_t **col, pastix_int_t **row, pastix_float_t **val, char **Type, char **RhsType, pastix_float_t **rhs)
{
  FILE *infile;
  pastix_int_t iterfile;
  char line[BUFSIZ];
  pastix_int_t rowlocal,rowglobal;
  pastix_int_t nzlocal,nzglobal;
  long filenamenumber;
  char **filenametab;
  pastix_int_t nbreltperline=6; /* nbr of elt per line */
  long tempint1,tempint2,tempint3,tempint4,tempint5,tempint6;
  long double tempfloat1, tempfloat2, tempfloat3;

  *Nnzero=0;
  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(4*sizeof(char));

  (*Type)[0] = 'R';
  (*Type)[1] = 'U';
  (*Type)[2] = 'A';
  (*Type)[3] = '\0';
  (*RhsType)[0] = 'A';
  (*RhsType)[2] = 'A';
  (*RhsType)[3] = '\0';

  /* Read rsaname */
  infile = fopen(filename, "r");
  if (infile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }
  fgets(line, BUFSIZ, infile);
  sscanf(line, "%ld", &filenamenumber); /* Read number of filename */
  filenametab = (char **) malloc(filenamenumber*sizeof(char *));
  for (iterfile=0; iterfile<filenamenumber; iterfile++)
    {
      filenametab[iterfile] = (char *) malloc(64*sizeof(char));
      fgets(line, BUFSIZ, infile);
      sscanf(line, "%s", filenametab[iterfile]);
    }
  fclose(infile);
  
  /* Calcul nnz global */
  for (iterfile=0; iterfile<filenamenumber; iterfile++)
    {
      infile = fopen(filenametab[iterfile], "r");
      if (infile==NULL)
	{
	  fprintf(stderr,"cannot load %s\n", filenametab[iterfile]);
	  EXIT(MOD_SI,FILE_ERR);
	}
      fgets(line, BUFSIZ, infile);
      sscanf(line, "%ld %ld %ld",&tempint1,&tempint2,&tempint3);
      *Nrow = (pastix_int_t)tempint1;
      rowlocal = (pastix_int_t)tempint2;
      nzlocal = (pastix_int_t)tempint3;
      printf("Nrow %ld rowlocal %ld nzlocal %ld\n", (long) *Nrow, (long) rowlocal, (long) nzlocal);
      fclose(infile);
      *Nnzero += nzlocal;
    }
  *Ncol = *Nrow;
  printf("Nnzero global %ld\n", (long)*Nnzero);

  /* memory alloc */
  (*col) = (pastix_int_t *) malloc((*Nrow+1)*sizeof(pastix_int_t));
  if ((*col) == NULL)
    fprintf(stderr, "peerRead : Not enough memory for *col\n");
  (*row) = (pastix_int_t *) malloc(*Nnzero*sizeof(pastix_int_t));
  if ((*row) == NULL)
    fprintf(stderr, "peerRead : Not enough memory for *row\n");
  (*val) = (pastix_float_t *) malloc(*Nnzero*sizeof(pastix_float_t));
  if ((*val) == NULL)
    fprintf(stderr, "peerRead : Not enough memory for *val\n");
  (*rhs) = (pastix_float_t *) malloc(*Nrow*sizeof(pastix_float_t));
  if ((*rhs) == NULL)
    fprintf(stderr, "peerRead : Not enough memory for *rhs\n");

  rowglobal=0;
  nzglobal=0;
  for (iterfile=0; iterfile<filenamenumber; iterfile++)
    {
      pastix_int_t iterelt;
      
      infile = fopen(filenametab[iterfile], "r");
      if (infile==NULL)
	{
	  fprintf(stderr,"cannot load %s\n", filenametab[iterfile]);
	  EXIT(MOD_SI,FILE_ERR);
	}
      fgets(line,BUFSIZ,infile);
      sscanf(line, "%ld %ld %ld", &tempint1, &tempint2, &tempint3);
      *Nrow = (pastix_int_t)tempint1;
      rowlocal = (pastix_int_t)tempint2;
      nzlocal = (pastix_int_t)tempint3;

      /* read col */
      for (iterelt=0; iterelt<rowlocal+1+1-nbreltperline;iterelt++)
	{
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld %ld %ld %ld", &tempint1, &tempint2, &tempint3, &tempint4, &tempint5, &tempint6);
	  (*col)[iterelt+rowglobal] = tempint1+nzglobal;
	  (*col)[iterelt+rowglobal+1] = tempint2+nzglobal;
	  (*col)[iterelt+rowglobal+2] = tempint3+nzglobal;
	  (*col)[iterelt+rowglobal+3] = tempint4+nzglobal;
	  (*col)[iterelt+rowglobal+4] = tempint5+nzglobal;
	  (*col)[iterelt+rowglobal+5] = tempint6+nzglobal;
	  iterelt+=5;
	}
      switch (rowlocal-iterelt+1)
	{
	case 1:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld",&tempint1);
	  (*col)[iterelt+rowglobal] += tempint1+nzglobal;
	  iterelt++;
	  break;
	case 2:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld", &tempint1, &tempint2);
	  (*col)[iterelt+rowglobal] = tempint1+nzglobal;
	  (*col)[iterelt+rowglobal+1] = tempint2+nzglobal;
	  iterelt+=2;
	  break;
	case 3:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld", &tempint1, &tempint2, &tempint3);
	  (*col)[iterelt+rowglobal] = tempint1+nzglobal;
	  (*col)[iterelt+rowglobal+1] = tempint2+nzglobal;
	  (*col)[iterelt+rowglobal+2] = tempint3+nzglobal;
	  iterelt+=3;
	  break;
	case 4:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld %ld", &tempint1, &tempint2, &tempint3, &tempint4);
	  (*col)[iterelt+rowglobal] = tempint1+nzglobal;
	  (*col)[iterelt+rowglobal+1] = tempint2+nzglobal;
	  (*col)[iterelt+rowglobal+2] = tempint3+nzglobal;
	  (*col)[iterelt+rowglobal+3] = tempint4+nzglobal;
	  iterelt+=4;
	  break;
	case 5:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld %ld %ld", &tempint1, &tempint2, &tempint3, &tempint4, &tempint5);
	  (*col)[iterelt+rowglobal] = tempint1+nzglobal;
	  (*col)[iterelt+rowglobal+1] = tempint2+nzglobal;
	  (*col)[iterelt+rowglobal+2] = tempint3+nzglobal;
	  (*col)[iterelt+rowglobal+3] = tempint4+nzglobal;
	  (*col)[iterelt+rowglobal+4] = tempint5+nzglobal;
	  iterelt+=5;
	  break;
	}

      /* read row */
      for (iterelt=0; iterelt<nzlocal+1-nbreltperline; iterelt++)
	{
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld %ld %ld %ld", &tempint1,&tempint2,&tempint3,&tempint4,&tempint5,&tempint6);
	  (*row)[iterelt+nzglobal]   = (pastix_int_t)tempint1;
	  (*row)[iterelt+nzglobal+1] = (pastix_int_t)tempint2;
	  (*row)[iterelt+nzglobal+2] = (pastix_int_t)tempint3;
	  (*row)[iterelt+nzglobal+3] = (pastix_int_t)tempint4;
	  (*row)[iterelt+nzglobal+4] = (pastix_int_t)tempint5;
	  (*row)[iterelt+nzglobal+5] = (pastix_int_t)tempint6;
	  iterelt+=5;
	}
      switch (nzlocal-iterelt)
	{
	case 1:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld",&tempint1);
	  (*row)[iterelt+nzglobal] = (pastix_int_t)tempint1;
	  iterelt++;
	  break;
	case 2:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld", &tempint1, &tempint2);
	  (*row)[iterelt+nzglobal]   = (pastix_int_t)tempint1;
	  (*row)[iterelt+nzglobal+1] = (pastix_int_t)tempint2;
	  iterelt+=2;
	  break;
	case 3:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld", &tempint1, &tempint2, &tempint3);
	  (*row)[iterelt+nzglobal]   = (pastix_int_t)tempint1;
	  (*row)[iterelt+nzglobal+1] = (pastix_int_t)tempint2;
	  (*row)[iterelt+nzglobal+2] = (pastix_int_t)tempint3;
	  iterelt+=3;
	  break;
	case 4:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld %ld", &tempint1, &tempint2, &tempint3, &tempint4);
	  (*row)[iterelt+nzglobal]   = (pastix_int_t)tempint1;
	  (*row)[iterelt+nzglobal+1] = (pastix_int_t)tempint2;
	  (*row)[iterelt+nzglobal+2] = (pastix_int_t)tempint3;
	  (*row)[iterelt+nzglobal+3] = (pastix_int_t)tempint4;
	  iterelt+=4;
	  break;
	case 5:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%ld %ld %ld %ld %ld", &tempint1, &tempint2, &tempint3, &tempint4, &tempint5);
	  (*row)[iterelt+nzglobal]   = (pastix_int_t)tempint1;
	  (*row)[iterelt+nzglobal+1] = (pastix_int_t)tempint2;
	  (*row)[iterelt+nzglobal+2] = (pastix_int_t)tempint3;
	  (*row)[iterelt+nzglobal+3] = (pastix_int_t)tempint4;
	  (*row)[iterelt+nzglobal+4] = (pastix_int_t)tempint5;
	  iterelt+=5;
	  break;
	}
      
      nbreltperline=3; /* nbr of elt per line */
      
      /* read val */
      for (iterelt=0; iterelt<nzlocal+1-nbreltperline; iterelt++)
	{
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3);
	  (*val)[iterelt+nzglobal]   = (pastix_float_t)tempfloat1;
	  (*val)[iterelt+nzglobal+1] = (pastix_float_t)tempfloat2;
	  (*val)[iterelt+nzglobal+2] = (pastix_float_t)tempfloat3;
	  iterelt+=2;
	}
      switch (nzlocal-iterelt)
	{
	case 1:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf",&tempfloat1);
	  (*val)[iterelt+nzglobal] = (pastix_float_t)tempfloat1;
	  iterelt++;
	  break;
	case 2:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf",&tempfloat1,&tempfloat2);
	  (*val)[iterelt+nzglobal]   = (pastix_float_t)tempfloat1;
	  (*val)[iterelt+nzglobal+1] = (pastix_float_t)tempfloat2;
	  iterelt+=2;
	  break;
	case 3:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3);
	  (*val)[iterelt+nzglobal]   = (pastix_float_t)tempfloat1;
	  (*val)[iterelt+nzglobal+1] = (pastix_float_t)tempfloat2;
	  (*val)[iterelt+nzglobal+2] = (pastix_float_t)tempfloat3;
	  iterelt+=3;
	  break;
	}
      nzglobal += nzlocal;

      /* read rhs */
      for (iterelt=0; iterelt<rowlocal+1-nbreltperline; iterelt++)
	{
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3);
	  (*rhs)[iterelt+rowglobal]   = (pastix_float_t)tempfloat1;
	  (*rhs)[iterelt+rowglobal+1] = (pastix_float_t)tempfloat2;
	  (*rhs)[iterelt+rowglobal+2] = (pastix_float_t)tempfloat3;
	  iterelt+=2;
	}

      switch (rowlocal-iterelt)
	{
	case 1:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf",&tempfloat1);
	  (*rhs)[iterelt+rowglobal] = (pastix_float_t)tempfloat1;
	  iterelt++;
	  break;
	case 2:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf",&tempfloat1,&tempfloat2);
	  (*rhs)[iterelt+rowglobal]   = (pastix_float_t)tempfloat1;
	  (*rhs)[iterelt+rowglobal+1] = (pastix_float_t)tempfloat2;
	  iterelt++;
	  break;
	case 3:
	  fgets(line,BUFSIZ,infile);
	  sscanf(line,"%Lf %Lf %Lf",&tempfloat1,&tempfloat2,&tempfloat3);
	  (*rhs)[iterelt+rowglobal]   = (pastix_float_t)tempfloat1;
	  (*rhs)[iterelt+rowglobal+1] = (pastix_float_t)tempfloat2;
	  (*rhs)[iterelt+rowglobal+2] = (pastix_float_t)tempfloat3;
	  iterelt++;
	  break;
	}
      rowglobal += rowlocal;

      fclose(infile);
    }
}

void mumpsReadHeader(FILE *infile, pastix_int_t *Nrow, pastix_int_t *Ncol, pastix_int_t *Nnzero, char *Type)
{
  char line[BUFSIZ];
  long temp1,temp2,temp3;
 
  Type[0] = 'R';
  Type[1] = 'U';
  Type[2] = 'A';
  Type[3] = '\0';

  /* ncol nrow nnzero */
  fgets(line,BUFSIZ,infile);

  sscanf(line, "%ld %ld %ld", &temp1,&temp2,&temp3);
  if (temp1!=temp2)
    {
      temp2=temp1;
      fgets(line,BUFSIZ,infile);
      sscanf(line, "%ld", &temp3);
    }

  *Nrow = (pastix_int_t)temp1;
  *Ncol = (pastix_int_t)temp2;
  *Nnzero = (pastix_int_t)temp3;
}

void mumpsRead(char const *dirname, 
	      pastix_int_t *Ncol, 
	      pastix_int_t *Nrow, 
	      pastix_int_t *Nnzero, 
	      pastix_int_t **col, 
	      pastix_int_t **row, 
	      pastix_float_t **val, 
	      char **Type, 
	      char **RhsType){

  FILE * iaFile;
  FILE * jaFile;
  FILE * raFile;
  FILE * headerFile;
  char * filename;
  pastix_int_t * tempcol;
  pastix_int_t * temprow;
  pastix_float_t * tempval;
  pastix_int_t iter,baseval;
  pastix_int_t tmp,total,pos,limit; 
  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(1*sizeof(char));
  (*RhsType)[0] = '\0';

  filename = malloc(strlen(dirname)+10);
  
  sprintf(filename,"%s/header",dirname);
  headerFile = fopen (filename,"r");
  if (headerFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }
  mumpsReadHeader(headerFile,Nrow,Ncol,Nnzero,*Type);
  fclose (headerFile);

  sprintf(filename,"%s/ia_mumps",dirname); 
  iaFile = fopen(filename,"r");  
  if (iaFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }

  sprintf(filename,"%s/ja_mumps",dirname);
  jaFile = fopen(filename,"r");
  if (jaFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }
  sprintf(filename,"%s/ra_mumps",dirname);
  raFile = fopen(filename,"r");
  if (raFile==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }

  /* Allocation memoire */
  tempcol = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  temprow = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  tempval = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  
  if ((tempcol==NULL) || (temprow == NULL) || (tempval == NULL))
    {
      fprintf(stderr, "mumpsRead : Not enough memory (Nnzero %ld)\n",(long)*Nnzero);
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  
  /* Remplissage */
  for (iter=0; iter<(*Nnzero); iter++)
    {
      long temp1,temp2;
      double tempv;
      fscanf(iaFile,"%ld\n", &temp1);
      temprow[iter]=(pastix_int_t)temp1;
      fscanf(jaFile,"%ld\n", &temp2);
      tempcol[iter]=(pastix_int_t)temp2;
      fscanf(raFile,"%le\n", &tempv);
      tempval[iter]= (pastix_float_t)tempv;
    }
  
  fclose (iaFile);
  fclose (jaFile);
  fclose (raFile);

  (*col) = (pastix_int_t *) malloc((*Nrow+1)*sizeof(pastix_int_t));  
  memset(*col,0,(*Nrow+1)*sizeof(pastix_int_t));  
  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  memset(*row,0,(*Nnzero)*sizeof(pastix_int_t));
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  if (((*col)==NULL) || ((*row) == NULL) || ((*val) == NULL))
    {
      fprintf(stderr, "mumpsRead : Not enough memory (Nnzero %ld)\n",(long)*Nnzero);
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  
  for (iter = 0; iter < (*Nnzero); iter ++)      
    {
      (*col)[tempcol[iter]-1]++;
    }

  baseval=1; /* Attention on base a 1 */
  total = baseval;
  
  for (iter = 0; iter < (*Ncol)+1; iter ++)      
    {
      tmp = (*col)[iter];
      (*col)[iter]=total;
      total+=tmp;
    }

  for (iter = 0; iter < (*Nnzero); iter ++)      
    {
      
      pos = (*col)[tempcol[iter]-1]-1;
      limit = (*col)[tempcol[iter]]-1;
      while((*row)[pos] != 0 && pos < limit)
	{
	  pos++;
	}
      if (pos == limit)
	fprintf(stderr, "Erreur de lecture\n");
      
      (*row)[pos] = temprow[iter];
      (*val)[pos] = tempval[iter];
    }      
  
  memFree_null(tempval);
  memFree_null(temprow);
  memFree_null(tempcol);
}


void MatrixMarketRead(char const *filename, 
		     pastix_int_t *Ncol, 
		     pastix_int_t *Nrow, 
		     pastix_int_t *Nnzero, 
		     pastix_int_t **col, 
		     pastix_int_t **row, 
		     pastix_float_t **val, 
		     char **Type, 
		     char **RhsType){

  FILE * file;
  pastix_int_t * tempcol;
  pastix_int_t iter,baseval;
  pastix_int_t * temprow;
  pastix_float_t * tempval;
  pastix_int_t total;
  pastix_int_t tmp;
  pastix_int_t pos;
  pastix_int_t limit;
  MM_typecode matcode;
  int tmpncol,tmpnrow,tmpnnzero;

  *Type = (char *) malloc(4*sizeof(char));
  *RhsType = (char *) malloc(1*sizeof(char));
  (*RhsType)[0] = '\0';

  file = fopen (filename,"r");
  if (file==NULL)
    {
      fprintf(stderr,"cannot load %s\n", filename);
      EXIT(MOD_SI,FILE_ERR);
    }
  
  if (mm_read_banner(file, &matcode) != 0)
    {
      fprintf(stderr,"Could not process Matrix Market banner.\n");
        exit(1);
    }
  
#ifdef    TYPE_COMPLEX
  (*Type)[0] = 'C';
  if (!mm_is_complex(matcode))
    {
      fprintf(stderr, "WARNING : Matrix should be complex.\n");
      (*Type)[0] = 'R';
    }
#else  /* TYPE_COMPLEX */
  (*Type)[0] = 'R';
  if (mm_is_complex(matcode))
    {
      fprintf(stderr, "WARNING : Matrix should not be complex. Only real part will be taken.\n");
    }
#endif /* TYPE_COMPLEX */

  (*Type)[1] = 'U';
  if (mm_is_symmetric(matcode))
    {
      (*Type)[1] = 'S';
    }
  (*Type)[2] = 'A';
  (*Type)[3] = '\0';
  /* find out size of sparse matrix .... */
  
  if (mm_read_mtx_crd_size(file, &tmpnrow, &tmpncol, &tmpnnzero) !=0)
    exit(1);

  *Ncol = tmpncol;
  *Nrow = tmpnrow;
  *Nnzero = tmpnnzero;

  /* Allocation memoire */
  tempcol = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  temprow = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  tempval = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  
  if ((tempcol==NULL) || (temprow == NULL) || (tempval == NULL))
    {
      fprintf(stderr, "MatrixMarketRead : Not enough memory (Nnzero %ld)\n",(long)*Nnzero);
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  
  /* Remplissage */
  {
    long temp1,temp2;
    double re,im;
    im = 0;
    for (iter=0; iter<(*Nnzero); iter++)
      {
	if (mm_is_complex(matcode)) 
	  {
	    fscanf(file,"%ld %ld %lg %lg\n", &temp1, &temp2, &re, &im);
	  }
	else
	  {
	    fscanf(file,"%ld %ld %lg\n", &temp1, &temp2, &re);
	  }
	temprow[iter]=(pastix_int_t)temp1;
	tempcol[iter]=(pastix_int_t)temp2;
#ifdef    TYPE_COMPLEX
	tempval[iter]=(pastix_float_t)(re+im*I);
#else  /* TYPE_COMPLEX */
	tempval[iter]=(pastix_float_t)(re);
#endif /* TYPE_COMPLEX */
      }
  }

  (*col) = (pastix_int_t *) malloc((*Nrow+1)*sizeof(pastix_int_t));  
  memset(*col,0,(*Nrow+1)*sizeof(pastix_int_t));  
  (*row) = (pastix_int_t *) malloc((*Nnzero)*sizeof(pastix_int_t));
  memset(*row,0,(*Nnzero)*sizeof(pastix_int_t));
  (*val) = (pastix_float_t *) malloc((*Nnzero)*sizeof(pastix_float_t));
  if (((*col)==NULL) || ((*row) == NULL) || ((*val) == NULL))
    {
      fprintf(stderr, "MatrixMarketRead : Not enough memory (Nnzero %ld)\n",(long)*Nnzero);
      EXIT(MOD_SI,OUTOFMEMORY_ERR);
    }
  
  for (iter = 0; iter < (*Nnzero); iter ++)      
    {
      (*col)[tempcol[iter]-1]++;
    }

  baseval=1; /* Attention on base a 1 */
  total = baseval;
  
  for (iter = 0; iter < (*Ncol)+1; iter ++)      
    {
      tmp = (*col)[iter];
      (*col)[iter]=total;
      total+=tmp;
    }

  for (iter = 0; iter < (*Nnzero); iter ++)      
    {
      
      pos = (*col)[tempcol[iter]-1]-1;
      limit = (*col)[tempcol[iter]]-1;
      while((*row)[pos] != 0 && pos < limit)
	{
	  pos++;
	}
      if (pos == limit)
	fprintf(stderr, "Erreur de lecture\n");
      
      (*row)[pos] = temprow[iter];
      (*val)[pos] = tempval[iter];
    }      
  
  memFree_null(tempval);
  memFree_null(temprow);
  memFree_null(tempcol);
}
