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
#include "raff_functions.h"

/*
 ** Section: Functions declarations
 */

/* Raffinement du second membre */
#define grad_smp          API_CALL(grad_smp)
#define grad_thread       API_CALL(grad_thread)

void* grad_smp         (void *arg);

/* Lancement d'une des fonctions seules */
void grad_thread (SolverMatrix *datacode, SopalinParam *sopaparam);
/*
 ** Section: Threads routines
 */


/*
 Function: API_CALL(grad_smp)

 Refine the solution using conjugate gradian method.

 Parameters:
 arg - Pointer to a <sopthread_data_t> structure containing
 the <Sopalin_Data_t> structure and the thread number ID.

 */
void* API_CALL(grad_smp)(void *arg)
{
  /* Choix du solveur */
  struct solver solveur = {NULL};
  Pastix_Solveur(&solveur);

  /* Variables */
  Clock  raff_clk;
  double t0      = 0;
  double t3      = 0;
  RAFF_FLOAT  tmp     = 0.0;
  RAFF_FLOAT  normr;
  RAFF_FLOAT  normb;
  RAFF_FLOAT  epsilon = solveur.Eps(arg);
  RAFF_INT    itermax = solveur.Itermax(arg);
  RAFF_INT    nb_iter = 0;
  RAFF_INT    n       = solveur.N(arg);

  RAFF_FLOAT *gradb = NULL;
  RAFF_FLOAT *gradr = NULL;
  RAFF_FLOAT *gradp = NULL;
  RAFF_FLOAT *gradz = NULL;
  RAFF_FLOAT *grad2 = NULL;
  RAFF_FLOAT *alpha = NULL;
  RAFF_FLOAT *beta  = NULL;
  RAFF_FLOAT *gradx = NULL;

  /* Initialisation des vecteurs */
  gradb = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  gradr = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  gradp = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  gradz = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  grad2 = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  alpha = solveur.Malloc(arg,     sizeof(RAFF_FLOAT));
  beta  = solveur.Malloc(arg,     sizeof(RAFF_FLOAT));
  gradx = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));

  gradb = solveur.Synchro(arg, (void*) gradb, 0);
  gradr = solveur.Synchro(arg, (void*) gradr, 1);
  gradp = solveur.Synchro(arg, (void*) gradp, 2);
  gradz = solveur.Synchro(arg, (void*) gradz, 3);
  grad2 = solveur.Synchro(arg, (void*) grad2, 4);
  gradx = solveur.Synchro(arg, (void*) gradx, 5);
  alpha = solveur.Synchro(arg, (void*) alpha, 6);
  beta  = solveur.Synchro(arg, (void*) beta,  7);

  RAFF_CLOCK_INIT;

  solveur.B(arg, gradb);
  solveur.X(arg, gradx);

  /* r=b-ax */
  solveur.bMAx(arg, gradb, gradx, gradr);
  normb = solveur.Norm(arg, gradb);
  normr = solveur.Norm(arg, gradr);
  tmp = normr / normb;

  /* z = M-1 r */
  solveur.Precond(arg, gradr, gradz, 1);

  solveur.Copy(arg, gradz, gradp, 1);

  while (((double)tmp > (double)epsilon) && (nb_iter < itermax))
    {
      RAFF_CLOCK_STOP;
      t0 = RAFF_CLOCK_GET;
      nb_iter++;

      /* grad2 = A * p */
      solveur.Ax(arg, gradp, grad2);

      /* alpha = <r, z> / <Ap, p> */
      solveur.Dotc(arg, gradr, gradz, beta, 0);
      solveur.Dotc(arg, grad2, gradp, alpha, 1);
      solveur.Div(arg, beta, alpha, alpha, 1);

      /* x = x + alpha * p */
      solveur.AXPY(arg, 1, alpha, gradx, gradp, 0);

      /* r = r - alpha * A * p */
      solveur.AXPY(arg, -1, alpha, gradr, grad2, 1);

      /* z = M-1 * r */
      solveur.Precond(arg, gradr, gradz, 1);

      /* beta = <r', z> / <r, z> */
      solveur.Dotc(arg, gradr, gradz, alpha, 0);
      solveur.Div(arg, alpha, beta, beta, 1);

      /* p = z + beta * p */
      solveur.BYPX(arg, beta, gradz, gradp, 1);

      normr = solveur.Norm(arg, gradr);
      tmp = normr / normb;

      RAFF_CLOCK_STOP;
      t3 = RAFF_CLOCK_GET;
      solveur.Verbose(arg, t0, t3, tmp, nb_iter);
      t0 = t3;
    }

  solveur.End(arg, tmp, nb_iter, t3, gradx);

  solveur.Free(arg, (void*) gradb);
  solveur.Free(arg, (void*) gradr);
  solveur.Free(arg, (void*) gradp);
  solveur.Free(arg, (void*) gradz);
  solveur.Free(arg, (void*) grad2);
  solveur.Free(arg, (void*) alpha);
  solveur.Free(arg, (void*) beta);
  solveur.Free(arg, (void*) gradx);
  return 0;
}

/*
 ** Section: Function creating threads
 */
void API_CALL(grad_thread)(SolverMatrix *datacode, SopalinParam *sopaparam)
{
  raff_thread(datacode, sopaparam, &API_CALL(grad_smp));
}
