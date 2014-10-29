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
#define gmres_smp         API_CALL(gmres_smp)
#define gmres_thread      API_CALL(gmres_thread)

void* gmres_smp        (void *arg);

/* Lancement d'une des fonctions seules */
void gmres_thread(SolverMatrix *datacode, SopalinParam *sopaparam);

/*
 ** Section: Threads routines
 */

/*
 Function: API_CALL(gmres_smp)

 Function computing GMRES iterative reffinement.

 Parameters:
 arg - Pointer to a <sopthread_data_t> structure containing
 the <Sopalin_Data_t> structure and the thread number ID.
 */

typedef struct gmres_s
{
  volatile PASTIX_INT gmresout_flag;     /*+ Flag for GMRES outter loop          +*/
  volatile PASTIX_INT gmresin_flag;      /*+ Flag for GMRES inner loop           +*/
  volatile double     gmresro;           /*+ Norm of GMRES residue               +*/
} gmres_t;

void* API_CALL(gmres_smp)(void *arg)
{
  struct solver solveur = {NULL};
  Pastix_Solveur(&solveur);

  RAFF_INT               n            = solveur.N(arg);
  Clock             raff_clk;
  double            t0           = 0;
  double            t3           = 0;
  RAFF_FLOAT          *  gmrestemp    = NULL;
  volatile RAFF_INT      gmresim      = 0;
  volatile RAFF_INT      gmresmaxits  = 0;
  RAFF_FLOAT          *  gmresb       = NULL;
  RAFF_FLOAT          ** gmresvv      = NULL;
  RAFF_FLOAT          ** gmreshh      = NULL;
  RAFF_FLOAT          *  gmresc       = NULL;
  RAFF_FLOAT          *  gmress       = NULL;
  RAFF_FLOAT          *  gmresrs      = NULL;
  RAFF_FLOAT          ** gmresw       = NULL;
  RAFF_FLOAT             gmresalpha;
  RAFF_FLOAT             gmrest;
  volatile RAFF_INT      gmresiters   = 0;
  RAFF_FLOAT          *  gmreswk1;
  RAFF_FLOAT          *  gmreswk2     = NULL;
  volatile double   gmreseps     = 0;
  volatile double   gmresnormb;
  volatile RAFF_INT      gmresi1      = 0;
  volatile RAFF_INT      i = 0;
  RAFF_INT               j, ii, k;
  RAFF_FLOAT             beta;
  RAFF_FLOAT          *  gmresx       = NULL;
  gmres_t        *  gmresdata;
  gmresim     = solveur.Krylov_Space(arg);
  gmresmaxits = solveur.Itermax(arg);
  gmreseps    = solveur.Eps(arg);

  gmrestemp = solveur.Malloc(arg, n           * sizeof(RAFF_FLOAT));
  gmresb    = solveur.Malloc(arg, n           * sizeof(RAFF_FLOAT));
  gmresc    = solveur.Malloc(arg, gmresim     * sizeof(RAFF_FLOAT));
  gmress    = solveur.Malloc(arg, gmresim     * sizeof(RAFF_FLOAT));
  gmresrs   = solveur.Malloc(arg, (gmresim+1) * sizeof(RAFF_FLOAT));
  gmresdata = solveur.Malloc(arg, 1           * sizeof(gmres_t));
  gmresx    = solveur.Malloc(arg, n           * sizeof(RAFF_FLOAT));

  MONO_BEGIN(arg);
  gmresvv = solveur.Malloc(arg, (gmresim+1) * sizeof(RAFF_FLOAT*));
  gmreshh = solveur.Malloc(arg, gmresim     * sizeof(RAFF_FLOAT*));
  gmresw  = solveur.Malloc(arg, gmresim     * sizeof(RAFF_FLOAT*));
  for (i=0; i<gmresim; i++)
    {
      gmresvv[i] = solveur.Malloc(arg, n           * sizeof(RAFF_FLOAT));
      gmreshh[i] = solveur.Malloc(arg, (gmresim+1) * sizeof(RAFF_FLOAT));
      gmresw[i]  = solveur.Malloc(arg, n           * sizeof(RAFF_FLOAT));
    }
  gmresvv[gmresim] = solveur.Malloc(arg, n * sizeof(RAFF_FLOAT));
  MONO_END(arg);
  SYNCHRO(arg);

  /* Synchronisations */
  gmrestemp  = (RAFF_FLOAT * )solveur.Synchro(arg, (void*) gmrestemp, 0);
  gmresb     = (RAFF_FLOAT * )solveur.Synchro(arg, (void*) gmresb,    1);
  gmresc     = (RAFF_FLOAT * )solveur.Synchro(arg, (void*) gmresc,    2);
  gmress     = (RAFF_FLOAT * )solveur.Synchro(arg, (void*) gmress,    3);
  gmresrs    = (RAFF_FLOAT * )solveur.Synchro(arg, (void*) gmresrs,   4);
  gmresvv    = (RAFF_FLOAT **)solveur.Synchro(arg, (void*) gmresvv,   6);
  gmreshh    = (RAFF_FLOAT **)solveur.Synchro(arg, (void*) gmreshh,   7);
  gmresw     = (RAFF_FLOAT **)solveur.Synchro(arg, (void*) gmresw,    8);
  gmresdata  = (gmres_t*)solveur.Synchro(arg, (void*) gmresdata, 9);

  gmresnormb = (double)(*((double*)solveur.Synchro(arg, (void*) &gmresnormb, 10)));
  gmresx     = (RAFF_FLOAT * )solveur.Synchro(arg, (void*) gmresx,    11);

  gmresdata->gmresro = 0.0;
  gmresdata->gmresout_flag = 1;

  solveur.B(arg, gmresb);
  gmresnormb = solveur.Norm(arg, gmresb);

  solveur.X(arg, gmresx);

  gmresalpha = -1.0;
  gmresiters = 0;

  RAFF_CLOCK_INIT;

  while (gmresdata->gmresout_flag)
    {
      gmreswk2 = gmresvv[0];

      /* gmresvv[0] = b - A * x */
      solveur.bMAx(arg, gmresb, gmresx, gmresvv[0]);

      /* ro = vv[0].vv[0] */
      solveur.Dotc_Gmres(arg, gmresvv[0], gmresvv[0], &beta, 0);

#ifdef CPLX
      gmresdata->gmresro = (RAFF_FLOAT)csqrt(beta);
#else
      gmresdata->gmresro = (RAFF_FLOAT)sqrt(beta);
#endif

      if ((double)ABS_FLOAT((RAFF_FLOAT)gmresdata->gmresro) <=
          gmreseps)
        {
          gmresdata->gmresout_flag = 0;
          break;
        }

      gmrest = (RAFF_FLOAT)(1.0/gmresdata->gmresro);

      solveur.Scal(arg, gmrest, gmresvv[0], 1);

      gmresrs[0] = (RAFF_FLOAT)gmresdata->gmresro;
      gmresdata->gmresin_flag = 1;
      i=-1;

      while(gmresdata->gmresin_flag)
        {
          RAFF_CLOCK_STOP;
          t0 = RAFF_CLOCK_GET;

          i++;
          gmresi1 = i+1;

          gmreswk1 = gmresvv[i];
          gmreswk2 = gmresw[i];

          SYNCHRO(arg);
          solveur.Precond(arg, gmreswk1, gmreswk2, 1);

          gmreswk1 = gmresvv[gmresi1];

          /* vv[i1] = A*wk2 */
          solveur.Ax(arg, gmreswk2, gmreswk1);

          /* classical gram - schmidt */
          for (j=0; j<=i; j++)
            {
              /* vv[j]*vv[i1] */
              solveur.Dotc_Gmres(arg,gmresvv[gmresi1], gmresvv[j], &beta, 0);

              gmreshh[i][j] = (RAFF_FLOAT)beta;
            }

          SYNCHRO(arg);

          for (j=0;j<=i;j++)
            {
              gmresalpha = -gmreshh[i][j];
              solveur.AXPY(arg, 1.0, &gmresalpha, gmresvv[gmresi1], gmresvv[j], 0);
            }

          SYNCHRO(arg);
          solveur.Dotc_Gmres(arg, gmresvv[gmresi1], gmresvv[gmresi1], &beta, 0);

#ifdef CPLX
      gmrest = (RAFF_FLOAT)csqrt(beta);
#else
      gmrest = (RAFF_FLOAT)sqrt(beta);
#endif

          gmreshh[i][gmresi1] = gmrest;

          if (ABS_FLOAT(gmrest) > 10e-50)
            {
              gmrest = fun / gmrest;
              solveur.Scal(arg, gmrest, gmresvv[gmresi1], 0);
            }

          SYNCHRO(arg);
          MONO_BEGIN(arg);

          if (i != 0)
            {
              for (j=1; j<=i;j++)
                {
                  gmrest = gmreshh[i][j-1];
#ifdef CPLX
                  gmreshh[i][j-1] = (RAFF_FLOAT)conj(gmresc[j-1])*gmrest +
                    (RAFF_FLOAT)conj(gmress[j-1])*gmreshh[i][j];
#else /* CPLX */
                  gmreshh[i][j-1] =  gmresc[j-1]*gmrest +
                    gmress[j-1]*gmreshh[i][j];
#endif /* CPLX */
                  gmreshh[i][j]   = -gmress[j-1]*gmrest +
                    gmresc[j-1]*gmreshh[i][j];
                }
            }
#ifdef CPLX
          gmrest = (RAFF_FLOAT)csqrt(ABS_FLOAT(gmreshh[i][i]*gmreshh[i][i])+
                                       gmreshh[i][gmresi1]*gmreshh[i][gmresi1]);
#else
          gmrest = (RAFF_FLOAT)sqrt(gmreshh[i][i]*gmreshh[i][i]+
                                      gmreshh[i][gmresi1]*gmreshh[i][gmresi1]);
#endif
          if (ABS_FLOAT(gmrest) <= gmreseps)
            gmrest = (RAFF_FLOAT)gmreseps;

          gmresc[i] = gmreshh[i][i]/gmrest;
          gmress[i] = gmreshh[i][gmresi1]/gmrest;
          gmresrs[gmresi1] = -gmress[i]*gmresrs[i];

#ifdef CPLX
          gmresrs[i] = (RAFF_FLOAT)conj(gmresc[i])*gmresrs[i];
          gmreshh[i][i] = (RAFF_FLOAT)conj(gmresc[i])*gmreshh[i][i] +
          gmress[i]*gmreshh[i][gmresi1];
#else
          gmresrs[i] = gmresc[i]*gmresrs[i];
          gmreshh[i][i] = gmresc[i]*gmreshh[i][i] +
          gmress[i]*gmreshh[i][gmresi1];
#endif
          gmresdata->gmresro = ABS_FLOAT(gmresrs[gmresi1]);

          MONO_END(arg);

          gmresiters++;

          MONO_BEGIN(arg);
          if ((i+1 >= gmresim) || (gmresdata->gmresro/gmresnormb <= gmreseps) || (gmresiters >= gmresmaxits))
            {
              gmresdata->gmresin_flag = 0;
            }
          MONO_END(arg);

          RAFF_CLOCK_STOP;
          t3 = RAFF_CLOCK_GET;
          solveur.Verbose(arg, t0, t3, gmresdata->gmresro/gmresnormb, gmresiters);
          SYNCHRO(arg);
        }

      MONO_BEGIN(arg);

      gmresrs[i] = gmresrs[i]/gmreshh[i][i];
      for (ii=2; ii<=i+1; ii++)
        {
          k = i-ii+1;
          gmrest = gmresrs[k];
          for (j=k+1; j<=i; j++)
            {
              gmrest = gmrest - gmreshh[j][k]*gmresrs[j];
            }
          gmresrs[k] = gmrest/gmreshh[k][k];
        }

      MONO_END(arg);
      SYNCHRO(arg);

      for (j=0; j<=i;j++)
        {
          gmrest = gmresrs[j];
          solveur.AXPY(arg, 1.0, &gmrest, gmresx, gmresw[j], 0);
        }
      SYNCHRO(arg);

      if ((gmresdata->gmresro/gmresnormb<= gmreseps) || (gmresiters >= gmresmaxits))
        {
          gmresdata->gmresout_flag = 0;
        }
    }

  RAFF_CLOCK_STOP;
  t3 = RAFF_CLOCK_GET;

  solveur.End(arg, gmresdata->gmresro/gmresnormb, gmresiters, t3, gmresx);

  solveur.Free(arg, (void*) gmrestemp);
  solveur.Free(arg, (void*) gmresb);
  solveur.Free(arg, (void*) gmresc);
  solveur.Free(arg, (void*) gmress);
  solveur.Free(arg, (void*) gmresrs);
  solveur.Free(arg, (void*) gmresdata);
  solveur.Free(arg, (void*) gmresx);

  MONO_BEGIN(arg);
  for (i=0; i<gmresim; i++)
    {
      solveur.Free(arg, gmresvv[i]);
      solveur.Free(arg, gmreshh[i]);
      solveur.Free(arg, gmresw[i]);
    }

  solveur.Free(arg, gmresvv[gmresim]);

  solveur.Free(arg, gmresvv);
  solveur.Free(arg, gmreshh);
  solveur.Free(arg, gmresw);
  MONO_END(arg);

  return 0;
}

/*
** Section: Function creating threads
*/
void API_CALL(gmres_thread)(SolverMatrix *datacode, SopalinParam *sopaparam)
{
  raff_thread(datacode, sopaparam, &API_CALL(gmres_smp));
}
