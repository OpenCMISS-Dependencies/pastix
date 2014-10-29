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

#include "common_pastix.h"
#include "param_blend.h"

/* we set default option */
PASTIX_INT blendParamInit(BlendParam *param)
{
  param->hpf_filename     = NULL;
  param->trace_filename   = "traceBlend.trf";
  param->ps_filename      = "matrix.ps";
  param->hpf              = 0;
  param->tracegen         = 0;
  param->ps               = 0;
  param->assembly         = 0;
  param->solvmtx_filename = "solvmtx.";
  param->sequentiel       = 0;
  param->count_ops        = 1;
  param->debug            = 0;
  param->timer            = 1;
  param->recover          = 1;
  param->blcolmin         = 60;
  param->blcolmax         = 120;
  param->blblokmin        = 90;
  param->blblokmax        = 140;
  param->leader           = 0;
  param->allcand          = 0;
  param->nocrossproc      = 0;
  param->forceE2          = 0;
  param->level2D          = 100000000;
  param->candcorrect      = 0;
  param->clusterprop      = 0;
  param->costlevel        = 1;
  param->autolevel        = 0;
  param->malt_limit       = -1;
  param->smpnbr           = 0;
  /* On suppose que tous les noeuds smp utilisés ont la même configuration */
  param->procnbr          = sysconf(_SC_NPROCESSORS_ONLN);
  param->ratiolimit       = 0.0;
  param->dense_endblock   = 0;
  param->ooc              = 0;
  param->oocmemlimit      = 4e7;
  param->abs              = 4;
  param->ricar            = 0;
  
  return 1;
}

void blendParamExit(BlendParam *param)
{
  memFree_null(param->hpf_filename);
  memFree_null(param->trace_filename);
  memFree_null(param->ps_filename);
  memFree_null(param->solvmtx_filename);
  memFree_null(param);
}


