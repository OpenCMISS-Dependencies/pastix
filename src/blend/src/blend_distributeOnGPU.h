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
#include "dof.h"
#include "cost.h"
#include "ftgt.h"
#include "symbol.h"
#include "queue.h"
#include "bulles.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "queue.h"
#include "sopalin_acces.h"

#ifndef GPU_MAX_FILL
#  define GPU_MAX_FILL 0.7
#endif
#ifndef GPU_MIN_NBPAGES
#  define GPU_MIN_NBPAGES 0
#endif
#ifndef GPU_MIN_UPDATES
#  define GPU_MIN_UPDATES 0
#endif
#ifndef GPU_MIN_FLOP
#  define GPU_MIN_FLOP 0
#endif


/* Alter threadid of tasks to put some of them on GPUs (threadid >= thread_nbr)
 * cblk are taken from the queue until all GPUs are alocated maxMem memory.
 */
#define blend_distributeOnGPU PASTIX_PREFIX(blend_distributeOnGPU)
int blend_distributeOnGPU(SolverMatrix  * solvmtr,
                          double          maxMem,
                          int             pageSize,
                          int             criterium,
                          enum API_GPU_CRITERIUM nGPUs,
                          enum API_FLOAT  floatType,
                          enum API_FACT   factType);
