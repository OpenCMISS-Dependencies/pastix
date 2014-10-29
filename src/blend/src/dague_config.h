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
#ifndef DAGUE_CONFIG_H_HAS_BEEN_INCLUDED
#define DAGUE_CONFIG_H_HAS_BEEN_INCLUDED

#define DAGUE_VERSION_MAJOR 1
#define DAGUE_VERSION_MINOR 0

/* Scheduling engine */
#define DAGUE_SCHED_DEPS_MASK
/* #undef DAGUE_SCHED_REPORT_STATISTICS */

/* Communication engine */
/* #undef DAGUE_DIST_WITH_MPI */
#define DAGUE_DIST_THREAD
#define DAGUE_DIST_PRIORITIES
/* #undef DAGUE_DIST_COLLECTIVES */
#define DAGUE_DIST_EAGER_LIMIT 0
#define DAGUE_DIST_SHORT_LIMIT 0

/* GPU Support */
#define DAGUE_GPU_WITH_CUDA
/* #undef DAGUE_GPU_CUDA_ALLOC_PER_TILE */
/* #undef DAGUE_GPU_WITH_OPENCL */
#define DAGUE_HAVE_PEER_DEVICE_MEMORY_ACCESS

/* debug */
#define DAGUE_DEBUG_VERBOSE 0
/* #undef DAGUE_DEBUG_HISTORY */
/* #undef DAGUE_DEBUG */
#define DAGUE_DEBUG_LIFO_USE_ATOMICS
/* #undef DAGUE_CALL_TRACE */

/* profiling */
/* #undef DAGUE_PROF_TRACE */
/* #undef DAGUE_PROF_STATS */
/* #undef DAGUE_PROF_PAPI */
/* #undef DAGUE_PROF_GRAPHER */
/* #undef DAGUE_PROF_DRY_RUN */
/* #undef DAGUE_PROF_DRY_BODY */
/* #undef DAGUE_PROF_DRY_DEP */

/* Simulating */
/* #undef DAGUE_SIM */

/* DSPARSE Options */
/* #undef DSPARSE_WITH_SOLVE */

/* system */
#define HAVE_PTHREAD
#define HAVE_SCHED_SETAFFINITY
#define HAVE_CLOCK_GETTIME
#define DAGUE_ATOMIC_USE_GCC_32_BUILTINS
#define DAGUE_ATOMIC_USE_GCC_64_BUILTINS
/* #undef DAGUE_ATOMIC_USE_XLC_32_BUILTINS */
/* #undef DAGUE_ATOMIC_USE_XLC_64_BUILTINS */
/* #undef DAGUE_ATOMIC_USE_MIPOSPRO_32_BUILTINS */
/* #undef DAGUE_ATOMIC_USE_MIPOSPRO_64_BUILTINS */
/* #undef DAGUE_ATOMIC_USE_SUN_32 */
/* #undef DAGUE_ATOMIC_USE_SUN_64 */
#define HAVE_ASPRINTF
#define HAVE_VASPRINTF
#define HAVE_STDARG_H
#define HAVE_UNISTD_H
#define HAVE_VA_COPY
/* #undef HAVE_UNDERSCORE_VA_COPY */
#define HAVE_GETOPT_LONG
#define HAVE_GETRUSAGE
#define HAVE_GETOPT_H
#define HAVE_ERRNO_H
#define HAVE_STDDEF_H
#define HAVE_LIMITS_H
#define HAVE_STRING_H
#define HAVE_COMPLEX_H
/* #undef ARCH_X86 */
#define ARCH_X86_64
/* #undef ARCH_PPC */
/* #undef MAC_OS_X */

/* Optional packages */
#define HAVE_HWLOC
#define HAVE_HWLOC_BITMAP
#define HAVE_HWLOC_PARENT_MEMBER
#define HAVE_HWLOC_CACHE_ATTR
#define HAVE_HWLOC_OBJ_PU

/* #undef HAVE_PAPI */
#define HAVE_CUDA
/* #undef HAVE_OPENCL */
/* #undef HAVE_MPI */
/* #undef HAVE_MPI_20 */

/* #undef HAVE_AYUDAME */

#if !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif  /* !defined(_GNU_SOURCE) */

#ifdef ARCH_PPC
#define inline __inline__
#define restrict 
#endif

#include "dague_config_bottom.h"

#endif  /*DAGUE_CONFIG_H_HAS_BEEN_INCLUDED */
