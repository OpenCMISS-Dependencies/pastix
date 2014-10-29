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
#if (defined WITH_STARPU)
#  include <hwloc.h>
#  include <starpu.h>
#endif /* WITH_STARPU */
#include "common_pastix.h"
#include "compute_context_nbr.h"

int compute_context_nbr(int * ctxnbr, int threadnbr, int verbose)
{
#if (defined WITH_STARPU && defined STARPU_CONTEXT)
  if (*ctxnbr == -1)
    {
      int depth;
      unsigned i, n;
      int ncpu_per_socket, nsocket, ncpu;
      hwloc_topology_t topology;
      hwloc_obj_t obj;
      /* Allocate and initialize topology object. */
      hwloc_topology_init(&topology);
      /* ... Optionally, put detection configuration here to ignore
         some objects types, define a synthetic topology, etc....
         The default is to detect all the objects of the machine that
         the caller is allowed to access.  See Configure Topology
         Detection. */
      /* Perform the topology detection. */
      hwloc_topology_load(topology);

      depth = hwloc_get_type_depth(topology, HWLOC_OBJ_SOCKET);
      if (depth == HWLOC_TYPE_DEPTH_UNKNOWN) {
        /* number of socket is unknow, let say we have quadcore... */
        ncpu_per_socket = 4;
      } else {
        ncpu    = hwloc_get_nbobjs_by_type(topology, HWLOC_OBJ_CORE);
        nsocket = hwloc_get_nbobjs_by_depth(topology, depth);
        ncpu_per_socket = ncpu/nsocket;
      }
      if (verbose > API_VERBOSE_NO)
        fprintf(stdout, "ncpu_per_socket %d\n", ncpu_per_socket);
      nsocket = threadnbr/ncpu_per_socket;
      if (threadnbr%ncpu_per_socket)
        nsocket ++;
      *ctxnbr = nsocket + 1;

      hwloc_topology_destroy(topology);
    }
  if (threadnbr + 1 < *ctxnbr)
    *ctxnbr = threadnbr+1;

  if (*ctxnbr == 1)
    *ctxnbr = 2;

  /* can't have more than STARPU_NMAX_SCHED_CTXS CTX */
  {
    PASTIX_INT nctx_bot = *ctxnbr - 1;
    while (*ctxnbr > STARPU_NMAX_SCHED_CTXS)
      {
        nctx_bot /= 2;
        *ctxnbr = nctx_bot + 1;
      }
  }
#endif /* WITH_STARPU */
  return 0;
}
