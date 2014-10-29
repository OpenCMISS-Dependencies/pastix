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
#include <sys/time.h>
#include <limits.h>

#ifndef MURGE_GETLOCALELEMENTLIST_C
#  define MURGE_GETLOCALELEMENTLIST_C
#else
#  ifndef SUFFIX
#    error murge_getLocalElementList.c can be included only once, or use SUFFIX(name)
#  endif /* SUFFIX */
#endif
#ifndef VERT_PER_ELEMENT
#  error VERT_PER_ELEMENT(MURGE_UserData_t d) need to be defined
#endif
#ifndef GET_VERTICES
#  error GET_VERTICES(INTS i_element, INTS * vertices, MURGE_UserData_t * d) must be defined
#endif
#ifdef SUFFIX
#  define MURGE_UserData_t          SUFFIX(MURGE_UserData_t)
#  define MURGE_GetLocalElementNbr  SUFFIX(MURGE_GetLocalElementNbr)
#  define MURGE_GetLocalElementList SUFFIX(MURGE_GetLocalElementList)
#  define MURGE_ElementListStorage  SUFFIX(MURGE_ElementListStorage)
#  define MURGE_ElementListStorageSize  \
  SUFFIX(MURGE_ElementListStorageSize)
#  define _element_  SUFFIX(_element_)
#  define _element_t SUFFIX(_element_t)
#endif

INTS * MURGE_ElementListStorage = NULL;
INTS   MURGE_ElementListStorageSize;

#ifdef MURGE_TIME
#  define MURGE_START_TIMER do {                                        \
    MPI_Barrier(comm);							\
    gettimeofday(&tv, NULL);                                            \
    t1 = ((double) tv.tv_sec + (double) tv.tv_usec * (double)1.0e-6L);  \
  } while(0)

#  define MURGE_STOP_TIMER(str) do {                                    \
    MPI_Barrier(comm);							\
    gettimeofday(&tv, NULL);                                            \
    t2 = ((double) tv.tv_sec + (double) tv.tv_usec * (double)1.0e-6L);  \
    if (rank == 0)                                                      \
      fprintf(stdout, " ++ time for '%-40s' %.2e s ++\n", str, t2-t1);  \
  } while(0)
#else /* not MURGE_TIME */
#  define MURGE_START_TIMER do {} while (0)
#  define MURGE_STOP_TIMER(str) do {} while (0)
#endif /* MURGE_TIME */

#define MPI_INTS (sizeof(long) == sizeof(INTS))?MPI_LONG:MPI_INT

#define MURGE_CHKERR(ierr) do {                 \
  if (ierr != MURGE_SUCCESS) return ierr;       \
  } while (0)

#define MURGE_MEMALLOC(ptr, size, type)                                 \
  do {                                                                  \
    if (NULL == ((ptr) = (type *) malloc((size) * sizeof(type))))       \
      {                                                                 \
        fprintf(stderr, "ERROR: Memory allocation error");              \
        return MURGE_ERR_ALLOCATE;                                      \
      }                                                                 \
  } while(0)

#define MURGE_FREE(ptr) do {                    \
    if (ptr) {                                  \
      free(ptr);                                \
      ptr = NULL;                               \
    }                                           \
  } while (0)

struct _element_ {
  INTS idx;
  INTS owner;
  INTS count;
};
typedef struct _element_ _element_t;

static inline
int cmp_element(const void *p1, const void *p2)
{
  _element_t * e1 = (_element_t *)p1;
  _element_t * e2 = (_element_t *)p2;


  if (e1->owner != e2->owner)
    return e1->owner - e2->owner;
  else
    if (e1->count != e2->count)
      return e2->count - e1->count;
    else
      return e1->idx - e2->idx;
}

#define INTSORTNAME                 sortElement
#define INTSORTSIZE                 (sizeof (_element_t))
#define INTSORTSWAP(p,q)                                                \
  do {                                                                  \
    _element_t t;                                                       \
    t = *((_element_t *) (p));                                          \
    *((_element_t *) (p)) = *((_element_t*) (q));                       \
    *((_element_t *) (q)) = t;                                          \
  } while (0)
#define INTSORTCMP(p,q)             (cmp_element((_element_t*)p,(_element_t*)q) < 0)


#ifndef MAX_THRESH

#define MAX_THRESH 6

#define max_thresh                  (MAX_THRESH * INTSORTSIZE) /* Variable turned into constant */

/* Stack node declarations used to store unfulfilled partition obligations. */
typedef struct
  {
    char *lo;
    char *hi;
  } stack_node;

/* The next 4 #defines implement a very fast in-line stack abstraction. */
/* The stack needs log (total_elements) entries (we could even subtract
   log(MAX_THRESH)).  Since total_elements has type size_t, we get as
   upper bound for log (total_elements):
   bits per unsigned char (CHAR_BIT) * sizeof(size_t).  */
#define STACK_SIZE	(CHAR_BIT * sizeof (size_t))
#define PUSH(low, high)	((void) ((top->lo = (low)), (top->hi = (high)), ++top))
#define	POP(low, high)	((void) (--top, (low = top->lo), (high = top->hi)))
#define	STACK_NOT_EMPTY	(stack < top)

#endif /* MAX_THRESH */

/* Order size using quicksort.  This implementation incorporates
   four optimizations discussed in Sedgewick:

   1. Non-recursive, using an explicit stack of pointer that store the
      next array partition to sort.  To save time, this maximum amount
      of space required to store an array of SIZE_MAX is allocated on the
      stack.  Assuming a 32-bit (64 bit) integer for size_t, this needs
      only 32 * sizeof(stack_node) == 256 bytes (for 64 bit: 1024 bytes).
      Pretty cheap, actually.

   2. Chose the pivot element using a median-of-three decision tree.
      This reduces the probability of selecting a bad pivot value and
      eliminates certain extraneous comparisons.

   3. Only quicksorts TOTAL_ELEMS / MAX_THRESH partitions, leaving
      insertion sort to order the MAX_THRESH items within each partition.
      This is a big win, since insertion sort is faster for small, mostly
      sorted array segments.

   4. The larger of the two sub-partitions is always pushed onto the
      stack first, with the algorithm then concentrating on the
      smaller partition.  This *guarantees* no more than log (total_elems)
      stack size is needed (actually O(1) in this case)!  */

/* To be defined :
** INTSORTNAME : Name of function
** INTSORTSIZE : Size of elements to sort
** INTSORTSWAP : Swapping macro
** INTSORTCMP  : Comparison function
*/

void
INTSORTNAME (
void * const                pbase,                /*+ Array to sort             +*/
const size_t                   total_elems)          /*+ Number of entries to sort +*/
{
  register char *base_ptr = (char *) pbase;

  if (total_elems == 0)
    /* Avoid lossage with unsigned arithmetic below.  */
    return;

  if (total_elems > MAX_THRESH)
    {
      char *lo = base_ptr;
      char *hi = &lo[INTSORTSIZE * (total_elems - 1)];
      stack_node stack[STACK_SIZE];
      stack_node *top = stack;

      PUSH (NULL, NULL);

      while (STACK_NOT_EMPTY)
        {
          char *left_ptr;
          char *right_ptr;

	  /* Select median value from among LO, MID, and HI. Rearrange
	     LO and HI so the three values are sorted. This lowers the
	     probability of picking a pathological pivot value and
	     skips a comparison for both the LEFT_PTR and RIGHT_PTR in
	     the while loops. */

	  char *mid = lo + INTSORTSIZE * ((hi - lo) / INTSORTSIZE >> 1);

	  if (INTSORTCMP ((void *) mid, (void *) lo))
	    INTSORTSWAP (mid, lo);
	  if (INTSORTCMP ((void *) hi, (void *) mid))
	    INTSORTSWAP (mid, hi);
	  else
	    goto jump_over;
	  if (INTSORTCMP ((void *) mid, (void *) lo))
	    INTSORTSWAP (mid, lo);
	jump_over:;

	  left_ptr  = lo + INTSORTSIZE;
	  right_ptr = hi - INTSORTSIZE;

	  /* Here's the famous ``collapse the walls'' section of quicksort.
	     Gotta like those tight inner loops!  They are the main reason
	     that this algorithm runs much faster than others. */
	  do
	    {
	      while (INTSORTCMP ((void *) left_ptr, (void *) mid))
		left_ptr += INTSORTSIZE;

	      while (INTSORTCMP ((void *) mid, (void *) right_ptr))
		right_ptr -= INTSORTSIZE;

	      if (left_ptr < right_ptr)
		{
		  INTSORTSWAP (left_ptr, right_ptr);
		  if (mid == left_ptr)
		    mid = right_ptr;
		  else if (mid == right_ptr)
		    mid = left_ptr;
		  left_ptr += INTSORTSIZE;
		  right_ptr -= INTSORTSIZE;
		}
	      else if (left_ptr == right_ptr)
		{
		  left_ptr += INTSORTSIZE;
		  right_ptr -= INTSORTSIZE;
		  break;
		}
	    }
	  while (left_ptr <= right_ptr);

          /* Set up pointers for next iteration.  First determine whether
             left and right partitions are below the threshold size.  If so,
             ignore one or both.  Otherwise, push the larger partition's
             bounds on the stack and continue sorting the smaller one. */

          if ((size_t) (right_ptr - lo) <= max_thresh)
            {
              if ((size_t) (hi - left_ptr) <= max_thresh)
		/* Ignore both small partitions. */
                POP (lo, hi);
              else
		/* Ignore small left partition. */
                lo = left_ptr;
            }
          else if ((size_t) (hi - left_ptr) <= max_thresh)
	    /* Ignore small right partition. */
            hi = right_ptr;
          else if ((right_ptr - lo) > (hi - left_ptr))
            {
	      /* Push larger left partition indices. */
              PUSH (lo, right_ptr);
              lo = left_ptr;
            }
          else
            {
	      /* Push larger right partition indices. */
              PUSH (left_ptr, hi);
              hi = right_ptr;
            }
        }
    }

  /* Once the BASE_PTR array is partially sorted by quicksort the rest
     is completely sorted using insertion sort, since this is efficient
     for partitions below MAX_THRESH size. BASE_PTR points to the beginning
     of the array to sort, and END_PTR points at the very last element in
     the array (*not* one beyond it!). */

#define min(x, y) ((x) < (y) ? (x) : (y))

  {
    char *const end_ptr = &base_ptr[INTSORTSIZE * (total_elems - 1)];
    char *tmp_ptr = base_ptr;
    char *thresh = min(end_ptr, base_ptr + max_thresh);
    register char *run_ptr;

    /* Find smallest element in first threshold and place it at the
       array's beginning.  This is the smallest array element,
       and the operation speeds up insertion sort's inner loop. */

    for (run_ptr = tmp_ptr + INTSORTSIZE; run_ptr <= thresh; run_ptr += INTSORTSIZE)
      if (INTSORTCMP ((void *) run_ptr, (void *) tmp_ptr))
        tmp_ptr = run_ptr;

    if (tmp_ptr != base_ptr)
      INTSORTSWAP (tmp_ptr, base_ptr);

    /* Insertion sort, running from left-hand-side up to right-hand-side.  */

    run_ptr = base_ptr + INTSORTSIZE;
    while ((run_ptr += INTSORTSIZE) <= end_ptr)
      {
	tmp_ptr = run_ptr - INTSORTSIZE;
	while (INTSORTCMP ((void *) run_ptr, (void *) tmp_ptr))
	  tmp_ptr -= INTSORTSIZE;

	tmp_ptr += INTSORTSIZE;
        if (tmp_ptr != run_ptr)
          {
            char *trav;

	    trav = run_ptr + INTSORTSIZE;
	    while (--trav >= run_ptr)
              {
                char c = *trav;
                char *hi, *lo;

                for (hi = lo = trav; (lo -= INTSORTSIZE) >= tmp_ptr; hi = lo)
                  *hi = *lo;
                *hi = c;
              }
          }
      }
  }
}

#undef INTSORTNAME
#undef INTSORTSIZE
#undef INTSORTSWAP
#undef INTSORTCMP

INTS MURGE_GetLocalElementNbr(INTS id,
                              INTS N,
                              INTS globalElementNbr,
                              INTS * localElementNbr,
                              INTS mode,
                              MURGE_UserData_t * d)
{
  MPI_Comm comm;
  int   rank;
  int   size;
  INTS  baseval;
  INTS  start;
  INTS  end;
  INTS  nodenbr;
  INTS *nodelist;
  INTS *owner;
  INTS  i,k;
  INTS  vertex_per_element = VERT_PER_ELEMENT(d);
  INTS *idx;
  INTS ierr;
  INTL nnz;
#ifdef MURGE_TIME
  struct timeval tv;
  double t1, t2;
#endif

  MURGE_FREE(MURGE_ElementListStorage);
  MURGE_MEMALLOC(idx, vertex_per_element, INTS);

  ierr = MURGE_GetCommRank(id, &rank); MURGE_CHKERR(ierr);
  ierr = MURGE_GetCommSize(id, &size); MURGE_CHKERR(ierr);
  ierr = MURGE_GetOptionINT(id, MURGE_IPARAM_BASEVAL, &baseval); MURGE_CHKERR(ierr);
  ierr = MURGE_GetComm(id, &comm); MURGE_CHKERR(ierr);

  start = rank*(globalElementNbr/size) +
    ((globalElementNbr%size) < rank ? (globalElementNbr%size) : rank);
  end = start + globalElementNbr/size +
    ((globalElementNbr%size) > rank);

  /* Build the graph */
  MURGE_START_TIMER;
  nnz = vertex_per_element*vertex_per_element*(end-start);
  ierr = MURGE_GraphBegin(id, N, nnz); MURGE_CHKERR(ierr);

  for (i=start; i<end; i++) {
    GET_VERTICES(i, idx, d);
    ierr = MURGE_GraphSetBlockEdges(id,
                                    vertex_per_element, idx,
                                    vertex_per_element, idx); MURGE_CHKERR(ierr);
  }
  MURGE_STOP_TIMER("Setting the graph");
  MURGE_START_TIMER;
  ierr = MURGE_GraphEnd(id); MURGE_CHKERR(ierr);
  MURGE_STOP_TIMER("GraphEnd");
  MURGE_START_TIMER;
  ierr = MURGE_GetLocalNodeNbr(id, &nodenbr);
  MURGE_STOP_TIMER("Getting local node number");
  MURGE_MEMALLOC(nodelist, nodenbr, INTS);
  MURGE_STOP_TIMER("Allocating node list");
  ierr = MURGE_GetLocalNodeList(id, nodelist);
  MURGE_STOP_TIMER("Getting local node list");

  /* owner contain -1 if the node is not local
   * and rank if it's local.
   */
  MURGE_MEMALLOC(owner, N, INTS);
  for (i =0; i < N; i++)
    owner[i] = -1;
  for (i =0; i < nodenbr; i++)
    owner[nodelist[i]-baseval] = rank;

  switch(mode)
    {
    case MURGE_DUPLICATE_ELEMENTS: {
      *localElementNbr = 0;
      for (i =0; i < globalElementNbr; i++)
        {
          GET_VERTICES(i, idx, d);
          for (k = 0; k < vertex_per_element; k++)
            if (owner[idx[k]-baseval] == rank) {
              (*localElementNbr)++;
              break;
            }
        }

      MURGE_ElementListStorageSize = (*localElementNbr);
      MURGE_MEMALLOC(MURGE_ElementListStorage, (*localElementNbr), INTS);
      (*localElementNbr) = 0;
      for (i =0; i < globalElementNbr; i++)
        {
          GET_VERTICES(i, idx, d);
          for (k = 0; k < vertex_per_element; k++)
            if (owner[idx[k]-baseval] == rank) {
              MURGE_ElementListStorage[(*localElementNbr)++] = i;
              break;
            }
        }
      break;
    }
    case MURGE_DISTRIBUTE_ELEMENTS: {
      {
        INTS * intarray;
        INTS * vertPerCPU;
        INTS * eltsPerCPU;
        INTS * eltsPerCPU_final;
        INTS * eltsPerCPU_curr;
        INTS * gowner;
        INTS   nb_inc;
        INTS   localElementNbr_avg;
        _element_t * elements;

        MURGE_MEMALLOC(intarray, 4*size, INTS);
        vertPerCPU = intarray; /* size : size */
        eltsPerCPU = intarray+size; /* size :
                                     *  - size + 1 durint creation,
                                     *  - size after
                                     */
        eltsPerCPU_final = intarray + 2*size;
        eltsPerCPU_curr  = intarray + 3*size;
        MURGE_MEMALLOC(elements, globalElementNbr, _element_t);
        /* For each node, gowner contain the rank of the processor
         * to which the node belong.
         */
        MURGE_MEMALLOC(gowner, N, INTS);
        ierr = MPI_Allreduce(owner, gowner, N, MPI_INTS, MPI_MAX, comm);

        if (ierr != MPI_SUCCESS) return MURGE_ERR_MPI;

        MURGE_START_TIMER;
        for (i =0; i < globalElementNbr; i++)
          {
            GET_VERTICES(i, idx, d);
            memset(vertPerCPU, 0, size*sizeof(INTS));
            /* count the number of vertices of i that belong
             * to each processor.
             */
            for (k = 0; k < vertex_per_element; k++)
              vertPerCPU[gowner[idx[k]-baseval]]++;
            /* we push to list of the processor that owns most
             * of the columns.
             */
            elements[i].idx   = i;
            elements[i].owner = 0;
            for (k = 1; k < size; k++)
              {
                if (vertPerCPU[elements[i].owner] < vertPerCPU[k])
                  elements[i].owner = k;
              }
            elements[i].count = vertPerCPU[elements[i].owner];
          }
        MURGE_STOP_TIMER("Preliminary distribution");
        MURGE_START_TIMER;
        sortElement(elements, globalElementNbr);
        MURGE_STOP_TIMER("qsort");

        /* count the number of elements per CPU */
        MURGE_START_TIMER;
        memset(eltsPerCPU, 0, size*sizeof(INTS));
        for (i=0; i < globalElementNbr; i++)
          {
            eltsPerCPU[elements[i].owner]++;
          }
        MURGE_STOP_TIMER("Counting element per CPU");

        /* Decide the number of element per CPU we want */
        MURGE_START_TIMER;
        localElementNbr_avg = globalElementNbr/size;
        nb_inc = 0;
        for (i = 0; i < size; i++)
          {
            eltsPerCPU_final[i] = localElementNbr_avg;
            if (eltsPerCPU[i] > localElementNbr_avg &&
                nb_inc < globalElementNbr%size)
              {
                nb_inc++;
                /* if (k == rank) */
                eltsPerCPU_final[i]++;
              }
          }
        MURGE_STOP_TIMER("Deciding number of element per CPU");

        /* Final distribution */
        {
          memset(eltsPerCPU_curr, 0, size*sizeof(INTS));
          MURGE_START_TIMER;
          for (i = 0; i < globalElementNbr; i++)
            {
              if (eltsPerCPU_curr[elements[i].owner] <
                  eltsPerCPU_final[elements[i].owner])
                eltsPerCPU_curr[elements[i].owner]++;
              else
                {
                  /* this elements must go to an other processor */
                  INTS orig_owner = elements[i].owner;
                  INTS orig_count = elements[i].count;
                  INTS max_count  = orig_count;

                  GET_VERTICES(i, idx, d);
                  memset(vertPerCPU, 0, size*sizeof(INTS));
                  for (k = 0; k < vertex_per_element; k++)
                    {
                      vertPerCPU[gowner[idx[k]-baseval]]++;
                    }
                  for (k = 0; k < size; k++)
                    {
                      /* only processor who have space available and
                       * possess less or equal number of vertices
                       * are availble */
                      if (max_count >= vertPerCPU[k] &&
                          eltsPerCPU[k] < localElementNbr_avg)
                        {
                          if (elements[i].owner == orig_owner &&
                              k != orig_owner)
                            {
                              /* we need to change owner,
                               * take the first we found */
                              elements[i].owner = k;
                              elements[i].count = vertPerCPU[k];
                            }
                          else
                            {
                              if (elements[i].count < vertPerCPU[k]) {
                                /* we found a better one */
                                elements[i].owner = k;
                                elements[i].count = vertPerCPU[k];
                              }
                            }
                        }
                      }
                  /* the element got a new owner */
                  eltsPerCPU[elements[i].owner]++;
                  eltsPerCPU[orig_owner]--;
                  eltsPerCPU_curr[elements[i].owner]++;
                }
            }
        }
        MURGE_STOP_TIMER("Computing final distribution");

        /* Now we should have the correct number of element on each processor */
        if (eltsPerCPU_curr[rank] != eltsPerCPU[rank])
          {
            fprintf(stdout, "INTERNAL_ERROR %ld != %ld\n",
                    (long)eltsPerCPU_curr[rank],
                    (long)eltsPerCPU[rank]);
            return MURGE_ERR_INTERNAL;
          }
        MURGE_START_TIMER;
        MURGE_ElementListStorageSize = eltsPerCPU_curr[rank];
        MURGE_MEMALLOC(MURGE_ElementListStorage,
                       eltsPerCPU_curr[rank], INTS);
        (*localElementNbr) = 0;
        for(i = 0; i < globalElementNbr; i++)
          if (elements[i].owner == rank) {
            MURGE_ElementListStorage[(*localElementNbr)++] = i;
          }
        MURGE_STOP_TIMER("Filling element list");
        MURGE_FREE(elements);
        MURGE_FREE(gowner);
        MURGE_FREE(intarray);
      }
      break;
    }
    default:
      fprintf(stderr, "Invalid mode");
      MURGE_FREE(nodelist);
      MURGE_FREE(owner);
      MURGE_FREE(idx);
      return MURGE_ERR_PARAMETER;
    }
  MURGE_FREE(nodelist);
  MURGE_FREE(owner);
  MURGE_FREE(idx);
  return MURGE_SUCCESS;
}

INTS MURGE_GetLocalElementList(INTS id, INTS * element_list) {
  INTS i, baseval;
  INTS ierr;
#ifdef MURGE_TIME
  struct timeval tv;
  double t1, t2;
  int rank;
  MPI_Comm comm;
  ierr = MURGE_GetComm(id, &comm); MURGE_CHKERR(ierr);
  ierr = MURGE_GetCommRank(id, &rank); MURGE_CHKERR(ierr);
#endif
  MURGE_START_TIMER;
  ierr = MURGE_GetOptionINT(id, MURGE_IPARAM_BASEVAL, &baseval); MURGE_CHKERR(ierr);
  if (MURGE_ElementListStorage == NULL)
    return MURGE_ERR_ORDER;
  for (i = 0; i < MURGE_ElementListStorageSize; i++)
    {
      element_list[i] = MURGE_ElementListStorage[i] + baseval;
    }
  MURGE_FREE(MURGE_ElementListStorage);
  MURGE_STOP_TIMER("Copying element list into user array");
  return MURGE_SUCCESS;
}

#ifdef SUFFIX
#  undef MURGE_UserData_t
#  undef MURGE_GetLocalElementNbr
#  undef MURGE_GetLocalElementList
#  undef MURGE_ElementListStorage
#  undef MURGE_ElementListStorageSize
#  undef _element_
#  undef _element_t
#endif

#undef MPI_INTS

#undef MURGE_CHKERR
#undef MURGE_MEMALLOC
#undef MURGE_FREE
#undef MURGE_START_TIMER
#undef MURGE_STOP_TIMER
