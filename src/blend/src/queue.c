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
/*
  File: Queue.c

  Operations on the queue structure.

*/
#include <stdio.h>

#include "common_pastix.h"
#include "queue.h"


/*
  Function: QueueInit

  Allocate the queue array
  
  Parameters:
    q    - The queue to initialize
    size - The initial size of the queue.

  Return:
    NO_ERR - If all goes well.
*/
int queueInit(Queue *q,
	      PASTIX_INT    size)
{
    q->size = size;
    q->used = 0;
    if (q->size != 0)
      {
	MALLOC_INTERN(q->elttab,  size, PASTIX_INT);
	MALLOC_INTERN(q->keytab,  size, double);
	MALLOC_INTERN(q->keytab2, size, PASTIX_INT);
      }
    else
      {
	q->elttab  = NULL;
	q->keytab  = NULL;
	q->keytab2 = NULL;
      }
    return NO_ERR;
}
/*
  Function: queueCopy

  Perform a copy of a queue.

  Parameters:
    dst - Destination of the copy.
    src - Source of the copy.

  Returns:
    The copy address or NULL if the destination
    or the source is NULL.
*/
Queue * queueCopy(Queue *dst, 
		  Queue *src)
{
  if(src == NULL || dst == NULL)
    return NULL;
  memcpy(dst, src, sizeof(Queue));
  MALLOC_INTERN(dst->elttab,  src->size, PASTIX_INT);
  memcpy(dst->elttab, src->elttab, src->size * sizeof(PASTIX_INT));

  MALLOC_INTERN(dst->keytab,  src->size, double);
  memcpy(dst->keytab, src->keytab, src->size * sizeof(double));

  MALLOC_INTERN(dst->keytab2, src->size, PASTIX_INT);
  memcpy(dst->keytab2, src->keytab2, src->size * sizeof(PASTIX_INT));

  return dst;
}

/*
  Function: queueExit

  Free a queue structure.

  Parameters:
    q - The queue to free.
*/
void queueExit(Queue *q)
{
  if(q->size != 0)
    {
      memFree_null(q->elttab);
      memFree_null(q->keytab);
      memFree_null(q->keytab2);
    }
  q->size = 0;
}

/*
  Function: queueAdd

  Add an element in a queue, following one
  double value as key.

  The element has to be a positive integer.

  Parameters:
    q   - the queue to fill.
    elt - the element to add.
    key - The key associated to the element.
*/
void queueAdd(Queue *q, 
	      PASTIX_INT    elt, 
	      double key)
{
  queueAdd2(q, elt, key, 0);
}

/*
  Function: queueAdd2

  Add an element in a queue, following one
  double value and an integer value as keys.

  The element has to be a positive integer.

  Parameters:
    q    - the queue to fill.
    elt  - the element to add.
    key  - The double key associated to the element.
    key2 - The integer key associated to the element.
*/
void queueAdd2(Queue *q, 
	       PASTIX_INT    elt,
	       double key, 
	       PASTIX_INT    key2)
{
    PASTIX_INT i;
    PASTIX_INT   swap_elt;
    double swap_key;
    PASTIX_INT   swap_key2;
    PASTIX_INT *tmp;
    double * tmp2;
    PASTIX_INT   * tmp3;


    /** Allocate more space if necessary **/
    if(q->size == q->used)
	{
	    tmp  = q->elttab;
	    tmp2 = q->keytab;
	    tmp3 = q->keytab2;
	    /* OIMBE Realloc ?? */
	    MALLOC_INTERN(q->elttab, q->size*2+1, PASTIX_INT);
	    memcpy(q->elttab, tmp, q->size * sizeof(PASTIX_INT));

	    MALLOC_INTERN(q->keytab, q->size*2+1, double);
	    memcpy(q->keytab, tmp2, q->size * sizeof(double));

	    MALLOC_INTERN(q->keytab2, q->size*2+1, PASTIX_INT);
	    memcpy(q->keytab2, tmp3, q->size * sizeof(PASTIX_INT));

	    q->size = q->size*2 +1;
	    if (tmp != NULL)
	      memFree_null(tmp);
	    if (tmp2 != NULL)
	      memFree_null(tmp2);
	    if (tmp3 != NULL)
	      memFree_null(tmp3);
	}

    q->elttab[q->used] = elt;
    q->keytab[q->used] = key;
    q->keytab2[q->used] = key2;
    q->used++;
    i = q->used;

    while( (i>1) 
	   &&  compWith2keys(q, i-1, i/2-1))
	{
	    swap_elt = q->elttab[i-1];
	    swap_key = q->keytab[i-1];
	    swap_key2 = q->keytab2[i-1];
	    q->elttab[i-1] = q->elttab[i/2-1];
	    q->keytab[i-1] = q->keytab[i/2-1];
	    q->keytab2[i-1] = q->keytab2[i/2-1];

	    q->elttab[i/2-1] = swap_elt;
	    q->keytab[i/2-1] = swap_key;
	    q->keytab2[i/2-1] = swap_key2;
	    i=i/2;
	}
}

/*
  Function: queueGet

  Get next element of the queue and 
  remove it from the queue.

  Parameters: 
    q - The queue from which user wants an element.
    
  Returns:
    The element if it was found, or -1 if it wasn't.
*/
PASTIX_INT queueGet(Queue *q)
{
    PASTIX_INT i, j;
    PASTIX_INT return_elt;
    PASTIX_INT swap_elt;
    double swap_key;
    PASTIX_INT   swap_key2;

    if (q->used == 0)
      return -1;

    return_elt = q->elttab[0];
    
    q->elttab[0]  = q->elttab[q->used-1];
    q->keytab[0]  = q->keytab[q->used-1];
    q->keytab2[0] = q->keytab2[q->used-1];
    q->used--;
    
    i=1;
    
    while(i <= (q->used/2))
	{
	    if( (2*i == q->used)
		|| compWith2keys(q, 2*i-1, 2*i))     /*(q->keytab[2*i-1] < q->keytab[2*i]))*/
		{
		    j = 2*i;
		}
	    else
		{
		    j = 2*i+1;
		}
	    if (!compWith2keys(q, i-1, j-1))         /*(q->keytab[i-1] >= q->keytab[j-1])*/
		{
		    swap_elt = q->elttab[i-1];
		    swap_key = q->keytab[i-1];
		    swap_key2 = q->keytab2[i-1];

		    q->elttab[i-1]  = q->elttab[j-1];
		    q->keytab[i-1]  = q->keytab[j-1];
		    q->keytab2[i-1] = q->keytab2[j-1];

		    q->elttab[j-1] = swap_elt;
		    q->keytab[j-1] = swap_key;
		    q->keytab2[j-1] = swap_key2;
		    
		    i=j;
		}
	    else
		break;
	}
    return return_elt;
}


/*
  Function: queueSize
 
  Compute the size of a queue.

  Parameters:
    q - the queue.

  Returns:
    The size of the queue.
*/
PASTIX_INT queueSize(Queue *q)
{
  return q->used;
}


void queueClear(Queue *q)
{
  q->used = 0;
}


/*
  Function: queueRead

  Read the next element that 'll be given by queueGet 
  but not suppress it from the queue 

  Parameters:
    q - The queue.

  Returns:
    The next element.
*/
PASTIX_INT queueRead(Queue *q)
{
  return q->elttab[0];
}

/*
  Function: compwith2keys

  Compare 2 elements following their two keys.

  Parameters:
    q    - compare two elements
    elt1 - index of the first element in the queue.
    elt2 - index of the second element in the queue.
*/
PASTIX_INT compWith2keys(Queue *q, 
		  PASTIX_INT    elt1, 
		  PASTIX_INT    elt2)
{
  /* if elt1 < elt2 return 1  */
  /* if elt1 = elt2 return 0  */
  /* if elt1 > elt2 return 0 */
  
  if(q->keytab[elt1] < q->keytab[elt2])
    return 1;
  if(q->keytab[elt1] > q->keytab[elt2]) 
    return 0;
  if(q->keytab2[elt1] < q->keytab2[elt2])
    return 1;
  return 0;
}
/*
  Function: queueGet2

  Get next element of the queue and remove it from the queue.
  
  Parameters:
    q    - The queue.
    key  - The first key (double) of the element.
    key2 - The second key (integer) of the element.

  Returns:
    The element, or -1 if not found.
*/
PASTIX_INT queueGet2(Queue  *q,  
	      double *key, 
	      PASTIX_INT    *key2)
{
    PASTIX_INT i, j;
    PASTIX_INT return_elt;
    PASTIX_INT swap_elt;
    double swap_key;
    PASTIX_INT   swap_key2;

    if (q->used == 0)
      return -1;

    return_elt = q->elttab[0];
    if (key  != NULL) (*key)  = q->keytab[0];
    if (key2 != NULL) (*key2) = q->keytab2[0];

    q->elttab[0]  = q->elttab[q->used-1];
    q->keytab[0]  = q->keytab[q->used-1];
    q->keytab2[0] = q->keytab2[q->used-1];
    q->used--;
    
    i=1;
    
    while(i <= (q->used/2))
	{
	    if( (2*i == q->used)
		|| compWith2keys(q, 2*i-1, 2*i))     /*(q->keytab[2*i-1] < q->keytab[2*i]))*/
		{
		    j = 2*i;
		}
	    else
 		{
		    j = 2*i+1;
		}
	    if (!compWith2keys(q, i-1, j-1))         /*(q->keytab[i-1] >= q->keytab[j-1])*/
		{
		    swap_elt = q->elttab[i-1];
		    swap_key = q->keytab[i-1];
		    swap_key2 = q->keytab2[i-1];

		    q->elttab[i-1]  = q->elttab[j-1];
		    q->keytab[i-1]  = q->keytab[j-1];
		    q->keytab2[i-1] = q->keytab2[j-1];

		    q->elttab[j-1] = swap_elt;
		    q->keytab[j-1] = swap_key;
		    q->keytab2[j-1] = swap_key2;
		    
		    i=j;
		}
	    else
		break;
	}

    return return_elt;
}
/*
  Function: queuePosses

  Check if an element belongs to a queue.

  Parameters:
    q   - The queue.
    elt - The searched element.

  Returns:
    API_YES - if the element was found.
    API_NO  - if the element was not found.
*/
int queuePossess(Queue * q, 
		 PASTIX_INT     elt ){
  PASTIX_INT i;
  for (i = 0; i < q->used; i++)
    if (q->elttab[i] == elt)
      return API_YES;

  return API_NO;
}

/*
  Function: queueInit
  
  Print a queue entries to standar error output.

  Parameters:
    q - The queue.
*/
void queuePrint(Queue *q)
{
  PASTIX_INT i;
  fprintf(stderr, "Queue :");
  for (i = 0; i < q->used; i++)
    fprintf(stderr, "(%d %f %d) ", 
	    (int)q->elttab[i], 
	    (double)q->keytab[i], 
	    (int)q->keytab2[i] );
  fprintf(stderr, "\n");
}
