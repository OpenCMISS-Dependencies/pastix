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
/************************************************************/
/**                                                        **/
/**   NAME       : queue.h                                 **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : queue of PASTIX_INT that sorts elements        **/
/**                in ascending way according to a         **/
/**                PASTIX_FLOAT key                               **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

#ifndef QUEUE_H
#define QUEUE_H

/*
**  The type and structure definitions.
*/

typedef struct Queue_ {
  PASTIX_INT        size;                  /*+ Allocated memory size             +*/ 
  PASTIX_INT        used;                  /*+ Number of element in the queue    +*/
  PASTIX_INT    *   elttab;                /*+ Array of the element              +*/
  double *   keytab;                /*+ Array of keys                     +*/
  PASTIX_INT    *   keytab2;               /*+ Another array of keys             +*/
} Queue;


#define static

int     queueInit       (Queue *, PASTIX_INT size);
void    queueExit       (Queue *);
Queue * queueCopy       (Queue *dst, Queue *src);
void    queueAdd        (Queue *, PASTIX_INT, double);
void    queueAdd2       (Queue *, PASTIX_INT, double, PASTIX_INT);
PASTIX_INT     queueGet        (Queue *);
PASTIX_INT     queueSize       (Queue *);
void    queueClear      (Queue *);
PASTIX_INT     queueRead       (Queue *);
PASTIX_INT     queueGet2       (Queue *, double *, PASTIX_INT *);
int     queuePossess    (Queue *, PASTIX_INT);
void    queuePrint      (Queue *);

static PASTIX_INT compWith2keys(Queue *, PASTIX_INT, PASTIX_INT);
#undef static
#endif
