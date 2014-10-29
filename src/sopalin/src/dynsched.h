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
  File: dynsched.h

  Main code to select next task to execute in dynamic scheduler.

  Authors:
    Mathieu Faverge - faverge@labri.fr

  Date:
    Version 0.0 - february 2003
*/

#if !defined(PASTIX_DYNSCHED_WITH_TREE)
/*
 * Here the stealing algorithm travel through the tree thanks to a
 * breadth-first search algorithm starting from the actual node and
 * not from the root of the tree.
 */
static inline PASTIX_INT
API_CALL(sopalin_dynsched_getNexTask)(Sopalin_Data_t *sopalin_data,
                                      SolverMatrix   *datacode,
                                      Thread_Data_t  *thread_data,
                                      PASTIX_INT *itaskptr,
                                      PASTIX_INT *itaskptr2,
                                      PASTIX_INT *bloknum,
                                      PASTIX_INT me)
{
  PASTIX_INT itasktab = *itaskptr;
  int position;
  PASTIX_INT itasktab2;
  PASTIX_INT i, restart;

 debd:

  /* On remonte dans l'arbre si on a rien a faire a ce niveau */
  restart = 0;
  for( position=0; position <  datacode->thrdnbr; position++ ) {
    itasktab2 = thread_data->tabtravel[position];

    if ( !(sopalin_data->tasktab_indice[itasktab2] < sopalin_data->datacode->ttsknbr[itasktab2]) )
      continue;

    MUTEX_LOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
    i = queueGet2(&(sopalin_data->taskqueue[itasktab2]), NULL, bloknum);
    MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));

    if ( i != -1 ) {
      *itaskptr  = itasktab;
      *itaskptr2 = itasktab2;
      return i;
    }

    restart = 1;
  }

  /* Crapy way to detect the end */
  /* || TASK_CTRBCNT(SOLV_TASKNBR-1) != 0 */
  if ( ( restart ||
	 ( SOLV_TASKNBR > 0 && TASK_CTRBCNT(SOLV_TASKNBR-1) != 0) ) &&
       ( sopalin_data->step_comm == COMMSTEP_INIT  ||
	 sopalin_data->step_comm == COMMSTEP_FACTO ||
	 sopalin_data->step_comm == COMMSTEP_DOWN  ||
	 sopalin_data->step_comm == COMMSTEP_UP    ) ) {
    MUTEX_LOCK(&(sopalin_data->tasktab_mutex[me]));
    COND_TIMEWAIT(&(sopalin_data->tasktab_cond[me]), &(sopalin_data->tasktab_mutex[me]));
    MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[me]));
    goto debd;
  }

  return -1;
}

#else

/*
 * Standard version, the stealing algorithm go up into the tree to find work
 */
static inline PASTIX_INT
API_CALL(sopalin_dynsched_getNexTask)(Sopalin_Data_t *sopalin_data,
                                      SolverMatrix   *datacode,
                                      Thread_Data_t  *thread_data,
                                      PASTIX_INT *itaskptr,
                                      PASTIX_INT *itaskptr2,
                                      PASTIX_INT *bloknum,
                                      PASTIX_INT me)
{
  PASTIX_INT itasktab = *itaskptr;
  PASTIX_INT itasktab2;
  PASTIX_INT i;

 deb:
  /* If there is nothing to do at our level, we just climb the tree */
  itasktab = me;
  while (itasktab != -1 &&
         !(sopalin_data->tasktab_indice[itasktab] < datacode->ttsknbr[itasktab]) )
    {
      itasktab =  BFATHER(datacode->btree, itasktab);
    }

  /* There is nothing more to do, so we exit the loop */
  if ((itasktab == -1) &&
      !(sopalin_data->tasktab_indice[me] < datacode->ttsknbr[me]))
    return -1;
  else if ( sopalin_data->tasktab_indice[me] < datacode->ttsknbr[me] ) {
    itasktab2 = me;
  }
  else {
    itasktab2 = itasktab;
  }

  MUTEX_LOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
  while(queueSize(&(sopalin_data->taskqueue[itasktab2])) == 0)
    {
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
      itasktab2 = BFATHER(datacode->btree, itasktab2);
      if (itasktab2 == -1)
        {
          itasktab2 = itasktab;
          MUTEX_LOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
          break;
        }
      MUTEX_LOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
    }

  i = queueGet2(&(sopalin_data->taskqueue[itasktab2]), NULL, bloknum);
#ifndef	NOSTEAL_BROTHER
  if (i == -1)
    {
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
      itasktab2 = itasktab;
      while (itasktab2 != -1)
        {
          int j;
          for (j=datacode->btree->nodetab[itasktab2].fsonnum;
               j < datacode->btree->nodetab[itasktab2].fsonnum+datacode->btree->nodetab[itasktab2].sonsnbr;
               j++)
            {
              int bubnum = datacode->btree->sonstab[j];
              MUTEX_LOCK(&(sopalin_data->tasktab_mutex[bubnum]));
              i = queueGet2(&(sopalin_data->taskqueue[bubnum]), NULL, bloknum);
              if (i != -1)
                {
                  itasktab2 = bubnum;
                  MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[bubnum]));
                  goto fin;
                }
              else
                MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[bubnum]));
            }
          itasktab2 = BFATHER(datacode->btree, itasktab2);
    }

      MUTEX_LOCK(&(sopalin_data->tasktab_mutex[itasktab]));
      COND_TIMEWAIT(&(sopalin_data->tasktab_cond[itasktab]), &(sopalin_data->tasktab_mutex[itasktab]));
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab]));
      goto deb;
    }
  else
    {
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
    }
#else /* NOSTEAL_BROTHER */
  if (i == -1)
    {
      COND_TIMEWAIT(&(sopalin_data->tasktab_cond[itasktab2]), &(sopalin_data->tasktab_mutex[itasktab2]));
      MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
      goto deb;
    }
  MUTEX_UNLOCK(&(sopalin_data->tasktab_mutex[itasktab2]));
#endif

 fin:
  *itaskptr  = itasktab;
  *itaskptr2 = itasktab2;

  return i;
}

#endif  /* !defined(PASTIX_DYNSCHED_WITH_TREE) */

