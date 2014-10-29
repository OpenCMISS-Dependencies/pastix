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
#ifndef OOC_H
#define OOC_H

/* Return values */
#define EXIT_FAILURE_CBLK_NOT_NULL      3
#define EXIT_FAILURE_SAVING_NULL_BUFFER 4
#define EXIT_FAILURE_OUT_OF_MEMORY      5
#define EXIT_FAILURE_FILE_OPENING       6
#define EXIT_FAILURE_FILE_TRUNCATED     7
#define EXIT_SUCCESS_HACK               8
#define EXIT_SUCCESS_ALL_LOADED         9
#define EXIT_FAILURE_CBLK_USED          10

/* OOC Step */
#define OOCSTEP_COEFINIT                1
#define OOCSTEP_SOPALIN                 2
#define OOCSTEP_DOWN                    3
#define OOCSTEP_DIAG                    4
#define OOCSTEP_UP                      5


#ifdef OOC

#define OOC_RECEIVING ooc_receiving(sopalin_data)
#define OOC_RECEIVED ooc_received(sopalin_data)
#define OOC_THREAD_NBR sopalin_data->sopar->iparm[IPARM_OOC_THREAD]

/* sets values for the global ooc structure 
 * sopalin_data : Sopalin_Data_t global structure
 * limit        : memory limit set by use
 */
void *ooc_thread(void * arg);

/* Init / Clean */
int ooc_init          (Sopalin_Data_t * sopalin_data, PASTIX_INT limit);
int ooc_exit          (Sopalin_Data_t * sopalin_data);

/* Step */
int ooc_stop_thread   (Sopalin_Data_t * sopalin_data);
int ooc_freeze        (Sopalin_Data_t * sopalin_data);
int ooc_defreeze      (Sopalin_Data_t * sopalin_data);
int ooc_set_step      (Sopalin_Data_t * sopalin_data, int step);

/* Cblk */
int ooc_wait_for_cblk (Sopalin_Data_t * sopalin_data, PASTIX_INT cblk, int me);
int ooc_hack_load     (Sopalin_Data_t * sopalin_data, PASTIX_INT cblk, int me);
int ooc_save_coef     (Sopalin_Data_t * sopalin_data, PASTIX_INT task, PASTIX_INT cblk, int me);

void ooc_receiving    (Sopalin_Data_t * sopalin_data);
void ooc_received     (Sopalin_Data_t * sopalin_data);
void ooc_wait_task    (Sopalin_Data_t * sopalin_data, PASTIX_INT task, int me);
#else /* OOC */

#define OOC_RECEIVING 
#define OOC_RECEIVED 
#define OOC_THREAD_NBR 0
#define ooc_thread     NULL

#define ooc_init(sopalin_data, limit)
#define ooc_exit(sopalin_data)

#define ooc_stop_thread(sopalin_data)
#define ooc_freeze(sopalin_data)
#define ooc_defreeze(sopalin_data)
#define ooc_set_step(sopalin_data, step)

#define ooc_wait_for_cblk(sopalin_data, cblk, me)
#define ooc_save_coef(sopalin_data, task, cblk, me)
#define ooc_hack_load(sopalin_data, cblknum, me)

#define ooc_wait_task(sopalin_data, task, me) 

#endif /* OOC */

#ifdef OOC_FTGT
/* Ftgt */
int ooc_wait_for_ftgt (Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum, int me);
int ooc_reset_ftgt    (Sopalin_Data_t * sopalin_data, PASTIX_INT ftgtnum, int me);
int ooc_save_ftgt     (Sopalin_Data_t * sopalin_data, PASTIX_INT tasknum, PASTIX_INT ftgtnum, int me);

#else

#define ooc_wait_for_ftgt(sopalin_data, ftgtnum, me)
#define ooc_reset_ftgt(sopalin_data, ftgtnum, me)
#define ooc_save_ftgt(sopalin_data, tasknum, ftgtnum, me)
#endif

#endif /* OOC_H */
