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
#ifndef STARPU_SUBMIT_TASKS_H
#define STARPU_SUBMIT_TASKS_H

#ifdef WITH_STARPU
#define starpu_loop_data_  API_CALL(starpu_loop_data_)
#define starpu_loop_data_t API_CALL(starpu_loop_data_t)
typedef struct starpu_loop_data_ starpu_loop_data_t;

#define starpu_submit_one_trf API_CALL(starpu_submit_one_trf)
int starpu_submit_one_trf (PASTIX_INT itertask, Sopalin_Data_t * sopalin_data);
#define starpu_submit_bunch_of_gemm API_CALL(starpu_submit_bunch_of_gemm)
int starpu_submit_bunch_of_gemm (PASTIX_INT itertask, Sopalin_Data_t * sopalin_data);

#define starpu_submit_tasks API_CALL(starpu_submit_tasks)
int starpu_submit_tasks(Sopalin_Data_t   * sopalin_data);
#endif

#endif /* STARPU_SUBMIT_TASKS_H */
