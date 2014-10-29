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
#ifndef UPDO_SENDRECV_H
#define UPDO_SENDRECV_H

#define updo_up_WaitCtrb_storage API_CALL(updo_up_WaitCtrb_storage)
void updo_up_WaitCtrb_storage  ( Sopalin_Data_t *sopalin_data,
                                 PASTIX_INT             updo_buffer_size,
                                 void           *updo_buffer,
                                 PASTIX_INT             me,
                                 PASTIX_INT             i);

#  define send_waitall      API_CALL(send_waitall)
#  define send_waitall_down API_CALL(send_waitall_down)
#  define send_waitall_up   API_CALL(send_waitall_up)
#  define send_waitone_down API_CALL(send_waitone_down)
#  define send_waitone_up   API_CALL(send_waitone_up)
#  define send_free_down    API_CALL(send_free_down)
#  define send_free_up      API_CALL(send_free_up)
#  define send_testall_down API_CALL(send_testall_down)
#  define send_testall_up   API_CALL(send_testall_up)

void send_free_down    ( Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT s_index );
void send_free_up      ( Sopalin_Data_t *sopalin_data, PASTIX_INT me, PASTIX_INT s_index );
void send_testall_down ( Sopalin_Data_t *sopalin_data, PASTIX_INT me );
void send_testall_up   ( Sopalin_Data_t *sopalin_data, PASTIX_INT me );
void send_waitall      ( Sopalin_Data_t *sopalin_data, PASTIX_INT me,
                         void (*funcfree)(Sopalin_Data_t*, PASTIX_INT, PASTIX_INT) );
void send_waitall_down ( Sopalin_Data_t *sopalin_data, PASTIX_INT me );
void send_waitall_up   ( Sopalin_Data_t *sopalin_data, PASTIX_INT me );
int  send_waitone_down ( Sopalin_Data_t *sopalin_data, PASTIX_INT me );
int  send_waitone_up   ( Sopalin_Data_t *sopalin_data, PASTIX_INT me );

#ifdef PASTIX_UPDO_ISEND

#  define test_all_downsend API_CALL(test_all_downsend)
#  define test_all_upsend   API_CALL(test_all_upsend)
#  define wait_all_downsend API_CALL(wait_all_downsend)
#  define wait_all_upsend   API_CALL(wait_all_upsend)

void test_all_downsend ( Sopalin_Data_t *sopalin_data, PASTIX_INT me, int tag );
void test_all_upsend   ( Sopalin_Data_t *sopalin_data, PASTIX_INT me, int tag );
void wait_all_downsend ( Sopalin_Data_t *sopalin_data, PASTIX_INT me );
void wait_all_upsend   ( Sopalin_Data_t *sopalin_data, PASTIX_INT me );

#endif

/* Réceptions des contributions pour la remontée */
#define probe_updown               API_CALL(probe_updown)
#define updo_up_WaitCtrb_nostorage API_CALL(updo_up_WaitCtrb_nostorage)
#define updo_down_send             API_CALL(updo_down_send)
#define updo_down_recv             API_CALL(updo_down_recv)
#define updo_up_send               API_CALL(updo_up_send)
#define updo_up_recv               API_CALL(updo_up_recv)
#ifndef FORCE_NOMPI

/* /\* Thread de communication *\/ */
/* void* API_CALL(updo_thread_comm)(void * arg){ return NULL; } */
int  probe_updown ( MPI_Comm, PASTIX_INT );

#  ifndef STORAGE
void updo_up_WaitCtrb_nostorage ( Sopalin_Data_t *sopalin_data,
                                  PASTIX_INT, void *, PASTIX_INT, PASTIX_INT );
#  endif

/* Réceptions des comms MPI */
static inline
void updo_down_send ( Sopalin_Data_t *sopalin_data, PASTIX_INT, PASTIX_INT, PASTIX_INT );
static inline
void updo_down_recv ( Sopalin_Data_t *sopalin_data, void *,
                      MPI_Status, PASTIX_INT);
static inline
void updo_up_send   ( Sopalin_Data_t *sopalin_data, PASTIX_INT, PASTIX_INT, PASTIX_INT);
static inline
void updo_up_recv   ( Sopalin_Data_t *sopalin_data, void *,
                      MPI_Status, PASTIX_INT);

/* Thread de communication */
#  define updo_thread_comm API_CALL(updo_thread_comm)
void* updo_thread_comm ( void * );


#else /* FORCE_NOMPI */
void updo_up_WaitCtrb_nostorage ( Sopalin_Data_t *sopalin_data,
                                  PASTIX_INT buf_size, void *buf,
                                  PASTIX_INT me, PASTIX_INT i)
{
  (void)sopalin_data; (void)buf_size; (void)buf; (void)me; (void)i;
}
int  probe_updown ( MPI_Comm comm, PASTIX_INT tag )
{
  (void)comm; (void)tag;
  return 0;
}

#endif /* FORCE_NOMPI */

#endif /* not UPDO_SENDRECV_H */
