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
  File: nompi.h

  Header redefining all MPI keywords in order to allow compilation without MPI.

  Used when compiling with -DFORCE_NOMPI

  Authors:
    Mathieu Faverge - faverge@labri.fr
    Pierre Ramet     - ramet@labri.fr
*/
#ifndef MPI_FALSE_H
#define MPI_FALSE_H

#define MPI_Datatype int
#define MPI_Op int
#define MPI_Request int
#define MPI_Aint int
#define MPI_Comm int
#define MPI_Fint int

#define MPI_CHAR           1
#define MPI_UNSIGNED_CHAR  2
#define MPI_BYTE           3
#define MPI_SHORT          4
#define MPI_UNSIGNED_SHORT 5
#define MPI_INT            6
#define MPI_UNSIGNED       7
#define MPI_LONG           8
#define MPI_UNSIGNED_LONG  9
#define MPI_FLOAT          10
#define MPI_DOUBLE         11
#define MPI_LONG_DOUBLE    12
#define MPI_LONG_LONG_INT  13
#define MPI_LONG_LONG      13
#define MPI_PACKED         14
#define MPI_LB             15
#define MPI_UB             16
#define MPI_FLOAT_INT         17
#define MPI_DOUBLE_INT        18
#define MPI_LONG_INT          19
#define MPI_SHORT_INT         20
#define MPI_2INT              21
#define MPI_LONG_DOUBLE_INT   22
#define MPI_COMPLEX           23
#define MPI_DOUBLE_COMPLEX    24
#define MPI_LOGICAL           25
#define MPI_REAL              26
#define MPI_DOUBLE_PRECISION  27
#define MPI_INTEGER           28
#define MPI_2INTEGER          29
#define MPI_2COMPLEX          30
#define MPI_2DOUBLE_COMPLEX   31
#define MPI_2REAL             32
#define MPI_2DOUBLE_PRECISION 33
#define MPI_INTEGER4          34
#define MPI_INTEGER8          35
#define MPI_CHARACTER         1

#define MPI_COMM_WORLD 0
#define MPI_COMM_SELF  1
#define MPI_ANY_SOURCE 2
#define MPI_ANY_TAG 3
#define MPI_REQUEST_NULL 4
#define MPI_UNDEFINED 5
#define MPI_SUM 8
#define MPI_MAX 9
#define MPI_BOTTOM NULL
#define MPI_THREAD_MULTIPLE 11
#define MPI_THREAD_SINGLE 12
#define MPI_THREAD_FUNNELED 13
#define MPI_THREAD_SERIALIZED 14

#define MPI_MAX_PROCESSOR_NAME 1

typedef void (MPI_User_function) ( void *, void *, int *, MPI_Datatype * );

typedef struct MPI_Status{
  int MPI_SOURCE;
  int MPI_TAG;
  int MPI_ERROR;
} MPI_Status;

#define MPI_Get_processor_name(__name, __len) { (void)(__name); *(__len) = 0; }
#define MPI_Send(buf, count, datatype, dest, tag, comm)
#define MPI_Recv(buf, count, datatype, source, tag, comm, status)
#define MPI_Irecv(buf, count, datatype, source, tag, comm, request);
#define MPI_Isend(buf, count, datatype, dest, tag, comm, request)
#define MPI_Wait(request, status)
#define MPI_Waitany(count, array_of_requests, index, status);
#define MPI_Cancel(request)
#define MPI_Test(request, flag, status)
#define MPI_Testany(count, array_of_requests, index, flag, array_of_statuses)
#define MPI_Iprobe(source, tag, comm, flag, status)
#define MPI_Recv_init(buf, count, datatype, dest, tag, comm, request)
#define MPI_Start(request)
#define MPI_Startall(count, array_of_requests)
#define MPI_Type_contiguous(count, oldtype, newtype)
#define MPI_Type_struct(count, array_of_blocklengths, array_of_displacement, \
                        oldtype, newtype)
#define MPI_Address(location, newtype)
#define MPI_Type_commit(datatype)
#define MPI_Type_free(datatype)
#define MPI_Request_free(request)
#define MPI_Barrier(comm)
#define MPI_Op_create(function, commute, op)
#define MPI_Init(argc, argv)
#define MPI_Finalize()
#define MPI_Comm_split(comm, color, id, new_comm)

#define MPI_Allgather(sendbuf, sendcount, sendtype, recvbuf,            \
                      recvcount, recvtype, comm)    recvbuf = sendbuf;

#define MPI_Get_count(status, datatype ,count) *(count) = 0;
#define MPI_Comm_size(comm, size) *(size)=1;
#define MPI_Comm_rank(comm, rank) *(rank)=0;

#define MPI_Gather(sendbuf, sendcount, sendtype, recvbuf,   \
                   recvcount, recvtype, root, comm)         \
  do {                                                      \
    switch (sendtype) {                                     \
    case MPI_INT:                                           \
      memcpy(recvbuf,sendbuf,sendcount*sizeof(int));        \
      break;                                                \
    case MPI_LONG:                                          \
      memcpy(recvbuf,sendbuf,sendcount*sizeof(long));       \
      break;                                                \
    case MPI_INTEGER4:                                      \
      memcpy(recvbuf,sendbuf,sendcount*sizeof(int32_t));    \
        break;                                              \
    case MPI_INTEGER8:                                      \
      memcpy(recvbuf,sendbuf,sendcount*sizeof(int64_t));    \
      break;                                                \
    case MPI_FLOAT:                                         \
      memcpy(recvbuf,sendbuf,sendcount*sizeof(float));      \
      break;                                                \
    case MPI_DOUBLE:                                        \
      memcpy(recvbuf,sendbuf,sendcount*sizeof(double));     \
      break;                                                \
    case MPI_COMPLEX:                                       \
      memcpy(recvbuf,sendbuf,sendcount*2*sizeof(float));    \
      break;                                                \
    case MPI_DOUBLE_COMPLEX:                                \
      memcpy(recvbuf,sendbuf,sendcount*2*sizeof(double));   \
      break;                                                \
    default:                                                \
      fprintf(stderr,"error in #define MPI_Gather\n");      \
    }                                                       \
  } while(0)

#define MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount,  \
                     recvtype, comm)                                    \
  MPI_Allreduce(sendbuf, recvbuf, sendcount, sendtype, 1, comm)

#define MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm)  \
  switch (datatype) {                                               \
  case MPI_INT:                                                     \
    memcpy(recvbuf,sendbuf,count*sizeof(int));                      \
    break;                                                          \
  case MPI_LONG:                                                    \
    memcpy(recvbuf,sendbuf,count*sizeof(long));                     \
    break;                                                          \
  case MPI_INTEGER4:                                                \
    memcpy(recvbuf,sendbuf,count*sizeof(int32_t));                  \
    break;                                                          \
  case MPI_INTEGER8:                                                \
    memcpy(recvbuf,sendbuf,count*sizeof(int64_t));                  \
    break;                                                          \
  case MPI_FLOAT:                                                   \
    memcpy(recvbuf,sendbuf,count*sizeof(float));                    \
    break;                                                          \
  case MPI_DOUBLE:                                                  \
    memcpy(recvbuf,sendbuf,count*sizeof(double));                   \
    break;                                                          \
  case MPI_COMPLEX:                                                 \
    memcpy(recvbuf,sendbuf,count*2*sizeof(float));                  \
    break;                                                          \
  case MPI_DOUBLE_COMPLEX:                                          \
    memcpy(recvbuf,sendbuf,count*2*sizeof(double));                 \
    break;                                                          \
  default:                                                          \
    fprintf(stderr,"error in #define MPI_Allreduce\n");             \
  }

#define MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm) \
  switch (datatype) {                                                 \
  case MPI_INT:                                                       \
    memcpy(recvbuf,sendbuf,count*sizeof(int));                        \
    break;                                                            \
  case MPI_LONG:                                                      \
    memcpy(recvbuf,sendbuf,count*sizeof(long));                       \
    break;                                                            \
  case MPI_FLOAT:                                                     \
    memcpy(recvbuf,sendbuf,count*sizeof(float));                      \
    break;                                                            \
  case MPI_INTEGER4:                                                  \
    memcpy(recvbuf,sendbuf,count*sizeof(int32_t));                    \
    break;                                                            \
  case MPI_INTEGER8:                                                  \
    memcpy(recvbuf,sendbuf,count*sizeof(int64_t));                    \
    break;                                                            \
  case MPI_DOUBLE:                                                    \
    memcpy(recvbuf,sendbuf,count*sizeof(double));                     \
    break;                                                            \
  case MPI_COMPLEX:                                                   \
    memcpy(recvbuf,sendbuf,count*2*sizeof(float));                    \
    break;                                                            \
  case MPI_DOUBLE_COMPLEX:                                            \
    memcpy(recvbuf,sendbuf,count*2*sizeof(double));                   \
    break;                                                            \
  default:                                                            \
    fprintf(stderr,"error in #define MPI_Reduce\n");                  \
  }


#define MPI_Bcast(buffer, count, datatype, root, comm)

#define MPI_Comm_f2c(comm) 0;
#define MPI_Init_thread(argc, argv,required,provided)

#endif /* MPI_FALSE_H */
