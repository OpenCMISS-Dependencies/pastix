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
/*---------------------------------------------------------------*/
/* This program times blocking send/receives, and reports the    */
/* latency and bandwidth of the communication system.  It is     */
/* designed to run with an even number of mpi tasks.             */
/*---------------------------------------------------------------*/
#include "commun.h"
#include <mpi.h>

#define MAXPTS 1250000
#define NCOUNTS 12

int main(int argc, char * argv[])
{
   double *sbuf, *rbuf;
   int iter, maxiter, repeats[NCOUNTS];
   int count[NCOUNTS];
   int nc, nbytes;
   int taskid, ntasks;
   int itag = 99;
   double etime;
   double latency, bw;

   sbuf = (double*) malloc(MAXPTS*sizeof(double));
   rbuf = (double*) malloc(MAXPTS*sizeof(double));

   MPI_Status mpi_status[2];
   MPI_Request mpi_request[2];

   /*----------------------------------------------*/
   /* define an array of counts for 8-unsigned char objects */
   /*----------------------------------------------*/
   count[0] = 0;
   count[1] = 1;
   count[2] = 4;
   count[3] = 12;
   count[4] = 40;
   count[5] = 125;
   count[6] = 400;
   count[7] = 1250;
   count[8] = 4000;
   count[9] = 12500;
   count[10] = 40000;
   count[11] = 125000;

   repeats[0] = 100;
   repeats[1] = 100;
   repeats[2] = 100;
   repeats[3] = 100;
   repeats[4] = 100;
   repeats[5] = 100;
   repeats[6] = 100;
   repeats[7] = 100;
   repeats[8] = 100;
   repeats[9] = 50;
   repeats[10] = 25;
   repeats[11] = 10;

   /*-----------------------------------------------------------*/
   /* set-up the parallel environment: assign ntasks and taskid */
   /*-----------------------------------------------------------*/
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
   MPI_Comm_rank(MPI_COMM_WORLD, &taskid);
   if ((ntasks % 2) != 0 && taskid == 0)
   {
     fprintf(stdout,"You must specify an even number of MPI tasks.\n");
     exit(-1);
   }

   /*----------------------------------------------------*/
   /* send/recv to ensure that the routines are loaded   */
   /*----------------------------------------------------*/
   nc = 1;
   if ((taskid % 2) == 0)
      MPI_Isend(sbuf, count[nc], MPI_DOUBLE, taskid+1, 
                itag, MPI_COMM_WORLD, &(mpi_request[0]));
   else
      MPI_Irecv(rbuf, count[nc], MPI_DOUBLE, taskid-1, 
                itag, MPI_COMM_WORLD, &(mpi_request[0]));

   if ((taskid % 2) == 1)
      MPI_Isend(sbuf, count[nc], MPI_DOUBLE, taskid-1, 
                itag, MPI_COMM_WORLD, &(mpi_request[1]));
   else
      MPI_Irecv(rbuf, count[nc], MPI_DOUBLE, taskid+1, 
                itag, MPI_COMM_WORLD, &(mpi_request[1]));
   MPI_Waitall(2,mpi_request,mpi_status);

   /*--------------------------------------------------------*/
   /* send or receive messages, and measure round-trip time. */
   /* even tasks send, odd tasks receive, then the reverse.  */
   /*--------------------------------------------------------*/
   for (nc=0; nc<NCOUNTS; nc++)
   {

      MPI_Barrier(MPI_COMM_WORLD); /* synchronize here */

      TEST_CLOCK_INIT
      maxiter = repeats[nc];
      for (iter=0; iter<maxiter; iter++)
      {
         /*--------------------------------------------*/
         /* send in one direction i->i+1               */
         /*--------------------------------------------*/
         if ((taskid % 2) == 0)
             MPI_Isend(sbuf, count[nc], MPI_DOUBLE, taskid+1, 
                       itag, MPI_COMM_WORLD, &(mpi_request[0]));
         else
             MPI_Irecv(rbuf, count[nc], MPI_DOUBLE, taskid-1, 
                       itag, MPI_COMM_WORLD, &(mpi_request[0]));

         /*--------------------------------------------*/
         /* send in the reverse direction i+1->i       */
         /*--------------------------------------------*/
         if ((taskid % 2) == 1)
             MPI_Isend(sbuf, count[nc], MPI_DOUBLE, taskid-1, 
                       itag, MPI_COMM_WORLD, &(mpi_request[1]));
         else
             MPI_Irecv(rbuf, count[nc], MPI_DOUBLE, taskid+1, 
                       itag, MPI_COMM_WORLD, &(mpi_request[1]));
         MPI_Waitall(2,mpi_request,mpi_status);

      }  /* end the repeat loop */
      TEST_CLOCK_STOP

      /*-----------------------------------------*/
      /* write timing data for each message size */
      /*-----------------------------------------*/
      nbytes = 8*count[nc]; /* 8 bytes per entry */
      etime = 0.5e3*(TEST_CLOCK_GET)/maxiter;
      if (taskid == 0)
      {
        fprintf(stdout,"msglen = %8d bytes,   elapsed time = %.4lf msec\n", 
                nbytes, etime);
      }
      if (nc == 0) latency = 1.0e3*etime;
      if (nc == (NCOUNTS-1))  bw = nbytes/(1.0e3*etime);

   }  /* end the loop over message sizes */

   /*--------------------------------------------------------*/
   /*report apporximate numbers for bandwidth and latency    */
   /*--------------------------------------------------------*/
   if (taskid == 0)
   {
     fprintf(stdout,"\nlatency = %.1lf microseconds\n", latency);
     fprintf(stdout,"bandwidth =  %.2lf MBytes/sec\n", bw);
     fprintf(stdout,"(approximate values for MPI_Isend/MPI_Irecv)\n");
   }

   MPI_Finalize();

   free(sbuf);
   free(rbuf);

   return(0);

}
