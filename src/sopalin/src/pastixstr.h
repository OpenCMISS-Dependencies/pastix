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
#ifndef PASTIX_STR_H
#define PASTIX_STR_H

#ifndef SOPALIN_3D_H
#error "sopalin3d.h must be included before pastixstr.h"
#endif

#ifndef ORDER_H
#error "order.h must be included before pastixstr.h"
#endif
/*
   struct: pastix_data_t

   Structure used to store datas for a step by step execution.
*/

struct pastix_data_t {
  SolverMatrix     solvmatr;         /*+ Matrix informations                                                 +*/
  CscMatrix	   cscmtx;	     /*+ Compress Sparse Column matrix                                       +*/
  SymbolMatrix    *symbmtx;          /*+ Symbol Matrix                                                       +*/
  SopalinParam     sopar;            /*+ Sopalin parameters                                                  +*/
  Order            ordemesh;         /*+ Order                                                               +*/
#ifdef WITH_SCOTCH
  SCOTCH_Graph     grafmesh;         /*+ Graph                                                               +*/
  int              malgrf;           /*+ boolean indicating if grafmesh has been allocated                   +*/
#endif /* WITH_SCOTCH */
#ifdef DISTRIBUTED
#ifdef WITH_SCOTCH
  SCOTCH_Dordering ordedat;          /*+ distributed scotch order                                            +*/
  SCOTCH_Dgraph    dgraph;
  PASTIX_INT             *PTS_permtab;
  PASTIX_INT             *PTS_peritab;
#endif /* WITH_SCOTCH */
  PASTIX_INT             *glob2loc;         /*+ local column number of global column, or -(owner+1) is not local    +*/
  PASTIX_INT              ncol_int;         /*+ Number of local columns in internal CSCD                            +*/
  PASTIX_INT             *l2g_int;          /*+ Local to global column numbers in internal CSCD                     +*/
  int              malrhsd_int;      /*+ Indicates if internal distributed rhs has been allocated            +*/
  int              mal_l2g_int;
  PASTIX_FLOAT           *b_int;            /*+ Local part of the right-hand-side                                   +*/
  PASTIX_INT             *loc2glob2;        /*+ local2global column number                                          +*/
#endif /* DISTRIBUTED */
  PASTIX_INT              gN;               /*+ global column number                                                +*/
  PASTIX_INT              n;                /*+ local column number                                                 +*/
  PASTIX_INT             *iparm;            /*+ Vecteur de parametres entiers                                       +*/
  double          *dparm;            /*+ Vecteur de parametres floattant                                     +*/
  PASTIX_INT              n2;               /*+ Number of local columns                                             +*/
  PASTIX_INT             *col2;             /*+ column tabular for the CSC matrix                                   +*/
				     /*+ (index of first element of each col in row and values tabulars)     +*/
  PASTIX_INT             *row2;             /*+ tabular containing row number of each element of                    +*/
				     /*+  the CSC matrix, ordered by column.                                 +*/
  int              bmalcolrow;       /*+ boolean indicating if col2 ans row2 have been allocated             +*/
  int              malord;           /*+ boolean indicating if ordemesh has been allocated                   +*/
  int              malcsc;           /*+ boolean indicating if solvmatr->cscmtx has beek allocated           +*/
  int              malsmx;           /*+ boolean indicating if solvmatr->updovct.sm2xtab has been allocated  +*/
  int              malslv;           /*+ boolean indicating if solvmatr has been allocated                   +*/
  int              malcof;           /*+ boolean indicating if coeficients tabular(s) has(ve) been allocated +*/
  MPI_Comm         pastix_comm;      /*+ PaStiX MPI communicator                                             +*/
  MPI_Comm         intra_node_comm;  /*+ PaStiX intra node MPI communicator                                  +*/
  MPI_Comm         inter_node_comm;  /*+ PaStiX inter node MPI communicator                                  +*/
  int              procnbr;          /*+ Number of MPI tasks                                                 +*/
  int              procnum;          /*+ Local MPI rank                                                      +*/
  int              intra_node_procnbr; /*+ Number of MPI tasks in node_comm                                  +*/
  int              intra_node_procnum; /*+ Local MPI rank in node_comm                                       +*/
  int              inter_node_procnbr; /*+ Number of MPI tasks in node_comm                                  +*/
  int              inter_node_procnum; /*+ Local MPI rank in node_comm                                       +*/
  int             *bindtab;          /*+ Tabular giving for each thread a CPU to bind it too                 +*/
  PASTIX_INT              nschur;           /*+ Number of entries for the Schur complement.                         +*/
  PASTIX_INT             *listschur;        /*+ List of entries for the schur complement.                           +*/
  PASTIX_FLOAT           *schur_tab;
  PASTIX_INT              schur_tab_set;
  int              scaling;          /*+ Indicates if the matrix has been scaled                             +*/
  PASTIX_FLOAT           *scalerowtab;      /*+ Describes how the matrix has been scaled                            +*/
  PASTIX_FLOAT           *iscalerowtab;
  PASTIX_FLOAT           *scalecoltab;
  PASTIX_FLOAT           *iscalecoltab;
#ifdef WITH_SEM_BARRIER
  sem_t           *sem_barrier;      /*+ Semaphore used for AUTOSPLIT_COMM barrier                           +*/
#endif
  PASTIX_INT              pastix_id;        /*+ Id of the pastix instance (PID of first MPI task)                   +*/
  int                     cscInternFilled;
};

#endif /* PASTIX_STR_H */
