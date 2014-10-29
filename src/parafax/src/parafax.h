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
#ifndef PARAFAX_H
#define PARAFAX_H

int pfax_vect_merge       (PASTIX_INT * v1, PASTIX_INT * v2);
int pfax_vect_build       (int com_ind, int * compo, CSCD * cscd_dgraph, PASTIX_INT ** collect);
int pfax_vect_send        (PASTIX_INT * vect, int dest, int comp_ind, MPI_Comm mpi_comm);
int pfax_receive_vect     (PASTIX_INT * vect, PASTIX_INT comp_ind, MPI_Comm mpi_comm);
int pfax_receive_merge    (PASTIX_INT * tree, PASTIX_INT ** collect, int cblk, int myrank, MPI_Comm mpi_comm);

int parafax               (CSCD * cscd_dgraph, int myrank, PASTIX_INT * owner, PASTIX_INT cblk, PASTIX_INT * tree);
int pfax_vect_compo_build (PASTIX_INT * beg_compo, PASTIX_INT * end_compo, PASTIX_INT * isLeaf, PASTIX_INT * supnode, PASTIX_INT cblk);
int pfax_isLeaf_build     (PASTIX_INT * isLeaf, PASTIX_INT * tree, PASTIX_INT cblk);

#endif /* PARAFAX_H */
