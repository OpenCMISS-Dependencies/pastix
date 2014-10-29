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
/**   NAME       : ps.h                                    **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                write postscript with draws             **/
/**                of symbol matrix and elimination tree   **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     09 sep 1998     **/
/**                                                        **/
/************************************************************/

#ifndef PS_H
#define static
#endif

FILE    *       ps_open(char *);
void            ps_close(FILE *);
void            ps_write_matrix(SymbolMatrix *, FILE *, PASTIX_INT *);
void            ps_write_tree(const CostMatrix *, const EliminTree *, FILE *, PASTIX_INT *);
static double    ps_rec_write_tree(PASTIX_INT , const CostMatrix *, const EliminTree *, FILE *,
			    void (*ps_draw_node)(FILE *,PASTIX_INT ,const CostMatrix *, const EliminTree *,double ,double ,double )
				  );
static void     ps_draw_node_num(FILE *, PASTIX_INT , const CostMatrix *,const EliminTree *, double , double , double );
void            ps_write_tree_owner(PASTIX_INT *,const CostMatrix *, const EliminTree *, FILE *, PASTIX_INT *);
static double    ps_rec_write_tree_owner(PASTIX_INT , PASTIX_INT *, const CostMatrix *, const EliminTree *, FILE *,
			    void (*ps_draw_node)(FILE *,PASTIX_INT, PASTIX_INT, const CostMatrix *, const EliminTree *,double ,double ,double )
				  );
static void     ps_draw_node_owner(FILE *, PASTIX_INT, PASTIX_INT, const CostMatrix *,const EliminTree *, double , double , double );




#undef static
