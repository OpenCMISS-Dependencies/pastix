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
#ifndef NOT_EMPTY_H
#define NOT_EMPTY_H
#include "redefine_functions.h"

#ifdef SOPALIN_LU
#ifndef CHOL_SOPALIN
#define CHOL_SOPALIN
#endif
#endif

#ifdef CHOL_SOPALIN
#  ifdef SOPALIN_LU
#    define NOT_EMPTY(filename)                                \
  void PASTIX_PREFIX_F(ge ## _not_empty_ ## filename) (void){  \
    return;                                                    \
  }
#  else /* not SOPALIN_LU */
#    define NOT_EMPTY(filename)                                \
  void PASTIX_PREFIX_F(po ## _not_empty_ ## filename) (void){  \
    return;                                                    \
  }
#  endif /* not SOPALIN_LU */
#else /* not CHOL_SOPALIN */
#  ifdef HERMITIAN
#    define NOT_EMPTY(filename)                                \
  void PASTIX_PREFIX_F(he ## _not_empty_ ## filename) (void){  \
    return;                                                    \
  }
#  else /* not HERMITIAN */
#    define NOT_EMPTY(filename)                                \
  void PASTIX_PREFIX_F(sy ## _not_empty_ ## filename) (void){  \
    return;                                                    \
  }
#  endif /* not HERMITIAN */
#endif /* not CHOL_SOPALIN */
#endif /* not NOT_EMPTY_H */
