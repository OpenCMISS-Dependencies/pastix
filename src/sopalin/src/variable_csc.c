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
#define CHUNCKSIZE 10

#ifndef NAN
#  define NAN (0.0 / 0.0)
#endif
#define OUT_INTL PASTIX_INT
#define OUT_INTS PASTIX_INT
#define TAG_VCSC_COUNT    128
#define TAG_VCSC_COLPTR   129
#define TAG_VCSC_ROWS     130
#define TAG_VCSC_VALUES  131

struct variable_csc_s {
  INTS n;
  INTL nz;
  INTS dof;
  INTS * colsizes;
  INTS ** rows;
  COEF ** values;
};

static inline
INTS vcsc_init(variable_csc_t * vcsc, INTS n, INTL nz, INTS dof, INTS id) ALWAYS_INLINE;

static inline
INTS vcsc_add_node(variable_csc_t vcsc, INTS COL, INTS ROW, COEF * VALUE,
		   COEF (*op)(COEF , COEF), INTS id) ALWAYS_INLINE;

static inline
INTS vcsc_add(variable_csc_t vcsc, INTS COL, INTS ROW, COEF VALUE,
              COEF (*op)(COEF , COEF), INTS id) ALWAYS_INLINE;

static inline
INTS vcsc_to_cscd(variable_csc_t   vcsc,
                  MPI_Comm         comm,
                  OUT_INTS      *  n,
                  OUT_INTL      ** colptr_o,
                  OUT_INTS      ** rows_o,
                  COEF          ** values_o,
                  OUT_INTS      ** l2g_o,
                  OUT_INTS      ** g2l_o,
                  COEF (*op)(COEF , COEF),
                  INTS             dropnonlocal,
                  INTS             id) ALWAYS_INLINE;
static inline
INTS vcsc_destroy(variable_csc_t vcsc, INTS id) ALWAYS_INLINE;

static inline
INTS vcsc_init(variable_csc_t * vcsc, INTS n, INTL nz, INTS dof, INTS id) {
  MURGE_MEMALLOC(*vcsc, 1, struct variable_csc_s);
  (*vcsc)->n  = n;
  (*vcsc)->nz = nz;
  (*vcsc)->dof = dof;
  MURGE_MEMALLOC((*vcsc)->colsizes, (*vcsc)->n, INTS);
  memset((*vcsc)->colsizes, 0, (*vcsc)->n*sizeof(INTS));

  MURGE_MEMALLOC((*vcsc)->rows, (*vcsc)->n, INTS*);
  memset((*vcsc)->rows, 0, (*vcsc)->n*sizeof(INTS*));
  if (dof > 0) {
    MURGE_MEMALLOC((*vcsc)->values, (*vcsc)->n, COEF*);
    memset((*vcsc)->values, 0, (*vcsc)->n*sizeof(COEF*));
  }
  return MURGE_SUCCESS;
}

static inline
INTS vcsc_add_node(variable_csc_t vcsc, INTS COL, INTS ROW, COEF * VALUE,
                   COEF (*op)(COEF , COEF), INTS id) {
  INTS i;
  INTS dof2 = vcsc->dof*vcsc->dof;

  if (vcsc->colsizes[COL-1] == 0) {
    MURGE_MEMALLOC(vcsc->rows[COL-1], CHUNCKSIZE, INTS);
    if (vcsc->dof > 0)
      MURGE_MEMALLOC(vcsc->values[COL-1],
		     dof2*(CHUNCKSIZE),
		     COEF);

    vcsc->colsizes[COL-1]++;
    vcsc->rows[COL-1][0] = ROW;
    if (vcsc->dof)
      memcpy(vcsc->values[(COL-1)], VALUE, dof2*sizeof(COEF));

  } else {
    for (i = 0; i < vcsc->colsizes[COL-1]; i++) {
      if (vcsc->rows[COL-1][i] > ROW)
        break;
      if (vcsc->rows[COL-1][i] == ROW) {
        if (vcsc->dof) {
          INTS ii;
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
          if ((COL-1)*vcsc->dof <= MURGE_TRACE_COL-1 &&
              (COL)*vcsc->dof > MURGE_TRACE_COL-1 &&
              (ROW-1)*vcsc->dof <= MURGE_TRACE_ROW-1 &&
              ROW*vcsc->dof > MURGE_TRACE_ROW-1) {

            fprintf(stdout, "before vcsc[%d-1][%d*%d*%d + %d *%d + %d] = %g (vcsc.rows[index] %d)\n",
                    COL, i, vcsc->dof, vcsc->dof, (MURGE_TRACE_COL-1-(COL-1)*vcsc->dof), vcsc->dof, 
                    (MURGE_TRACE_ROW-1 - (ROW-1)*vcsc->dof),
                    vcsc->values[COL-1][i*vcsc->dof*vcsc->dof + 
                                        (MURGE_TRACE_COL-1-(COL-1)*vcsc->dof)*vcsc->dof +
                                        (MURGE_TRACE_ROW-1 - (ROW-1)*vcsc->dof)],
                    vcsc->rows[COL-1][i]);
          }
#endif

          for (ii = 0; ii < dof2; ii++) {
            if (isnan(vcsc->values[COL-1][i*dof2 + ii])) {
              vcsc->values[COL-1][i*dof2 + ii] = VALUE[ii];
            } else {
              vcsc->values[COL-1][i*dof2 + ii] = op(vcsc->values[COL-1][i*dof2 + ii],
                                                    VALUE[ii]);
            }
          }

#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
          if ((COL-1)*vcsc->dof <= MURGE_TRACE_COL-1 &&
              (COL)*vcsc->dof > MURGE_TRACE_COL-1 &&
              (ROW-1)*vcsc->dof <= MURGE_TRACE_ROW-1 &&
              ROW*vcsc->dof > MURGE_TRACE_ROW-1) {
            fprintf(stdout, "%s:%d after vcsc[%d-1][%d*%d*%d + %d *%d + %d] = %g"
		    " (vcsc.rows[index] %d)\n", __FILE__, __LINE__,
                    COL, i, vcsc->dof, vcsc->dof,
                    (MURGE_TRACE_COL-1-(COL-1)*vcsc->dof), vcsc->dof,
                    (MURGE_TRACE_ROW-1 - (ROW-1)*vcsc->dof),
                    vcsc->values[COL-1][i*vcsc->dof*vcsc->dof +
                                        (MURGE_TRACE_COL-1-(COL-1)*vcsc->dof)*vcsc->dof +
                                        (MURGE_TRACE_ROW-1 - (ROW-1)*vcsc->dof)],
                    vcsc->rows[COL-1][i]);
          }
#endif

        }
        return MURGE_SUCCESS;
      }
    }

    {
      /* not found */
      INTS ii;
      if (vcsc->colsizes[COL-1]%CHUNCKSIZE == 0) {
        /* no more room in the column add CHUNKSIZE elements */
        MURGE_REALLOC(vcsc->rows[COL-1], vcsc->colsizes[COL-1]+CHUNCKSIZE, INTS);
        if (vcsc->dof)
          MURGE_REALLOC(vcsc->values[COL-1],
                        dof2*(vcsc->colsizes[COL-1]+CHUNCKSIZE),
                        COEF);
      }

      for (ii = vcsc->colsizes[COL-1]-1; ii >= i; ii--) {
        vcsc->rows[COL-1][ii+1] = vcsc->rows[COL-1][ii];
      }

      if (vcsc->dof)
        for (ii = vcsc->colsizes[COL-1]-1; ii >= i; ii--)
          memcpy(&(vcsc->values[COL-1][dof2*(ii+1)]),
                 &(vcsc->values[COL-1][dof2*ii]),
                 dof2*sizeof(COEF));

      vcsc->colsizes[COL-1]++;
      vcsc->rows[COL-1][i] = ROW;
      if (vcsc->dof) {
        memcpy(&(vcsc->values[COL-1][dof2*i]),
               VALUE,
               dof2*sizeof(COEF));
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
        if ((COL-1)*vcsc->dof <= MURGE_TRACE_COL-1 &&
            (COL)*vcsc->dof > MURGE_TRACE_COL-1 &&
            (ROW-1)*vcsc->dof <= MURGE_TRACE_ROW-1 &&
            ROW*vcsc->dof > MURGE_TRACE_ROW-1) {

          fprintf(stdout, "%s:%d after vcsc[%d-1][%d*%d*%d + %d *%d + %d] = %g (vcsc.rows[index] %d)\n",
		  __FILE__, __LINE__,
                  COL, i, vcsc->dof, vcsc->dof, (MURGE_TRACE_COL-1-(COL-1)*vcsc->dof), vcsc->dof,
                  (MURGE_TRACE_ROW-1 - (ROW-1)*vcsc->dof),
                  vcsc->values[COL-1][i*vcsc->dof*vcsc->dof +
                                      (MURGE_TRACE_COL-1-(COL-1)*vcsc->dof)*vcsc->dof +
                                      (MURGE_TRACE_ROW-1 - (ROW-1)*vcsc->dof)],
                  vcsc->rows[COL-1][i]);
        }
#endif

      }

    }

  } 
  return MURGE_SUCCESS;
}

static inline
INTS vcsc_add(variable_csc_t vcsc, INTS COL, INTS ROW, COEF VALUE,
              COEF (*op)(COEF , COEF), INTS id) {
  INTS i;
  INTS dof2 = vcsc->dof*vcsc->dof;

  if (vcsc->dof < 2) {
    return vcsc_add_node(vcsc, COL, ROW, &VALUE, op, id);
  } else {
    INTS COL_NODE = (COL - 1)/vcsc->dof;
    INTS ROW_NODE = (ROW - 1)/vcsc->dof;
    INTS COL_INNODE = (COL-1)%vcsc->dof;
    INTS ROW_INNODE = (ROW-1)%vcsc->dof;
    COEF * V;
    if (!(vcsc->colsizes[COL_NODE] == 0)) {
      /* Column is not empty */
      for (i = 0; i < vcsc->colsizes[COL_NODE]; i++) {
        if (vcsc->rows[COL_NODE][i] > ROW_NODE+1)
          /* row not found */
          break;
        if (vcsc->rows[COL_NODE][i] == ROW_NODE+1) {
          /* row found */
          INTS idx = i*dof2 +
            COL_INNODE*vcsc->dof +
            ROW_INNODE;
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
          if (COL == MURGE_TRACE_COL &&
              ROW == MURGE_TRACE_ROW ) {
            fprintf(stdout, "before vcsc[%d][%d] = %g (vcsc.rows[%d] %d)\n",
                    COL_NODE, idx,
                    vcsc->values[COL_NODE][idx], i,
                    vcsc->rows[COL_NODE][i]);
          }
#endif
          if (isnan(vcsc->values[COL_NODE][idx])) {
            vcsc->values[COL_NODE][idx] = VALUE;
          } else {
            vcsc->values[COL_NODE][idx] = op(vcsc->values[COL_NODE][idx],
                                             VALUE);
          }
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
          if (COL == MURGE_TRACE_COL &&
              ROW == MURGE_TRACE_ROW ) {
            fprintf(stdout, "after vcsc[%d][%d] = %g (vcsc.rows[index] %d)\n",
                    COL_NODE, idx,
                    vcsc->values[COL_NODE][idx],
                    vcsc->rows[COL_NODE][i]);
          } else {
            PASTIX_INT MY_COL, MY_ROW, index;
            INTS MY_COL_INNODE, MY_ROW_INNODE;;
            MY_COL = (MURGE_TRACE_COL-1)/vcsc->dof + 1;
            MY_ROW = (MURGE_TRACE_ROW-1)/vcsc->dof + 1;
            MY_COL_INNODE = (MY_COL-1)%vcsc->dof;
            MY_ROW_INNODE = (MY_ROW-1)%vcsc->dof;
            for (i = 0; i < vcsc->colsizes[MY_COL]; i++) {
              if (vcsc->rows[MY_COL][i] > ROW_NODE+1)
                /* row not found */
                break;
              if (vcsc->rows[MY_COL][i] == ROW_NODE+1) {
                INTS idx = i*dof2 +
                  MY_COL_INNODE*vcsc->dof +
                  MY_ROW_INNODE;

                fprintf(stdout, "after vcsc[%d][%d] = %g (vcsc.rows[index] %d)\n",
                        MY_COL, idx,
                        vcsc->values[MY_COL][idx],
                        vcsc->rows[MY_COL][i]);
              }
            }
          }
#endif
          return MURGE_SUCCESS;
        }
      }
    }
    /* The node does not exists so we add a new node.
     * The new node is filled with NaN so that
     * we can ignore unset values.
     */
    MURGE_MEMALLOC(V, dof2, COEF);
    //memset(V, 0, dof2*sizeof(COEF));
    for (i = 0; i < dof2; i++)
      V[i] = (COEF)NAN;
    V[COL_INNODE*vcsc->dof+ROW_INNODE] = VALUE;
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
    if (COL == MURGE_TRACE_COL &&
        ROW == MURGE_TRACE_ROW ) {
      fprintf(stdout, " vcsc_add set COL %d ROW %d = %g\n",
              COL, ROW, VALUE);
    }
#endif

    vcsc_add_node(vcsc, COL_NODE+1, ROW_NODE+1, V, op, id);
    MURGE_FREE(V);
  }
  return MURGE_SUCCESS;
}

static inline
INTS vcsc_to_cscd(variable_csc_t   vcsc,
                  MPI_Comm         comm,
                  OUT_INTS      *  n,
                  OUT_INTL      ** colptr_o,
                  OUT_INTS      ** rows_o,
                  COEF          ** values_o,
                  OUT_INTS      ** l2g_o,
                  OUT_INTS      ** g2l_o,
                  COEF (*op)(COEF , COEF),
                  INTS             dropnonlocal,
                  INTS             id) {
  int procnbr,me;
  INTS iter;
#ifdef MURGE_TIME
  Clock            clock;
#endif
  CLOCK_INIT;
  MPI_Comm_size(comm, &procnbr);
  MPI_Comm_rank(comm, &me);

  if (*l2g_o == NULL) {
    INTS *sizecols_max;
    INTS *owners;

    /* Count the number of entries per column */
    MURGE_MEMALLOC(sizecols_max,
                   vcsc->n*procnbr,
                   INTS);

    MPI_Allreduce (vcsc->colsizes,
                   sizecols_max, vcsc->n, MPI_INTS, MPI_MAX,
                   comm);

    CLOCK_PRINT("- After Allreduce");
    /* The processor with more entries is the owner */
    MURGE_MEMALLOC(owners, 2*vcsc->n, INTS);
    for (iter = 0; iter < vcsc->n; iter++) {
      if (sizecols_max[iter] == vcsc->colsizes[iter])
        owners[iter] = me;
      else
        owners[iter] = -1;
    }
    /* If two processeurs have the same numbers of entries,
     * keep higher rank => less work on proc 0, more on proc commSize-1
     */
    MPI_Allreduce (owners, owners+vcsc->n,
                   vcsc->n, MPI_INTS, MPI_MAX, comm);
    CLOCK_PRINT("- After Allreduce 2");

    MURGE_MEMALLOC(*l2g_o, vcsc->n, OUT_INTS);
    *n=0;
    for (iter = 0; iter < vcsc->n; iter++) {
      if (owners[vcsc->n+iter] == me)
        (*l2g_o)[(*n)++] = iter+1;
    }
    MURGE_REALLOC(*l2g_o, *n, OUT_INTS);
    CLOCK_PRINT("computing l2g");
    {
      PASTIX_INT vcsc_n = (PASTIX_INT)vcsc->n;
      cscd_build_g2l(*n,
                     *l2g_o,
                     comm,
                     &vcsc_n,
                     g2l_o);
      MURGE_TRACE_MALLOC(PTR_MEMSIZE(*g2l_o), char);
    }
    CLOCK_PRINT("computing g2l");
    MURGE_FREE(owners);
    MURGE_FREE(sizecols_max);
  }

#undef LOC_NODE_COUNT
#undef NODE_OWNER
#undef NODE_COUNT
#undef NODE_IDX
  if (dropnonlocal == 1) {
    /* if we drop non local columns we don't need to send it...
     * Just delete unnecessary columns.
     */
    for (iter = 0; iter < vcsc->n; iter++) {
      if (!((*g2l_o)[iter] > 0)) {
        MURGE_FREE(vcsc->rows[iter]);
        if (vcsc->dof > 0) {
          MURGE_FREE(vcsc->values[iter]);
        }
        vcsc->colsizes[iter] = 0;
      }
    }

  } else {
    /* send non local data */

    INTS p;
    INTS dof2 = vcsc->dof*vcsc->dof;
    INTS * count     = NULL;
    INTS ** colptr   = NULL;
    INTS ** rows     = NULL;
    COEF ** values   = NULL;
    INTS * loccolnum = NULL;
    MPI_Status    status;
    MPI_Request * send_sizes_req  = NULL;
    MPI_Request * send_colptr_req = NULL;
    MPI_Request * send_rows_req   = NULL;
    MPI_Request * send_values_req = NULL;

    MURGE_MEMALLOC(count, 2*procnbr, INTS);
    MURGE_MEMALLOC(send_sizes_req, procnbr,  MPI_Request);
    MURGE_MEMALLOC(send_colptr_req, procnbr, MPI_Request);
    MURGE_MEMALLOC(send_rows_req, procnbr,   MPI_Request);
    if (vcsc->dof > 0) {
      MURGE_MEMALLOC(send_values_req, procnbr, MPI_Request);
    }
    memset(count, 0, 2*procnbr*sizeof(INTS));

    /* count number of entries to send to all processors */
    for (iter = 0; iter < vcsc->n; iter++) {
      if (!((*g2l_o)[iter] > 0)) {
        /* data is not local */
        INTS owner = -(*g2l_o)[iter];
        /* if a column does not belong to anyone */
        if (owner == procnbr) continue;
        count[2*owner]++;
        count[2*owner+1] += vcsc->colsizes[iter];
      }
    }

    MURGE_MEMALLOC(colptr, procnbr, INTS*);
    MURGE_MEMALLOC(rows,   procnbr, INTS*);
    if (vcsc->dof > 0) {
      MURGE_MEMALLOC(values, procnbr, COEF*);
    }
    MURGE_MEMALLOC(loccolnum, procnbr, INTS);

    for (p = 0; p < procnbr; p++) {

      if (p == me) continue;
      /* create new csc to send to processor p*/

      MURGE_MEMALLOC(colptr[p], count[2*p]+1, INTS);
      MURGE_MEMALLOC(rows[p],   count[2*p+1], INTS);
      if (vcsc->dof > 0) {
        MURGE_MEMALLOC(values[p], dof2*count[2*p+1], COEF);
      }
      colptr[p][0] = 1;
      loccolnum[p] = 0;
    }
    for (iter = 0; iter < vcsc->n; iter++) {
      if (!((*g2l_o)[iter] > 0)) {
        p = -(*g2l_o)[iter];
        colptr[p][loccolnum[p]+1] = colptr[p][loccolnum[p]] + vcsc->colsizes[iter];
        memcpy(&(rows[p][colptr[p][loccolnum[p]]-1]),
               vcsc->rows[iter],
               vcsc->colsizes[iter]*sizeof(INTS));
        MURGE_FREE(vcsc->rows[iter]);
        if (vcsc->dof > 0) {
          memcpy(&(values[p][(colptr[p][loccolnum[p]]-1)*dof2]),
                 vcsc->values[iter],
                 vcsc->colsizes[iter]*sizeof(COEF)*dof2);

#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
          if ((iter)*vcsc->dof <= MURGE_TRACE_COL-1 &&
              (iter+1)*vcsc->dof > MURGE_TRACE_COL-1) {
            int index = colptr[p][loccolnum[p]]-1;
            while (index < colptr[p][loccolnum[p]+1]-1 &&
                   rows[p][index] !=  (MURGE_TRACE_ROW-1)/vcsc->dof+1)
              index++;
            if ( (rows[p][index]-1)*vcsc->dof <=  (MURGE_TRACE_ROW-1) &&
                 (rows[p][index])*vcsc->dof > MURGE_TRACE_ROW) {
              fprintf(stdout, "send values values[p][%d*%d*%d + %d *%d + %d] = %g"
                      " (rows[p][index] %d)\n",
                      index, vcsc->dof, vcsc->dof,
                      (MURGE_TRACE_COL-1- (iter)*vcsc->dof), vcsc->dof,
                      (MURGE_TRACE_ROW-1 -(rows[p][index]-1)*vcsc->dof),
                      values[p][index*vcsc->dof*vcsc->dof +
                                (MURGE_TRACE_COL-1- (iter)*vcsc->dof)* vcsc->dof +
                                (MURGE_TRACE_ROW-1 -(rows[p][index]-1)*vcsc->dof)],
                      rows[p][index]);
            }
          }
#endif
          MURGE_FREE(vcsc->values[iter]);
        }
        vcsc->colsizes[iter] = 0;
        loccolnum[p]++;
      }
    }

    for (p = 0; p < procnbr; p++) {
      if (p == me) continue;

      /* isend the count + csc to p*/
      MPI_Isend(&(count[2*p]), 2,            MPI_INTS, p, TAG_VCSC_COUNT,  comm,
                &(send_sizes_req[p]));
      MPI_Isend(colptr[p],     count[2*p]+1, MPI_INTS, p, TAG_VCSC_COLPTR, comm,
                &(send_colptr_req[p]));
      MPI_Isend(rows[p],       count[2*p+1], MPI_INTS, p, TAG_VCSC_ROWS,   comm,
                &(send_rows_req[p]));
      if (vcsc->dof > 0) {
        MPI_Isend(values[p], dof2*count[2*p+1], MURGE_MPI_COEF, p, TAG_VCSC_VALUES, comm,
                  &(send_values_req[p]));

      }
    }
    CLOCK_PRINT("+ Sending data");
    for (p = 0; p < procnbr; p++) {

      if (p == me) continue;
      /* receive data */
      MPI_Recv(&(count[2*me]), 2, MPI_INTS, p, TAG_VCSC_COUNT, comm, &status);
      MURGE_MEMALLOC(colptr[me], count[2*me]+1, INTS);
      MURGE_MEMALLOC(rows[me],   count[2*me+1], INTS);
      if (vcsc->dof > 0) {
        MURGE_MEMALLOC(values[me], dof2*count[2*me+1], COEF);
      }
      MPI_Recv(colptr[me], count[2*me]+1, MPI_INTS, p, TAG_VCSC_COLPTR, comm, &status);
      MPI_Recv(rows[me], count[2*me+1],   MPI_INTS, p, TAG_VCSC_ROWS, comm, &status);
      if (vcsc->dof > 0) {
        MPI_Recv(values[me], dof2*count[2*me+1], MURGE_MPI_COEF, p, TAG_VCSC_VALUES, comm,
                 &status);
      }

      /* Add received CSC to vcsc */
      if (vcsc->dof > 0) {
        for (iter = 0; iter < count[2*me]; iter++) {
          INTS iter2;
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
          if (((*l2g_o)[iter]-1)*vcsc->dof <= MURGE_TRACE_COL-1 &&
              ((*l2g_o)[iter])*vcsc->dof > MURGE_TRACE_COL-1) {
            int index = colptr[me][iter]-1;
            while (index < colptr[me][iter+1]-1 &&
                   rows[me][index] !=  (MURGE_TRACE_ROW-1)/vcsc->dof+1)
              index++;
            if (rows[p][index] ==  (MURGE_TRACE_ROW-1)/vcsc->dof+1) {

              fprintf(stdout, "receved values values[p][%d*%d*%d + %d *%d + %d] = %g"
                      " (rows[p][index] %d)\n",
                      index, vcsc->dof, vcsc->dof,
                      (MURGE_TRACE_COL-1- ((*l2g_o)[iter]-1)*vcsc->dof), vcsc->dof,
                      (MURGE_TRACE_ROW-1 -(rows[me][index]-1)*vcsc->dof),
                      values[me][index*vcsc->dof*vcsc->dof +
                                 (MURGE_TRACE_COL-1- ((*l2g_o)[iter]-1)*vcsc->dof)* vcsc->dof +
                                 (MURGE_TRACE_ROW-1 -(rows[me][index]-1)*vcsc->dof)],
                      rows[me][index]);
            }
          }
#endif

          for (iter2 = colptr[me][iter]-1; iter2 < colptr[me][iter+1]-1; iter2++) {
            vcsc_add_node(vcsc,
                          (*l2g_o)[iter], rows[me][iter2],
                          &(values[me][iter2*dof2]), op, id);
          }
        }
      } else {
        for (iter = 0; iter < count[2*me]; iter++) {
          INTS iter2;
          COEF f = 0.0;
          for (iter2 = colptr[me][iter]-1; iter2 < colptr[me][iter+1]-1; iter2++) {
            vcsc_add_node(vcsc, (*l2g_o)[iter], rows[me][iter2],
                          &f, op, id);
          }
        }
      }
      MURGE_FREE(colptr[me]);
      MURGE_FREE(rows[me]);
      if (vcsc->dof > 0)
        MURGE_FREE(values[me]);
    }
    CLOCK_PRINT("+ Receiving and merging data");

    /* wait for reception */
    for (p = 0; p < procnbr; p++) {
      if (p == me) continue;

      MPI_Wait(&(send_sizes_req[p]), &status);
      MPI_Wait(&(send_colptr_req[p]), &status);
      MURGE_FREE(colptr[p]);
      MPI_Wait(&(send_rows_req[p]), &status);
      MURGE_FREE(rows[p]);
      if (vcsc->dof > 0) {
        MPI_Wait(&(send_values_req[p]), &status);
        MURGE_FREE(values[p]);
      }
    }

    MURGE_FREE(count);
    MURGE_FREE(colptr);
    MURGE_FREE(rows);
    MURGE_FREE(loccolnum);
    MURGE_FREE(send_sizes_req);
    MURGE_FREE(send_colptr_req);
    MURGE_FREE(send_rows_req);
    if (vcsc->dof > 0) {
      MURGE_FREE(values);
      MURGE_FREE(send_values_req);
    }
    CLOCK_PRINT("+ Completing sends");
  }
  /* Build the CSC from the vcsc */
  {
    INTS dof2 = vcsc->dof*vcsc->dof;
    if (*colptr_o == NULL) {
      /* build a new CSCd */
      MURGE_MEMALLOC(*colptr_o, (*n)+1, OUT_INTL);
      (*colptr_o)[0] = 1;
      for (iter = 0; iter < *n; iter++) {
        (*colptr_o)[iter + 1] = (*colptr_o)[iter] + vcsc->colsizes[(*l2g_o)[iter]-1];
      }
      MURGE_MEMALLOC(*rows_o, (*colptr_o)[*n]-1, OUT_INTS);
      if (vcsc->dof > 0 )
        MURGE_MEMALLOC(*values_o, dof2*(*colptr_o)[*n]-1, COEF);
      for (iter = 0; iter < *n; iter++) {
        INTS iter2;
        for (iter2 = (*colptr_o)[iter]-1; iter2 < (*colptr_o)[iter+1]-1; iter2++) {
          INTS vcsc_index = iter2-((*colptr_o)[iter]-1);
          (*rows_o)[iter2] = vcsc->rows[ (*l2g_o)[iter]-1 ][vcsc_index];
          if (vcsc->dof > 0 ) {
            INTS iterdof;
            for (iterdof = 0; iterdof < dof2; iterdof++)
              (*values_o)[dof2*iter2+iterdof] =
                vcsc->values[ (*l2g_o)[iter]-1 ][dof2*vcsc_index+iterdof];
          }
        }
      }
    } else {
      /* fill existing CSCd, the vcsc must fit in the given CSCd */
      for (iter = 0; iter < *n; iter++) {
        INTS iter2;
        INTS iter3 = (*colptr_o)[iter]-1;
        INTS MYCOL = (*l2g_o)[iter];
        for (iter2 = 0; iter2 < vcsc->colsizes[MYCOL-1]; iter2++) {
          INTS MYROW = vcsc->rows[MYCOL-1][iter2];
          for (;iter3 < (*colptr_o)[iter+1]-1; iter3++) {
            if (MYROW == (*rows_o)[iter3])
              break;
          }
          if (iter3 == (*colptr_o)[iter+1]-1) {
            errorPrint("%d %d is not existing in the previous built CSCd\n",
                       MYROW, MYCOL);
            return MURGE_ERR_PARAMETER;
          } else {
            INTS iterdof;
            for (iterdof = 0; iterdof < dof2; iterdof++) {
              if (!isnan(vcsc->values[MYCOL-1][dof2*iter2+iterdof])) {
                /* ignore NaN values from VCSC */
                (*values_o)[dof2*iter3+iterdof] = op(
                  (*values_o)[dof2*iter3+iterdof],
                  vcsc->values[MYCOL-1][dof2*iter2+iterdof]);
              }
            }
          }
        }
#if (defined MURGE_TRACE_COL && defined MURGE_TRACE_ROW)
        if (((*l2g_o)[iter]-1)*vcsc->dof <= MURGE_TRACE_COL-1 &&
            ((*l2g_o)[iter])*vcsc->dof > MURGE_TRACE_COL-1) {
          int index = (*colptr_o)[iter]-1;
          while (index < (*colptr_o)[iter+1]-1 &&
                 (*rows_o)[index] !=  (MURGE_TRACE_ROW-1)/vcsc->dof+1)
            index++;

          fprintf(stdout, "final value (*values_o)[%d*%d*%d + %d *%d + %d] = %g (rows[p][index] %d (*l2g_o)[iter] %d (%d))\n",
                  index, vcsc->dof, vcsc->dof,
                  (MURGE_TRACE_COL-1- ((*l2g_o)[iter]-1)*vcsc->dof), vcsc->dof,
                  (MURGE_TRACE_ROW-1 -((*rows_o)[index]-1)*vcsc->dof),
                  (*values_o)[index*vcsc->dof*vcsc->dof +
                              (MURGE_TRACE_COL-1- ((*l2g_o)[iter]-1)*vcsc->dof)* vcsc->dof +
                              (MURGE_TRACE_ROW-1 -((*rows_o)[index]-1)*vcsc->dof)],
                  (*rows_o)[index], (*l2g_o)[iter], iter);
        }
#endif
      }
    }
  }
  CLOCK_PRINT("+ Bulding final CSCd");
  return MURGE_SUCCESS;
}
#undef OUT_INTL
#undef OUT_INTS
#undef TAG_VCSC_COUNT
#undef TAG_VCSC_COLPTR
#undef TAG_VCSC_ROWS
#undef TAG_VCSC_VALUES

static inline
INTS vcsc_destroy(variable_csc_t vcsc, INTS id) {
  INTS iter;
  MURGE_FREE(vcsc->colsizes);
  for (iter = 0; iter < vcsc->n; iter++) {
    if (vcsc->rows[iter] != NULL) {
      MURGE_FREE(vcsc->rows[iter]);
      if (vcsc->dof > 0)
        MURGE_FREE(vcsc->values[iter]);
    }
  }
  MURGE_FREE(vcsc->rows);
  if (vcsc->dof > 0)
    MURGE_FREE(vcsc->values);
  MURGE_FREE(vcsc);
  return MURGE_SUCCESS;
}

#undef CHUNCKSIZE
