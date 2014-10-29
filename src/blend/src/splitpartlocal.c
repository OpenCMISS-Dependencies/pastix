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
#ifdef PASTIX_DYNSCHED
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <assert.h>
#include <sys/stat.h>

#include "common_pastix.h"
#include "queue.h"
#include "extendVector.h"
#include "cand.h"
#include "ftgt.h"
#include "symbol.h"
#include "simu.h"
#include "cost.h"
#include "elimin.h"
#include "bulles.h"
#include "param_blend.h"
#include "blendctrl.h"
#include "dof.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "costfunc.h"
#include "splitpartlocal.h"

/* Check before to activate, a problem occured to the merge and I'm not sure about the result of updatetreelocal */
//#define SPLITPARTLOCAL_BALANCING

extern void setTreeLevel    (Cand *, const EliminTree *);
extern void setTreeCostLevel(Cand *, const EliminTree *, const CostMatrix *);
int VerifTaskRepartition(BlendCtrl *ctrl, SimuCtrl * simuctrl);
void splitPartLocalUpdate(const BlendCtrl  *ctrl,
                          const SimuCtrl   *simuctrl,
                          BubbleTree       *btree,
                          double           *bcost,
                          const double      totalcost);

void splitPartLocal(BlendCtrl *ctrl, SimuCtrl * simuctrl,
                    SymbolMatrix *symbmtx, const Dof *dofptr)
{

  double *bubble_cost = NULL;
  PASTIX_INT     i, j, size;

  /* Reset the cost of non local computation and update cost of local ones
   * In the resulting tree, the sons are sorted by decreasing
   * cost of their subtree */
  subtreeUpdateCostLocal(ROOT(ctrl->etree), ctrl, symbmtx, simuctrl,
                         dofptr, ctrl->clustnum);

  setTreeLevel(ctrl->candtab, ctrl->etree);
  if(ctrl->option->costlevel)
    setTreeCostLevel(ctrl->candtab, ctrl->etree, ctrl->costmtx);

  /* Initialisation de la structure de d'arbre de bulles */
  Bubble_InitTree(ctrl->btree, ctrl->thrdlocnbr);

  MALLOC_INTERN(bubble_cost, ctrl->btree->nodemax, double);
  for (i=0; i< ctrl->btree->nodemax; i++)
    bubble_cost[i] = 0.0;

  /* Proportionnal mapping of the tree on the available threads
     with no crossing */
  propMappTreeLocal(ctrl, symbmtx, simuctrl, bubble_cost);

  /* construction de l'arbre de bulles */
  Bubble_BuildTree(ctrl->btree);
  ctrl->bublnbr = ctrl->btree->nodenbr;

  if (ctrl->option->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    {
      FILE *out;
      char  filename[18];
      sprintf(filename,"BubbleTree.%02d.dot", (int)ctrl->clustnum);
      filename[17] = '\0';
      OUT_OPENFILEINDIR(ctrl->option->iparm, out, filename, "w");
      Bubble_Print(ctrl->btree, bubble_cost,
                   ctrl->costmtx->cblktab[ROOT(ctrl->etree)].subtree, out);
      OUT_CLOSEFILEINDIR(out);
    }

  simuctrl->clustab[ctrl->clustnum].fprocnum = 0;
  simuRealloc(simuctrl, ctrl->procnbr, ctrl->bublnbr);

#if defined(SPLITPARTLOCAL_BALANCING)
  splitPartLocalUpdate(ctrl, simuctrl, ctrl->btree, bubble_cost,
                       ctrl->costmtx->cblktab[ROOT(ctrl->etree)].subtree);

  if (ctrl->option->iparm[IPARM_VERBOSE] > API_VERBOSE_NO)
    {
      FILE *out;
      char  filename[19];
      sprintf(filename,"BubbleTree2.%02d.dot", (int)ctrl->clustnum);
      filename[18] = '\0';
      OUT_OPENFILEINDIR(ctrl->option->iparm, out, filename, "w");
      Bubble_Print(ctrl->btree, bubble_cost,
                   ctrl->costmtx->cblktab[ROOT(ctrl->etree)].subtree, out);
      OUT_CLOSEFILEINDIR(out);
    }
#endif

  for (i=0; i<ctrl->bublnbr; i++) {
    size = queueSize(ctrl->btree->nodetab[i].taskheap);
    extendint_Init(simuctrl->proctab[i].tasktab, size);

    /* Move the queue to an extendint array */
    for(j=0; j<size; j++) {
      extendint_Add(simuctrl->proctab[i].tasktab,
                    queueGet(ctrl->btree->nodetab[i].taskheap));
    }
  }

  memFree_null(bubble_cost);
}


/*+ Repartition symbolic matrix using proportionnal mapping method +*/
void propMappTreeLocal(BlendCtrl *ctrl, const SymbolMatrix *symbmtx,
                       const SimuCtrl * simuctrl, double * bubble_cost)
{
  propMappSubtreeLocalNC(ctrl, symbmtx, simuctrl, ROOT(ctrl->etree),
                         0, ctrl->thrdlocnbr-1, bubble_cost);
}

static inline void
assignNodeToBubble(BlendCtrl *ctrl,
                   const SymbolMatrix *symbmtx,
                   const SimuCtrl * simuctrl,
                   PASTIX_INT rootnum,
                   PASTIX_INT bubblenum)
{
  PASTIX_INT i, j;
  PASTIX_INT tasknum;
  tasknum = simuctrl->bloktab[symbmtx->cblktab[rootnum].bloknum].tasknum;
  ASSERTDBG(tasknum != -1, MOD_BLEND);

  if (ctrl->candtab[rootnum].distrib == D1)
    {
      ASSERTDBG(simuctrl->tasktab[tasknum].taskid == COMP_1D, MOD_BLEND);
      queueAdd(ctrl->btree->nodetab[bubblenum].taskheap, tasknum,
               simuctrl->tasktab[tasknum].prionum);
    }
  /* D2 */
  else
    {
      /* Diag */

      if (ctrl->proc2clust[simuctrl->blprtab[symbmtx->cblktab[rootnum].bloknum]] == ctrl->clustnum)
        {
          tasknum = simuctrl->bloktab[symbmtx->cblktab[rootnum].bloknum].tasknum;

          ASSERTDBG(tasknum != -1, MOD_BLEND);
          ASSERTDBG(simuctrl->tasktab[tasknum].taskid == DIAG, MOD_BLEND);
          ASSERTDBG(ctrl->proc2clust[simuctrl->blprtab[simuctrl->tasktab[tasknum].bloknum]] ==
                    ctrl->clustnum, MOD_BLEND);

          queueAdd(ctrl->btree->nodetab[bubblenum].taskheap,
                   tasknum, simuctrl->tasktab[tasknum].prionum);
        }

      /* E1 tasks */
      for(i=symbmtx->cblktab[rootnum].bloknum+1;
          i<symbmtx->cblktab[rootnum+1].bloknum; i++)
        {
          tasknum = simuctrl->bloktab[i].tasknum;
          if (ctrl->proc2clust[simuctrl->blprtab[i]] == ctrl->clustnum)
            {
              ASSERTDBG(tasknum != -1, MOD_BLEND);
              ASSERTDBG(simuctrl->tasktab[tasknum].taskid == E1, MOD_BLEND);
              ASSERTDBG(ctrl->proc2clust[simuctrl->blprtab[simuctrl->tasktab[tasknum].bloknum]] ==
                        ctrl->clustnum, MOD_BLEND);

              queueAdd(ctrl->btree->nodetab[bubblenum].taskheap, tasknum,
                       simuctrl->tasktab[tasknum].prionum);
            }
          else
            {
              ASSERTDBG(ctrl->proc2clust[simuctrl->blprtab[simuctrl->tasktab[tasknum].bloknum]] !=
                        ctrl->clustnum, MOD_BLEND);
            }

          /* E2 tasks */
          for(j=i;j<symbmtx->cblktab[rootnum+1].bloknum;j++)
            {
              tasknum++;
              if (ctrl->proc2clust[simuctrl->blprtab[simuctrl->tasktab[tasknum].bloknum]] ==
                  ctrl->clustnum)
                queueAdd(ctrl->btree->nodetab[bubblenum].taskheap, tasknum,
                         simuctrl->tasktab[tasknum].prionum);

              ASSERTDBG(simuctrl->tasktab[tasknum].taskid == E2, MOD_BLEND);
            }
        }
    }
}

void propMappSubtreeLocalNC(BlendCtrl *ctrl, const SymbolMatrix *symbmtx,
                            const SimuCtrl * simuctrl, PASTIX_INT rootnum,
                            PASTIX_INT fcandnum, PASTIX_INT lcandnum, double *bubble_cost)
{
  PASTIX_INT i, ison, son;
  PASTIX_INT candnbr;
  PASTIX_INT nbmtsons = 0;
  PASTIX_INT fcand = 0;
  PASTIX_INT sonsnbr;
  PASTIX_INT bubblenum;
  double sonscost, cost, isocost;

  candnbr = lcandnum - fcandnum + 1;

  /*
   * If the subtree is handled by another process,
   * the cost is null, so we skip it
   */
  if (ctrl->costmtx->cblktab[rootnum].subtree == 0.0)
    return;

  /* Get the index of the bubble requested.
     Generate the new bubble if needed.  */
  bubblenum = Bubble_Add(ctrl->btree, fcandnum, lcandnum,
                         ctrl->candtab[rootnum].costlevel,
                         ctrl->candtab[rootnum].treelevel);

  /* If we have only one candidate, let's do it simpler and faster */
  if ( fcandnum == lcandnum ) {
    propMappSubtreeLocalOn1P(ctrl, symbmtx, simuctrl, rootnum,
                             bubblenum);

    /* Update given computation costs */
    cost = ctrl->costmtx->cblktab[rootnum].subtree;
    bubble_cost[bubblenum] += cost;
    return;
  }

  /* Substract the local cost to the subtrees */
  sonscost = ctrl->costmtx->cblktab[rootnum].subtree -
    ctrl->costmtx->cblktab[rootnum].total;
  sonsnbr  = ctrl->etree->nodetab[rootnum].sonsnbr;

  /* If this is a leaf, let's return */
  if(ctrl->etree->nodetab[rootnum].sonsnbr == 0)
      goto end;

  isocost = sonscost / candnbr;

  /*
   * Let's found how many sons require several threads and we put them
   * in queue ordered by the inverse of the cost
   * On thread is attached to each subtree, and the remaining threads
   * are distributed dynnamically to the less important ratio.
   */
  {
    int *nbcand = NULL;
    int nbthreads;
    Queue *queue_tree;

    MALLOC_INTERN(queue_tree, 1, Queue);
    queueInit(queue_tree, sonsnbr);
    for(i=0; i<sonsnbr; i++ )
      {
        son = TSON(ctrl->etree, rootnum, i);

        cost = ctrl->costmtx->cblktab[son].subtree;

        /*  We asume that the sons are sorted by decreasing order
            of the cost */
        if ( cost <= isocost )
          break;

        queueAdd2(queue_tree, i, 1. / cost, -cost );
        nbmtsons++;
      }

    if (nbmtsons > 0) {
      MALLOC_INTERN(nbcand, nbmtsons, int);

      for(i=0; i<nbmtsons; i++) {
        nbcand[i] = 1;
      }

      /* We already affected one thread per son */
      nbthreads = candnbr - nbmtsons;
      ASSERTDBG( nbthreads >= 0, MOD_BLEND );
      assert( nbthreads >= 0 );

      /* Let's dispatch fairly the remaining threads */
      while ( nbthreads > 0 ) {
        i = queueGet(queue_tree);

        son = TSON(ctrl->etree, rootnum, i);
        cost = ctrl->costmtx->cblktab[son].subtree;

        nbcand[i]++;
        nbthreads--;

        queueAdd2(queue_tree, i,
                  (double)nbcand[i] / cost, -cost );
      }

      fcand = 0;
      for( i=0; i<nbmtsons; i++ ) {
        son = TSON(ctrl->etree, rootnum, i);

        /* Compute the number of threads to give to this subtree */
        propMappSubtreeLocalNC(ctrl, symbmtx, simuctrl, son,
                               fcandnum+fcand,
                               fcandnum+fcand+nbcand[i]-1,
                               bubble_cost);
        fcand += nbcand[i];
      }
    }
    queueExit(queue_tree);
    memFree(queue_tree);

    if( nbmtsons > 0 )
      memFree(nbcand);
  }

end:
  {
    Queue *queue_proc;

    /* allocate queue of proc  */
    MALLOC_INTERN(queue_proc, 1, Queue);
    queueInit(queue_proc, candnbr);

    /* Fill queue proc order by remain cost descending */
    for (i=fcandnum; i<=lcandnum; i++) {
      queueAdd(queue_proc, i, bubble_cost[i]);
    }

    /*
     * Let's take care of the remaining subtrees that can be distributed
     * on only one candidate
     */
#ifdef TOP_BUBBLES
        PASTIX_INT bubblenum_save = bubblenum;
#endif
    if ( nbmtsons < sonsnbr &&
         ctrl->costmtx->cblktab[TSON(ctrl->etree, rootnum, nbmtsons)].subtree > 0.0 )
      {
        /* Distribute all subtree executed by only one thread */
        for(ison=nbmtsons; ison<sonsnbr; ison++ )
          {
            son = TSON(ctrl->etree, rootnum, ison);
            bubblenum = queueGet(queue_proc);

            /** Cost in the current subtree to be mapped **/
            cost = ctrl->costmtx->cblktab[son].subtree;

            if ( cost <= 0.0 )
              break;

            propMappSubtreeLocalOn1P(ctrl, symbmtx, simuctrl,
                                     son, bubblenum);

            bubble_cost[bubblenum] += cost;

            queueAdd(queue_proc, bubblenum, bubble_cost[bubblenum]);
          }

      }
#ifdef TOP_BUBBLES
    bubblenum = bubblenum_save;
#endif
    /* Let's add the work associated to this node to the current bubble */
    if (ctrl->costmtx->cblktab[rootnum].total != 0.0 )
      {
#ifndef TOP_BUBBLES
        bubblenum = queueGet(queue_proc);
#endif
        assignNodeToBubble( ctrl, symbmtx, simuctrl, rootnum, bubblenum );

        ctrl->candtab[rootnum].cluster = bubblenum;
        bubble_cost[bubblenum] += ctrl->costmtx->cblktab[rootnum].total;
      }

    queueExit(queue_proc);
    memFree(queue_proc);
  }
  return;
}

void propMappSubtreeLocalOn1P(BlendCtrl *ctrl, const SymbolMatrix *symbmtx,
                              const SimuCtrl *simuctrl,
                              PASTIX_INT rootnum, PASTIX_INT bubblenum)
{
  PASTIX_INT i;
  PASTIX_INT sonsnbr;

  /* Si le sous-arbre est traité par un autre processus MPI */
  if (ctrl->costmtx->cblktab[rootnum].subtree == 0.0)
    return;

  /* work that each processor is intended to get from this treenode */
  if (ctrl->costmtx->cblktab[rootnum].total != 0.0 )
    {
      assignNodeToBubble( ctrl, symbmtx, simuctrl, rootnum, bubblenum );
      ctrl->candtab[rootnum].cluster = bubblenum;
    }

  sonsnbr = ctrl->etree->nodetab[rootnum].sonsnbr;

  /* Add all sons to the candidate */
  for(i=0;i<sonsnbr;i++)
    propMappSubtreeLocalOn1P(ctrl, symbmtx, simuctrl,
                             TSON(ctrl->etree, rootnum, i),
                             bubblenum);

  return;
}


int VerifTaskRepartition(BlendCtrl *ctrl, SimuCtrl * simuctrl){

  int i, j, clust;

  for (clust=0; clust<ctrl->clustnbr; clust++){
    if (ctrl->clustnum == clust)
      {
        int *indice = NULL;
        int *value  = NULL;
        int *fin    = NULL;
        int nbfin=0, nberror;
        Queue ** copyqueue = NULL;

        MALLOC_INTERN(indice,    ctrl->proclocnbr, int);
        MALLOC_INTERN(fin,       ctrl->thrdlocnbr, int);
        MALLOC_INTERN(value,     ctrl->thrdlocnbr, int);
        MALLOC_INTERN(copyqueue, ctrl->thrdlocnbr, Queue*);

        for (i=0; i<ctrl->procnbr; i++)
          indice[i] = 0;

        for (i=0; i<ctrl->thrdlocnbr; i++){
          fin[i] = 1;
          value[i] = -1;
          MALLOC_INTERN(copyqueue[i], 1, Queue);
          copyqueue[i] = queueCopy(copyqueue[i],
                                   ctrl->btree->nodetab[i].taskheap);
        }

        while (nbfin < ctrl->thrdlocnbr ){

          for (i=0; i<ctrl->thrdlocnbr ; i++){
            if (fin[i] && value[i] < 0)
              {
                if (queueSize(copyqueue[i]) > 0)
                  value[i] = queueGet(copyqueue[i]);
                else{
                  value[i] = -1;
                  fin[i]   = 0;
                  nbfin++;
                }
              }
          }

          if (!(nbfin<ctrl->thrdlocnbr))
            break;

          printf("tasktab ");
          for (i=simuctrl->clustab[ctrl->clustnum].fprocnum;
               i<simuctrl->clustab[ctrl->clustnum].lprocnum+1; i++)
            printf("%ld ",
                   (indice[i]<extendint_Size(simuctrl->proctab[i].tasktab))?
                   extendint_Read(simuctrl->proctab[i].tasktab, indice[i]):-1);
          printf("\n");

          printf("queue ");
          for (i=0; i<ctrl->thrdlocnbr; i++)
            printf("%d ", value[i]);
          printf("\n");

          nberror = 0;
          for(i=simuctrl->clustab[ctrl->clustnum].fprocnum;
              i<simuctrl->clustab[ctrl->clustnum].lprocnum+1; i++)
            {
              for(j=0; j<ctrl->thrdlocnbr; j++)
                {
                  if (value[j] != -1 &&
                      indice[i] < extendint_Size(simuctrl->proctab[i].tasktab) &&
                      value[j] == extendint_Read(simuctrl->proctab[i].tasktab,
                                                 indice[i]))
                    {
                      indice[i]++;
                      value[j] = -2;
                      break;
                    }
                  else
                    nberror++;
                }
            }

          if (!(nberror < ctrl->proclocnbr*ctrl->thrdlocnbr)){
            printf("tasktab ");
            for(i=simuctrl->clustab[ctrl->clustnum].fprocnum;
                i<simuctrl->clustab[ctrl->clustnum].lprocnum+1; i++)
              printf("%d ", extendint_Read(simuctrl->proctab[i].tasktab,
                                           indice[i]));
            printf("\n");

            printf("queue ");
            for (i=0; i<ctrl->thrdlocnbr; i++)
              printf("%d ", value[i]);
            printf("\n");

          }
          ASSERT(nberror < ctrl->procnbr*ctrl->thrdlocnbr, MOD_BLEND);

        }
      }
    else
      sleep(10);
  }
  return 1;
}


static inline void updateCosts(const BubbleTree * btree,
                               const double     * bcost,
                               double *subtrees_costs)
{
  PASTIX_INT i;

  memset( subtrees_costs, 0, btree->nodenbr * sizeof(double) );
  for (i=0; i<btree->nodenbr; i++) {
    int father = btree->nodetab[i].fathnum;
    subtrees_costs[i] += bcost[i];

    while ( father != -1 ) {
      subtrees_costs[father] += bcost[i];
      father = btree->nodetab[father].fathnum;
    }
  }
}

void splitPartLocalUpdate(const BlendCtrl *ctrl,
                          const SimuCtrl  *simuctrl,
                          BubbleTree      *btree,
                          double          *bcost,
                          const double     totalcost)
{
  PASTIX_INT i, father, size;
  PASTIX_INT fprocnum, lprocnum;
  double percent, percenttgt;
  double *subtrees_costs = (double*)malloc( btree->nodenbr * sizeof(double) );
  updateCosts( btree, bcost, subtrees_costs );

  for (i=0; i<btree->nodenbr; i++)
    {
      father = btree->nodetab[i].fathnum;
      if (father == -1)
        continue;

      fprocnum = btree->nodetab[i].fcandnum;
      lprocnum = btree->nodetab[i].lcandnum;

      percent    = subtrees_costs[i] / totalcost;
      percenttgt = (double)( lprocnum-fprocnum+1 ) / (double)(btree->leavesnbr);

      size = queueSize( btree->nodetab[i].taskheap );
      if ( percent > percenttgt && size > 0) {
        PASTIX_INT j;
        double key1, cost;
        PASTIX_INT value;
        Queue *tmpqueue;

        MALLOC_INTERN(tmpqueue, 1, Queue);
        queueInit(tmpqueue, size );

        for(j=0; j<size; j++) {
          value = queueGet2( btree->nodetab[i].taskheap, &key1, NULL );
          queueAdd( tmpqueue, value, -key1 );
        }

        while( percent > percenttgt ) {
          value = queueGet2( tmpqueue, &key1, NULL );
          queueAdd(  btree->nodetab[father].taskheap, value, -key1 );

          if ( simuctrl->tasktab[value].cost == -1. ) {
            cost = ctrl->costmtx->cblktab[ simuctrl->tasktab[value].cblknum ].total;
          }
          else {
            cost = simuctrl->tasktab[value].cost;
          }

          print_debug(DBG_BUBBLESPLIT,
                      "percentages in bubble %ld (%lf):\n"
                      "\t\tcost= %lf, percent=%lf, taskcost= %lf\n"
                      "\t\tfather( %ld, %lf)\n",
                      (long)i, percenttgt, subtrees_costs[i], percent, cost,
                      (long)father, subtrees_costs[father]);

          bcost[i]          -= cost;
          bcost[father]     += cost;
          subtrees_costs[i] -= cost;
          percent = subtrees_costs[i] / totalcost;

          print_debug(DBG_BUBBLESPLIT,
                      "\t\tnew: cost=%lf, percent=%lf, fcost= %lf\n",
                      subtrees_costs[i], percent, subtrees_costs[father]);
        }

        size = queueSize(tmpqueue);
        for(j=0; j<size; j++) {
          value = queueGet2( tmpqueue, &key1, NULL );
          queueAdd( btree->nodetab[i].taskheap, value, -key1 );
        }

        queueExit(tmpqueue);
        memFree(tmpqueue);
        updateCosts( btree, bcost, subtrees_costs );
      }
    }
  free(subtrees_costs);
}
#else
/* ISO C forbids an empty source file */
#include "not_empty.h"
NOT_EMPTY(splitpartlocal)
#endif
