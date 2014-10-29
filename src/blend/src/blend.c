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
  File: blend.c

  Main blend source code.

  Re-split column blocs and distribute them on processors.

  Authors:
    - Pascal  Henon
    - Pierre  Ramet
    - Mathieu Faverge
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <memory.h>
#include <sys/stat.h>
#include <unistd.h>

#include "common_pastix.h"
#include "out.h"
#include "dof.h"
#include "ftgt.h"
#include "cost.h"
#include "symbol.h"
#include "queue.h"
#include "bulles.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "solverRealloc.h"
#include "elimin.h"
#include "extrastruct.h"
#include "extendVector.h"
#include "cand.h"
#include "param_blend.h"
#include "blendctrl.h"
#include "eliminfunc.h"
#include "splitpart.h"
#include "write_ps.h"
#include "simu.h"
#include "costfunc.h"
#include "splitpartlocal.h"
#include "distribPart.h"
#include "assembly.h"
#include "assemblyGener.h"
#include "solverMatrixGen.h"
#include "solver_check.h"
#include "symbol_cost.h"
#include "task.h"
#include "solver_check.h"
#include "fanboth2.h"
#include "blend.h"
#include "order.h"
#include "bordi.h"

/* #include "assert.h" */

/*
 * Function: solverBlend
 *
 * Main blend function
 *
 * Build the elimination graph from the symbolic partition.
 *
 * Build the cost matrix from the symbolic partition.
 *
 * Build the elimination tree from the symbolic partition.
 *
 * Distribute each column bloc on candidate processors.
 *
 * Build a new symbol matrix...
 *
 * Parameters:
 *   solvmtx    - Solver matrix structure.
 *   assemb1D   -
 *   assemb2D   -
 *   clustnbr   - Number of MPI processes.
 *   thrdlocnbr - Number of threads.
 *   cudanbr    - Number of cuda devices.
 *   clustnum   - Processor ID number.
 *   option     - Blend parameters.
 *   dofptr     -
 */
void solverBlend(SolverMatrix *solvmtx,
                 SymbolMatrix *symbmtx,
                 Assembly1D   *assemb1D,
                 Assembly2D   *assemb2D,
                 int           clustnbr,
                 int           thrdlocnbr,
                 int           cudanbr,
                 int           clustnum,
                 BlendParam   *option,
                 const Dof    *dofptr)
{
    BlendCtrl *ctrl;
    SimuCtrl  *simuctrl;
    FILE      *ps_file       = NULL;
    Clock      timer_all;
    Clock      timer_current;
    PASTIX_INT       *bcofind       = NULL;
    PASTIX_INT        page          = 0;
#ifdef DEBUG_BLEND
    PASTIX_INT        leader;
#endif
    /* initialisation of the control structure */
    MALLOC_INTERN(ctrl, 1, BlendCtrl);
    blendCtrlInit(ctrl, clustnbr, thrdlocnbr, cudanbr, clustnum, option);
#ifdef DEBUG_BLEND
    leader = ctrl->option->leader;
#endif
    if(clustnum >= clustnbr)
      {
        errorPrint("solverBlend parameter clustnum(%ld) is greater"
                   " than clustnbr(%ld).",
                   (long) clustnum, (long) clustnbr);
        EXIT(MOD_BLEND,INTERNAL_ERR);
      }

    clockInit(&timer_all);
    clockInit(&timer_current);
    if(ctrl->option->timer)
        {
            clockStart(&timer_all);
        }

    if(ctrl->option->ps)
        {
            ps_file = ps_open(ctrl->option->ps_filename);
            page = 1;
        }

    if ((ctrl->option->leader == clustnum) &&
        (ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO))
      {
        fprintf(stdout, OUT_CLUSTNBR, (long)clustnbr);
        fprintf(stdout, OUT_PROCNBR,  (long)ctrl->proclocnbr);
        fprintf(stdout, OUT_THRDNBR,  (long)ctrl->thrdlocnbr);
      }

    /* Blend */
    symbolRealloc(symbmtx);
    symbolBase(symbmtx,0);
    symbolCheck(symbmtx);

    /* Rustine */
#define RUSTINE
#ifdef RUSTINE
    printf("Debut Rustine\n");
    symbolRustine(symbmtx, symbmtx);
    printf("Fin Rustine\n");
#endif


#ifdef DEBUG_BLEND
    /** Verify the coherence of the initial symbol matrix **/
    if(ctrl->option->debug)
      {
        if ((ctrl->option->leader == clustnum) &&
            (ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO))
          fprintf(stdout, OUT_BLEND_CHKSMBMTX);
        symbolCheck(symbmtx);
      }
#endif

    /** OOC works only with 1D structures **/
    if(ctrl->option->ooc)
      {
        if ((ctrl->option->leader == clustnum) &&
            (ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO))
          fprintf(stdout, "Force 1D distribution because of OOC \n");
        ctrl->option->ratiolimit = INTVALMAX;
        ctrl->option->malt_limit = -1;
      }

    if(ctrl->option->count_ops && ctrl->option->leader == clustnum)
      symbCost(option->iparm, option->dparm, symbmtx, dofptr);

    if(ctrl->option->ps && symbmtx->nodenbr <= 961)
      {
            if(ctrl->option->leader == clustnum &&
               ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
                fprintf(stdout, "Genering a symbolic matrix draw in Post-Script format \n");
            ps_write_matrix(symbmtx, ps_file, &page);
      }


    if(ctrl->option->timer)
        {
            clockInit(&timer_current);
            clockStart(&timer_current);
        }

    if(ctrl->option->leader == clustnum &&
       ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, OUT_BLEND_ELIMGRAPH);
    MALLOC_INTERN(ctrl->egraph, 1, EliminGraph);

    egraphInit(ctrl->egraph);

    /** build the elimination graph from the symbolic partition **/
    eliminGraphBuild(symbmtx, ctrl->egraph);

    if(ctrl->option->timer)
        {
            clockStop(&timer_current);
            printf("--Graph build at time: %g --\n", clockVal(&timer_current));
        }


    MALLOC_INTERN(ctrl->costmtx, 1, CostMatrix);
    costInit(ctrl->costmtx);
    if(ctrl->option->leader == clustnum &&
       ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, OUT_BLEND_COSTMATRIX);
    if(ctrl->option->timer)
        {
            clockInit(&timer_current);
            clockStart(&timer_current);
        }

    /** Build the cost matrix from the symbolic partition **/
    costMatrixBuild(ctrl->costmtx, symbmtx, dofptr);

    if(ctrl->option->timer)
        {
            clockStop(&timer_current);
            printf("--Cost Matrix build at time: %g --\n", clockVal(&timer_current));
        }

    MALLOC_INTERN(ctrl->etree, 1, EliminTree);
    treeInit(ctrl->etree);
    if(ctrl->option->leader == clustnum &&
       ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, OUT_BLEND_ELIMTREE);
    if(ctrl->option->timer)
        {
            clockInit(&timer_current);
            clockStart(&timer_current);
        }

    /** Build the elimination tree from the symbolic partition **/
    eliminTreeBuild(symbmtx, ctrl);

    if(ctrl->option->timer)
        {
            clockStop(&timer_current);
            printf("--Tree build at time: %g --\n", clockVal(&timer_current));
        }


    if(ctrl->option->ps)
        {
            if(ctrl->option->leader == clustnum &&
               ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
                fprintf(stdout, "Genering an elimination tree draw in Post-Script format \n");
            ps_write_tree(ctrl->costmtx, ctrl->etree, ps_file, &page);
        }



    if(ctrl->option->leader == clustnum &&
       ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, "Spliting initial partition \n");
    if(ctrl->option->timer)
        {
            clockInit(&timer_current);
            clockStart(&timer_current);
        }

    /** repartitioning of the initial symbolic factorization
      and processing of candidate processors group for
      each colum bloc **/
    splitPart(symbmtx, ctrl, dofptr);

    if ( (ctrl->option->leader == clustnum) && (ctrl->option->tracegen == 1))
      {
        FILE *out;
        OUT_OPENFILEINDIR(ctrl->option->iparm, out, "elimintree.dot", "w");
        treePlot(ctrl->etree, out);
        OUT_CLOSEFILEINDIR(out);
      }

    if(ctrl->option->timer)
      {
        clockStop(&timer_current);
        printf("--Split build at time: %g --\n", clockVal(&timer_current));
      }

#ifdef DEBUG_BLEND
    /** Verify the coherence of the new symbol matrix **/
    if(ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO &&
       ctrl->option->debug)
      {
        if (ctrl->option->leader == clustnum)
          fprintf(stdout, "\n Checking the new symbol matrix \n");
        symbolCheck(symbmtx);
      }
#endif

    if(ctrl->option->count_ops && ctrl->option->leader == clustnum)
      symbCost(option->iparm, option->dparm, symbmtx, dofptr);

    if ((ctrl->option->leader == clustnum) &&
        (ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO))
      fprintf(stdout, "** New Partition: cblknbr=  %ld     bloknbr=  %ld     ratio=%f ** \n",
              (long)symbmtx->cblknbr, (long)symbmtx->bloknbr,
              (float)symbmtx->bloknbr/(float)symbmtx->cblknbr);

    if(ctrl->option->ps && (symbmtx->nodenbr <= 961))
      {
        if(ctrl->option->leader == clustnum &&
           ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
          fprintf(stdout, "Genering a symbolic matrix draw in Post-Script format \n");
        ps_write_matrix(symbmtx, ps_file, &page);
      }

    if(ctrl->option->count_ops &&
       ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO &&
       ctrl->option->leader == clustnum)
      fprintf(stdout, "Factorization of the new symbol matrix by Crout blok algo takes : %g \n",
              recursive_sum(0, symbmtx->cblknbr-1, crout_blok, symbmtx, dofptr ));

    if(ctrl->option->timer)
        {
            clockInit(&timer_current);
            clockStart(&timer_current);
        }
    /** the former graph can't be used any more **/
    egraphExit(ctrl->egraph);
    /** Build a new one on the new symbol matrix **/

    if( ctrl->option->leader == clustnum &&
        ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, OUT_BLEND_ELIMGRAPH2);
    MALLOC_INTERN(ctrl->egraph, 1, EliminGraph);

    egraphInit(ctrl->egraph);

    /** build the elimination graph from the symbolic partition **/
    eliminGraphBuild(symbmtx, ctrl->egraph);

    if(ctrl->option->timer)
        {
            clockStop(&timer_current);
            printf("--Graph build at time: %g --\n", clockVal(&timer_current));
        }

    if(ctrl->option->timer)
        {
            clockInit(&timer_current);
            clockStart(&timer_current);
        }

    /* initialize simu structure control */
    MALLOC_INTERN(simuctrl, 1, SimuCtrl);
    simuInit(simuctrl, symbmtx, ctrl->clustnbr, ctrl->procnbr,
             symbmtx->cblknbr, symbmtx->bloknbr, ctrl->candtab);

    /** Build tasks **/
    if( ctrl->option->leader == clustnum &&
        ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, OUT_BLEND_TASKGRAPH);
    taskBuild(simuctrl, symbmtx, ctrl->candtab, dofptr, ctrl->egraph, ctrl);

    if( ctrl->option->leader == clustnum &&
        ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
      fprintf(stdout, OUT_BLEND_NBTASK, (long)simuctrl->tasknbr);
    if(ctrl->option->timer)
      {
        clockStop(&timer_current);
        printf("--Task built at time: %g --\n", clockVal(&timer_current));
      }

    /** Distribution Phase **/
    if(ctrl->option->timer)
      {
        clockInit(&timer_current);
        clockStart(&timer_current);
      }

#ifdef DEBUG_BLEND
    ASSERT(check_candidat(symbmtx, ctrl)>=0,MOD_BLEND);
#endif

    if((ctrl->option->leader == clustnum) &&
       (ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO))
        fprintf(stdout, OUT_BLEND_DISTPART);
    distribPart(symbmtx, simuctrl, ctrl, dofptr);

    if(ctrl->option->timer)
        {
            clockStop(&timer_current);
            printf("--Distribution computed at time: %g --\n", clockVal(&timer_current));
        }

    if(ctrl->option->assembly)
      {
        /** Gener the Assembly structures **/
        if(ctrl->option->timer)
          {
            clockInit(&timer_current);
            clockStart(&timer_current);
          }

        if(ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
          {
            fprintf(stdout, "%ld : Genering assembly structure \n", (long)clustnum);
          }

#ifdef OLD_ASSEMBLY
        assemblyGener(assemb, ctrl->procnbr, symbmtx, simuctrl->ownetab);
#else
        assemblyGener(clustnum, assemb1D, assemb2D, clustnbr, symbmtx, simuctrl->blprtab, ctrl, dofptr);
#endif

        if(ctrl->option->timer)
          {
            clockStop(&timer_current);
            printf("--Assembly computed at time: %g --\n", clockVal(&timer_current));
          }
      }

#ifdef PASTIX_DYNSCHED /* 2 eme passe de splitpart */

    if(ctrl->option->timer)
        {
            clockInit(&timer_current);
            clockStart(&timer_current);
        }

    /** repartitioning of the initial symbolic factorization
        and processing of candidate processors group for
        each colum bloc **/

    splitPartLocal(ctrl, simuctrl, symbmtx, dofptr);

    if(ctrl->option->timer)
      {
        clockStop(&timer_current);
        printf("--Split build at time: %g --\n", clockVal(&timer_current));
      }

#endif

    /** Free some memory **/
    costExit(ctrl->costmtx);
    treeExit(ctrl->etree);

    if(ctrl->option->ps)
        ps_close(ps_file);

    /** gener the final solverMarix for this processor
      i.e. relative bloc numbering **/
    if(ctrl->option->timer)
        {
            clockInit(&timer_current);
            clockStart(&timer_current);
        }

    if(ctrl->option->sequentiel)
      {
        if(ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
          fprintf(stdout, "%ld : Genering all SolverMatrix files \n", (long)clustnum);
        allSolverMatrixSave(ctrl->option->solvmtx_filename,  symbmtx, simuctrl, ctrl, dofptr);
      }
    else
      {
        if(ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
          fprintf(stdout, "%ld : Genering final SolverMatrix \n", (long)clustnum);
        if (ctrl->option->iparm[IPARM_DOF_COST] != 0) {
          Dof dofstr;
          dofInit(&dofstr);
          dofConstant(&dofstr, 0, symbmtx->nodenbr, ctrl->option->iparm[IPARM_DOF_NBR]);
          bcofind = solverMatrixGen(ctrl->clustnum, solvmtx, symbmtx, simuctrl, ctrl, &dofstr);
          dofExit(&dofstr);
        }
        else{
          bcofind = solverMatrixGen(ctrl->clustnum, solvmtx, symbmtx, simuctrl, ctrl, dofptr);
        }
      }

    if(ctrl->option->timer)
        {
            clockStop(&timer_current);
            printf("--SolverMatrix computed at time: %g --\n", clockVal(&timer_current));
        }

    /*if(ctrl->option->count_ops)
      {
        printSolverInfo(stderr, solvmtx, dofptr);
      }*/

    if( ctrl->option->leader == clustnum &&
        ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
      fprintf(stdout, "** End of Partition & Distribution phase ** \n");

    /** Time end **/
    if(ctrl->option->timer)
      {
        clockStop(&timer_all);
        printf("---- Total execution at time: %g ----\n",clockVal(&timer_all));
        set_dparm(option->dparm, DPARM_ANALYZE_TIME, clockVal(&timer_all));
      }

    /** Free allocated memory **/
    simuExit(simuctrl, ctrl->clustnbr, ctrl->procnbr, ctrl->bublnbr);

    /*costExit(ctrl->costmtx);*/
    egraphExit(ctrl->egraph);

    if(ctrl->option->debug)
      {
        if(!ctrl->option->sequentiel)
          setBcofPtr(solvmtx, bcofind);
        if( ctrl->option->leader == clustnum &&
            ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
          fprintf(stdout, OUT_BLEND_CHKSOLVER);
        solverCheck(solvmtx);
      }

    if(!ctrl->option->sequentiel)
      {
        /***************************************
         * Malt : Moderate AmaLgamaTion        *
         ***************************************/
        if( ctrl->option->iparm[IPARM_VERBOSE]>API_VERBOSE_NO &&
            ctrl->option->malt_limit >= 0)
          {
            fprintf(stdout, "** Malt  phase ** \n");
            fprintf(stdout, "** Attemp to reduce  AUB memmory to a limit of %ld percents ** \n", (long)ctrl->option->malt_limit);
          }
        if(ctrl->option->malt_limit >= 0)
          {
            PASTIX_INT maxalloc;
            /*	maxalloc = Malt(solvmtx, INTVALMAX);
                fprintf(stderr, "Maxalloc for AUB = %ld\n", (long)maxalloc);*/
            maxalloc = Malt2(solvmtx, (float)ctrl->option->malt_limit);
            fprintf(stdout, "Max of Memory allocation %ld Bytes\n", (long)maxalloc);
            if(ctrl->option->debug)
              solverCheck(solvmtx);
          }
      }


    /***************************************
     * Realloc Memory in a contiguous way  *
     ***************************************/
    if(!ctrl->option->sequentiel)
      {
        blendCtrlExit(ctrl);
        printf("Contiguous reallocation of the solverMatrix ...\n");
        solverRealloc(solvmtx, bcofind);
        printf("Done \n");
      }
    else
      {
        /** Set the bcofptr in the block targets **/
        /** It is already done in solverRealloc if option sequentiel is active **/
        setBcofPtr(solvmtx, bcofind);

        /** Set the local btagptr **/
        setLocalBtagPtr(solvmtx);
      }

#ifdef DEBUG_BLEND
    if (leader == clustnum)
      fprintf(stdout, OUT_BLEND_CHKSOLVER);
    if (option->ricar) {
      if (leader == clustnum)
        errorPrintW("No solverMatrix checking in incomplete factorisation.");
    }else {
      solverCheck(solvmtx);
    }
#endif
    if (bcofind != NULL)
      memFree_null(bcofind);
}
