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
  File: redefine_functions.h

  This file redefine function names, adding _PASTIX_
  in prefix to all internal function.

  Authors:
     Xavier Lacoste - lacoste@labri.fr
*/
#ifndef REDEFINE_FUNCTIONS_H
#define REDEFINE_FUNCTIONS_H

#ifndef WITH_TYPE_PREFIX
#define NO_TYPE_PREFIX
#endif

#ifdef PASTIX_RENAME
#ifdef NO_TYPE_PREFIX
#define PASTIX_PREFIX_F(x)   _PASTIX_ ## x
#define PASTIX_EXTERN_F(x)   x
#else
#ifdef PREC_DOUBLE
#ifdef TYPE_COMPLEX
#define PASTIX_PREFIX_F(x)   _PASTIX_Z_ ## x
#define PASTIX_EXTERN_F(x)   z_ ## x
#else
#define PASTIX_PREFIX_F(x)   _PASTIX_D_ ## x
#define PASTIX_EXTERN_F(x)   d_ ## x
#endif /* TYPE_COMPLEX */
#else  /* PREC_DOUBLE */
#ifdef TYPE_COMPLEX
#define PASTIX_PREFIX_F(x)   _PASTIX_C_ ## x
#define PASTIX_EXTERN_F(x)   c_ ## x
#else
#define PASTIX_PREFIX_F(x)   _PASTIX_S_ ## x
#define PASTIX_EXTERN_F(x)   s_ ## x
#endif /* TYPE_COMPLEX */
#endif /* PREC_DOUBLE  */
#endif
#define PASTIX_PREFIX(x)     _PASTIX_ ## x
#define PASTIX_EXTERN(x)     x

#else  /* PASTIX_RENAME */

#ifdef NO_TYPE_PREFIX
#define PASTIX_PREFIX_F(x)   x
#define PASTIX_EXTERN_F(x)   x
#else
#ifdef PREC_DOUBLE
#ifdef TYPE_COMPLEX
#define PASTIX_PREFIX_F(x)   Z_ ## x
#define PASTIX_EXTERN_F(x)   z_ ## x
#else
#define PASTIX_PREFIX_F(x)   D_ ## x
#define PASTIX_EXTERN_F(x)   d_ ## x
#endif /* TYPE_COMPLEX */
#else  /* PREC_DOUBLE */
#ifdef TYPE_COMPLEX
#define PASTIX_PREFIX_F(x)   C_ ## x
#define PASTIX_EXTERN_F(x)   c_ ## x
#else
#define PASTIX_PREFIX_F(x)   S_ ## x
#define PASTIX_EXTERN_F(x)   s_ ## x
#endif /* TYPE_COMPLEX */
#endif /* PREC_DOUBLE  */
#endif
#define PASTIX_PREFIX(x)     x
#define PASTIX_EXTERN(x)     x


#endif /* PASTIX_RENAME */

/* Common functions */
#define api_dumparm        PASTIX_PREFIX_F(api_dumparm)
#define clockGet           PASTIX_PREFIX_F(clockGet)
#define errorPrint         PASTIX_PREFIX(errorPrint)
#define errorPrintW        PASTIX_PREFIX(errorPrintW)
#define errorProg          PASTIX_PREFIX(errorProg)
#define intAscn            PASTIX_PREFIX(intAscn)
#define intLoad            PASTIX_PREFIX(intLoad)
#define intPerm            PASTIX_PREFIX(intPerm)
#define intRandInit        PASTIX_PREFIX(intRandInit)
#define intRandReset       PASTIX_PREFIX(intRandReset)
#define intSave            PASTIX_PREFIX(intSave)
#define intSort1asc1       PASTIX_PREFIX(intSort1asc1)
#define intSort2asc1       PASTIX_PREFIX(intSort2asc1)
#define intSort2asc2       PASTIX_PREFIX(intSort2asc2)
#define intSort3asc1       PASTIX_PREFIX(intSort3asc1)
#define memAllocGroup      PASTIX_PREFIX(memAllocGroup)
#define memOffset          PASTIX_PREFIX(memOffset)
#define memReallocGroup    PASTIX_PREFIX(memReallocGroup)
#define memAllocGetCurrent PASTIX_PREFIX(memAllocGetCurrent)
#define memAllocGetMax     PASTIX_PREFIX(memAllocGetMax)
#define memAllocTraceReset PASTIX_PREFIX(memAllocTraceReset)
#define memAlloc_func      PASTIX_PREFIX(memAlloc_func)
#ifdef MEMORY_USAGE
#define memFree            PASTIX_PREFIX(memFree)
#define memRealloc_func    PASTIX_PREFIX(memRealloc_func)
#endif
#define qsortIntFloatAsc   PASTIX_PREFIX_F(qsortIntFloatAsc)
#define qsort2IntFloatAsc  PASTIX_PREFIX_F(qsort2IntFloatAsc)
#define qsort2IntAsc       PASTIX_PREFIX_F(qsort2IntAsc)
#define qsort2SmallIntAsc  PASTIX_PREFIX_F(qsort2SmallIntAsc)
#define set_dparm          PASTIX_PREFIX_F(set_dparm)
#define set_iparm          PASTIX_PREFIX_F(set_iparm)
#define usagePrint         PASTIX_PREFIX_F(usagePrint)
#define sortColRowValueAsc PASTIX_PREFIX(sortColRowValueAsc)
#define sortColRowAsc      PASTIX_PREFIX(sortColRowAsc)

/* Blend funcions */
#define allSolverMatrixSave     PASTIX_PREFIX(allSolverMatrixSave)
#define assemblyGener           PASTIX_PREFIX(assemblyGener)
#define blendCtrlExit           PASTIX_PREFIX(blendCtrlExit)
#define blendCtrlInit           PASTIX_PREFIX(blendCtrlInit)
#define blendParamExit          PASTIX_PREFIX(blendParamExit)
#define blendParamInit          PASTIX_PREFIX(blendParamInit)
#define blokUpdateCost          PASTIX_PREFIX(blokUpdateCost)
#define Bubble_Add              PASTIX_PREFIX(Bubble_Add)
#define Bubble_buildtree        PASTIX_PREFIX(Bubble_buildtree)
#define Bubble_Free             PASTIX_PREFIX(Bubble_Free)
#define Bubble_InitTree         PASTIX_PREFIX(Bubble_InitTree)
#define build_cblk              PASTIX_PREFIX(build_cblk)
#define build_smx               PASTIX_PREFIX(build_smx)
#define cblkComputeCost         PASTIX_PREFIX(cblkComputeCost)
#define cblkComputeCost2D       PASTIX_PREFIX(cblkComputeCost2D)
#define cblkCost                PASTIX_PREFIX(cblkCost)
#define cblkMaxCost             PASTIX_PREFIX(cblkMaxCost)
#define cblkNbr                 PASTIX_PREFIX(cblkNbr)
#define check_candidat          PASTIX_PREFIX(check_candidat)
#define cholesky                PASTIX_PREFIX(cholesky)
#define comp_int                PASTIX_PREFIX(comp_int)
#define compTimer               PASTIX_PREFIX(compTimer)
#define computeBlockCtrbNbr     PASTIX_PREFIX(computeBlockCtrbNbr)
#define computeCost             PASTIX_PREFIX(computeCost)
#define computeTaskReceiveTime  PASTIX_PREFIX(computeTaskReceiveTime)
#define compWith2keys           PASTIX_PREFIX(compWith2keys)
#define contribAddCost          PASTIX_PREFIX(contribAddCost)
#define contribCompCost         PASTIX_PREFIX(contribCompCost)
#define costExit                PASTIX_PREFIX(costExit)
#define costFtgtAdd             PASTIX_PREFIX(costFtgtAdd)
#define costFtgtSend            PASTIX_PREFIX(costFtgtSend)
#define costInit                PASTIX_PREFIX(costInit)
#define costMatrixBuild         PASTIX_PREFIX(costMatrixBuild)
#define costMatrixCorrect       PASTIX_PREFIX(costMatrixCorrect)
#define countBlok               PASTIX_PREFIX(countBlok)
#define crout_2t                PASTIX_PREFIX(crout_2t)
#define crout_3t                PASTIX_PREFIX(crout_3t)
#define crout_blok              PASTIX_PREFIX(crout_blok)
#define crout_hyb               PASTIX_PREFIX(crout_hyb)
#define DIAGCost                PASTIX_PREFIX(DIAGCost)
#define distribPart             PASTIX_PREFIX(distribPart)
#define E1Cost                  PASTIX_PREFIX(E1Cost)
#define E2Cost                  PASTIX_PREFIX(E2Cost)
#define egraphExit              PASTIX_PREFIX(egraphExit)
#define egraphInit              PASTIX_PREFIX(egraphInit)
#define eliminGraphBuild        PASTIX_PREFIX(eliminGraphBuild)
#define eliminTreeBuild         PASTIX_PREFIX(eliminTreeBuild)
#define etreeBuild              PASTIX_PREFIX(etreeBuild)
#define extendint_Add           PASTIX_PREFIX(extendint_Add)
#define extendint_Clear         PASTIX_PREFIX(extendint_Clear)
#define extendint_Exit          PASTIX_PREFIX(extendint_Exit)
#define extendint_incr          PASTIX_PREFIX(extendint_incr)
#define extendint_Init          PASTIX_PREFIX(extendint_Init)
#define extendint_Read          PASTIX_PREFIX(extendint_Read)
#define extendint_Size          PASTIX_PREFIX(extendint_Size)
#define extendint_ToSize        PASTIX_PREFIX(extendint_ToSize)
#define extracostExit           PASTIX_PREFIX(extracostExit)
#define extracostInit           PASTIX_PREFIX(extracostInit)
#define extra_inc_blok          PASTIX_PREFIX(extra_inc_blok)
#define extra_inc_cblk          PASTIX_PREFIX(extra_inc_cblk)
#define extrasymbolExit         PASTIX_PREFIX(extrasymbolExit)
#define extrasymbolInit         PASTIX_PREFIX(extrasymbolInit)
#define getFaceBlockE2          PASTIX_PREFIX(getFaceBlockE2)
#define getFtgtInd2             PASTIX_PREFIX(getFtgtInd2)
#define getFtgtNextAccess       PASTIX_PREFIX(getFtgtNextAccess)
#define getNextProc             PASTIX_PREFIX(getNextProc)
#define getNextTaskNextProc     PASTIX_PREFIX(getNextTaskNextProc)
#define getTaskUnmapped         PASTIX_PREFIX(getTaskUnmapped)
#define hazard                  PASTIX_PREFIX(hazard)
#define Malt2                   PASTIX_PREFIX(Malt2)
#define maxProcCost             PASTIX_PREFIX(maxProcCost)
#define memorySpaceCost         PASTIX_PREFIX(memorySpaceCost)
#define nnz                     PASTIX_PREFIX(nnz)
#define nodeTreeLevel           PASTIX_PREFIX(nodeTreeLevel)
#define P1D                     PASTIX_PREFIX(P1D)
#define P2D                     PASTIX_PREFIX(P2D)
#define partBuild               PASTIX_PREFIX(partBuild)
#define perfcluster2            PASTIX_PREFIX(perfcluster2)
#define PERF_limit_2D           PASTIX_PREFIX(PERF_limit_2D)
#define printSolverInfo         PASTIX_PREFIX(printSolverInfo)
#define printSymbolMatrix       PASTIX_PREFIX(printSymbolMatrix)
#define printTree               PASTIX_PREFIX(printTree)
#define propMappSubtree         PASTIX_PREFIX(propMappSubtree)
#define propMappSubtreeNC       PASTIX_PREFIX(propMappSubtreeNC)
#define propMappSubtreeNoSplit  PASTIX_PREFIX(propMappSubtreeNoSplit)
#define propMappSubtreeOn1P     PASTIX_PREFIX(propMappSubtreeOn1P)
#define propMappTree            PASTIX_PREFIX(propMappTree)
#define propMappTreeNoSplit     PASTIX_PREFIX(propMappTreeNoSplit)
#define ps_close                PASTIX_PREFIX(ps_close)
#define ps_draw_node_num        PASTIX_PREFIX(ps_draw_node_num)
#define ps_draw_node_owner      PASTIX_PREFIX(ps_draw_node_owner)
#define ps_open                 PASTIX_PREFIX(ps_open)
#define ps_rec_write_tree       PASTIX_PREFIX(ps_rec_write_tree)
#define ps_rec_write_tree_owner PASTIX_PREFIX(ps_rec_write_tree_owner)
#define ps_write_matrix         PASTIX_PREFIX(ps_write_matrix)
#define ps_write_tree           PASTIX_PREFIX(ps_write_tree)
#define ps_write_tree_owner     PASTIX_PREFIX(ps_write_tree_owner)
#define putInReadyQueue         PASTIX_PREFIX(putInReadyQueue)
#define queueAdd                PASTIX_PREFIX(queueAdd)
#define queueAdd2               PASTIX_PREFIX(queueAdd2)
#define queueClear              PASTIX_PREFIX(queueClear)
#define queueCopy               PASTIX_PREFIX(queueCopy)
#define queueExit               PASTIX_PREFIX(queueExit)
#define queueGet                PASTIX_PREFIX(queueGet)
#define queueGet2               PASTIX_PREFIX(queueGet2)
#define queueInit               PASTIX_PREFIX(queueInit)
#define queuePossess            PASTIX_PREFIX(queuePossess)
#define queuePrint              PASTIX_PREFIX(queuePrint)
#define queueRead               PASTIX_PREFIX(queueRead)
#define queueReorder            PASTIX_PREFIX(queueReorder)
#define queueSize               PASTIX_PREFIX(queueSize)
#define recursive_sum           PASTIX_PREFIX(recursive_sum)
#define setBcofPtr              PASTIX_PREFIX(setBcofPtr)
#define setDistribType          PASTIX_PREFIX(setDistribType)
#define setSubtreeBlokNbr       PASTIX_PREFIX(setSubtreeBlokNbr)
#define setSubtreeCostLevel     PASTIX_PREFIX(setSubtreeCostLevel)
#define setSubtreeDistribType   PASTIX_PREFIX(setSubtreeDistribType)
#define setSubtreeLevel         PASTIX_PREFIX(setSubtreeLevel)
#define setTreeCostLevel        PASTIX_PREFIX(setTreeCostLevel)
#define setTreeLevel            PASTIX_PREFIX(setTreeLevel)
#define simuExit                PASTIX_PREFIX(simuExit)
#define simuInit                PASTIX_PREFIX(simuInit)
#define simuRealloc             PASTIX_PREFIX(simuRealloc)
#define solverBlend             PASTIX_PREFIX(solverBlend)
#define solverCheck             PASTIX_PREFIX(solverCheck)
#define solverExit              PASTIX_PREFIX(solverExit)
#define solverInit              PASTIX_PREFIX(solverInit)
#define solverLoad              PASTIX_PREFIX(solverLoad)
#define solverMatrixGen         PASTIX_PREFIX(solverMatrixGen)
#define solverRealloc           PASTIX_PREFIX(solverRealloc)
#define solverSave              PASTIX_PREFIX(solverSave)
#define solverSpaceCost         PASTIX_PREFIX(solverSpaceCost)
#define splitCblk               PASTIX_PREFIX(splitCblk)
#define splitOnProcs            PASTIX_PREFIX(splitOnProcs)
#define splitPart               PASTIX_PREFIX(splitPart)
#define splitSeqCblk2D          PASTIX_PREFIX(splitSeqCblk2D)
#define subtreeSetCand          PASTIX_PREFIX(subtreeSetCand)
#define subtreeUpdateCost       PASTIX_PREFIX(subtreeUpdateCost)
#define subtreeUpdateCostLocal  PASTIX_PREFIX(subtreeUpdateCostLocal)
#define symbCost                PASTIX_PREFIX(symbCost)
#define symbolGener             PASTIX_PREFIX(symbolGener)
#define symbolRand              PASTIX_PREFIX(symbolRand)
#define symbolSpaceCost         PASTIX_PREFIX(symbolSpaceCost)
#define taskBuild               PASTIX_PREFIX(taskBuild)
#define taskExec_COMP1D         PASTIX_PREFIX(taskExec_COMP1D)
#define taskExec_DIAG           PASTIX_PREFIX(taskExec_DIAG)
#define taskExec_E1             PASTIX_PREFIX(taskExec_E1)
#define taskExec_E2             PASTIX_PREFIX(taskExec_E2)
#define taskSendCost            PASTIX_PREFIX(taskSendCost)
#define timerAdd                PASTIX_PREFIX(timerAdd)
#define timerSet                PASTIX_PREFIX(timerSet)
#define timerVal                PASTIX_PREFIX(timerVal)
#define totalCost               PASTIX_PREFIX(totalCost)
#define treeExit                PASTIX_PREFIX(treeExit)
#define treeInit                PASTIX_PREFIX(treeInit)
#define treeLeaveNbr            PASTIX_PREFIX(treeLeaveNbr)
#define treeLevel               PASTIX_PREFIX(treeLevel)
#define updateFtgtStruct        PASTIX_PREFIX(updateFtgtStruct)
#define virtualSplit            PASTIX_PREFIX(virtualSplit)
#define dofConstant             PASTIX_PREFIX(dofConstant)
#define dofExit                 PASTIX_PREFIX(dofExit)
#define dofGraph                PASTIX_PREFIX(dofGraph)
#define dofInit                 PASTIX_PREFIX(dofInit)
#define dofSave                 PASTIX_PREFIX(dofSave)
#define treePlot                PASTIX_PREFIX(treePlot)

/* Fax functions */
#define symbolCompact   PASTIX_PREFIX(symbolCompact)
#define symbolCosti     PASTIX_PREFIX(symbolCosti)
#define symbolFax       PASTIX_PREFIX(symbolFax)
#define symbolFaxDgraph PASTIX_PREFIX(symbolFaxDgraph)
#define symbolFaxGraph  PASTIX_PREFIX(symbolFaxGraph)
#define symbolFaxi      PASTIX_PREFIX(symbolFaxi)
#define symbolFaxiGraph PASTIX_PREFIX(symbolFaxiGraph)

/* kass functions */
#define amalgamate               PASTIX_PREFIX(amalgamate)
#define Build_SymbolMatrix       PASTIX_PREFIX(Build_SymbolMatrix)
#define cblk_time_fact           PASTIX_PREFIX(cblk_time_fact)
#define cblk_time_solve          PASTIX_PREFIX(cblk_time_solve)
#define cleanCS                  PASTIX_PREFIX(cleanCS)
#define col_match                PASTIX_PREFIX(col_match)
#define col_merge                PASTIX_PREFIX(col_merge)
#define compact_graph            PASTIX_PREFIX(compact_graph)
#define compute_elimination_tree PASTIX_PREFIX(compute_elimination_tree)
#define compute_subtree_size     PASTIX_PREFIX(compute_subtree_size)
#define CSnnz                    PASTIX_PREFIX(CSnnz)
#define CS_Perm                  PASTIX_PREFIX(CS_Perm)
#define CS_RowPerm               PASTIX_PREFIX(CS_RowPerm)
#define find_supernodes          PASTIX_PREFIX(find_supernodes)
#define get_son                  PASTIX_PREFIX(get_son)
#define ifax                     PASTIX_PREFIX(ifax)
#define initCS                   PASTIX_PREFIX(initCS)
#define is_col_match             PASTIX_PREFIX(is_col_match)
#define kass                     PASTIX_PREFIX(kass)
#define kass_symbol              PASTIX_PREFIX(kass_symbol)
#define KSupernodes              PASTIX_PREFIX(KSupernodes)
#define merge_col                PASTIX_PREFIX(merge_col)
#define merge_cost               PASTIX_PREFIX(merge_cost)
#define Merge_Cost               PASTIX_PREFIX(Merge_Cost)
#define merge_gain               PASTIX_PREFIX(merge_gain)
#define Patch_SymbolMatrix       PASTIX_PREFIX(Patch_SymbolMatrix)
#define post_order               PASTIX_PREFIX(post_order)
#define SF_Direct                PASTIX_PREFIX(SF_Direct)
#define SF_level                 PASTIX_PREFIX(SF_level)
#define sort_row                 PASTIX_PREFIX_F(sort_row)
#define UnionSet                 PASTIX_PREFIX(UnionSet)


/* sopalin functions */
#define Matrix_Unscale_Unsym       PASTIX_PREFIX_F(Matrix_Unscale_Unsym)
#define Matrix_Unscale_Sym         PASTIX_PREFIX_F(Matrix_Unscale_Sym)
#define SolverMatrix_Unscale_Unsym PASTIX_PREFIX_F(SolverMatrix_Unscale_Unsym)
#define SolverMatrix_Unscale_Sym   PASTIX_PREFIX_F(SolverMatrix_Unscale_Sym)
#define CscMatrix_Unscale_Unsym    PASTIX_PREFIX_F(CscMatrix_Unscale_Unsym)
#define CscMatrix_Unscale_Sym      PASTIX_PREFIX_F(CscMatrix_Unscale_Sym)
#define SolverMatrix_DiagMult      PASTIX_PREFIX_F(SolverMatrix_DiagMult)
#define SolverMatrix_ColMult       PASTIX_PREFIX_F(SolverMatrix_ColMult)
#define SolverMatrix_RowMult       PASTIX_PREFIX_F(SolverMatrix_RowMult)
#define CscMatrix_ColMult          PASTIX_PREFIX_F(CscMatrix_ColMult)
#define CscMatrix_RowMult          PASTIX_PREFIX_F(CscMatrix_RowMult)


#define CSC_isolate                     PASTIX_PREFIX_F(CSC_isolate)
#define CSC_colPerm                     PASTIX_PREFIX_F(CSC_colPerm)
#define CSC_colScale                    PASTIX_PREFIX_F(CSC_colScale)
#define CSC_Fnum2Cnum                   PASTIX_PREFIX_F(CSC_Fnum2Cnum)
#define CSC_Cnum2Fnum                   PASTIX_PREFIX_F(CSC_Cnum2Fnum)
#define CSC_buildZerosAndNonZerosGraphs PASTIX_PREFIX_F(CSC_buildZerosAndNonZerosGraphs)
#define CSC_rowScale                    PASTIX_PREFIX_F(CSC_rowScale)
#define CSC_sort                        PASTIX_PREFIX_F(CSC_sort)
#define add_two_floats                  PASTIX_PREFIX_F(add_two_floats)
#define cmp_colrow                      PASTIX_PREFIX_F(cmp_colrow)
#define csc_checksym                    PASTIX_PREFIX_F(csc_checksym)
#define csc_cyclic_distribution         PASTIX_PREFIX_F(csc_cyclic_distribution)
#define cscd_addlocal                   PASTIX_EXTERN_F(cscd_addlocal)
#define cscd_addlocal_int               PASTIX_PREFIX_F(cscd_addlocal_int)
#define cscd_checksym                   PASTIX_PREFIX_F(cscd_checksym)
#define cscd_load                       PASTIX_EXTERN_F(cscd_load)
#define cscd_noDiag                     PASTIX_PREFIX_F(cscd_noDiag)
#define cscd_redispatch_int             PASTIX_PREFIX_F(cscd_redispatch_int)
#define cscd_redispatch_scotch          PASTIX_PREFIX_F(cscd_redispatch_scotch)
#define cscd_save                       PASTIX_EXTERN_F(cscd_save)
#define csc_save                        PASTIX_EXTERN_F(csc_save)
#define csc_load                        PASTIX_EXTERN_F(csc_load)
#define cscd_symgraph                   PASTIX_PREFIX_F(cscd_symgraph)
#define cscd_symgraph_int               PASTIX_PREFIX_F(cscd_symgraph_int)
#define csc_build_g2l                   PASTIX_PREFIX_F(csc_build_g2l)
#define csc_noDiag                      PASTIX_PREFIX_F(csc_noDiag)
#define csc_symgraph                    PASTIX_PREFIX_F(csc_symgraph)
#define csc_check_doubles               PASTIX_PREFIX_F(csc_check_doubles)
#define csc_simple_distribution         PASTIX_PREFIX_F(csc_simple_distribution)
#define csc_symgraph_int                PASTIX_PREFIX_F(csc_symgraph_int)
#define get_max                         PASTIX_PREFIX_F(get_max)
#define get_min                         PASTIX_PREFIX_F(get_min)
#define keep_first                      PASTIX_PREFIX_F(keep_first)
#define keep_last                       PASTIX_PREFIX_F(keep_last)
#define dim_dgeam                       PASTIX_PREFIX_F(dim_dgeam)
#define bordi                           PASTIX_PREFIX(bordi)
#define buildGlob2loc                   PASTIX_PREFIX_F(buildGlob2loc)
#define buildUpdoVect                   PASTIX_PREFIX_F(buildUpdoVect)
#define cmpint                          PASTIX_PREFIX_F(cmpint)
#define CoefMatrix_Allocate             PASTIX_PREFIX_F(CoefMatrix_Allocate)
#define CoefMatrix_Free                 PASTIX_PREFIX_F(CoefMatrix_Free)
#define CoefMatrix_Init                 PASTIX_PREFIX_F(CoefMatrix_Init)
#define correct2                        PASTIX_PREFIX_F(correct2)
#define Csc2solv_cblk                   PASTIX_PREFIX_F(Csc2solv_cblk)
#define Csc2updown                      PASTIX_PREFIX_F(Csc2updown)
#define Csc2updown_X0                   PASTIX_PREFIX_F(Csc2updown_X0)
#define CscBLoad                        PASTIX_PREFIX_F(CscBLoad)
#define CscBSave                        PASTIX_PREFIX_F(CscBSave)
#define CscdOrdistrib                   PASTIX_PREFIX_F(CscdOrdistrib)
#define CscdRhsUpdown                   PASTIX_PREFIX_F(CscdRhsUpdown)
#define CscdUpdownRhs                   PASTIX_PREFIX_F(CscdUpdownRhs)
#define CscExit                         PASTIX_PREFIX_F(CscExit)
#define CscLoad                         PASTIX_PREFIX_F(CscLoad)
#define CscOrdistrib                    PASTIX_PREFIX_F(CscOrdistrib)
#define CscRhsUpdown                    PASTIX_PREFIX_F(CscRhsUpdown)
#define CscSave                         PASTIX_PREFIX_F(CscSave)
#define CscUpdownRhs                    PASTIX_PREFIX_F(CscUpdownRhs)
#define dump1                           PASTIX_PREFIX_F(dump1)
#define dump2                           PASTIX_PREFIX_F(dump2)
#define dump3                           PASTIX_PREFIX_F(dump3)
#define dump3_LU                        PASTIX_PREFIX_F(dump3_LU)
#define dump4                           PASTIX_PREFIX_F(dump4)
#define dump5                           PASTIX_PREFIX_F(dump5)
#define dump6                           PASTIX_PREFIX_F(dump6)
#define dump7                           PASTIX_PREFIX_F(dump7)

#define GetMpiSum        PASTIX_PREFIX_F(GetMpiSum)
#define FreeMpiSum       PASTIX_PREFIX_F(FreeMpiSum)
#define GetMpiType       PASTIX_PREFIX_F(GetMpiType)
#define FreeMpiType      PASTIX_PREFIX_F(FreeMpiType)
#define global2localperm PASTIX_PREFIX_F(global2localperm)
#define global2localrhs  PASTIX_PREFIX_F(global2localrhs)
#define mysum            PASTIX_PREFIX_F(mysum)
#define orderSplit       PASTIX_PREFIX(orderSplit)
#define orderSplit2      PASTIX_PREFIX(orderSplit2)
#define orderSplit3      PASTIX_PREFIX(orderSplit3)

#define redispatch_rhs        PASTIX_PREFIX_F(redispatch_rhs)
#define sizeofsolver          PASTIX_PREFIX_F(sizeofsolver)
#define solve_backup          PASTIX_PREFIX_F(solve_backup)
#define solve_restore         PASTIX_PREFIX_F(solve_restore)
#define sopalin_backup        PASTIX_PREFIX_F(sopalin_backup)
#define sopalin_bindthread    PASTIX_PREFIX(sopalin_bindthread)
#define sopalin_check_param   PASTIX_PREFIX_F(sopalin_check_param)
#define sopalin_clean         PASTIX_PREFIX_F(sopalin_clean)
#define sopalin_clean_smp     PASTIX_PREFIX_F(sopalin_clean_smp)
#define sopalin_init          PASTIX_PREFIX_F(sopalin_init)
#define sopalin_init_smp      PASTIX_PREFIX_F(sopalin_init_smp)
#define sopalin_launch_comm   PASTIX_PREFIX(sopalin_launch_comm)
#define sopalin_launch_thread PASTIX_PREFIX(sopalin_launch_thread)
#define sopalin_option        PASTIX_PREFIX_F(sopalin_option)
#define sopalin_restore       PASTIX_PREFIX_F(sopalin_restore)
#define symbolCostn           PASTIX_PREFIX(symbolCostn)
#define symbolRustine         PASTIX_PREFIX(symbolRustine)
#define symbolSplit           PASTIX_PREFIX(symbolSplit)

#define Usopalin_launch            PASTIX_PREFIX_F(Usopalin_launch)
#define Usopalin_thread            PASTIX_PREFIX_F(Usopalin_thread)
#define Usopalin_updo_thread       PASTIX_PREFIX_F(Usopalin_updo_thread)
#define Usopalin_updo_gmres_thread PASTIX_PREFIX_F(Usopalin_updo_gmres_thread)
#define Usopalin_updo_grad_thread  PASTIX_PREFIX_F(Usopalin_updo_grad_thread)
#define Usopalin_updo_pivot_thread PASTIX_PREFIX_F(Usopalin_updo_pivot_thread)
#define Uupdo_thread               PASTIX_PREFIX_F(Uupdo_thread)
#define Upivot_thread              PASTIX_PREFIX_F(Upivot_thread)
#define Ugmres_thread              PASTIX_PREFIX_F(Ugmres_thread)
#define Ugrad_thread               PASTIX_PREFIX_F(Ugrad_thread)

#define Dsopalin_launch            PASTIX_PREFIX_F(Dsopalin_launch)
#define Dsopalin_thread            PASTIX_PREFIX_F(Dsopalin_thread)
#define Dsopalin_updo_thread       PASTIX_PREFIX_F(Dsopalin_updo_thread)
#define Dsopalin_updo_gmres_thread PASTIX_PREFIX_F(Dsopalin_updo_gmres_thread)
#define Dsopalin_updo_grad_thread  PASTIX_PREFIX_F(Dsopalin_updo_grad_thread)
#define Dsopalin_updo_pivot_thread PASTIX_PREFIX_F(Dsopalin_updo_pivot_thread)
#define Dupdo_thread               PASTIX_PREFIX_F(Dupdo_thread)
#define Dpivot_thread              PASTIX_PREFIX_F(Dpivot_thread)
#define Dgmres_thread              PASTIX_PREFIX_F(Dgmres_thread)
#define Dgrad_thread               PASTIX_PREFIX_F(Dgrad_thread)

#define Lsopalin_launch            PASTIX_PREFIX_F(Lsopalin_launch)
#define Lsopalin_thread            PASTIX_PREFIX_F(Lsopalin_thread)
#define Lsopalin_updo_thread       PASTIX_PREFIX_F(Lsopalin_updo_thread)
#define Lsopalin_updo_gmres_thread PASTIX_PREFIX_F(Lsopalin_updo_gmres_thread)
#define Lsopalin_updo_grad_thread  PASTIX_PREFIX_F(Lsopalin_updo_grad_thread)
#define Lsopalin_updo_pivot_thread PASTIX_PREFIX_F(Lsopalin_updo_pivot_thread)
#define Lupdo_thread               PASTIX_PREFIX_F(Lupdo_thread)
#define Lpivot_thread              PASTIX_PREFIX_F(Lpivot_thread)
#define Lgmres_thread              PASTIX_PREFIX_F(Lgmres_thread)
#define Lgrad_thread               PASTIX_PREFIX_F(Lgrad_thread)

/* Symbol functions */

#define symbolBase        PASTIX_PREFIX(symbolBase)
#define symbolCheck       PASTIX_PREFIX(symbolCheck)
#define symbolCost        PASTIX_PREFIX(symbolCost)
#define symbolDraw        PASTIX_PREFIX(symbolDraw)
#define symbolDrawColor   PASTIX_PREFIX(symbolDrawColor)
#define symbolDrawFunc    PASTIX_PREFIX(symbolDrawFunc)
#define symbolExit        PASTIX_PREFIX(symbolExit)
#define symbolInit        PASTIX_PREFIX(symbolInit)
#define symbolKeepAdd     PASTIX_PREFIX(symbolKeepAdd)
#define symbolKeepCompute PASTIX_PREFIX(symbolKeepCompute)
#define symbolKeepDel     PASTIX_PREFIX(symbolKeepDel)
#define symbolKeepExit    PASTIX_PREFIX(symbolKeepExit)
#define symbolKeepHisto   PASTIX_PREFIX(symbolKeepHisto)
#define symbolKeepInit    PASTIX_PREFIX(symbolKeepInit)
#define symbolKeepPurge   PASTIX_PREFIX(symbolKeepPurge)
#define symbolKeepView    PASTIX_PREFIX(symbolKeepView)
#define symbolLevf        PASTIX_PREFIX(symbolLevf)
#define symbolLoad        PASTIX_PREFIX(symbolLoad)
#define symbolNonzeros    PASTIX_PREFIX(symbolNonzeros)
#define symbolRealloc     PASTIX_PREFIX(symbolRealloc)
#define symbolSave        PASTIX_PREFIX(symbolSave)
#define symbolTree        PASTIX_PREFIX(symbolTree)


#define csc2cscd                        PASTIX_EXTERN_F(csc2cscd)
#define csc_dispatch                    PASTIX_EXTERN_F(csc_dispatch)
#define cscd2csc                        PASTIX_EXTERN_F(cscd2csc)
#define cscd2csc_int                    PASTIX_EXTERN_F(cscd2csc_int)
#define cscd_build_g2l                  PASTIX_PREFIX_F(cscd_build_g2l)
#define cscd_redispatch                 PASTIX_EXTERN_F(cscd_redispatch)
#define dpastix                         PASTIX_EXTERN_F(dpastix)
#define dpastix_task_fax                PASTIX_PREFIX_F(dpastix_task_fax)
#define pastix                          PASTIX_EXTERN_F(pastix)
#define pastix_setSchurUnknownList      PASTIX_EXTERN_F(pastix_setSchurUnknownList)
#define pastix_getSchur                 PASTIX_EXTERN_F(pastix_getSchur)
#define pastix_getSchurLocalNodeNbr     PASTIX_EXTERN_F(pastix_getSchurLocalNodeNbr)
#define pastix_getSchurLocalUnkownNbr   PASTIX_EXTERN_F(pastix_getSchurLocalUnkownNbr)
#define pastix_getSchurLocalNodeList    PASTIX_EXTERN_F(pastix_getSchurLocalNodeList)
#define pastix_getSchurLocalUnknownList PASTIX_EXTERN_F(pastix_getSchurLocalUnknownList)
#define pastix_setSchurArray            PASTIX_EXTERN_F(pastix_setSchurArray)
#define pastix_checkMatrix              PASTIX_EXTERN_F(pastix_checkMatrix)
#define pastix_checkMatrix_int          PASTIX_PREFIX_F(pastix_checkMatrix_int)
#define pastix_fillin_csc               PASTIX_PREFIX_F(pastix_fillin_csc)
#define pastix_getLocalNodeLst          PASTIX_EXTERN_F(pastix_getLocalNodeLst)
#define pastix_getLocalNodeNbr          PASTIX_EXTERN_F(pastix_getLocalNodeNbr)
#define pastix_getLocalUnknownLst       PASTIX_EXTERN_F(pastix_getLocalUnknownLst)
#define pastix_getLocalUnknownNbr       PASTIX_EXTERN_F(pastix_getLocalUnknownNbr)
#define pastix_initParam                PASTIX_EXTERN_F(pastix_initParam)
#define pastix_print_memory_usage       PASTIX_PREFIX_F(pastix_print_memory_usage)
#define pastix_bindThreads              PASTIX_EXTERN_F(pastix_bindThreads)
#define pastix_unscale                  PASTIX_EXTERN_F(pastix_unscale)
#define pastix_task_blend               PASTIX_PREFIX_F(pastix_task_blend)
#define pastix_task_clean               PASTIX_PREFIX_F(pastix_task_clean)
#define pastix_task_fax                 PASTIX_PREFIX_F(pastix_task_fax)
#define pastix_task_init                PASTIX_PREFIX_F(pastix_task_init)
#define pastix_task_raff                PASTIX_PREFIX_F(pastix_task_raff)
#define pastix_task_scotch              PASTIX_PREFIX_F(pastix_task_scotch)
#define dpastix_task_scotch             PASTIX_PREFIX_F(dpastix_task_scotch)
#define pastix_task_sopalin             PASTIX_PREFIX_F(pastix_task_sopalin)
#define pastix_task_updown              PASTIX_PREFIX_F(pastix_task_updown)
#define pastix_welcome_print            PASTIX_PREFIX_F(pastix_welcome_print)
#define api_dparmreader                 PASTIX_EXTERN_F(api_dparmreader)
#define api_iparmreader                 PASTIX_EXTERN_F(api_iparmreader)
#define pastix_order_prepare_csc        PASTIX_PREFIX_F(pastix_order_prepare_csc)
#define dpastix_order_prepare_cscd      PASTIX_PREFIX_F(dpastix_order_prepare_cscd)
#define pastix_order_load               PASTIX_PREFIX_F(pastix_order_load)
#define pastix_order_save               PASTIX_PREFIX_F(pastix_order_save)
/* OOC */

#ifdef OOC
#define get_ooc_thread       PASTIX_EXTERN_F(get_ooc_thread)
#define ooc_init_thread      PASTIX_EXTERN_F(ooc_init_thread)
#define ooc_clock_get        PASTIX_EXTERN_F(ooc_clock_get)
#define ooc_do_save_coef     PASTIX_EXTERN_F(ooc_do_save_coef)
#define ooc_do_load_coef     PASTIX_EXTERN_F(ooc_do_load_coef)
#define ooc_allocate         PASTIX_EXTERN_F(ooc_allocate)
#ifdef OOC_FTGT
#define ooc_do_load_ftgt     PASTIX_EXTERN_F(ooc_do_load_ftgt)
#define ooc_remove_ftgt      PASTIX_EXTERN_F(ooc_remove_ftgt)
#define ooc_allocate_ftgt    PASTIX_EXTERN_F(ooc_allocate_ftgt)
#define ooc_do_save_ftgt     PASTIX_EXTERN_F(ooc_do_save_ftgt)
#endif /* OOC_FTGT */
#define reduceMem            PASTIX_EXTERN_F(reduceMem)
#define cblkNextAccess       PASTIX_EXTERN_F(cblkNextAccess)
#define cblkAndContribSize   PASTIX_EXTERN_F(cblkAndContribSize)
#define updo_init_mem        PASTIX_EXTERN_F(updo_init_mem)
#define updo_init_smp_mem    PASTIX_EXTERN_F(updo_init_smp_mem)
#define sopalin_init_smp_mem PASTIX_EXTERN_F(sopalin_init_smp_mem)
#define raff_mem             PASTIX_EXTERN_F(raff_mem)
#define recursive_mkdir      PASTIX_EXTERN_F(recursive_mkdir)
#define ooc_thread           PASTIX_EXTERN_F(ooc_thread)
#define ooc_init             PASTIX_EXTERN_F(ooc_init)
#define ooc_exit             PASTIX_EXTERN_F(ooc_exit)
#define ooc_stop_thread      PASTIX_EXTERN_F(ooc_stop_thread)
#define ooc_freeze           PASTIX_EXTERN_F(ooc_freeze)
#define ooc_defreeze         PASTIX_EXTERN_F(ooc_defreeze)
#define ooc_set_step         PASTIX_EXTERN_F(ooc_set_step)
#define ooc_wait_for_cblk    PASTIX_EXTERN_F(ooc_wait_for_cblk)
#define ooc_hack_load        PASTIX_EXTERN_F(ooc_hack_load)
#define ooc_save_coef        PASTIX_EXTERN_F(ooc_save_coef)
#define ooc_receiving        PASTIX_EXTERN_F(ooc_receiving)
#define ooc_received         PASTIX_EXTERN_F(ooc_received)
#define ooc_wait_task        PASTIX_EXTERN_F(ooc_wait_task)
#ifdef OOC_FTGT
#define ooc_wait_for_ftgt    PASTIX_EXTERN_F(ooc_wait_for_ftgt)
#define ooc_reset_ftgt       PASTIX_EXTERN_F(ooc_reset_ftgt)
#define ooc_save_ftgt        PASTIX_EXTERN_F(ooc_save_ftgt)
#endif /* OOC_FTGT */
#endif                          /* OOC */


#endif /* REDEFINE_FUNCTIONS_H */
