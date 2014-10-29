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
/**   NAME       : param_blend.h                           **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                These lines are the data declarations   **/
/**                for the parameters of blend.            **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 22 jul 1998     **/
/**                                 to     08 sep 1998     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ The parameters structure definition +*/

typedef struct BlendParam_ {
  char * hpf_filename;          /*+ file name for HPF distribution  +*/
  char * trace_filename;        /*+ file name for Paragraph traces  +*/
  char * ps_filename;           /*+ file name for matrix postscript +*/
  
  PASTIX_INT    hpf;             /*+ gener an HPF distribution file                                    +*/
  PASTIX_INT    tracegen ;       /*+ gener a simulated Paragraph execution trace file                  +*/
  PASTIX_INT    ps ;             /*+ gener a post-script of symbol matrix and elimination tree         +*/
  PASTIX_INT    assembly;        /*+ Gener the info structure needed to assemble                       +*/
  char * solvmtx_filename;/*+ name of solver matrix files (in sequential mode                   +*/ 
  PASTIX_INT    sequentiel;      /*+ Exec blend in sequentiel mode -> all solver matrix files generated+*/
  PASTIX_INT    count_ops ;      /*+ print costs in term of number of elementary operations            +*/
  PASTIX_INT    debug ;          /*+ make some check at certains execution points                      +*/
  PASTIX_INT    timer;           /*+ print execution time                                              +*/
  PASTIX_INT    recover;         /*+ take acount of a recover time estimation for ftgt                 +*/  
  PASTIX_INT    blcolmin ;       /*+ minimun number of column for a good use of BLAS primitives        +*/
  PASTIX_INT    blcolmax;
  PASTIX_INT    blblokmin ;      /*+ size of blockage for a good use of BLAS primitives  in 2D distribution +*/
  PASTIX_INT    blblokmax;
  PASTIX_INT    abs;             /*+ adaptative block size: := 0 all block are cut to blcolmin else try to make (ncand*abs) column +*/ 
  PASTIX_INT    leader;          /*+ Processor leader for not parallele task (ex: gener assembly1D     +*/
  PASTIX_INT    allcand;         /*+ All processor are candidat for a splitted cblk                    +*/
  PASTIX_INT    nocrossproc;     /*+ Crossing processor forbiden in the splitting phase                +*/
  PASTIX_INT    forceE2;
  PASTIX_INT    level2D;         /*+ number of level to treat with a 2D distribution                   +*/
  PASTIX_INT    candcorrect;   
  PASTIX_INT    clusterprop;     /*+ Proportionnal mapping with clustering for upper layers            +*/
  PASTIX_INT    costlevel;       /*+ Calcul du cout de chaque sous arbre dans candtab                  +*/
  PASTIX_INT    autolevel;       /*+ Level to shift 1D to 2D is automaticly computed                   +*/
  PASTIX_INT    malt_limit;      /*+ Limit for AUB memory allocations    (in octets)                   +*/
  PASTIX_INT    smpnbr;          /*+ Number of smp node                                                +*/
  PASTIX_INT    procnbr;         /*+ Number of physical processors in a smp node                       +*/
  double ratiolimit;  
  PASTIX_INT    dense_endblock;   /*+ Treat the square right lower part of the matrix as a dense matrix+*/  
  PASTIX_INT    ooc;              /*+ To use the out of core version of Pastix                         +*/
  PASTIX_INT    ricar;            /*+ If set to 1 then use blend dedicated to ricar                    +*/
  double oocmemlimit;      /*+ limit of physical memory for ooc                                 +*/
  PASTIX_INT   *iparm;            /*+ In/Out Integer parameters +*/
  double *dparm;           /*+ In/Out Float parameters   +*/
  PASTIX_INT     n;               /*+ Size of the matrix        +*/  
} BlendParam;




PASTIX_INT      blendParamInit(BlendParam *);
void     blendParamExit(BlendParam *);

