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
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <assert.h>

#include "common_pastix.h"
#include "dof.h"
#include "cost.h"
#include "symbol.h"
#include "elimin.h"
#include "extendVector.h"
#include "cand.h"
#include "param_blend.h"
#include "queue.h"
#include "bulles.h"
/* #include "extrastruct.h" */
#include "blendctrl.h"
#include "eliminfunc.h"
#include "ftgt.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "simu.h"
#include "costfunc.h"
#include "write_ps.h"


FILE *ps_open(char *filepath)
{
  FILE *out= fopen(filepath, "w");

  fprintf(out, "%%!PS-Adobe-2.0\n");
  fprintf(out, "%%%%Creator: rnd \n");
  fprintf(out, "%%%%Title: rnd \n");
  fprintf(out, "%%%%Pages: 6 \n");
  fprintf(out, "%%%%PageOrder: Ascend\n");
  fprintf(out, "%%%%EndComments\n");

  /* define some usefull things */
  fprintf(out, "/m {moveto} bind def\n");
  fprintf(out, "/rm {rmoveto} bind def\n");
  fprintf(out, "/sf {/Times-Roman findfont 2 1 roll scalefont setfont} def\n");
  fprintf(out, "/s {show} bind def\n");
  fprintf(out, "/l {lineto} bind def\n");
  fprintf(out, "/f {fill} bind def\n");
  fprintf(out, "/rl {rlineto} bind def\n");
  fprintf(out, "/np {newpath} bind def\n");
  fprintf(out, "/cp {closepath} bind def\n");
  fprintf(out, "/gs {gsave} bind def\n");
  fprintf(out, "/gr {grestore} bind def\n");
  fprintf(out, "/srgb {setrgbcolor} bind def\n");
  fprintf(out, "/sg {setgray} bind def\n");
  fprintf(out, "/sk {stroke} bind def\n");
  return (out);
}

void ps_close(FILE *out)
{
  fprintf(out, "%%%%Trailer\n");
  fprintf(out, "%%%%EOF\n");
  fclose(out);
}



/* imprime la matrice creuse en post-script */
void ps_write_matrix(SymbolMatrix *symbmtx, FILE *out, PASTIX_INT *page)
{
    PASTIX_INT p;
    PASTIX_INT i;
    double alpha, beta;
    double a, b;
    double s;
    double gray_owner;
    PASTIX_INT ncol;

    ncol= symbmtx->nodenbr;

    s= 500.0/ncol;

    fprintf(out, "%%%%Page: %ld %ld \n", (long)*page, (long)*page); (*page)++;
    fprintf(out, "%% Sparse Matrix \n");

    fprintf(out, "40 40 translate 1 1 scale\n");

    fprintf(out, "40 sf\n");
    fprintf(out, "%f %f m (%ldx%ld) s\n",
            (double)500.0*0.75, (double)500.0*0.75,
            (long)sqrt((double)(symbmtx->nodenbr)),
            (long)sqrt((double)(symbmtx->nodenbr)));


    for(i=0;i<symbmtx->cblknbr;i++)
        {
            alpha = (double)(symbmtx->cblktab[i].fcolnum);
            beta  = (double)(symbmtx->cblktab[i].lcolnum + 1);

            gray_owner= 0.5;


            fprintf(out, "np ");
            fprintf(out, "%f %f m ",alpha*s, 0.0);
            fprintf(out, "%f %f l ",beta*s, 0.0);
            fprintf(out, "%f %f l ",beta*s, (double)(ncol-beta)*s);
            fprintf(out, "%f %f l ",alpha*s, (double)(ncol-alpha)*s);

            fprintf(out, "cp gs 1 sg f gr sk\n");

            /* diag blok */
            fprintf(out, "np ");
            fprintf(out, "%f %f m ",beta*s, (ncol-beta)*s);
            fprintf(out, "%f %f l ",alpha*s, (ncol-alpha)*s);
            fprintf(out, "%f %f l ",alpha*s, (ncol-beta)*s);
            fprintf(out, "cp gs %f sg f gr sk\n",
                    gray_owner);

            /* blok num */
            fprintf(out, "%f sf\n", (beta-alpha)*s);
            fprintf(out, "%f %f m (%ld) s\n",(alpha+beta)*s/2, (ncol-(alpha+beta)/2)*s, (long)i);

            for (p = symbmtx->cblktab[i].bloknum+1 ; p < symbmtx->cblktab[i+1].bloknum ; p++)
                {
                    a= (double)(symbmtx->bloktab[p].frownum);
                    b= (double)(symbmtx->bloktab[p].lrownum + 1);

                    fprintf(out, "np ");
                    fprintf(out, "%f %f m ", alpha*s, (ncol-a)*s);
                    fprintf(out, "%f %f l ", beta*s, (ncol-a)*s);
                    fprintf(out, "%f %f l ", beta*s, (ncol-b)*s);
                    fprintf(out, "%f %f l ", alpha*s, (ncol-b)*s);
                    fprintf(out, "cp gs %f sg f gr sk\n", gray_owner);
                }

        }

    /* column number */
    fprintf(out, "%f sf\n", s*0.9);
    for (p= 1; p<=symbmtx->nodenbr; p++)
        {
            fprintf(out, "%f %f m (%ld) s\n",
                    p*s, -s, (long)p);
        }
    fprintf(out, "showpage\n");
    fflush(out);
}



/* imprime l'arbre d'elimination en post-script */
void ps_write_tree(const CostMatrix *costmtx, const EliminTree *etree, FILE *out, PASTIX_INT *page)
     /* PASTIX_INT (*ps_draw_node)(FILE *out, cbl c, cbls allcbls, double s, double x, double y) */
{
  fprintf(out, "%%%%Page: %ld %ld \n", (long)*page, (long)*page); (*page)++;
  fprintf(out, "%% Elimination Tree (color by Layer)\n");

  fprintf(out, "40 40 translate 1 1 scale\n");
  ps_rec_write_tree(ROOT(etree), costmtx, etree, out, ps_draw_node_num);
  fprintf(out, "showpage\n");

  fflush(out);

}


double ps_rec_write_tree(PASTIX_INT nodenum, const CostMatrix *costmtx, const EliminTree *etree, FILE *out,
                            void (*ps_draw_node)(FILE *out, PASTIX_INT nodenum, const CostMatrix *costmtx,
                                                const EliminTree *etree, double s,
                                                double x, double y))
{
  static PASTIX_INT cx;
  PASTIX_INT nodelevel;
  double pos=0, lpos=0, rpos=0;
  double s, sy;

  nodelevel = nodeTreeLevel(nodenum, etree);

  sy= 1.4*500.0/treeLevel(etree);
  s= 500.0/treeLeaveNbr(etree);

  if (nodenum == ROOT(etree))
    {
      cx= 0;
    }

  pos= 0;

  /* draw sons */
  switch (etree->nodetab[nodenum].sonsnbr)
    {
    case 2:
      {
        rpos= ps_rec_write_tree(TSON(etree, nodenum, 1), costmtx, etree, out, ps_draw_node);
        pos+= rpos;
      }
    case 1:
      {
        lpos= ps_rec_write_tree(TSON(etree, nodenum, 0), costmtx, etree, out, ps_draw_node);
        pos+= lpos;
        pos/= etree->nodetab[nodenum].sonsnbr;
        break;
      }
    case 0:
      {
        pos= (double)cx;
        cx++;
        break;
      }
    default:
      {
          errorPrint("Erreur dans l'arbre d'elimination \n\t Cblk %ld sonsnbr %ld \n\tAttention on ne peut imprimer l'arbre d'elimination que pour les grilles !!", (long)nodenum, (long)(etree->nodetab[nodenum].sonsnbr));
          EXIT(MOD_BLEND,INTERNAL_ERR);
      }
    }

  /* draw links to sons */
  switch (etree->nodetab[nodenum].sonsnbr)
    {
    case 2:
      {
        fprintf(out, "np %f %f m %f %f l %f %f l sk\n",
                lpos*s, (nodelevel+1)*sy,
                pos*s, nodelevel*sy,
                rpos*s, (nodelevel+1)*sy);
        break;
      }
    case 1:
      {
        fprintf(out, "np %f %f m %f %f l sk\n",
                lpos*s, (nodelevel+1)*sy,
                pos*s, (nodelevel)*sy);
        break;
      }
    }


  ps_draw_node(out, nodenum, costmtx, etree, sy, pos*s,
               sy*nodelevel);

  return pos;
}

void ps_draw_node_num(FILE *out, PASTIX_INT nodenum, const CostMatrix *costmtx, const EliminTree *etree,
                        double s, double x, double y)
{
  double ss;

  ss= s * costmtx->cblktab[nodenum].total/cblkMaxCost(etree->nodenbr, costmtx);

  fprintf(out, "gs\n");

  fprintf(out, "np %f %f m "
          "%f %f rl %f %f rl %f %f rl cp ",
          x-ss/2, y-ss/2,
          ss, (double)0, (double)0, ss, -ss, (double)0);
  fprintf(out, "sk \n");
  fprintf(out, "%f sf\n", ss);
  fprintf(out, "%f %f m (%ld) s \n",
          x-0.9*ss/2, y-0.9*ss/2, (long)nodenum);
  fprintf(out, "gr\n");
}





/* imprime l'arbre d'elimination en post-script */
void ps_write_tree_owner(PASTIX_INT *ownertab, const CostMatrix *costmtx, const EliminTree *etree, FILE *out, PASTIX_INT *page)
     /* PASTIX_INT (*ps_draw_node)(FILE *out, cbl c, cbls allcbls, double s, double x, double y) */
{
  fprintf(out, "%%%%Page: %ld %ld \n", (long)*page, (long)*page); (*page)++;
  fprintf(out, "%% Elimination Tree (color by Layer)\n");

  fprintf(out, "40 40 translate 1 1 scale\n");
  ps_rec_write_tree_owner(ROOT(etree), ownertab, costmtx, etree, out, ps_draw_node_owner);
  fprintf(out, "showpage\n");

  fflush(out);

}


double ps_rec_write_tree_owner(PASTIX_INT nodenum, PASTIX_INT *ownertab, const CostMatrix *costmtx, const EliminTree *etree, FILE *out,
                            void (*ps_draw_node)(FILE *out, PASTIX_INT nodenum, PASTIX_INT procnum, const CostMatrix *costmtx,
                                                const EliminTree *etree, double s,
                                                double x, double y))
{
  static PASTIX_INT cx;
  PASTIX_INT nodelevel;
  double pos=0, lpos=0, rpos=0;
  double s, sy;

  nodelevel = nodeTreeLevel(nodenum, etree);

  sy= 1.4*500.0/treeLevel(etree);
  s= 500.0/treeLeaveNbr(etree);

  if (nodenum == ROOT(etree))
    {
      cx= 0;
    }

  pos= 0;

  /* draw sons */
  switch (etree->nodetab[nodenum].sonsnbr)
    {
    case 2:
      {
        rpos= ps_rec_write_tree_owner(TSON(etree, nodenum, 1), ownertab, costmtx, etree, out, ps_draw_node);
        pos+= rpos;
      }
    case 1:
      {
        lpos= ps_rec_write_tree_owner(TSON(etree, nodenum, 0), ownertab, costmtx, etree, out, ps_draw_node);
        pos+= lpos;
        pos/= etree->nodetab[nodenum].sonsnbr;
        break;
      }
    case 0:
      {
        pos= (double)cx;
        cx++;
        break;
      }
    default:
      {
          errorPrint("Erreur dans l'arbre d'elimination \n\tCblk %ld sonsnbr %ld \n\tAttention on ne peut imprimer l'arbre d'elimination que pour les grilles !!", (long)nodenum, (long)(etree->nodetab[nodenum].sonsnbr));
          EXIT(MOD_BLEND,INTERNAL_ERR);
      }
    }

  /* draw links to sons */
  switch (etree->nodetab[nodenum].sonsnbr)
    {
    case 2:
      {
        fprintf(out, "np %f %f m %f %f l %f %f l sk\n",
                lpos*s, (nodelevel+1)*sy,
                pos*s, nodelevel*sy,
                rpos*s, (nodelevel+1)*sy);
        break;
      }
    case 1:
      {
        fprintf(out, "np %f %f m %f %f l sk\n",
                lpos*s, (nodelevel+1)*sy,
                pos*s, (nodelevel)*sy);
        break;
      }
    }


  ps_draw_node_owner(out, nodenum, ownertab[nodenum], costmtx, etree, sy, pos*s,
               sy*nodelevel);

  return pos;
}













void ps_draw_node_owner(FILE *out, PASTIX_INT nodenum, PASTIX_INT procnum, const CostMatrix *costmtx, const EliminTree *etree,
                        double s, double x, double y)
{
  double ss;

  ss= s * costmtx->cblktab[nodenum].total/cblkMaxCost(etree->nodenbr, costmtx);

  fprintf(out, "gs\n");

  fprintf(out, "np %f %f m "
          "%f %f rl %f %f rl %f %f rl cp ",
          x-ss/2, y-ss/2,
          ss, (double)0, (double)0, ss, -ss, (double)0);
  fprintf(out, "sk \n");
  fprintf(out, "%f sf\n", ss);
  fprintf(out, "%f %f m (%ld) s \n",
          x-0.9*ss/2, y-0.9*ss/2, (long)procnum);
  fprintf(out, "gr\n");
}
