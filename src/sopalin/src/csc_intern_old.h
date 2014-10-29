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
/*== Creation/Destruction de CSC ==*/

/******************************************************************************/
/* void CscHbRead(CscMatrix *thecsc, char *Type, char *RhsType, PASTIX_FLOAT **rhs,  */
/*                PASTIX_FLOAT **rhs2, const Order *ord, const char *rsaname)        */
/*                                                                            */
/* Lecture de la csc a partir d'un fichier au format HB et permutation a      */
/* partir du vecteur permutation fournit par Scotch                           */
/*                                                                            */
/* thecsc : La csc                                                            */
/* Type : type du HB (RUA, RSA ....)                                          */
/* RhsType : type pour les seconds membres                                    */
/* rhs : Vecteur second membre                                                */
/* rhs2 : Vecteur solution                                                    */
/* ord : la permutation                                                       */
/* rsaname : nom du fichier HB                                                */
/*                                                                            */
/* Type doit etre alloue avant l'appel, char Type[4]                          */
/* RhsType doit etre alloue avant l'appel, char RhsType[4]                    */
/* rhs et rhs2 sont alloue si besoin est.                                     */
/* RhsType[0] == '\0' si pas de second membre dans le fichier                 */
/******************************************************************************/

/* !!!!!!!!!!!!!FONCTION INUTILISEE !!!!!!!! */
void CscOrder(CscMatrix *thecsc,
	      char *Type, char *RhsType, PASTIX_FLOAT **transcsc,
	      PASTIX_FLOAT **rhs, PASTIX_FLOAT **rhs2,
	      const Order *ord,
	      PASTIX_INT Nrow, PASTIX_INT Ncol, PASTIX_INT Nnzero, 
	      PASTIX_INT *colptr, PASTIX_INT *rowind, PASTIX_FLOAT *val, PASTIX_INT forcetrans);

/*== Distribution/Remplissage ==*/
/******************************************************************************/
/* void CscDistrib(SymbolMatrix *symbmtx, CscMatrix *cscmtx,                  */
/*                 CscMatrix *cscdist)                                        */
/*                                                                            */
/* Distribution de la csc                                                     */
/*                                                                            */
/* symbmtx : symbol matrix locale                                             */
/* cscmtx : csc globale                                                       */
/* cscdist : csc locale                                                       */
/******************************************************************************/
void CscDistrib(const SymbolMatrix *symbmtx, const CscMatrix *cscmtx,
		CscMatrix *cscdist, const PASTIX_FLOAT *transcsc, PASTIX_FLOAT **trandcsc);

/*== Divers ==*/
/******************************************************************************/
/* void CscTrans(CscMatrix *cscmtx, CscMatrix *csctrp)                        */
/*                                                                            */
/* Transpose une csc                                                          */
/*                                                                            */
/* cscmtx : csc                                                               */
/* csctrp : csc transposee                                                    */
/******************************************************************************/
void CscTrans(const CscMatrix *cscmtx, CscMatrix *csctrp);

void CscScaling2(char *Type, PASTIX_INT Ncol, PASTIX_INT *col, PASTIX_INT *row, PASTIX_FLOAT *val, PASTIX_FLOAT *rhs, PASTIX_FLOAT *rhs2);
void CscScaling(CscMatrix *cscmtx, PASTIX_FLOAT *transcsc, PASTIX_FLOAT *rhs1, PASTIX_FLOAT *rhs2);

/******************************************************************************/
/* void CscVerifUpdown(const UpDownVector *updovct,                           */
/*                     const SymbolMatrix *symbmtx; const PASTIX_FLOAT *rhs2)        */
/*                                                                            */
/* Verification entre le second membre fournit dans le fichier HB et le second*/
/* membre calcule.                                                            */
/*                                                                            */
/* updovct : vecteur second membre calcule                                    */
/* symbmtx : Symbol matrix                                                    */
/* rhs2 : vecteur second membre solution fournit dans le fichier HB           */
/******************************************************************************/
void CscVerifUpdown(const UpDownVector *updovct, const SymbolMatrix *symbmtx,
		    const PASTIX_FLOAT *rhs2);

/******************************************************************************/
/* void CscUpDown(UpDownVector *updovct, const SymbolMatrix *symbmtx,         */
/*                const PASTIX_FLOAT *rhs)                                           */
/*                                                                            */
/* Remplissage du vector second membre a partir de celui fournit dans le      */
/* fichier HB.                                                                */
/*                                                                            */
/******************************************************************************/
void CscUpdown(UpDownVector *updovct, /*const*/ SymbolMatrix *symbmtx,
	       const PASTIX_FLOAT *rhs);
void CscUpdown2(UpDownVector *updovct, /*const*/ SymbolMatrix *symbmtx,
		const PASTIX_FLOAT *rhs);



/******************************************************************************/
/* void CscDiagDom(CscMatrix *cscmtx)                                         */
/*                                                                            */
/* Transforme la csc en csc a diagonale dominante                             */
/*                                                                            */
/* cscmtx : Csc                                                               */
/******************************************************************************/
void CscDiagDom(CscMatrix *cscmtx);

PASTIX_INT CscStrucSym(CscMatrix *cscmtx);
