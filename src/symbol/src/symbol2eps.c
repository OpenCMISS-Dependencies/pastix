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
/**   NAME       : symbol2eps.c                            **/
/**                                                        **/
/**   AUTHORS    : Pascal HENON                            **/
/**                Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : Part of a parallel direct block solver. **/
/**                This module draws symbolic matrices in  **/
/**                PostScript (tm) format .                **/
/**                                                        **/
/**   DATES      : # Version 0.0  : from : 29 sep 1998     **/
/**                                 to     29 sep 1998     **/
/**                # Version 1.0  : from : 26 jun 2002     **/
/**                                 to     26 jun 2002     **/
/**                # Version 1.3  : from : 10 apr 2003     **/
/**                                 to     10 jun 2003     **/
/**                                                        **/
/************************************************************/

/*
**  The defines and includes.
*/

#include "common_pastix.h"
#include "symbol.h"
#include "symbol_draw.h"

#define SYMBOL2EPS_EXENAME "symbol2eps"

int main(int argc, char **argv) {
  int i;
  int nb_param_option = 0;
  int last_nb_param_option = 0;
  int lu_option = 0;
  int nb_infile = 0;
  SymbolMatrix symbmtx;
  SymbolMatrix * symbmtxtab;
  char * filename_pattern;
  char filename[256];
  FILE *temp_stream;

  fprintf(stderr, "[SYMBOL2EPS] Generation d'une representation de la SymbolMatrix\n");

  /* Le nom de l'exe doit etre valide */
  if (strcmp(argv[0],SYMBOL2EPS_EXENAME) != 0) {
    fprintf(stderr, "Non d'executable non supporte ...\n");
    exit(0);
  }

  if ((argc < 2) || (argc > 4)) {
    fprintf(stderr, "Usage : %s [-lu] symbmtx_filename [nb_files]\n",SYMBOL2EPS_EXENAME);
    exit(0);  
  }
  
  /* on recherche le nombre d'option */
  for (i=0; i < argc; i++) {
    
    /* test du prefix de l'argument */
    if (strncmp(argv[1+nb_param_option],"-",1) == 0) {
      last_nb_param_option = nb_param_option;
      
      /* on doit etre en mesure de trouver une option */
      if (strcmp(argv[1+nb_param_option],"-lu") == 0) {
	lu_option = 1;
	nb_param_option++;
      }

      if ((last_nb_param_option == nb_param_option ) || (argc < 2) || (argc > 4)) {
	fprintf(stderr, "Invalid option %s \n",argv[nb_param_option +1]);
	fprintf(stderr, "Usage : %s [-lu] symbmtx_filename [nb_files]\n",SYMBOL2EPS_EXENAME);
	exit(0);

      }
    }
  }

  if (lu_option == 1)
    symbolDrawSetLU();

  if (nb_param_option + 3 == argc) {
    /* il s'agit d'un cas distribue */
    nb_infile = atoi(argv[argc - 1]);
    fprintf(stderr, "Visualisation de la distribution on %ld processors\n",nb_infile);
  }

  /* pattern */
  filename_pattern = argv[nb_param_option +1];
  fprintf(stderr,"Pattern : %s\n",filename_pattern);


  if (nb_infile == 0) {

    /* chargement ... */
    symbolInit(&symbmtx);
    sprintf(filename,"%s",filename_pattern);
    temp_stream = fopen(filename, "r");
    if (symbolLoad(&symbmtx,temp_stream) != 0)
      {
	errorPrint ("test: cannot load %s",filename_pattern);
	exit (1);
      } 
    fclose(temp_stream);
    fprintf(stderr,"Chargement ok !\n");
    
    sprintf(filename,"%s.eps",filename_pattern);
    temp_stream = fopen(filename, "w");
    if (symbolDrawFunc (&symbmtx, NULL, NULL, NULL, temp_stream) != 0)
      {
	errorPrint ("test: cannot generate %s",filename);
	exit (1);
      } 
    fclose(temp_stream); 
    fprintf(stderr,"Generation ok !\n");

  } else {

    /* allocation du tableau recevant les symbols matrix */
    MALLOC_INTERN(symbmtxtab, nb_infile, void* );
    
    /* remplissage de ce tableau ... */
    for (i=0; i < nb_infile; i++) {
      symbolInit(&symbmtxtab[i]);
      sprintf(filename,"%s%ld",filename_pattern,i);
      /*fprintf(stderr,"tentative d'ouvertue %s",filename);*/
      fprintf(stderr,".");
      temp_stream = fopen(filename, "r");
      if (symbolLoad(&symbmtxtab[i],temp_stream) != 0)
	{
	  errorPrint ("test: cannot load %s%ld",filename_pattern,i);
	  exit (1);
	} 
      fclose(temp_stream);       
    }
    fprintf(stderr,"\nChargement ok !\n");
    
    sprintf(filename,"%s.eps",filename_pattern);
    temp_stream = fopen(filename, "w");
    if (dsymbolDrawFunc (nb_infile,symbmtxtab, temp_stream) != 0)
      {
	errorPrint ("test: cannot generate %s",filename);
	exit (1);
      } 
    fclose(temp_stream);
    fprintf(stderr,"Generation ok !\n");

  } 
}
