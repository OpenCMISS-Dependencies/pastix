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
#include "blend_distributeOnGPU.h"
#include "flops.h"

#define symbol_get_cblk_stride(  _datacode_, _cblknum_ ) (SOLV_STRIDE((_cblknum_)))
#define symbol_get_cblk_width(   _datacode_, _cblknum_ ) (SYMB_LCOLNUM((_cblknum_))-SYMB_FCOLNUM((_cblknum_))+1)
#define symbol_get_blok_coefind( _datacode_, _bloknum_ ) (SOLV_COEFIND((_bloknum_)))
#define symbol_get_blok_height(  _datacode_, _bloknum_ ) (BLOK_ROWNBR((_bloknum_)))
#include "fifo.h"

typedef struct _cblk_elem_ {
  dague_list_item_t item;
  /* for the sort */
  int               criterium;
  SolverCblk       *cblk;
  int cblk_size;
  int updates;
  float flop; /* in Mflop */
}cblk_elem;

typedef struct _gpu_elem_ {
  dague_list_item_t item;
  //gpu_device_t* gpu_device;
  int current_avail_mem;
  int gpu_total_mem;
  int id;
}gpu_elem;


int blend_distributeOnGPU(SolverMatrix * datacode,
                          double         maxMem,
                          int            pageSize,
                          int            criterium,
                          enum API_GPU_CRITERIUM nGPUs,
                          enum API_FLOAT floatType,
                          enum API_FACT  factType) {
  SolverCblk   *cblktab = datacode->cblktab;
  UpDownVector * updovect = &(datacode->updovct);
  PASTIX_INT   cblknbr, cblknum;
  dague_list_t cblklist;
  dague_list_t gpulist;
  PASTIX_INT cblk_offset  = -1;
  PASTIX_INT gpu_offset   = -1;
  PASTIX_INT color_offset = -1;
  PASTIX_INT gcblk2list   = -1;
  PASTIX_INT ndevices;
  int i;
  PASTIX_INT j, bloknum, m, k, n;
  float flop;
  int color;
  gpu_elem *new_gpu;
  cblk_elem *new_cblk;
  gpu_elem *gpu_cursor;
  cblk_elem *cblk_cursor;
  int type_sze;
  size_t unit_size;
  ndevices = nGPUs;

  switch (floatType) {
  case API_REALSINGLE:
    type_sze=sizeof(float);
    break;
  case API_REALDOUBLE:
    type_sze=sizeof(double);
    break;
  case API_COMPLEXSINGLE:
    type_sze=2*sizeof(float);
    break;
  case API_COMPLEXDOUBLE:
    type_sze=2*sizeof(double);
    break;
  default:
    errorPrint("Unkwnown type");
    return FLOAT_TYPE_ERR;
  }
  unit_size = pageSize/type_sze;

  int* devices_cblk = malloc((1+ndevices)*sizeof(int));
  int* devices_gemm = malloc((1+ndevices)*sizeof(int));
  memset(devices_cblk,0, (1+ndevices)*sizeof(int));
  memset(devices_gemm,0, (1+ndevices)*sizeof(int));

  /*
   * Compute the GPU distribution
   */
  dague_list_construct( &cblklist );
  dague_list_construct( &gpulist );

  /* Sort the cblk according to criterium */
  cblknbr = SYMB_CBLKNBR;
  //fprintf(stdout,"start loop on  cblk \n");
  for(cblknum = 0; cblknum < cblknbr; cblknum++){
    gcblk2list = UPDOWN_GCBLK2LIST(UPDOWN_LOC2GLOB(cblknum));

    /* if it's a leaf, color = cpu */
    if ( gcblk2list == -1 || ndevices <=0) {
      cblktab[cblknum].color = -1;
      devices_cblk[0]++;
    } else {
      size_t cblksize;
      int    nbpages;
      PASTIX_INT updates = UPDOWN_LISTPTR(gcblk2list+1)
        -                 UPDOWN_LISTPTR(gcblk2list);

      cblktab[cblknum].color = -2;
      /* If we are here, we have at least one update */
      assert(updates > 0);
      new_cblk = (cblk_elem*)malloc(sizeof(cblk_elem));
      new_cblk->updates = updates;
      dague_list_item_construct( (dague_list_item_t*)new_cblk );
      new_cblk->cblk = cblktab + cblknum;

      /* FLOP */
      int cblknum2;
      flop = 0.0;
      //fprintf(stdout,"start flops computation for cblk %ld \n", cblknum);
      for(j=UPDOWN_LISTPTR(gcblk2list); j<UPDOWN_LISTPTR(gcblk2list+1);j++){
        bloknum = UPDOWN_LISTBLOK( j );
        cblknum2 = updovect->listcblk[j];
        //cblknum2 = sparse_matrix_get_lcblknum(datacode, bloknum);
        m = symbol_get_cblk_stride( datacode, cblknum2 ) - symbol_get_blok_coefind(datacode, bloknum);
        k = symbol_get_cblk_width( datacode, cblknum2 );
        n = symbol_get_blok_height( datacode, bloknum );
        flop += (float)FLOPS_SGEMM(m,n,k)/(float)(2*1e6);
      }
      //fprintf(stdout,"stop flops computation for cblk %ld \n", cblknum);

      /* if (floatType == API_COMPLEXSINGLE || floatType == API_COMPLEXDOUBLE) { */
      /*   flop *= 4.0; */
      /* } */

      /* if (factType == API_FACT_LU){ *\/ */
      /*   flop *= 2.0; */
      /* } */

      new_cblk->flop = flop;
      /* check for int overflow */
      assert(flop < 4294967296.);

      /* Amount of memory to push for the cblk */
      cblksize = symbol_get_cblk_width( datacode, cblknum ) * symbol_get_cblk_stride( datacode, cblknum );
      nbpages = ( cblksize + unit_size - 1) / unit_size;

      if (factType == API_FACT_LU){
          nbpages *= 2.0;
      }

      new_cblk->cblk_size = nbpages;
      /* Sort criterium */
      switch(criterium){
      case API_GPU_CRITERION_UPDATES:
        new_cblk->criterium = new_cblk->updates;
        break;
      case API_GPU_CRITERION_CBLKSIZE:
        new_cblk->criterium = new_cblk->cblk_size;
        break;
      case API_GPU_CRITERION_FLOPS:
        new_cblk->criterium = (int)new_cblk->flop;
        break;
      default :
        new_cblk->criterium = -cblknum;
      }

      if( cblk_offset < 0 ) {
        cblk_offset = (void*)&(new_cblk->criterium) - (void*)new_cblk;
        assert( cblk_offset > 0 );
      }

      /* Sort */
      dague_list_nolock_push_sorted(&cblklist, (dague_list_item_t*)new_cblk, cblk_offset);
    }
  }
  // fprintf(stdout,"end loop on  cblk \n");

  /* Sort the GPUs according to available memory */
  for(i = 0; i < ndevices; i++) {
    new_gpu = (gpu_elem*)malloc(sizeof(gpu_elem));
    dague_list_item_construct( (dague_list_item_t*)new_gpu );
    //new_gpu->gpu_device = gpu_enabled_devices[i];
    new_gpu->id = (PASTIX_INT)i;
    if( gpu_offset < 0 ) {
      gpu_offset = (void *)&(new_gpu->current_avail_mem) - (void*)new_gpu;
      assert( gpu_offset > 0 );
    }
    new_gpu->gpu_total_mem     = (int)(GPU_MAX_FILL * maxMem /pageSize);
    //new_gpu->gpu_total_mem     =  new_gpu->gpu_device->memory->max_segment);
    new_gpu->current_avail_mem = new_gpu->gpu_total_mem;
    dague_list_nolock_push_sorted(&gpulist, (dague_list_item_t*)new_gpu, gpu_offset);

    fprintf(stdout, "Push init: gpu(%d) available mem=%d\n",
            new_gpu->id, new_gpu->current_avail_mem);
    }

    /* Create association between cblks and gpus */
  if(ndevices > 0){
    while(! dague_list_nolock_is_empty( &cblklist ) ){
      gpu_cursor = NULL;

      /* Get the cblk with the highest criterium */
      cblk_cursor = (cblk_elem*)dague_list_nolock_pop_front(&cblklist);

      /* The cblk has a predicted GPU */
      color = cblk_cursor->cblk->color;

      if (color != -2) {
        assert( color != -1 );
        if(color_offset == -1){
          gpu_elem *elem = (gpu_elem*)(&gpulist);
          color_offset = (uintptr_t)(&(elem->id)) - (uintptr_t)(&gpulist);
        }
        gpu_cursor = (gpu_elem*)dague_list_extract( &gpulist, color, (int)color_offset );
        if(gpu_cursor != NULL){
          if(gpu_cursor->current_avail_mem > cblk_cursor->cblk_size){
            dague_list_nolock_push_sorted( &gpulist, (dague_list_item_t*)gpu_cursor, gpu_offset );
            gpu_cursor = NULL;

          }
        }
      }

      /* Get the gpu with the highest criterium */
      if (gpu_cursor == NULL)
        gpu_cursor  = (gpu_elem*) dague_list_nolock_pop_front(&gpulist);

      //fprintf(stdout, "Pop: cblk(%ld) size=%ld, nbpages=%d, updates=%d\n",
      //(cblk_cursor->cblk) - cblktab,
      //(cblk_cursor->cblk->lcolnum - cblk_cursor->cblk->fcolnum + 1) * cblk_cursor->cblk->stride,
      //cblk_cursor->cblk_size, cblk_cursor->updates );
      //fprintf(stdout, "Pop: gpu(%d) available mem=%d\n",
      //gpu_cursor->id, gpu_cursor->current_avail_mem);

      assert(gpu_cursor != NULL);
      if( (gpu_cursor->current_avail_mem >= cblk_cursor->cblk_size)
          && (gpu_cursor->current_avail_mem > 0)
          && (cblk_cursor->updates >= GPU_MIN_UPDATES)
          && (cblk_cursor->cblk_size >= GPU_MIN_NBPAGES)
          && (cblk_cursor->flop >= GPU_MIN_FLOP)
          ) {
        /* association */
        gpu_cursor->current_avail_mem -= cblk_cursor->cblk_size;
        cblk_cursor->cblk->color = gpu_cursor->id;

        /* Add prediction to contibutors */


        for(j=UPDOWN_LISTPTR(gcblk2list); j<UPDOWN_LISTPTR(gcblk2list+1);j++){
          bloknum = UPDOWN_LISTBLOK( j );
          cblknum = updovect->listcblk[j];
          //cblknum = sparse_matrix_get_lcblknum(datacode, bloknum);
          if(cblktab[cblknum].color == -2 && gcblk2list != -1){
            cblktab[cblknum].color = gpu_cursor->id;
          }
        }

      } else {
        /* send to CPU */
        cblk_cursor->cblk->color = -1;
      }

      devices_cblk[cblk_cursor->cblk->color+1]++;
      devices_gemm[cblk_cursor->cblk->color+1]+=cblk_cursor->updates;

      free(cblk_cursor);

      //fprintf(stdout, "Push: gpu(%d) available mem=%d\n",
      //gpu_cursor->id, gpu_cursor->current_avail_mem);

      dague_list_nolock_push_sorted( &gpulist, (dague_list_item_t*)gpu_cursor, gpu_offset );
    }
  }

  fprintf(stdout,"cpu :  %d cblk, %d gemm\n", devices_cblk[0], devices_gemm[0]);

  while(! dague_list_nolock_is_empty( &gpulist ) ){
    gpu_cursor  = (gpu_elem*) dague_list_nolock_pop_front(&gpulist);

    i = gpu_cursor->id;
    fprintf(stdout,"gpu %d:  %d cblk, %d gemm, memory : %d / %d \n",
            i, devices_cblk[i+1], devices_gemm[i+1],
            gpu_cursor->current_avail_mem, gpu_cursor->gpu_total_mem );
    free(gpu_cursor);
  }

  dague_list_destruct( &cblklist );
  dague_list_destruct( &gpulist );

  return cblknbr;
}
