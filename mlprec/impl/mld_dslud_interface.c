/*
 * 
 *                            MLD2P4  version 2.1
 *   MultiLevel Domain Decomposition Parallel Preconditioners Package
 *              based on PSBLAS (Parallel Sparse BLAS version 3.5)
 *   
 *   (C) Copyright 2008-2018
 * 
 *      Salvatore Filippone                             
 *      Ambra Abdullahi Hassan                     
 *      Alfredo Buttari        CNRS-IRIT, Toulouse, FR                            
 *      Pasqua D'Ambra                                  
 *      Daniela di Serafino   
 * 
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *    2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions, and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *    3. The name of the MLD2P4 group or the names of its contributors may
 *       not be used to endorse or promote products derived from this
 *       software without specific written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE MLD2P4 GROUP OR ITS CONTRIBUTORS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 * 
 *
 * File: mld_dslud_interface.c
 *
 * Functions: mld_dsludist_fact, mld_dsludist_solve, mld_dsludist_free.
 *
 * This file is an interface to the SuperLU_dist routines for sparse factorization and
 * solve. It was obtained by modifying the c_fortran_dgssv.c file from the SuperLU_dist
 * source distribution; original copyright terms are reproduced below.
 * 
 */

/*		=====================

Copyright (c) 2003, The Regents of the University of California, through
Lawrence Berkeley National Laboratory (subject to receipt of any required 
approvals from U.S. Dept. of Energy) 

All rights reserved. 

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

(1) Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer. 
(2) Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution. 
(3) Neither the name of Lawrence Berkeley National Laboratory, U.S. Dept. of
Energy nor the names of its contributors may be used to endorse or promote
products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
  
*/

/*
 * -- Distributed SuperLU routine (version 2.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * March 15, 2003
 *
 */

#ifdef Have_SLUDist_
#include <math.h>
#include "superlu_ddefs.h"

#define HANDLE_SIZE  8

#if defined(SLUD_VERSION_6)
typedef struct {
  SuperMatrix *A;
  dLUstruct_t *LUstruct;
  gridinfo_t *grid;
  dScalePermstruct_t *ScalePermstruct;
} factors_t;
#else
typedef struct {
  SuperMatrix *A;
  LUstruct_t *LUstruct;
  gridinfo_t *grid;
  ScalePermstruct_t *ScalePermstruct;
} factors_t;
#endif

#else

#include <stdio.h>

#endif


int mld_dsludist_fact(int n, int nl, int nnzl, int ffstr,
		  double *values, int *rowptr, int *colind,
		  void **f_factors, int nprow, int npcol)
{
/* 
 * This routine can be called from Fortran.
 *  performs LU decomposition.
 *
 * f_factors (input/output) void** 
 *      On  output contains the pointer pointing to
 *       the structure of the factored matrices.
 *
 */
 
#ifdef Have_SLUDist_
    SuperMatrix *A;
    NRformat_loc *Astore;

#if defined(SLUD_VERSION_6)
    dScalePermstruct_t *ScalePermstruct;
    dLUstruct_t *LUstruct;
    dSOLVEstruct_t SOLVEstruct;
#else
    ScalePermstruct_t *ScalePermstruct;
    LUstruct_t *LUstruct;
    SOLVEstruct_t SOLVEstruct;
#endif
    gridinfo_t *grid;
    int      i, panel_size, permc_spec, relax, info;
    trans_t  trans;
    double   drop_tol = 0.0, b[1], berr[1];
#if defined(SLUD_VERSION_6)||defined(SLUD_VERSION_5)
    superlu_dist_options_t options;
#elif defined(SLUD_VERSION_4)||defined(SLUD_VERSION_3)
    superlu_options_t options;
#else
    choke_on_me;
#endif
    SuperLUStat_t stat;
    factors_t *LUfactors;
    int fst_row;
    int *icol,*irpt;
    double *ival;

    trans = NOTRANS;
    grid = (gridinfo_t *) SUPERLU_MALLOC(sizeof(gridinfo_t));
    superlu_gridinit(MPI_COMM_WORLD, nprow, npcol, grid);
    /* Initialize the statistics variables. */
    PStatInit(&stat);
    fst_row = (ffstr);
    
    A  = (SuperMatrix *) malloc(sizeof(SuperMatrix));
    dCreate_CompRowLoc_Matrix_dist(A, n, n, nnzl, nl, fst_row,
				   values, colind, rowptr,
				   SLU_NR_loc, SLU_D, SLU_GE);
    
    /* Initialize ScalePermstruct and LUstruct. */
#if defined(SLUD_VERSION_6)
    ScalePermstruct = (dScalePermstruct_t *) SUPERLU_MALLOC(sizeof(dScalePermstruct_t));
    LUstruct = (dLUstruct_t *) SUPERLU_MALLOC(sizeof(dLUstruct_t));
    dScalePermstructInit(n,n, ScalePermstruct);
#else
    ScalePermstruct = (ScalePermstruct_t *) SUPERLU_MALLOC(sizeof(ScalePermstruct_t));
    LUstruct = (LUstruct_t *) SUPERLU_MALLOC(sizeof(LUstruct_t));
    ScalePermstructInit(n,n, ScalePermstruct);
#endif
#if defined(SLUD_VERSION_6)
    dLUstructInit(n, LUstruct);
#elif defined(SLUD_VERSION_4) || defined(SLUD_VERSION_5)
    LUstructInit(n, LUstruct);
#elif defined(SLUD_VERSION_3) 
    LUstructInit(n,n, LUstruct);
#else
    choke_on_me;
#endif

    /* Set the default input options. */
    set_default_options_dist(&options);
    options.IterRefine=NO;
    options.PrintStat=NO;

    pdgssvx(&options, A, ScalePermstruct, b, nl, 0,
	    grid, LUstruct, &SOLVEstruct, berr, &stat, &info);
    
    if ( info == 0 ) {
      ;
    } else {
      printf("pdgssvx() error returns INFO= %d\n", info);
      if ( info <= n ) { /* factorization completes */
	; 
      }
    }
    if (options.SolveInitialized) {
      dSolveFinalize(&options,&SOLVEstruct);
    }
    
    
    /* Save the LU factors in the factors handle */
    LUfactors = (factors_t *) SUPERLU_MALLOC(sizeof(factors_t));
    LUfactors->LUstruct = LUstruct;
    LUfactors->grid     = grid;
    LUfactors->A        = A;
    LUfactors->ScalePermstruct = ScalePermstruct;  
/*     fprintf(stderr,"slud factor: LUFactors %p \n",LUfactors);  */
/*     fprintf(stderr,"slud factor: A %p %p\n",A,LUfactors->A);  */
/*     fprintf(stderr,"slud factor: grid %p %p\n",grid,LUfactors->grid);  */
/*     fprintf(stderr,"slud factor: LUstruct %p %p\n",LUstruct,LUfactors->LUstruct);  */
    *f_factors = (void *) LUfactors;
    PStatFree(&stat);
    return(info);
#else
    fprintf(stderr," SLUDist Not Configured, fix make.inc and recompile\n");
    return(-1);
#endif
}


int mld_dsludist_solve(int itrans, int n, int nrhs, 
		       double *b, int ldb, void *f_factors)

{
/* 
 * This routine can be called from Fortran.
 *      performs triangular solve
 *
 */
#ifdef Have_SLUDist_ 
    SuperMatrix *A;
#if defined(SLUD_VERSION_6)
    dScalePermstruct_t *ScalePermstruct;
    dLUstruct_t *LUstruct;
    dSOLVEstruct_t SOLVEstruct;
#else
    ScalePermstruct_t *ScalePermstruct;
    LUstruct_t *LUstruct;
    SOLVEstruct_t SOLVEstruct;
#endif
    gridinfo_t *grid;
    int      i, panel_size, permc_spec, relax, info;
    trans_t  trans;
    double   drop_tol = 0.0;
    double *berr;
#if defined(SLUD_VERSION_6)||defined(SLUD_VERSION_5)
    superlu_dist_options_t options;
#elif defined(SLUD_VERSION_4)||defined(SLUD_VERSION_3)
    superlu_options_t options;
#else
    choke_on_me;
#endif
    SuperLUStat_t stat;
    factors_t *LUfactors;

    LUfactors       = (factors_t *) f_factors   ;
    A               = LUfactors->A              ;
    LUstruct        = LUfactors->LUstruct       ;
    grid            = LUfactors->grid           ;

    ScalePermstruct = LUfactors->ScalePermstruct;
/*     fprintf(stderr,"slud solve: LUFactors %p \n",LUfactors);  */
/*     fprintf(stderr,"slud solve: A %p %p\n",A,LUfactors->A);  */
/*     fprintf(stderr,"slud solve: grid %p %p\n",grid,LUfactors->grid);  */
/*     fprintf(stderr,"slud solve: LUstruct %p %p\n",LUstruct,LUfactors->LUstruct);  */


    if (itrans == 0) {
      trans = NOTRANS;
    } else if (itrans ==1) {
      trans = TRANS;
    } else if (itrans ==2) {
      trans = CONJ;
    } else {
      trans = NOTRANS;
    }

/*     fprintf(stderr,"Entry to sludist_solve\n"); */
    berr = (double *) malloc((nrhs) *sizeof(double));

    /* Initialize the statistics variables. */
    PStatInit(&stat);
    
    /* Set the default input options. */
    set_default_options_dist(&options);
    options.IterRefine = NO;
    options.Fact       = FACTORED;
    options.PrintStat  = NO;

    pdgssvx(&options, A, ScalePermstruct, b, ldb, nrhs, 
	    grid, LUstruct, &SOLVEstruct, berr, &stat, &info);
    
/*     fprintf(stderr,"Double check: after solve %d %lf\n",*info,berr[0]); */
    if (options.SolveInitialized) {
      dSolveFinalize(&options,&SOLVEstruct);
    }
    PStatFree(&stat);
    free(berr);
    return(info);
#else
    fprintf(stderr," SLUDist Not Configured, fix make.inc and recompile\n");
    return(-1);
#endif
    
}


int mld_dsludist_free(void *f_factors)
{
/* 
 * This routine can be called from Fortran.
 *
 *      free all storage in the end
*
 */
#ifdef Have_SLUDist_ 
    SuperMatrix *A;
#if defined(SLUD_VERSION_6)
    dScalePermstruct_t *ScalePermstruct;
    dLUstruct_t *LUstruct;
    dSOLVEstruct_t SOLVEstruct;
#else
    ScalePermstruct_t *ScalePermstruct;
    LUstruct_t *LUstruct;
    SOLVEstruct_t SOLVEstruct;
#endif
    gridinfo_t *grid;
    int      i, panel_size, permc_spec, relax;
    trans_t  trans;
    double   drop_tol = 0.0;
    double *berr;
#if defined(SLUD_VERSION_6)||defined(SLUD_VERSION_5)
    superlu_dist_options_t options;
#elif defined(SLUD_VERSION_4)||defined(SLUD_VERSION_3)
    superlu_options_t options;
#else
    choke_on_me;
#endif
    SuperLUStat_t stat;
    factors_t *LUfactors;


    if (f_factors == NULL) 
      return(0);
    LUfactors       = (factors_t *)  f_factors  ;
    A               = LUfactors->A              ;
    LUstruct        = LUfactors->LUstruct       ;
    grid            = LUfactors->grid           ;
    ScalePermstruct = LUfactors->ScalePermstruct;

    // Memory leak: with SuperLU_Dist 3.3
    // we either have a leak or a segfault here.
    // To be investigated further. 
    //Destroy_CompRowLoc_Matrix_dist(A);
#if defined(SLUD_VERSION_6)
    dScalePermstructFree(ScalePermstruct);
    dLUstructFree(LUstruct);
#else
    ScalePermstructFree(ScalePermstruct);
    LUstructFree(LUstruct);
#endif
    superlu_gridexit(grid);

    free(grid);
    free(LUstruct);
    free(LUfactors);
    return(0);

#else
    fprintf(stderr," SLUDist Not Configured, fix make.inc and recompile\n");
    return(-1);
#endif
}


