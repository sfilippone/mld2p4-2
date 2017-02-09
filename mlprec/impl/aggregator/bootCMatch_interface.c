#include <string.h>
#include <stdio.h>

#include "bcm.h"

bcm_CSRMatrix bootCMatch(bcm_CSRMatrix *C, int *match_algorithm, int *n_sweeps, int *max_nlevels, int *max_csize, bcm_Vector *w);
bcm_CSRMatrix bootCMatch(bcm_CSRMatrix *C, int *match_algorithm, int *n_sweeps, int *max_nlevels, int *max_csize, bcm_Vector *w){
   bcm_Vector *w_temp;
   int info;
   //double *w_inp;
   //w_inp=bcm_VectorData(w); 

   bcm_CSRMatrix *P;
   bcm_CSRMatrix *Ac;
   int ftcoarse=1;
   // Here I am building Ac but I won't use it. 
   Ac=bcm_CSRMatchingAgg(C, &w, &P, *match_algorithm, *n_sweeps, *max_nlevels,*max_csize , &ftcoarse);
   //w_inp=bcm_VectorData(w); 
 return *P;
}
