// This program checks the solution
// ----------------------------------------------------------------------------
#include "ApplyA_functions.h"
#include <stdio.h>
// Fortran interface routines
void zero(int *iStart,int *iEnd,int *iSize,double *VectorToInitialize) {
Zero(iStart,iEnd,iSize,VectorToInitialize);
}
void zero_(int *iStart,int *iEnd,int *iSize,double *VectorToInitialize) {
Zero(iStart,iEnd,iSize,VectorToInitialize);
}
// Main subroutine header
void Zero(int *iStart,int *iEnd,int *iSize,double *VectorToInitialize) {
//PROTEUS_Int,  INTENT(IN)    :: iStart,iEnd,iSize
//PROTEUS_Real, INTENT(INOUT) :: VectorToInitialize(iSize)
int I;

   for (I = *iStart; I < *iEnd+1; I++) { // DO I = iStart,iEnd
      VectorToInitialize[I-1] = 0.0;
   } // END DO

} // END SUBROUTINE Zero

void zero_threaded(int *iSize,double *VectorToInitialize) {
Zero_Threaded(iSize,VectorToInitialize);
}
void zero_threaded_(int *iSize,double *VectorToInitialize) {
Zero_Threaded(iSize,VectorToInitialize);
}
// Main subroutine header
void Zero_Threaded(int *iSize,double *VectorToInitialize) {
#ifdef WITHOMP
#include <omp.h>
#endif
//IMPLICIT NONE
//PROTEUS_Int,  INTENT(IN)    :: iSize
//PROTEUS_Real, INTENT(INOUT) :: VectorToInitialize(iSize)
// Local
int I;

#ifdef WITHOMP
#pragma omp parallel for schedule(static)
#endif
   for (I = 1; I < *iSize+1; I++) { // DO I = 1,iSize
      VectorToInitialize[I-1] = 0.0;
   } // END DO

} // END SUBROUTINE Zero_Threaded
