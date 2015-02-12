//---------------------------------------------------------------------------------------------------------------------------------
// This subroutine sets up the spatial and angular matrices
// It also forms the assembled matrix-vector system
//---------------------------------------------------------------------------------------------------------------------------------
#include "ApplyA_functions.h"
#include <stdio.h>
// Fortran interface routines
void buildanglecubature(int *NumAngles,double *Omega,double *OmegaOmega,double *AngleWeights) {
BuildAngleCubature(NumAngles,Omega,OmegaOmega,AngleWeights);
}
void buildanglecubature_(int *NumAngles,double *Omega,double *OmegaOmega,double *AngleWeights) {
BuildAngleCubature(NumAngles,Omega,OmegaOmega,AngleWeights);
}
// Main subroutine header
void BuildAngleCubature(int *NumAngles,double *Omega,double *OmegaOmega,double *AngleWeights) {
// Passed in
//int  NumAngles
//PROTEUS_Real Omega(NumAngles,3),OmegaOmega(NumAngles,6),AngleWeights(NumAngles)

// Local variables
int I;
double Om1,Om2,Om3,SizeOm;
double ConstCoef;

// Initialize Omega with made-up numbers
ConstCoef = *NumAngles;
ConstCoef = 1.0 / ConstCoef;
for (I = 1;I<*NumAngles+1;I++) {
   Om1 = (I-1)*(I-1);
   Om2 = (I+1)*(I+1);
   Om3 = -2*(I*I)-1;
   SizeOm = sqrt( Om1*Om1 + Om2*Om2 + Om3*Om3);
   Omega[             I-1] = Om1 / SizeOm;
   Omega[*NumAngles  +I-1] = Om2 / SizeOm;
   Omega[*NumAngles*2+I-1] = Om3 / SizeOm;
   AngleWeights[I-1] = ConstCoef;
   //WRITE(*,'(I,I4,3F10.3)') 'Direction',I,Omega(I,1), Omega(I,2), Omega(I,3)  
   OmegaOmega[             I-1] = Omega[             I-1] * Omega[             I-1];  // 1 is 1*1
   OmegaOmega[*NumAngles  +I-1] = Omega[*NumAngles  +I-1] * Omega[*NumAngles  +I-1];  // 2 is 2*2
   OmegaOmega[*NumAngles*2+I-1] = Omega[*NumAngles*2+I-1] * Omega[*NumAngles*2+I-1];  // 3 is 3*3
   OmegaOmega[*NumAngles*3+I-1] = Omega[             I-1] * Omega[*NumAngles  +I-1];  // 4 is 1*2
   OmegaOmega[*NumAngles*4+I-1] = Omega[             I-1] * Omega[*NumAngles*2+I-1];  // 5 is 1*3
   OmegaOmega[*NumAngles*5+I-1] = Omega[*NumAngles  +I-1] * Omega[*NumAngles*2+I-1];  // 6 is 2*3
} // END DO

} // END SUBROUTINE BuildAngleCubature
