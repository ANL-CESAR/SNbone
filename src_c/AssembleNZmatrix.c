//---------------------------------------------------------------------------------------------------------------------------------
// This subroutine sets up the spatial and angular matrices
// It also forms the assembled matrix-vector system
//---------------------------------------------------------------------------------------------------------------------------------
#include <stdlib.h>
#include "ApplyA_common.h"
#include "ApplyA_functions.h"
//#define Local_Debug
//#define Debug_DumpAssembledMatrix
//#define Local_DumpDebugNZS
//#define Local_DebugNZS
void assemblenzmatrix(int *Input_Scheme) {
AssembleNZmatrix(Input_Scheme);
}
void assemblenzmatrix_(int *Input_Scheme) {
AssembleNZmatrix(Input_Scheme);
}

// Main subroutine header
void AssembleNZmatrix(int *Input_Scheme) {
// Local variables
int I,J,K;
int I1,I2,I3,I4;
double ConstCoef;

if (*Input_Scheme == 0 || *Input_Scheme == 3)  {
// Threadable section eventually
//   DO iColor = 1,AS_NumColors
//      DO iTask = 1,TasksPerThread
//         iThread = ThreadID*TasksPerThread + iTask
//         DO K = AS_ThreadWiseWork(1,iThread,iColor),AS_ThreadWiseWork(2,iThread,iColor)

   // Form and store the matrix
   for (K = 1;K<NumElements+1;K++) {
      // Row 1 ---------------------------------------------------------------
      I = Conn[(K-1)*FEVertices+1-1]; // The row vertex
      for (J = NZS_RowLoc[I-1];J<NZS_RowLoc[I+1-1];J++) {
         if (Conn[(K-1)*FEVertices+1-1] == NZS_ColNum[J-1]) {I1 = J;}
         if (Conn[(K-1)*FEVertices+2-1] == NZS_ColNum[J-1]) {I2 = J;}
         if (Conn[(K-1)*FEVertices+3-1] == NZS_ColNum[J-1]) {I3 = J;}
         if (Conn[(K-1)*FEVertices+4-1] == NZS_ColNum[J-1]) {I4 = J;}
      }
      // 1,1
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = 2.0*fcoef[K-1]*ConstF[K-1]                                                            
         - ConstU[K-1]*((Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+1-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+2-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+3-1])   
                 +(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+4-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+5-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+6-1])   
                 +(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+7-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+8-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+9-1]))  
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+1-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+2-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+3-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+4-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+5-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+6-1])     
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+7-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+8-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+9-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+10-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+11-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+12-1])     
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+13-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+14-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+15-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+16-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+17-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+18-1])     
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+19-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+20-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+21-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+22-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+23-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+24-1])     
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+25-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+26-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+27-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+28-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+29-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+30-1])     
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+31-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+32-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+33-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+34-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+35-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+36-1])    ;
         NZS_Data[(I1-1)*NumAngles+J-1] = NZS_Data[(I1-1)*NumAngles+J-1] + ConstCoef;
      }
      // 1,2
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                                                  
         + ConstU[K-1]* (Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+1-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+2-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+3-1])   
         - ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+1-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+2-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+3-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+4-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+5-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+6-1]);
         NZS_Data[(I2-1)*NumAngles+J-1] = NZS_Data[(I2-1)*NumAngles+J-1] + ConstCoef;
      }
      // 1,3
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                                                  
         + ConstU[K-1]* (Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+4-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+5-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+6-1])   
         - ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+7-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+8-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+9-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+10-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+11-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+12-1])     
         - ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+19-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+20-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+21-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+22-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+23-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+24-1]);
         NZS_Data[(I3-1)*NumAngles+J-1] = NZS_Data[(I3-1)*NumAngles+J-1] + ConstCoef;
      }
      // 1,4
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                                                  
         + ConstU[K-1]* (Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+7-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+8-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+9-1])   
         - ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+13-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+14-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+15-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+16-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+17-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+18-1])     
         - ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+19-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+20-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+21-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+22-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+23-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+24-1])     
         - ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+25-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+26-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+27-1]     
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+28-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+29-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+30-1])    
         - ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+31-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+32-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+33-1]     
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+34-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+35-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+36-1])    ;
         NZS_Data[(I4-1)*NumAngles+J-1] = NZS_Data[(I4-1)*NumAngles+J-1] + ConstCoef;
      }
      // Row 2 ---------------------------------------------------------------
      I = Conn[(K-1)*FEVertices+2-1]; // The row vertex
      for (J = NZS_RowLoc[I-1];J<NZS_RowLoc[I+1-1];J++) {
         if (Conn[(K-1)*FEVertices+1-1] == NZS_ColNum[J-1]) {I1 = J;}
         if (Conn[(K-1)*FEVertices+2-1] == NZS_ColNum[J-1]) {I2 = J;}
         if (Conn[(K-1)*FEVertices+3-1] == NZS_ColNum[J-1]) {I3 = J;}
         if (Conn[(K-1)*FEVertices+4-1] == NZS_ColNum[J-1]) {I4 = J;}
      }
      // 2,1
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                           
         - ConstU[K-1]*((Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+1-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+2-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+3-1])    
                 +(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+4-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+5-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+6-1])    
                 +(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+7-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+8-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+9-1]))   
         - ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+1-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+2-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+3-1]       
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+4-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+5-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+6-1])      
         - ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+19-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+20-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+21-1]       
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+22-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+23-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+24-1])      
         - ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+25-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+26-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+27-1]       
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+28-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+29-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+30-1]);
         NZS_Data[(I1-1)*NumAngles+J-1] = NZS_Data[(I1-1)*NumAngles+J-1] + ConstCoef;
      }
      // 2,2
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = 2.0*fcoef[K-1]*ConstF[K-1]                                                                  
         + ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+1-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+2-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+3-1])   
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+1-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+2-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+3-1]     
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+4-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+5-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+6-1]);
         NZS_Data[(I2-1)*NumAngles+J-1] = NZS_Data[(I2-1)*NumAngles+J-1] + ConstCoef;
      }
      // 2,3
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                                                  
         + ConstU[K-1]* (Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+4-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+5-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+6-1])   
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+19-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+20-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+21-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+22-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+23-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+24-1]);
         NZS_Data[(I3-1)*NumAngles+J-1] = NZS_Data[(I3-1)*NumAngles+J-1] + ConstCoef;
      }
      // 2,4
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                                                  
         + ConstU[K-1]* (Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+7-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+8-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+9-1])   
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+25-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+26-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+27-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+28-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+29-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+30-1]);
         NZS_Data[(I4-1)*NumAngles+J-1] = NZS_Data[(I4-1)*NumAngles+J-1] + ConstCoef;
      }
      // Row 3 ---------------------------------------------------------------
      I = Conn[(K-1)*FEVertices+3-1]; // The row vertex
      for (J = NZS_RowLoc[I-1];J<NZS_RowLoc[I+1-1];J++) {
         if (Conn[(K-1)*FEVertices+1-1] == NZS_ColNum[J-1]) {I1 = J;}
         if (Conn[(K-1)*FEVertices+2-1] == NZS_ColNum[J-1]) {I2 = J;}
         if (Conn[(K-1)*FEVertices+3-1] == NZS_ColNum[J-1]) {I3 = J;}
         if (Conn[(K-1)*FEVertices+4-1] == NZS_ColNum[J-1]) {I4 = J;}
      }
      // 3,1
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                           
         - (ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+1-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+2-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+3-1])    
           +ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+4-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+5-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+6-1])    
           +ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+7-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+8-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+9-1]))   
         -  ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+7-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+8-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+9-1]      
                  +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+10-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+11-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+12-1])     
         -  ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+31-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+32-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+33-1]      
                  +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+34-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+35-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+36-1]);
         NZS_Data[(I1-1)*NumAngles+J-1] = NZS_Data[(I1-1)*NumAngles+J-1] + ConstCoef;
      }
      // 3,2
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                          
         + ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+1-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+2-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+3-1]);
         NZS_Data[(I2-1)*NumAngles+J-1] = NZS_Data[(I2-1)*NumAngles+J-1] + ConstCoef;
      }
      // 3,3
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = 2.0*fcoef[K-1]*ConstF[K-1]                                    
         + ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+4-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+5-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+6-1])    
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+7-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+8-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+9-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+10-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+11-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+12-1]) ;
         NZS_Data[(I3-1)*NumAngles+J-1] = NZS_Data[(I3-1)*NumAngles+J-1] + ConstCoef;
      }
      // 3,4
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                          
         + ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+7-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+8-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+9-1])    
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+31-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+32-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+33-1]      
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+34-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+35-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+36-1]);
         NZS_Data[(I4-1)*NumAngles+J-1] = NZS_Data[(I4-1)*NumAngles+J-1] + ConstCoef;
      }
      // Row 4 ---------------------------------------------------------------
      I = Conn[(K-1)*FEVertices+4-1]; // The row vertex
      for (J = NZS_RowLoc[I-1];J<NZS_RowLoc[I+1-1];J++) {
         if (Conn[(K-1)*FEVertices+1-1] == NZS_ColNum[J-1]) {I1 = J;}
         if (Conn[(K-1)*FEVertices+2-1] == NZS_ColNum[J-1]) {I2 = J;}
         if (Conn[(K-1)*FEVertices+3-1] == NZS_ColNum[J-1]) {I3 = J;}
         if (Conn[(K-1)*FEVertices+4-1] == NZS_ColNum[J-1]) {I4 = J;}
      }
      //WRITE(6,'("I1-I4",10I6)') I,I1,I2,I3,I4
      // 4,1
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                           
         - (ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+1-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+2-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+3-1])    
           +ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+4-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+5-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+6-1])    
           +ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+7-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+8-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+9-1]))   
         -  ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+13-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+14-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+15-1]   
                  +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+16-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+17-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+18-1]);
         NZS_Data[(I1-1)*NumAngles+J-1] = NZS_Data[(I1-1)*NumAngles+J-1] + ConstCoef;
      }
      // 4,2
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                          
         + ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+1-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+2-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+3-1]);
         NZS_Data[(I2-1)*NumAngles+J-1] = NZS_Data[(I2-1)*NumAngles+J-1] + ConstCoef;
      }
      // 4,3
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = fcoef[K-1]*ConstF[K-1]                                                       
         + ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+4-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+5-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+6-1]);
         NZS_Data[(I3-1)*NumAngles+J-1] = NZS_Data[(I3-1)*NumAngles+J-1] + ConstCoef;
      }
      // 4,4
      for (J = 1;J<NumAngles+1;J++) {
         ConstCoef = 2.0*fcoef[K-1]*ConstF[K-1]                                                                      
         + ConstU[K-1]*(Omega[NumAngles*(1-1)+J-1]*ucoef[(K-1)*9+7-1] + Omega[NumAngles*(2-1)+J-1]*ucoef[(K-1)*9+8-1] + Omega[NumAngles*(3-1)+J-1]*ucoef[(K-1)*9+9-1])                 
         + ConstTau[K-1]*(OmegaOmega[NumAngles*(1-1)+J-1]*pcoef[(K-1)*36+13-1]+OmegaOmega[NumAngles*(2-1)+J-1]*pcoef[(K-1)*36+14-1]+OmegaOmega[NumAngles*(3-1)+J-1]*pcoef[(K-1)*36+15-1]  
                 +OmegaOmega[NumAngles*(4-1)+J-1]*pcoef[(K-1)*36+16-1]+OmegaOmega[NumAngles*(5-1)+J-1]*pcoef[(K-1)*36+17-1]+OmegaOmega[NumAngles*(6-1)+J-1]*pcoef[(K-1)*36+18-1])    ;
         NZS_Data[(I4-1)*NumAngles+J-1] = NZS_Data[(I4-1)*NumAngles+J-1] + ConstCoef;
      }
   } // Element loop

#ifdef Local_DumpDebugNZS
   {DO I = 1,NumVertices
      {DO K = NZS_RowLoc[I-1],NZS_RowLoc[I+1-1]-1
         for (J = 1;J<NumAngles+1;J++) {
            WRITE(23,777) 0,I,NZS_ColNum[K-1],J,NZS_Data(J,K)
         }
      }
   }
#endif

#ifdef Local_DebugNZS
   for (I = 1;I<NumVertices+1;I++) {
      printf("[SN-KERNEL] Row %8d stores data between %6d:%6d \n",I,NZS_RowLoc[I-1],NZS_RowLoc[I]-1);
      for (K=1;K<NZS_RowLoc[I]-NZS_RowLoc[I-1]+1;K++) {printf("    %3d    ",K);}
      printf("\n");
      for (K=NZS_RowLoc[I-1];K<NZS_RowLoc[I]-1+1;K++) {printf("    %3d    ",NZS_ColNum[K-1]);}
      printf("\n");
      for (J=1;J<NumAngles+1;J++) {
         for (K=NZS_RowLoc[I-1];K<NZS_RowLoc[I]-1+1;K++) {printf(" %9.3f ",NZS_Data[(K-1)*NumAngles+J-1]);}
         printf("\n");
      }
   }
#endif

   //{DO I = 1,NumVertices
   //   WRITE(Output_Unit,'("[SN-KERNEL] Row ",I8," stores data between ",I6,":",I6)') I,NZS_RowLoc[I-1],NZS_RowLoc[I+1-1]-1
   //}
} // Scheme 0 and 3 

} // END SUBROUTINE AssembleNZmatrix
