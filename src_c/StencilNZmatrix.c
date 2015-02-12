//---------------------------------------------------------------------------------------------------------------------------------
// This subroutine stencils the spatially connected FE system
//---------------------------------------------------------------------------------------------------------------------------------
#include "ApplyA_common.h"
#include "ApplyA_functions.h"
//#define Local_Debug
//#define Debug_DumpAssembledMatrix
//#define Local_DumpDebugNZS
//#define Local_DebugNZS
#ifdef WITHOMP
#include <omp.h>
#endif
// Fortran interface routines
void stencilnzmatrix(int *Input_Scheme) {
StencilNZmatrix(Input_Scheme);
}
void stencilnzmatrix_(int *Input_Scheme) {
StencilNZmatrix(Input_Scheme);
}
// Main subroutine header
void StencilNZmatrix(int *Input_Scheme) {
// Local variables
int I,J,K,L,M;
int Istart,IBatches,II,FixedBatchSize;
// Arrays needed to evaluate the derivatives
int     *Local_VertexFlag;
int     *Temp_NZS_ColNum;

if ((*Input_Scheme == 0) || (*Input_Scheme == 3)) {
   FixedBatchSize = 30;
   //printf("initial status of Local_VertexFlag %10d \n",Local_VertexFlag);
   //printf("initial status of NZS_RowLoc       %10d \n",NZS_RowLoc);
   Local_VertexFlag = (int*) malloc(NumVertices*FixedBatchSize*sizeof(int));
   NZS_RowLoc       = (int*) malloc((NumVertices+1)*sizeof(int));
   NZS_RowLoc[1-1]  = 1;

   IBatches = NumVertices/FixedBatchSize;
   if (IBatches * FixedBatchSize != NumVertices) {IBatches = IBatches + 1;}
   Istart    = 0;

   // Determine the maximum number of vertex connections per vertex

   // We need to stencil the spatial matrix and for (so by sweeping the connectivity matrix V times
   for (I = 1;I<IBatches+1;I++) {
      if (Istart+FixedBatchSize > NumVertices) {FixedBatchSize = NumVertices - Istart;}
      for (II = 1;II<FixedBatchSize+1;II++) {
         for (J = 1;J<NumVertices+1;J++) {Local_VertexFlag[(II-1)*NumVertices+J-1] = 0;} // initialize to 0
      }
      // Flag all of the vertices that are connected
      for (K = 1;K<NumElements+1;K++) {
         for (L = 1;L<FEVertices+1;L++) {
            for (II = 1;II<FixedBatchSize+1;II++) {
               if (Conn[(K-1)*4+L-1] == Istart+II) {
                  for (M = 1;M<FEVertices+1;M++) {
                     Local_VertexFlag[(II-1)*NumVertices+Conn[(K-1)*4+M-1]-1] = 1; // flag the vertices we touched
                  }
               }
            }
         }
      }

      // Update NZS_NonZeros and store the rowloc information
      K = NZS_NonZeros;
      for (II = 1;II<FixedBatchSize+1;II++) {
         for (J = 1;J<NumVertices+1;J++) {
            if (Local_VertexFlag[(II-1)*NumVertices+J-1] == 1) {K = K + 1;}
         }
         NZS_RowLoc[Istart+II+1-1] = K+1; // Start of the next row
      }
      // Resize ColNum
      if (NZS_NonZeros == 0) 
         {NZS_ColNum = (int*) malloc(K*sizeof(int));}
      else {
         Temp_NZS_ColNum = (int*) malloc(NZS_NonZeros*sizeof(int)); 
         for (J = 1;J<NZS_NonZeros+1;J++) {Temp_NZS_ColNum[J-1] = NZS_ColNum[J-1];}
         free(NZS_ColNum);
         NZS_ColNum = (int*) malloc(K*sizeof(int));
         for (J = 1;J<NZS_NonZeros+1;J++) {NZS_ColNum[J-1] = Temp_NZS_ColNum[J-1];}
         free(Temp_NZS_ColNum);
      }
      // Store ColNum for this set of data
      K = NZS_NonZeros;
      for (II = 1;II<FixedBatchSize+1;II++) {
         for (J = 1;J<NumVertices+1;J++) {
            if (Local_VertexFlag[(II-1)*NumVertices+J-1] == 1) {
               K = K + 1;
               NZS_ColNum[K-1] = J;
            }
         }
      }
      NZS_NonZeros = K;
      Istart = Istart + FixedBatchSize;
   } // End of Batch of vertices

   NZS_Data = (double *) malloc(NumAngles*NZS_NonZeros*sizeof(double));

   for (I = 1;I<NZS_NonZeros*NumAngles+1;I++) {
      NZS_Data[I-1] = 0.0;
      }
   free(Local_VertexFlag);
   }
else {
   NZS_RowLoc = (int*) malloc(1*sizeof(int));
   NZS_ColNum = (int*) malloc(1*sizeof(int));
   NZS_Data   = (double *) malloc(1*sizeof(double));
   } // Scheme 0 and 3 
} // END SUBROUTINE StencilNZmatrix
