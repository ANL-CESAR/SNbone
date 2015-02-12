// This program generates a solution for the solver to obtain
// A * x = b   where x = LHS and b = RHS
// ----------------------------------------------------------------------------
#include <stdlib.h>
#include "ApplyA_common.h"
#include "ApplyA_functions.h"
//#define Local_Debug

// Fortran interface routines
void generatexb(int *Input_Scheme,double *LHS_C, double *LHS_Answer, double *RHS_C) {
GenerateXb(Input_Scheme,LHS_C,LHS_Answer,RHS_C);
}
void generatexb_(int *Input_Scheme,double *LHS_C, double *LHS_Answer, double *RHS_C) {
GenerateXb(Input_Scheme,LHS_C,LHS_Answer,RHS_C);
}

void GenerateXb(int *Input_Scheme,double *LHS_C, double *LHS_Answer, double *RHS_C) {
// LHS_C(NumAngles,NumVertices),LHS_Answer(NumAngles,NumVertices),RHS_C(NumAngles,NumVertices)
#ifdef WITHOMP
#include <omp.h>
#endif

// Local
int I,J,K,MyThreadID;

// We need to define the desired solution
for (I = 1; I < NumVertices+1; I++) {
   K = (I-1)*NumAngles-1;
   for (J = 1; J < NumAngles+1; J++) {
         LHS_C[K+J]=0.0;
         RHS_C[K+J]=0.0;
         LHS_Answer[K+J] = 0.2;
      }
   }
// Construct the source that goes with this answer
#ifdef WITHOMP
#pragma omp parallel shared(Conn,AS_ThreadWiseWork,LHS_Answer,RHS_C) private(MyThreadID)
   {
   MyThreadID = omp_get_thread_num() + 1;
#else
   MyThreadID = 1;
#endif
   if ((*Input_Scheme == 0) || (*Input_Scheme == 3)) {SolveWGS_PassThrough_AVE3_NoHPM(LHS_Answer,RHS_C);}
   else if (*Input_Scheme == 2) {SolveWGS_PassThrough_AVE2_NoHPM(LHS_Answer,RHS_C);}
   else {SolveWGS_PassThrough_AVE1_NoHPM(LHS_Answer,RHS_C);}
#ifdef WITHOMP
   }
#endif

#ifdef Local_Debug
   printf("Answer and Source \n");
   for (J = 1;J<NumAngles+1;J++) {
      for (I = 1;I<NumVertices+1;I++) {
         printf("%6d  %6d   %13.6f %13.6f \n",I,J,LHS_Answer[(I-1)*NumAngles+J-1],RHS_C[(I-1)*NumAngles+J-1]);
      }
   }
#endif

} 

