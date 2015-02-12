//---------------------------------------------------------------------------------------------------------------------------------
// This subroutine is a pass through point for GMRES
//---------------------------------------------------------------------------------------------------------------------------------
#include "ApplyA_common.h"
#include "ApplyA_functions.h"
#ifdef WITHOMP
#include <omp.h>
#endif
// Fortran interface routines
void solvewgs_passthrough_ave3(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE3(RHS_C, LHS_C);
}
void solvewgs_passthrough_ave3_(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE3(RHS_C, LHS_C);
}
void solvewgs_passthrough_ave3_nohpm(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE3_NoHPM(RHS_C, LHS_C);
}
void solvewgs_passthrough_ave3_nohpm_(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE3_NoHPM(RHS_C, LHS_C);
}

// Main subroutine header
void SolveWGS_PassThrough_AVE3(double *RHS_C, double *LHS_C) {
// double LHS_C(NumAngles*NumVertices),RHS_C(NumAngles*NumVertices)
int MyThreadID;
int I,J,K,iRowStart,iRowEnd,iV_A,iVV_A,iOff;

#ifdef WITHBGQHPM
   if (MyThreadID == 1) {hpm_start('AVE3_ApplyA');}
#endif

#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1;
   I = NumVertices/NumThreads;
   iRowStart = (MyThreadID-1)*I + 1;
   iRowEnd   = MyThreadID*I;
   if (MyThreadID == NumThreads) {iRowEnd = NumVertices;}
#else
   MyThreadID = 1;
   iRowStart = 1;
   iRowEnd   = NumVertices;
#endif
// This barrier ensures that the incoming vector is fully defined by all threads
#ifdef WITHOMP
#pragma omp barrier
#endif
for (I = iRowStart; I < iRowEnd+1; I++) {
   for (J = NZS_RowLoc[I-1]; J < NZS_RowLoc[I]; J++) {
      iV_A = (I-1)*NumAngles-1;
      iVV_A = (NZS_ColNum[J-1]-1)*NumAngles-1;
      iOff  = (J-1)*NumAngles-1;
      for (K = 1; K < NumAngles+1; K++) {
         LHS_C[iV_A+K] = LHS_C[iV_A+K] + NZS_Data[iOff+K]*RHS_C[iVV_A+K];
      }
   }
}
// This second barrier is needed because the FGMRES threadwise splitting is currently different than the above. Likely can change that
#ifdef WITHOMP
#pragma omp barrier
#endif
#ifdef WITHBGQHPM
   IF (MyThreadID .EQ. 1) call hpm_stop('AVE3_ApplyA') ! Stops the hardware counters
#endif
} 

void SolveWGS_PassThrough_AVE3_NoHPM(double *RHS_C, double *LHS_C) {
// double LHS_C(NumAngles*NumVertices),RHS_C(NumAngles*NumVertices)
int MyThreadID;
int I,J,K,iRowStart,iRowEnd,iV_A,iVV_A,iOff;

#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1;
   I = NumVertices/NumThreads;
   iRowStart = (MyThreadID-1)*I + 1;
   iRowEnd   = MyThreadID*I;
   if (MyThreadID == NumThreads) {iRowEnd = NumVertices;}
#else
   MyThreadID = 1;
   iRowStart = 1;
   iRowEnd   = NumVertices;
#endif
// This barrier ensures that the incoming vector is fully defined by all threads
#ifdef WITHOMP
#pragma omp barrier
#endif
for (I = iRowStart; I < iRowEnd+1; I++) {
   for (J = NZS_RowLoc[I-1]; J < NZS_RowLoc[I]; J++) {
      iV_A = (I-1)*NumAngles-1;
      iVV_A = (NZS_ColNum[J-1]-1)*NumAngles-1;
      iOff  = (J-1)*NumAngles-1;
      for (K = 1; K < NumAngles+1; K++) {
         LHS_C[iV_A+K] = LHS_C[iV_A+K] + NZS_Data[iOff+K]*RHS_C[iVV_A+K];
      }
   }
}
// This second barrier is needed because the FGMRES threadwise splitting is currently different than the above. Likely can change that
#ifdef WITHOMP
#pragma omp barrier
#endif
}
