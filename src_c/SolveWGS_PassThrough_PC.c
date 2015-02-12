//---------------------------------------------------------------------------------------------------------------------------------
// This subroutine is a pass through point for GMRES
//---------------------------------------------------------------------------------------------------------------------------------
#include "ApplyA_common.h"
#include "ApplyA_functions.h"
#ifdef WITHOMP
#include <omp.h>
#endif
// Fortran interface routines
void solvewgs_passthrough_pc(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_PC(RHS_C, LHS_C);
}
void solvewgs_passthrough_pc_(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_PC(RHS_C, LHS_C);
}

// Main subroutine header
void SolveWGS_PassThrough_PC(double *RHS_C, double *LHS_C) {
// double LHS_C(NumAngles*NumVertices),RHS_C(NumAngles*NumVertices)
int I,J;
int Istart,Iend;

#ifdef WITHOMP
   I = omp_get_thread_num() + 1;
   J = NumAngles*NumVertices/NumThreads;
   Istart = (I-1)*J + 1;
   Iend   =  I*J;
   if (I == NumThreads) {Iend = NumAngles*NumVertices;}
#else
   Istart = 1;
   Iend   = NumAngles*NumVertices;
#endif
// No incoming barrier is needed as this preconditioner does not have any overlap, but we are by rule required to impose a barrier
for (I = Istart; I < Iend+1; I++) {
   LHS_C[I-1] = RHS_C[I-1]; // The identity operation
}
#ifdef WITHOMP
#pragma omp barrier
#endif

} 
