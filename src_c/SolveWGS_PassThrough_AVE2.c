//---------------------------------------------------------------------------------------------------------------------------------
// This subroutine is a pass through point for GMRES
//---------------------------------------------------------------------------------------------------------------------------------
#include "ApplyA_common.h"
#include "ApplyA_functions.h"
#ifdef WITHOMP
#include <omp.h>
#endif
// Fortran interface routines
void solvewgs_passthrough_ave2(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE2(RHS_C, LHS_C);
}
void solvewgs_passthrough_ave2_(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE2(RHS_C, LHS_C);
}
void solvewgs_passthrough_ave2_nohpm(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE2_NoHPM(RHS_C, LHS_C);
}
void solvewgs_passthrough_ave2_nohpm_(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE2_NoHPM(RHS_C, LHS_C);
}

// Main subroutine header
void SolveWGS_PassThrough_AVE2(double *RHS_C, double *LHS_C) {
// double LHS_C(NumAngles*NumVertices),RHS_C(NumAngles*NumVertices)

int MyThreadID,K,I;

#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1;
#else
   MyThreadID = 1;
#endif
#ifdef WITHBGQHPM
   if (MyThreadID == 1) {hpm_start('AVE2_ApplyA');}
#endif

I = (MyThreadID-1)*NumAngles*FEVertices+1-1;
K = (MyThreadID-1)*AS_NumColors*2 + 1 - 1;
// printf("[SN-KERNEL]...MyThreadID= %6d I= %6d ",MyThreadID,I);
// This barrier ensures that the incoming vector is fully defined by all threads
#ifdef WITHOMP
#pragma omp barrier
#endif
ApplyA_AVE2_Tet_SUPG(&NumElements,&NumAngles,&NumVertices,Conn,
                     &AS_NumColors,&AS_NumThreads,&TasksPerThread,&MyThreadID,AS_ThreadWiseWork,
                              ConstF, ConstU, ConstUT, FEShapeFunctions,FEDerivatives,FEDetJacandWgt, Omega,OmegaOmega,
                              LHS_C,RHS_C,&Scratch_V1[I],&Scratch_V2[I]);

} 

// Main subroutine header
void SolveWGS_PassThrough_AVE2_NoHPM(double *RHS_C, double *LHS_C) {
// double LHS_C(NumAngles*NumVertices),RHS_C(NumAngles*NumVertices)

int MyThreadID,K,I;

#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1;
#else
   MyThreadID = 1;
#endif

I = (MyThreadID-1)*NumAngles*FEVertices+1-1;
K = (MyThreadID-1)*AS_NumColors*2 + 1 - 1;

//printf("[SN-KERNEL]...MyThreadID= %6d I= %6d ",MyThreadID,I);

// This barrier ensures that the incoming vector is fully defined by all threads
#ifdef WITHOMP
#pragma omp barrier
#endif
ApplyA_AVE2_Tet_SUPG(&NumElements,&NumAngles,&NumVertices,Conn,
                     &AS_NumColors,&AS_NumThreads,&TasksPerThread,&MyThreadID,AS_ThreadWiseWork,
                              ConstF, ConstU, ConstUT, FEShapeFunctions,FEDerivatives,FEDetJacandWgt, Omega,OmegaOmega,
                              LHS_C,RHS_C,&Scratch_V1[I],&Scratch_V2[I]);

} 
