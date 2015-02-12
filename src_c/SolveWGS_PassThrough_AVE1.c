//---------------------------------------------------------------------------------------------------------------------------------
// This subroutine is a pass through point for GMRES
//---------------------------------------------------------------------------------------------------------------------------------
#include "ApplyA_common.h"
#include "ApplyA_functions.h"
#ifdef WITHOMP
#include <omp.h>
#endif
// Fortran interface routines
void solvewgs_passthrough_ave1(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE1(RHS_C, LHS_C);
}
void solvewgs_passthrough_ave1_(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE1(RHS_C, LHS_C);
}
void solvewgs_passthrough_ave1_nohpm(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE1_NoHPM(RHS_C, LHS_C);
}
void solvewgs_passthrough_ave1_nohpm_(double *RHS_C, double *LHS_C) {
SolveWGS_PassThrough_AVE1_NoHPM(RHS_C, LHS_C);
}

// Main subroutine header
void SolveWGS_PassThrough_AVE1(double *RHS_C, double *LHS_C) {
// double LHS_C(NumAngles*NumVertices),RHS_C(NumAngles*NumVertices)

int MyThreadID,K;

#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1;
#else
   MyThreadID = 1;
#endif
#ifdef WITHBGQHPM
   if (MyThreadID == 1) {hpm_start('AVE1_ApplyA');}
#endif

K = (MyThreadID-1)*AS_NumColors*2 + 1 - 1;

// This barrier ensures that the incoming vector is fully defined by all threads
#ifdef WITHOMP
#pragma omp barrier
#endif
ApplyA_AVE1_Tet_SUPG(&NumElements,&NumAngles,&NumVertices,Conn,
                     &AS_NumColors,&AS_NumThreads,&TasksPerThread,&MyThreadID,AS_ThreadWiseWork,
                              ConstF, ConstU, ConstUT, FEShapeFunctions,FEDerivatives,FEDetJacandWgt, Omega,OmegaOmega,
                              LHS_C,RHS_C);

} 


// Main subroutine header
void SolveWGS_PassThrough_AVE1_NoHPM(double *RHS_C, double *LHS_C) {
// double LHS_C(NumAngles*NumVertices),RHS_C(NumAngles*NumVertices)

int MyThreadID,K;

#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1;
#else
   MyThreadID = 1;
#endif

K = (MyThreadID-1)*AS_NumColors*2 + 1 - 1;

// This barrier ensures that the incoming vector is fully defined by all threads
#ifdef WITHOMP
#pragma omp barrier
#endif
ApplyA_AVE1_Tet_SUPG(&NumElements,&NumAngles,&NumVertices,Conn,
                     &AS_NumColors,&AS_NumThreads,&TasksPerThread,&MyThreadID,AS_ThreadWiseWork,
                              ConstF, ConstU, ConstUT, FEShapeFunctions,FEDerivatives,FEDetJacandWgt, Omega,OmegaOmega,
                              LHS_C,RHS_C);

}
