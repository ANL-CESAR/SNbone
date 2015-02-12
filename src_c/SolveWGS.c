//---------------------------------------------------------------------------------------------------------------------------------
// This subroutine drives the solution of the SN equations via GMRES
//---------------------------------------------------------------------------------------------------------------------------------
#include "ApplyA_functions.h"
#include "ApplyA_common.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef WITHOMP
#include <omp.h>
#endif

// Fortran interface routines
void solvewgs(int *Input_Iterations,int *IterationCount,int *iMethod,double *LHS_C,double *RHS_C) {
SolveWGS(Input_Iterations,IterationCount,iMethod,LHS_C,RHS_C);
}
void solvewgs_(int *Input_Iterations,int *IterationCount,int *iMethod,double *LHS_C,double *RHS_C) {
SolveWGS(Input_Iterations,IterationCount,iMethod,LHS_C,RHS_C);
}

void SolveWGS(int *Input_Iterations,int *IterationCount,int *iMethod,double *LHS_C,double *RHS_C) {
//PROTEUS_Int Input_Iterations,IterationCount,iMethod
//PROTEUS_Real LHS_C(NumAngles,NumVertices),RHS_C(NumAngles,NumVertices)

// Local
int  I,J,I_SizeVec,Output_Unit;
int  MyThreadID,iStart,iEnd;
int  ParallelComm,ParallelRank;
int  GuessIsNonZero;
int  ReasonForConvergence;       // Divergence (<0), MaxIterationCount (=0), Convergence (>0)
double ResidualNorm;               // Norm of the residual
// Additional Threaded stuff
double VectorNorm;
double *VectorNorm_Local;

I_SizeVec = NumAngles*NumVertices;
ParallelRank = 0;
ParallelComm = 0;

//printf("[SN-KERNEL]...Incoming variables 1 %5d \n",*IterationCount);
//printf("[SN-KERNEL]...Incoming variables 2 %5d \n",*iMethod);

GuessIsNonZero = 0;
Krylov_Maximum_Iterations = *Input_Iterations;
Krylov_Absolute_Tolerance = 1.0e-12;
Krylov_Relative_Tolerance = 1.0e-8;
if (*iMethod == 1) {printf("[SN-KERNEL]...Calling FGMRES solver for AVE1 \n");}
if (*iMethod == 2) {printf("[SN-KERNEL]...Calling FGMRES solver for AVE2 \n");}
if (*iMethod == 3) {printf("[SN-KERNEL]...Calling FGMRES solver for AVE3 \n");}
if (*iMethod == 4) {printf("[SN-KERNEL]...Calling FGMRES solver for AVE4 \n");}
if (*iMethod == 5) {printf("[SN-KERNEL]...Calling FGMRES solver for AVE5 \n");}

VectorNorm_Local = (double*) malloc(NumThreads*Krylov_BackVectors*sizeof(double));

#ifdef WITHOMP
#pragma omp parallel shared(ResidualNorm,VectorNorm,VectorNorm_Local,Conn,AS_ThreadWiseWork, \
                           DA_ThreadWiseWork,MM_ThreadWiseWork,LHS_C,RHS_C,Scratch_V1,Scratch_V2,NumThreads, \
                           GuessIsNonZero,ReasonForConvergence,IterationCount,ParallelComm,ParallelRank) \
                           private(MyThreadID,I,iStart,iEnd)
   {
   MyThreadID = omp_get_thread_num() + 1;
   I = (NumVertices*NumAngles)/NumThreads;
   iStart = (MyThreadID-1)*I + 1;
   if (MyThreadID == NumThreads) {iEnd = NumVertices*NumAngles;}
     else {iEnd = MyThreadID*I;}
#else
   MyThreadID = 1;
   iStart = 1;
   iEnd = NumVertices*NumAngles;
#endif
   //printf("[SN-KERNEL]...My thread %5d \n",MyThreadID);
   //printf("[SN-KERNEL]...iStart    %5d \n",iStart);
   //printf("[SN-KERNEL]...iEnd      %5d \n",iEnd);
       FGMRES_Threaded(&Output_Unit,
       &Krylov_Local_Owned,&Krylov_BackVectors,&Krylov_Maximum_Iterations,&Krylov_Iterations,
       &Krylov_Absolute_Tolerance,&Krylov_Relative_Tolerance,&Krylov_Divergence_Tolerance,
       Krylov_Basis,Krylov_Hessenberg,Krylov_Givens,Krylov_PC_Basis,Krylov_Modified_RHS,
       LHS_C,RHS_C,
       &MyThreadID,&iStart,&iEnd,
       &NumThreads,&GuessIsNonZero,&ReasonForConvergence,IterationCount,&ParallelComm,&ParallelRank,
       &ResidualNorm,&VectorNorm,VectorNorm_Local,VectorNorm_Local,
       iMethod);
#ifdef WITHOMP
   }
#endif
free(VectorNorm_Local);

printf("[SN-KERNEL]...FGMRES returned with an error of %13.6e after %6d iterations \n",ResidualNorm,*IterationCount);

} // END SUBROUTINE SolveWGS
