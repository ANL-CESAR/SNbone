//----------------------------------------------------------------------------------------------------
// Method_FGMRES.F90
// Driver of the FGMRES Routine to solve a linear system Ax = b
// To implement GMRES I followed the algorithm proposed by Y. Saad and M Schultz in GMRES: A generalized minimal residual 
// algorithm for solving nonsymmetric linear systems, SIAM J. Sci. Stat. Comput., Vol 7, No 3, 
// July 1986
// Then I edited it to add a varying right preconditioner (Flexible GMRES) following the paper 
// 'A set of flexible GMRES routines for real and complex arithmetics on high Performance computers', 
// V. Fraysse, L. Giraud, S. Gratton, CERFACS Technical Report TR/PA/06/09
//----------------------------------------------------------------------------------------------------
//  Copyright(c) 2005 Argonne National Laboratory
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
// Externally threaded version of the subroutine
//----------------------------------------------------------------------------------------------------
// Barriers if needed must be erected around all user subroutines
// Incoming vectors are shared and users must put barriers in place to prevent issues
// This routine supports combined MPI and OpenMP so long as MyThreadID=0 matches the MPI rank
// 
// Because this subroutine is called by each thread, all local variables will be unique to each thread
//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------
#include <math.h> 
#include "ApplyA_functions.h"
#include <stdio.h>
#ifdef WITHOMP
#include <omp.h>
#endif
//#define Local_Debug_FGMRES_Driver
// Fortran interface routines
void fgmres_threaded(int *Output_Unit,
  int *Krylov_Local_Owned,int *Krylov_BackVectors,int *Krylov_Maximum_Iterations,int *Krylov_Iterations,
  double *Krylov_Absolute_Tolerance,double *Krylov_Relative_Tolerance,double *Krylov_Divergence_Tolerance,
  double *Krylov_Basis,double *Krylov_Hessenberg,double *Krylov_Givens,double *Krylov_PC_Basis, double *Krylov_Modified_RHS,
  double *Solution,double *RightHandSide,
  int *MyThreadID,int *MyStart,int *MyEnd,
  int *NumThreads,int *GuessIsNonZero,int *ReasonForConvergence,int *IterationCount,int *ParallelComm,int *ParallelRank,
  double *ResidualNorm,double *VectorNorm,double *VectorNorm_Local,double *HessenNorm_Local,
  int *iMethod) {
FGMRES_Threaded(Output_Unit,
  Krylov_Local_Owned,Krylov_BackVectors,Krylov_Maximum_Iterations,Krylov_Iterations,
  Krylov_Absolute_Tolerance,Krylov_Relative_Tolerance,Krylov_Divergence_Tolerance,
  Krylov_Basis,Krylov_Hessenberg,Krylov_Givens,Krylov_PC_Basis,Krylov_Modified_RHS,
  Solution,RightHandSide,
  MyThreadID,MyStart,MyEnd,
  NumThreads,GuessIsNonZero,ReasonForConvergence,IterationCount,ParallelComm,ParallelRank,
  ResidualNorm,VectorNorm,VectorNorm_Local,HessenNorm_Local,
  iMethod);
}
void fgmres_threaded_(int *Output_Unit,
  int *Krylov_Local_Owned,int *Krylov_BackVectors,int *Krylov_Maximum_Iterations,int *Krylov_Iterations,
  double *Krylov_Absolute_Tolerance,double *Krylov_Relative_Tolerance,double *Krylov_Divergence_Tolerance,
  double *Krylov_Basis,double *Krylov_Hessenberg,double *Krylov_Givens,double *Krylov_PC_Basis, double *Krylov_Modified_RHS,
  double *Solution,double *RightHandSide,
  int *MyThreadID,int *MyStart,int *MyEnd,
  int *NumThreads,int *GuessIsNonZero,int *ReasonForConvergence,int *IterationCount,int *ParallelComm,int *ParallelRank,
  double *ResidualNorm,double *VectorNorm,double *VectorNorm_Local,double *HessenNorm_Local,
  int *iMethod) {
FGMRES_Threaded(Output_Unit,
  Krylov_Local_Owned,Krylov_BackVectors,Krylov_Maximum_Iterations,Krylov_Iterations,
  Krylov_Absolute_Tolerance,Krylov_Relative_Tolerance,Krylov_Divergence_Tolerance,
  Krylov_Basis,Krylov_Hessenberg,Krylov_Givens,Krylov_PC_Basis,Krylov_Modified_RHS,
  Solution,RightHandSide,
  MyThreadID,MyStart,MyEnd,
  NumThreads,GuessIsNonZero,ReasonForConvergence,IterationCount,ParallelComm,ParallelRank,
  ResidualNorm,VectorNorm,VectorNorm_Local,HessenNorm_Local,
  iMethod);
}
// Main subroutine header
void FGMRES_Threaded(int *Output_Unit,
  int *Krylov_Local_Owned,int *Krylov_BackVectors,int *Krylov_Maximum_Iterations,int *Krylov_Iterations,
  double *Krylov_Absolute_Tolerance,double *Krylov_Relative_Tolerance,double *Krylov_Divergence_Tolerance,
  double *Krylov_Basis,double *Krylov_Hessenberg,double *Krylov_Givens,double *Krylov_PC_Basis, double *Krylov_Modified_RHS,
  double *Solution,double *RightHandSide,
  int *MyThreadID,int *MyStart,int *MyEnd,
  int *NumThreads,int *GuessIsNonZero,int *ReasonForConvergence,int *IterationCount,int *ParallelComm,int *ParallelRank,
  double *ResidualNorm,double *VectorNorm,double *VectorNorm_Local,double *HessenNorm_Local,
  int *iMethod) {
#ifdef USEMPI
#include "mpif.h"
#endif
//#define Local_Debug_FGMRES_Driver
//INTEGER Output_Unit
//INTEGER Krylov_Local_Owned       //The number of dof assigned to this processor
//INTEGER Krylov_BackVectors
//INTEGER Krylov_Maximum_Iterations
//INTEGER Krylov_Iterations
//REAL*8  Krylov_Absolute_Tolerance
//REAL*8  Krylov_Relative_Tolerance
//REAL*8  Krylov_Divergence_Tolerance
//REAL*8  Krylov_Basis(Krylov_Local_Owned*Krylov_BackVectors)
//REAL*8  Krylov_Hessenberg(Krylov_BackVectors*Krylov_BackVectors+Krylov_BackVectors)
//REAL*8  Krylov_Givens(Krylov_BackVectors*2)
//REAL*8  Krylov_PC_Basis(Krylov_Local_Owned*Krylov_BackVectors)
//REAL*8  Krylov_Modified_RHS(Krylov_BackVectors)
//REAL*8 Solution(Krylov_Local_Owned)      //In: an initial guess/ Out: the updated approximation of the solution
//REAL*8 RightHandSide(Krylov_Local_Owned) //Right hand side of the equation (= b vector)
//Thread specific variables
//INTEGER  MyThreadID                 //The ID of the thread 1:NumThreads
//INTEGER  MyStart,MyEnd              //The starting/ending points of the thread
//Thread shared variables
//INTEGER  NumThreads                 //The number of threads
//LOGICAL  GuessIsNonZero             //If true the vector stored in Solution is used as the initial guess, otherwise a null vector is used
//INTEGER  ReasonForConvergence       //Divergence (<0), MaxIterationCount (=0), Convergence (>0)
//INTEGER  IterationCount             //Count the total number of inner iteration (number of time we apply A and orthogonalize)
//INTEGER  ParallelComm               //ParallelCommunicator, needed to compute vector norm via all reduce
//INTEGER  ParallelRank
//REAL*8 ResidualNorm,VectorNorm      //Norm of the residual that is returned (a shared variable between the threads)
//REAL*8 VectorNorm_Local(NumThreads) //An array used to store threadwise copies of the vector sum
//REAL*8 HessenNorm_Local(NumThreads*Krylov_BackVectors)
//INTEGER  iMethod
//Local Thread specific variables
int     Maximum_Outer;             //Maximum of outer iteration to stop even if we for (not get convergence 
double  LocalConst,Relative_Stop,Divergence_Stop; //These could be thread shared
double  Cosinus,Sinus,aconst,bconst; //These are only needed on the root thread
double  stupidnum;
int     I,J,K,L,Outer,Inner;            //Thes must be thread specific to avoid unnecessary barriers
int     BackVectorsp1;

//Grab the value from the structure
Maximum_Outer = *Krylov_Maximum_Iterations / *Krylov_BackVectors;
if (Maximum_Outer* *Krylov_BackVectors != *Krylov_Maximum_Iterations) {Maximum_Outer = Maximum_Outer + 1;}

#ifdef Local_Debug_FGMRES_Driver
   if (*MyThreadID == 1) {
      printf("[GMRES]...Initial guess %9d \n",*Krylov_Local_Owned);
      for (I = 1;I<*Krylov_Local_Owned+1;I++) {
         printf("[GMRES]...solution(%4d) = %13.6e RHS=%13.6e \n",I,Solution[I-1],RightHandSide[I-1]);
      }
   }
#endif

BackVectorsp1 = *Krylov_BackVectors+1;

if (*GuessIsNonZero == 1) { //Apply A to the initial guess
   for (I = *MyStart;I<*MyEnd+1;I++) {
      Krylov_Basis[I-1] = 0.0;
   }
   //Start: Compute the initial residual r0 = Ax0 - b and the first basis vector v1 = r0/||r0||
   LookDownBelow(iMethod,Solution,Krylov_Basis); //(1,1)
   for (I = *MyStart;I<*MyEnd+1;I++) {
      Krylov_Basis[I-1] = RightHandSide[I-1] - Krylov_Basis[I-1];
   } }
else { //Zero solution if GuessIsNonZero == False
   for (I = *MyStart;I<*MyEnd+1;I++) {
      Solution[I-1] = 0.0;
      Krylov_Basis[I-1] = RightHandSide[I-1];
   }
}

#ifdef Local_Debug_FGMRES_Driver
   if (*MyThreadID == 1) {
      printf("[GMRES]...Initial residual  \n");
      for (I = 1;I<*Krylov_Local_Owned+1;I++) {
         printf("[GMRES]...Basis(%4d) = %13.6e \n",I,Krylov_Basis[I-1]);
      }
   }
#endif

VectorNorm_Local[*MyThreadID-1] = 0.0;
for (I = *MyStart;I<*MyEnd+1;I++) {
   VectorNorm_Local[*MyThreadID-1] = VectorNorm_Local[*MyThreadID-1] + Krylov_Basis[I-1] * Krylov_Basis[I-1];
}
//Reduction over all threads. Two barriers are needed to ensure data is consistent at both points
#ifdef WITHOMP
#pragma omp barrier
#endif

if (*MyThreadID == 1) {
   *ReasonForConvergence = 0;
   *IterationCount = 0;
   stupidnum = 0.0;
   for (I = 1;I<*NumThreads+1;I++) {
      //printf("[GMRES]...vector norm local %13.6e \n",VectorNorm_Local[I-1]);
      stupidnum = stupidnum + VectorNorm_Local[I-1];
   }
#ifdef USEMPI
   LocalConst = stupidnum;
   CALL MPI_ALLREDUCE(LocalConst,ResidualNorm,1,MPI_DOUBLE_PRECISION,MPI_SUM,ParallelComm,J)
   CALL Basic_CheckError(Output_Unit,J,"MPI_ALLREDUCE for ResidualNorm in (Method_FGMRES)")
#endif
   //printf("[GMRES]...before sqrt residual norm %13.6e \n",stupidnum);

   *ResidualNorm = sqrt(stupidnum);
#ifdef Local_Debug_FGMRES_Driver
   printf("[GMRES]...Initial residual norm %13.6e \n",*ResidualNorm);
#endif
}
#ifdef WITHOMP
#pragma omp barrier
#endif

//Setup the relative tolerance limit
Relative_Stop = *ResidualNorm * *Krylov_Relative_Tolerance;
//Setup the divergence tolerance limit
Divergence_Stop = *ResidualNorm *(1.0 + *Krylov_Divergence_Tolerance);

#ifdef Local_Debug_FGMRES_Driver
   if (*MyThreadID == 1) {
      printf("[GMRES]...Max Outer        = %14d \n",Maximum_Outer);
      printf("[GMRES]...Max Inner        = %14d \n",*Krylov_BackVectors);
      printf("[GMRES]...Total Max Iter   = %14d \n",*Krylov_Maximum_Iterations);
      printf("[GMRES]...Absolute target  = %13.6e \n",*Krylov_Absolute_Tolerance);
      printf("[GMRES]...Relative target  = %13.6e \n",Relative_Stop);
      printf("[GMRES]...Divergence Stop  = %13.6e \n",Divergence_Stop);
   }
#endif

if (*ResidualNorm == 0.0) {goto theend;} //Yes, there is nothing to for (because we have the exact answer

for (Outer = 1;Outer<Maximum_Outer+1;Outer++) {
   //Serial work
   if (*MyThreadID == 1) {
      for (J = 1;J<*Krylov_BackVectors+1;J++) {
         for (I = 1;I<*Krylov_BackVectors+1;I++) {
            Krylov_Hessenberg[(J-1)*BackVectorsp1+I-1] = 0.0;
         }
      }
      for (J = 1;J<*Krylov_BackVectors+1;J++) {
         Krylov_Modified_RHS[J-1] = 0.0;
      }
      Krylov_Modified_RHS[1-1] = *ResidualNorm; //No barrier is needed for this as I assume we need one just before and after Apply_PC and Apply_A
   }

   //For the initialization the Residual is in Basis(*,1) and its norm in ResidualNorm    
   LocalConst = 1.0 / *ResidualNorm;
   for (I = *MyStart;I<*MyEnd+1;I++) {
      Krylov_Basis[I-1] =  Krylov_Basis[I-1] * LocalConst;
      Krylov_Basis[*Krylov_Local_Owned+I-1] = 0.0;
      Krylov_PC_Basis[I-1] = 0.0;
   }

   //Iterate: For j =1,2,...m do:
   //  h(i,j) = (Avj,vi), i = 1,...,j
   //  v*(j+1) = Av(j) - sum(h(i,j)*v(i),i=1,j)
   //  h(j+1,j) = ||v*(j+1)||
   //  v(j+1) = v*(j+1)/ h(j+1,j)
   for (Inner = 1;Inner<*Krylov_BackVectors+1;Inner++) {
      //Apply the preconditionner to the basis vector and store it in a intermediate vector
      SolveWGS_PassThrough_PC(&Krylov_Basis[(Inner-1)* *Krylov_Local_Owned+1-1],
                                           &Krylov_PC_Basis[(Inner-1)* *Krylov_Local_Owned+1-1] );
#ifdef Local_Debug_FGMRES_Driver
      if (*MyThreadID == 1) {
         printf("[GMRES]...At Outer %4d and Inner %4d some vectors: \n",Outer,Inner);
         for (I = 1;I<*Krylov_Local_Owned+1;I++) {
            printf("[GMRES]...residual/norm=V(%4d) = %13.6e PC*V(%4d) = %13.6e \n",
                 I,Krylov_Basis[(Inner-1)* *Krylov_Local_Owned+I-1],I,Krylov_PC_Basis[(Inner-1)* *Krylov_Local_Owned+I-1]);
         }
         for (I = 1;I<*Krylov_Local_Owned+1;I++) {
            printf("[GMRES]...InitV=(%4d) = %13.6e \n",I, Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1]);
         }

      }
#endif
      //Apply A to the intermediate vector and store it in the next basis point before orthonormalization
      LookDownBelow(iMethod,&Krylov_PC_Basis[(Inner-1)* *Krylov_Local_Owned+1-1],
                            &Krylov_Basis[(Inner+1-1)* *Krylov_Local_Owned+1-1]);
#ifdef Local_Debug_FGMRES_Driver
      if (*MyThreadID == 1) {
         printf("[GMRES]...At Outer %4d and Inner %4d intermediate vectors: \n",Outer,Inner);
         for (I = 1;I<*Krylov_Local_Owned+1;I++) {
            printf("[GMRES]...A*Z=V(%4d) = %13.6e %13.6e \n",I, Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1],
                                                                 Krylov_PC_Basis[(Inner-1)* *Krylov_Local_Owned+I-1]);
         }
      }
#endif
      //Update the jth column of the Hessenberg matrix
      //We need to reduce on Hessenberg(:,K) across all threads
      for (I = 1;I<Inner+1;I++) {
         HessenNorm_Local[(I-1)* *NumThreads+*MyThreadID-1] = 0.0;
         for (J = *MyStart;J<*MyEnd+1;J++) {
            HessenNorm_Local[(I-1)* *NumThreads+*MyThreadID-1] = 
            HessenNorm_Local[(I-1)* *NumThreads+*MyThreadID-1] 
            + Krylov_Basis[(Inner)* *Krylov_Local_Owned+J-1] * Krylov_Basis[(I-1)* *Krylov_Local_Owned+J-1];
         }
      }
#ifdef WITHOMP
#pragma omp barrier
#endif
      if (*MyThreadID == 1) {
         K = (BackVectorsp1-1) * BackVectorsp1;
         L = (Inner-1) * BackVectorsp1;
         for (I = 1;I<Inner+1;I++) {
            LocalConst = 0.0;
            for (J = 1;J<*NumThreads+1;J++) {
               LocalConst = LocalConst + HessenNorm_Local[(I-1)* *NumThreads+J-1];
            }
#ifdef USEMPI
            Krylov_Hessenberg[K+I-1] = LocalConst;
#else
            Krylov_Hessenberg[L+I-1] = LocalConst;
#endif
         }
#ifdef USEMPI
         //Reduce the dot product on the whole communicator space
         CALL MPI_ALLREDUCE(Krylov_Hessenberg(K+1),Krylov_Hessenberg(L+1),  &
                            *Krylov_Local_Owned,MPI_DOUBLE_PRECISION,MPI_SUM,ParallelComm,J)
         CALL Basic_CheckError(Output_Unit,J,'MPI_ALLREDUCE Hessenberg_Column in (Basic_FillHessenberg)')
#endif
      }
#ifdef WITHOMP
#pragma omp barrier
#endif

#ifdef Local_Debug_FGMRES_Driver
      if (*MyThreadID == 1) {
         printf("[GMRES]...Hessenberg matrix K=%3d \n",K);
         for (I = 1;I<Inner+1;I++) {
            printf("[GMRES]...Hessenberg(%3d) = ",I);
            for (J=1;J<Inner+1;J++) {
               printf("%13.6e ",Krylov_Hessenberg[(J-1)*BackVectorsp1+I-1]);
            }
            printf("\n");
         }
      }
#endif
      //printf("error is here? 0 \n");
      //Compute an new orthogonal vector
      for (J = 1;J<Inner+1;J++) {
         for (I = *MyStart;I<*MyEnd+1;I++) {
            Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1] = 
            Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1] 
            - Krylov_Hessenberg[(Inner-1)*BackVectorsp1+J-1]
            * Krylov_Basis[(J-1)* *Krylov_Local_Owned+I-1];
         }
      }
      //printf("error is here? 1 \n");

      //Compute its norm
      VectorNorm_Local[*MyThreadID-1] = 0.0;
      for (I = *MyStart;I<*MyEnd+1;I++) {
         VectorNorm_Local[*MyThreadID-1] = 
         VectorNorm_Local[*MyThreadID-1] 
         + Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1]
         * Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1];
      }
      //Reduction over all threads. Two barriers are needed to ensure data is consistent at both points
#ifdef WITHOMP
#pragma omp barrier
#endif
      //printf("error is here? 2 \n");
      if (*MyThreadID == 1) {
         *VectorNorm = 0.0;
         for (I = 1;I<*NumThreads+1;I++) {
            *VectorNorm = *VectorNorm + VectorNorm_Local[I-1];
         }
#ifdef USEMPI
         LocalConst = *VectorNorm;
         CALL MPI_ALLREDUCE(LocalConst,VectorNorm,1,MPI_DOUBLE_PRECISION,MPI_SUM,ParallelComm,J)
         CALL Basic_CheckError(Output_Unit,J,"MPI_ALLREDUCE for ResidualNorm in (Method_FGMRES)")
#endif
         *VectorNorm = sqrt(*VectorNorm);
         Krylov_Hessenberg[(Inner-1)*BackVectorsp1+Inner+1-1] = *VectorNorm;
         if (*VectorNorm == 0.0) {
            *ResidualNorm = *VectorNorm;} //This was moved from below and the check by all below will cause an exit
         else {
            *VectorNorm = 1.0 / *VectorNorm;
            //An efficient way to compute the norm of the residual is to use a QR decomposition of the Hessenberg matrix
            //We also solve the minimization problem
            //Apply previous rotation to the last column of Hessenberg
            for (I = 1;I<Inner-1+1;I++) {
               aconst = Krylov_Hessenberg[(Inner-1)*BackVectorsp1+I-1];
               bconst = Krylov_Hessenberg[(Inner-1)*BackVectorsp1+I+1-1]; //We cannot remove this because of the redefinition in the next two lines
               Krylov_Hessenberg[(Inner-1)*BackVectorsp1+I-1] = 
               Krylov_Givens[I-1]*aconst - Krylov_Givens[BackVectorsp1+I-1]*bconst;
               Krylov_Hessenberg[(Inner-1)*BackVectorsp1+I+1-1] = 
               Krylov_Givens[BackVectorsp1+I-1]*aconst + Krylov_Givens[I-1]*bconst;
            }
            //Compute the givens coefficients for this iteration
            stupidnum = Krylov_Hessenberg[(Inner-1)*BackVectorsp1+Inner-1];
            LocalConst = Krylov_Hessenberg[(Inner-1)*BackVectorsp1+Inner+1-1];
            *ResidualNorm = sqrt(stupidnum* stupidnum + LocalConst*LocalConst);
            LocalConst = 1.0 / *ResidualNorm;
            Krylov_Givens[Inner-1] =  Krylov_Hessenberg[(Inner-1)*BackVectorsp1+Inner-1] * LocalConst;
            Krylov_Givens[BackVectorsp1+Inner-1] = -Krylov_Hessenberg[(Inner-1)*BackVectorsp1+Inner+1-1] * LocalConst;

            //Apply this rotation to the Hessenberg matrix and to the modified right hand side
            aconst  = Krylov_Hessenberg[(Inner-1)*BackVectorsp1+Inner-1];
            bconst  = Krylov_Hessenberg[(Inner-1)*BackVectorsp1+Inner+1-1]; //We cannot remove this because of the redefinition in the next two lines
            Krylov_Hessenberg[(Inner-1)*BackVectorsp1+Inner-1] = 
            Krylov_Givens[Inner-1]*aconst - Krylov_Givens[BackVectorsp1+Inner-1]*bconst;
            Krylov_Hessenberg[(Inner-1)*BackVectorsp1+Inner+1-1] = 
            Krylov_Givens[BackVectorsp1+Inner-1]*aconst + Krylov_Givens[Inner-1]*bconst;
            aconst  = Krylov_Modified_RHS[Inner-1];
            bconst  = Krylov_Modified_RHS[Inner+1-1]; //We cannot remove this because of the redefinition in the next two lines
            Krylov_Modified_RHS[Inner-1]   = Krylov_Givens[Inner-1]*aconst - Krylov_Givens[BackVectorsp1+Inner-1]*bconst;
            Krylov_Modified_RHS[Inner+1-1] = Krylov_Givens[BackVectorsp1+Inner-1]*aconst + Krylov_Givens[Inner-1]*bconst;
            //Now we have the residual norm at no extra cost
            *ResidualNorm = Krylov_Modified_RHS[Inner+1-1];
            if (*ResidualNorm < 0.0) {*ResidualNorm = - *ResidualNorm;}
            //update the number of iteration before exit
            *IterationCount = *IterationCount + 1;
            Krylov_Iterations = Krylov_Iterations + 1;
         } //Vector norm = 0
      }
#ifdef WITHOMP
#pragma omp barrier
#endif
      //printf("error is here? 3 \n");

#ifdef Local_Debug_FGMRES_Driver
      if (*MyThreadID == 1) {
         printf("[GMRES]...VectorNorm= %13.6e \n",*VectorNorm);
         for (I = 1;I<*Krylov_Local_Owned+1;I++) {
            printf("[GMRES]...New XX(%4d) = %13.6e \n",I,Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1]* *VectorNorm);
         }
         for (I = 1;I<Inner+1+1;I++) {
            printf("[GMRES]...Rotated Hessenberg(%3d) = ",I);
            for (J=1;J<Inner+1;J++) {
               printf("%13.6e ",Krylov_Hessenberg[(J-1)*BackVectorsp1+I-1]);
            }
            printf("\n");
         }
         for (I = 1;I<Inner+1;I++) {
            printf("[GMRES]...Givens(%3d) = %13.6e  %13.6e \n",I,Krylov_Givens[I-1],Krylov_Givens[BackVectorsp1+I-1]);
         }
   
         for (I = 1;I<Inner+1;I++) {
            printf("[GMRES]......RHS(%4d) = %13.6e \n",I,Krylov_Modified_RHS[I-1]);
         }
         printf("[GMRES]...At Outer %4d and Inner %4d residual = %13.6e \n",Outer,Inner,*ResidualNorm);
      }
#ifdef WITHOMP
#pragma omp barrier
#endif
#endif
      //if Vector Norm == 0, we reach the happy break down.
      if (*VectorNorm == 0.0) {goto ExitInnerLoop;} //Inner iteration loop

      if (Inner == *Krylov_BackVectors) { //We have nothing left to initialize so...
         for (I = *MyStart;I<*MyEnd+1;I++) {
            Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1] = 
            Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1] * *VectorNorm;
         } }
      else {
         for (I = *MyStart;I<*MyEnd+1;I++) {
            Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1] = 
            Krylov_Basis[(Inner)* *Krylov_Local_Owned+I-1] * *VectorNorm;
            Krylov_Basis[(Inner+1)* *Krylov_Local_Owned+I-1]  = 0.0;
            Krylov_PC_Basis[(Inner)* *Krylov_Local_Owned+I-1] = 0.0;
         }
      }

      //Test if the Residual meets the stopping criterion
      if ( (*ResidualNorm < *Krylov_Absolute_Tolerance) || 
           (*ResidualNorm < Relative_Stop) ) {goto ExitInnerLoop;}

      //Test if we diverged
      if (*ResidualNorm > Divergence_Stop) {goto ExitInnerLoop;}

      //Test if the iteration limit was hit
      if (*IterationCount >= *Krylov_Maximum_Iterations) {goto ExitInnerLoop;}
   } //Inner

   ExitInnerLoop:

   if (Inner > *Krylov_BackVectors) {Inner = *Krylov_BackVectors;}
   //Form the approximate solution x_m = x_0 + V_m*y_m
   //where y_m minimizes ||Beta.e_1 - H_m.y||, y in R^m
   //All we need is the matrix H_m and the vectors (v1,...vm)
   //We only need to compute h(i,m) = (A*vm, vi)
#ifdef Local_Debug_FGMRES_Driver
   if (*MyThreadID == 1) {
      printf("[GMRES]...Exit inner loop with inner = %4d \n",Inner);
   }
#endif
   if (*ResidualNorm > Divergence_Stop) {
      *ReasonForConvergence = -4; //KSP_DIVERGED_DTOL
      goto theend;
   }

   //Compute Ym by solving Ym = R-1 * G
   //R is an upper triangular matrix, the modified Hessenberg matrix
   //OMP BARRIER  //This barrier is not needed as both instances of exiting out of the inner loop will have consistent data
   if (*MyThreadID == 1) {
      for (K = 1;K<Inner+1;K++) {
         I = Inner+1-K;
         for (J = I+1;J<Inner+1;J++) {
            Krylov_Modified_RHS[I-1] = 
            Krylov_Modified_RHS[I-1] 
            - Krylov_Hessenberg[(J-1)*BackVectorsp1+I-1]*Krylov_Modified_RHS[J-1];
         }
         Krylov_Modified_RHS[I-1] = 
         Krylov_Modified_RHS[I-1] / Krylov_Hessenberg[(I-1)*BackVectorsp1+I-1]; //This cannot be moved
      }
#ifdef Local_Debug_FGMRES_Driver
      for (I = 1;I<*Krylov_Local_Owned+1;I++) {
         printf("[GMRES]...Exiting PC Basis(%3d) = ",I);
         for (J=1;J<Inner+1;J++) {
            printf("%13.6e ",Krylov_PC_Basis[(J-1)* *Krylov_Local_Owned+I-1]);
         }
         printf("\n");
      }
      for (I = 1;I<*Krylov_Local_Owned+1;I++) {
         printf("[GMRES]...Exiting Basis(%3d) = ",I);
         for (J=1;J<Inner+1;J++) {
            printf("%13.6e ",Krylov_Basis[(J-1)* *Krylov_Local_Owned+I-1]);
         }
         printf("\n");
      }
      for (I = 1;I<Inner+1+1;I++) {
         printf("[GMRES]...Exiting RHS(%4d) = %13.6e \n",I,Krylov_Modified_RHS[I-1]);
      }
#endif
   }
#ifdef WITHOMP
#pragma omp barrier
#endif
   //Update the solution Xm = X0 + Zm*Y
   for (J = 1;J<Inner+1;J++) {
      for (I = *MyStart;I<*MyEnd+1;I++) {
         Solution[I-1] = Solution[I-1] 
         + Krylov_PC_Basis[(J-1)* *Krylov_Local_Owned+I-1]*Krylov_Modified_RHS[J-1];
      }
   }

#ifdef Local_Debug_FGMRES_Driver
   if (*MyThreadID == 1) {
      printf("[GMRES]...Solution after outer iteration %4d \n",Outer);
      for (I = 1;I<*Krylov_Local_Owned+1;I++) {
         printf("[GMRES]...X(%4d) = %13.6e \n",I,Solution[I-1]);
      }
   }
#endif

   //Initialize the basis so the user does not have to for (it
   for (I = *MyStart;I<*MyEnd+1;I++) {
      Krylov_Basis[I-1] = 0.0;
   }

   //Compute r_m = f- A x_m and store it in Scratch_Basis 
#ifdef WITHOMP
#pragma omp barrier
#endif
   LookDownBelow(iMethod,Solution,Krylov_Basis);
#ifdef WITHOMP
#pragma omp barrier
#endif
   
   VectorNorm_Local[*MyThreadID-1] = 0.0;
   for (I = *MyStart;I<*MyEnd+1;I++) {
      Krylov_Basis[I-1] = RightHandSide[I-1] - Krylov_Basis[I-1];
      VectorNorm_Local[*MyThreadID-1] = VectorNorm_Local[*MyThreadID-1] + Krylov_Basis[I-1]*Krylov_Basis[I-1];
   }
   //Reduction over all threads. Two barriers are needed to ensure data is consistent at both points
#ifdef WITHOMP
#pragma omp barrier
#endif
   if (*MyThreadID == 1) {
      stupidnum = 0.0;
      for (I = 1;I<*NumThreads+1;I++) {
         stupidnum = stupidnum + VectorNorm_Local[I-1];
      }
#ifdef USEMPI
      LocalConst = stupidnum
      CALL MPI_ALLREDUCE(LocalConst,VectorNorm,1,MPI_DOUBLE_PRECISION,MPI_SUM,ParallelComm,J)
      CALL Basic_CheckError(Output_Unit,J,"MPI_ALLREDUCE for ResidualNorm in (Method_FGMRES)")
#endif
      *ResidualNorm = sqrt(stupidnum);
   }
#ifdef WITHOMP
#pragma omp barrier
#endif

#ifdef Local_Debug_FGMRES_Driver
   if (*MyThreadID == 1) {printf("[GMRES]...At FinalCheck residual = %13.6e \n",*ResidualNorm);}
#endif
       //Test if the ResidualNorm meet the stopping criterion
    if (*ResidualNorm <= 0.0) { //Happy breakdown convergence
       *ReasonForConvergence = 5; //KSP_CONVERGED_BREAKDOWN
       goto theend; 
    }
    if (*ResidualNorm < *Krylov_Absolute_Tolerance) {
       *ReasonForConvergence = 3; //KSP_CONVERGED_ATOL 
       goto theend;
    }
    if (*ResidualNorm < Relative_Stop) {
       *ReasonForConvergence = 2; //KSP_CONVERGED_RTOL
       goto theend;
    }
    if (*IterationCount >= *Krylov_Maximum_Iterations) {
       *ReasonForConvergence = 1; //KSP_HIT_ITERATION_LIMIT
       goto theend;
    }
    if (*MyThreadID == 1) {
         printf("[SN-KERNEL]...Outer %6d after %6d inners has residual = %13.6e \n",Outer,Inner,*ResidualNorm);
    }
} //Outer

theend:
I=0;

} // END SUBROUTINE Method_FGMRES_Threaded

void LookDownBelow(int *iMethod,double *RHS_C,double *LHS_C) {
//INTEGER iMethod
//REAL*8  RHS_C,LHS_C
//printf("iMethod = %4d",*iMethod);
if (*iMethod == 1) {SolveWGS_PassThrough_AVE1(RHS_C,LHS_C);}
if (*iMethod == 2) {SolveWGS_PassThrough_AVE2(RHS_C,LHS_C);}
if (*iMethod == 3) {SolveWGS_PassThrough_AVE3(RHS_C,LHS_C);}

} // END SUBROUTINE LookDownBelow
