#include <stdlib.h>
#include "ApplyA_functions.h"
#ifdef WITHOMP
#include <omp.h>
#endif
// This subroutine performs the combined F,UT,U matrix application for a bunch of elements and angles for a set of tets
// The matrix operation can be viewed as: LHS(A,E) = C(A,E)*RHS(A,E)
// ------------
// FLOP SUMMARY
// There is significant pointer math which I ignore
//  832 additions/subtractions /E/A
// 4096 multiplications /E/A
// 4928 flops per element per angle
// ----
// Fortran interface routines
void applya_ave2_tet_supg(int *NE, int *NA, int *NV, int *C,
                          int *AS_NumColors,int *AS_NumThreads,int *TasksPerThread,int *MyThreadID,int *ThreadWiseWork,
                              double *CF, double *CU, double *CP, double *FES, double *FED, double *FEW, double *OM, double *OO,
                              double *LHS, double *RHS,double *SRHS,double *SLHS) {
ApplyA_AVE2_Tet_SUPG(NE,NA,NV,C,
                     AS_NumColors,AS_NumThreads,TasksPerThread,MyThreadID,ThreadWiseWork,
                              CF, CU, CP, FES, FED, FEW, OM, OO,
                              LHS, RHS,SRHS,SLHS);
}
void applya_ave2_tet_supg_(int *NE, int *NA, int *NV, int *C,
                           int *AS_NumColors,int *AS_NumThreads,int *TasksPerThread,int *MyThreadID,int *ThreadWiseWork,
                              double *CF, double *CU, double *CP, double *FES, double *FED, double *FEW, double *OM, double *OO,
                              double *LHS, double *RHS,double *SRHS,double *SLHS) {
ApplyA_AVE2_Tet_SUPG(NE,NA,NV,C,
                     AS_NumColors,AS_NumThreads,TasksPerThread,MyThreadID,ThreadWiseWork,
                              CF, CU, CP, FES, FED, FEW, OM, OO,
                              LHS, RHS,SRHS,SLHS);
}

// Main subroutine header
void ApplyA_AVE2_Tet_SUPG(int *NE, int *NA, int *NV, int *C,
                           int *AS_NumColors,int *AS_NumThreads,int *TasksPerThread,int *MyThreadID,int *ThreadWiseWork,
                              double *CF, double *CU, double *CP, double *FES, double *FED, double *FEW, double *OM, double *OO,
                              double *LHS, double *RHS,double *SRHS,double *SLHS) {
//IMPLICIT NONE 
//#include "PROTEUS_Preprocess.h"
//#define Local_TryLHStoo
// Problem size arguments and mesh data
//PROTEUS_Int,     INTENT(IN) :: NE            // Number of Es
//PROTEUS_Int,     INTENT(IN) :: NA              // Number of angles
//PROTEUS_Int,     INTENT(IN) :: NV                // Number of vertices in the mesh
//PROTEUS_Int,     INTENT(IN) :: C(4,NE)       // The connectivity matrix
//// Thread size arguments and work seperation
//PROTEUS_Int,     INTENT(IN) :: AS_NumColors,AS_NumThreads,TasksPerThread,MyThreadID
//PROTEUS_Int,     INTENT(IN) :: ThreadWiseWork(2,AS_NumColors,AS_NumThreads) // The element start/stop for each thread in each color
//// Cross section vectors
//PROTEUS_Real,    INTENT(IN) :: CF(NE)
//PROTEUS_Real,    INTENT(IN) :: CU(NE)
//PROTEUS_Real,    INTENT(IN) :: CP(NE)
//// FE shape functions
//PROTEUS_FE_Real, INTENT(IN) :: FES(4,4)               // Vertices,GP     // Shape
//PROTEUS_FE_Real, INTENT(IN) :: FED(4,3,4,NE) // Vertices,Dim,GP // Derivatives in real space
//PROTEUS_FE_Real, INTENT(IN) :: FEW(4,NE)     // GP              // Weight * |Jac|
//// Angular cubature
//PROTEUS_FE_Real, INTENT(IN) :: OM(NA,3)                   // Angular cubature evaluation
//PROTEUS_FE_Real, INTENT(IN) :: OO(NA,6)                   // Basically it is OM(A,*) * OM(A,*)
//// Vectors
//PROTEUS_Real, INTENT(INOUT) :: LHS(NA,NV) // LHS (to be modified)
//PROTEUS_Real,    INTENT(IN) :: RHS(NA,NV) // RHS (to be multiplied)
//// Scratch vectors
//PROTEUS_Real, INTENT(INOUT) :: SRHS(NA,4),SLHS(NA,4)

// Local
int E,A,V,VV,GP,iNA;
int iCV_E,iCVV_E,iGP,iEGP,iEGP1,iEGP2,iEGP3;
int iColor,iStart,iEnd;
int AS_Thread,iTask,ifun;
double LocalConst1,LocalConst2,LocalConst3,Products[11];

iNA = *NA;

//printf("NA=%d NE=%d NV=%d \n",*NA,*NE,*NV);

for (iColor = 1; iColor < *AS_NumColors+1; iColor++) {
#ifdef WITHOMP
 for (iTask = 1; iTask < *TasksPerThread+1; iTask++) {
   AS_Thread = *TasksPerThread * (*MyThreadID-1) + iTask;
   ifun = *AS_NumColors * (AS_Thread-1)*2 + 2*(iColor-1);
   iStart   = ThreadWiseWork[ifun+1-1];
   iEnd     = ThreadWiseWork[ifun+2-1];
//   iStart   = ThreadWiseWork[(iColor-1)*2+1-1];
//   iEnd     = ThreadWiseWork[(iColor-1)*2+2-1];
#else
   iStart = 1;
   iEnd   = *NE;
#endif
   for (E = iStart; E < iEnd+1; E++) {
      for (V = 1; V < 5; V++) {
         // Pull a copy of the vectors
         iCV_E = (C[(E-1)*4+V-1]-1)*iNA-1;
         for (A = 1; A < iNA+1; A++) {
            SRHS[(V-1)*iNA+A-1] = RHS[iCV_E + A];
#ifdef Local_TryLHStoo
            SLHS[(V-1)*iNA+A-1] = LHS[iCV_E + A];
#endif
         }
      }
      for (V = 1; V < 5; V++) {
#ifndef Local_TryLHStoo
         iCV_E = (C[(E-1)*4+V-1]-1)*iNA-1; // II = C(V,E)
#endif
         for (GP = 1; GP < 5; GP++) {
            iGP  = (GP-1)*4-1;
            iEGP  = (E-1)*4+GP-1;
            iEGP1 = (E-1)*4*3*4 + (GP-1)*3*4 + (1-1)*4-1;
            iEGP2 = (E-1)*4*3*4 + (GP-1)*3*4 + (2-1)*4-1;
            iEGP3 = (E-1)*4*3*4 + (GP-1)*3*4 + (3-1)*4-1;

            LocalConst1 = FEW[iEGP]*CF[E-1];
            LocalConst2 = FEW[iEGP]*CU[E-1];
            LocalConst3 = FEW[iEGP]*CP[E-1];
            for (VV = 1; VV < 5; VV++) {
               iCVV_E = (VV-1)*iNA-1; //(C[(E-1)*4+VV-1]-1)*iNA-1;
               Products[ 1] = LocalConst1*FES[iGP+V]  *FES[iGP+VV];
               Products[ 2] = LocalConst2*FED[iEGP1+V]*FES[iGP+VV];
               Products[ 3] = LocalConst2*FED[iEGP2+V]*FES[iGP+VV];
               Products[ 4] = LocalConst2*FED[iEGP3+V]*FES[iGP+VV];
               Products[ 5] = LocalConst3*FED[iEGP1+V]*FED[iEGP1+VV];
               Products[ 6] = LocalConst3*FED[iEGP2+V]*FED[iEGP2+VV];
               Products[ 7] = LocalConst3*FED[iEGP3+V]*FED[iEGP3+VV];
               Products[ 8] = LocalConst3*(FED[iEGP1+V]*FED[iEGP2+VV]+FED[iEGP2+V]*FED[iEGP1+VV]);
               Products[ 9] = LocalConst3*(FED[iEGP1+V]*FED[iEGP3+VV]+FED[iEGP3+V]*FED[iEGP1+VV]);
               Products[10] = LocalConst3*(FED[iEGP2+V]*FED[iEGP3+VV]+FED[iEGP3+V]*FED[iEGP2+VV]);
               for (A = 1; A < iNA+1; A++) {
#ifdef Local_TryLHStoo
                  SLHS[iCVV_E+A] = SLHS[iCVV_E+A] 
                                        + Products[ 1]              *SRHS[iCVV_E+A]    // F
                                        + Products[ 2]*OM[A-1]      *SRHS[iCVV_E+A]    // U^T(1)
                                        + Products[ 3]*OM[iNA+A-1]  *SRHS[iCVV_E+A]    // U^T(2)    
                                        + Products[ 4]*OM[2*iNA+A-1]*SRHS[iCVV_E+A]    // U^T(3)
                                        + Products[ 5]*OO[A-1]      *SRHS[iCVV_E+A]    // P(1,1)
                                        + Products[ 6]*OO[iNA+A-1]  *SRHS[iCVV_E+A]    // P(2,2)
                                        + Products[ 7]*OO[2*iNA+A-1]*SRHS[iCVV_E+A]    // P(3,3)
                                        + Products[ 8]*OO[3*iNA+A-1]*SRHS[iCVV_E+A]    // P(1,2)
                                        + Products[ 9]*OO[4*iNA+A-1]*SRHS[iCVV_E+A]    // P(1,3)
                                        + Products[10]*OO[5*iNA+A-1]*SRHS[iCVV_E+A];   // P(2,3)
#else
                  LHS[iCV_E+A] = LHS[iCV_E+A] 
                                        + Products[ 1]              *SRHS[iCVV_E+A]    // F
                                        + Products[ 2]*OM[A-1]      *SRHS[iCVV_E+A]    // U^T(1)
                                        + Products[ 3]*OM[iNA+A-1]  *SRHS[iCVV_E+A]    // U^T(2)    
                                        + Products[ 4]*OM[2*iNA+A-1]*SRHS[iCVV_E+A]    // U^T(3)
                                        + Products[ 5]*OO[A-1]      *SRHS[iCVV_E+A]    // P(1,1)
                                        + Products[ 6]*OO[iNA+A-1]  *SRHS[iCVV_E+A]    // P(2,2)
                                        + Products[ 7]*OO[2*iNA+A-1]*SRHS[iCVV_E+A]    // P(3,3)
                                        + Products[ 8]*OO[3*iNA+A-1]*SRHS[iCVV_E+A]    // P(1,2)
                                        + Products[ 9]*OO[4*iNA+A-1]*SRHS[iCVV_E+A]    // P(1,3)
                                        + Products[10]*OO[5*iNA+A-1]*SRHS[iCVV_E+A];   // P(2,3)
#endif
                  //SLHS(A,V) = SLHS(A,V) + Products[ 1)        *SRHS[iCVV_E+A]   // F
                  //SLHS(A,V) = SLHS(A,V) + Products[ 2)*OM(A,1)*SRHS[iCVV_E+A]   // U^T(1)
                  //SLHS(A,V) = SLHS(A,V) + Products[ 3)*OM(A,2)*SRHS[iCVV_E+A]   // U^T(2)    
                  //SLHS(A,V) = SLHS(A,V) + Products[ 4)*OM(A,3)*SRHS[iCVV_E+A]   // U^T(3)
                  //SLHS(A,V) = SLHS(A,V) + Products[ 5)*OO(A,1)*SRHS[iCVV_E+A]   // P(1,1)
                  //SLHS(A,V) = SLHS(A,V) + Products[ 6)*OO(A,2)*SRHS[iCVV_E+A]   // P(2,2)
                  //SLHS(A,V) = SLHS(A,V) + Products[ 7)*OO(A,3)*SRHS[iCVV_E+A]   // P(3,3)
                  //SLHS(A,V) = SLHS(A,V) + Products[ 8)*OO(A,4)*SRHS[iCVV_E+A]   // P(1,2)
                  //SLHS(A,V) = SLHS(A,V) + Products[ 9)*OO(A,5)*SRHS[iCVV_E+A]   // P(1,3)
                  //SLHS(A,V) = SLHS(A,V) + Products[10)*OO(A,6)*SRHS[iCVV_E+A]   // P(2,3)
               } // A
            } // VV
         } // GP
      } // V
      // Store the LHS result
#ifdef Local_TryLHStoo
      for (V = 1; V < 5; V++) {
         iCV_E = C[(E-1)*4+V-1]-1; // VV = C(V,E)
               for (A = 1; A < iNA+1; A++) {
            LHS[iCV_E*iNA + A-1] = SLHS[(V-1)*iNA+A-1];
         }
      }
#endif
   } // E
#ifdef WITHOMP
 } // iTask
#pragma omp barrier
#else
   break;
#endif
} // iColor

} // END SUBROUTINE ApplyA_AVE2_Tet_SUPG
