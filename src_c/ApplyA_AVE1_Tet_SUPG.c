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
void applya_ave1_tet_supg(int *NE, int *NA, int *NV, int *C,
                          int *AS_NumColors,int *AS_NumThreads,int *TasksPerThread,int *MyThreadID,int *ThreadWiseWork,
                              double *CF, double *CU, double *CP, double *FES, double *FED, double *FEW, double *OM, double *OO,
                              double *LHS, double *RHS) {
ApplyA_AVE1_Tet_SUPG(NE,NA,NV,C,
                     AS_NumColors,AS_NumThreads,TasksPerThread,MyThreadID,ThreadWiseWork,
                              CF, CU, CP, FES, FED, FEW, OM, OO,
                              LHS, RHS);
}
void applya_ave1_tet_supg_(int *NE, int *NA, int *NV, int *C,
                           int *AS_NumColors,int *AS_NumThreads,int *TasksPerThread,int *MyThreadID,int *ThreadWiseWork,
                              double *CF, double *CU, double *CP, double *FES, double *FED, double *FEW, double *OM, double *OO,
                              double *LHS, double *RHS) {
ApplyA_AVE1_Tet_SUPG(NE,NA,NV,C,
                     AS_NumColors,AS_NumThreads,TasksPerThread,MyThreadID,ThreadWiseWork,
                              CF, CU, CP, FES, FED, FEW, OM, OO,
                              LHS, RHS);
}

// Main subroutine header
void ApplyA_AVE1_Tet_SUPG(int *NE, int *NA, int *NV, int *C,
                           int *AS_NumColors,int *AS_NumThreads,int *TasksPerThread,int *MyThreadID,int *ThreadWiseWork,
                              double *CF, double *CU, double *CP, double *FES, double *FED, double *FEW, double *OM, double *OO,
                              double *LHS, double *RHS) {
// Problem size arguments and mesh data
//INTEGER,     INTENT(IN) :: NE           // Number of Es
//INTEGER,     INTENT(IN) :: NA                    // Number of angles
//INTEGER,     INTENT(IN) :: NV           // Number of vertices in the mesh
//INTEGER,     INTENT(IN) :: C(4*NE)      // The connectivity matrix
// Thread size arguments and work seperation
//INTEGER,     INTENT(IN) :: AS_NumColors,AS_NumThreads,TasksPerThread,MyThreadID
//INTEGER,     INTENT(IN) :: ThreadWiseWork(2,AS_NumColors,AS_NumThreads) ! The element start/stop for each thread in each color
// Cross section vectors
//REAL*8,    INTENT(IN) :: CF(NE)
//REAL*8,    INTENT(IN) :: CU(NE)
//REAL*8,    INTENT(IN) :: CP(NE)
// FE shape functions
//REAL*8, INTENT(IN) :: FES(4*4)               // Vertices,GP     // Shape
//REAL*8, INTENT(IN) :: FED(4*3*4*NE) // Vertices,Dim,GP // Derivatives in real space
//REAL*8, INTENT(IN) :: FEW(4*NE)     // GP              // Weight * |Jac|
// Angular cubature
//REAL*8, INTENT(IN) :: OM(NA*3)                   // Angular cubature evaluation
//REAL*8, INTENT(IN) :: OO(NA*6)                   // Basically it is OM(A,*) * OM(A,*)
// Vectors
//REAL*8, INTENT(INOUT) :: LHS(NA*NV) // LHS (to be modified)
//REAL*8,    INTENT(IN) :: RHS(NA*NV) // RHS (to be multiplied)

// Local
int E,A,V,VV,GP,iNA;
int iCV_E,iCVV_E,iGP,iEGP,iEGP1,iEGP2,iEGP3;
int iColor,iStart,iEnd;
int AS_Thread,iTask,ifun;

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
         //printf("E=%d V=%d ii=%d \n",E,V,(E-1)*4+V-1);
         iCV_E = (C[(E-1)*4+V-1]-1)*iNA-1;
         for (GP = 1; GP < 5; GP++) {
            iGP  = (GP-1)*4-1;
            iEGP  = (E-1)*4+GP-1;
            iEGP1 = (E-1)*4*3*4 + (GP-1)*3*4 + (1-1)*4-1;
            iEGP2 = (E-1)*4*3*4 + (GP-1)*3*4 + (2-1)*4-1;
            iEGP3 = (E-1)*4*3*4 + (GP-1)*3*4 + (3-1)*4-1;
            for (VV = 1; VV < 5; VV++) {
               iCVV_E = (C[(E-1)*4+VV-1]-1)*iNA-1;
               for (A = 1; A < iNA+1; A++) {
//               LHS[iCV_E+A] = LHS[iCV_E+A] &
//+ FEW[iEGP]*RHS[iCVV_E+A]*(CF[E-1]*FES(iGP+V)*FES(iGP+VV) &  // F
//                        +CU[E-1]*(OM[A]*FED(iEGP1+V)*FES(iGP+VV)+OM(iNA+A)*FED(iEGP2+V)*FES(iGP+VV)+OM(2*iNA+A)*FED(iEGP3+V)*FES(iGP+VV))
//+CP[E-1]*(OO(A)*FED(iEGP1+V)*FED(iEGP1+VV)+OO(A)*FED(iEGP2+V)*FED(iEGP2+VV)+OO(2*iNA+A)*FED(iEGP3+V)*FED(iEGP3+VV)
//       +OO(3*iNA+A)*FED(iEGP1+V)*FED(iEGP2+VV)+OO(3*iNA+A)*FED(iEGP2+V)*FED(iEGP1+VV)+OO(4*iNA+A)*FED(iEGP1+V)*FED(iEGP3+VV) &
//       +OO(4*iNA+A)*FED(iEGP3+V)*FED(iEGP1+VV)+OO(5*iNA+A)*FED(iEGP2+V)*FED(iEGP3+VV)+OO(5*iNA+A)*FED(iEGP3+V)*FED(iEGP2+VV))  )
               // There are 4+15+45*GP*V*V=64*4*4*4=4096 mults per element-angle and 13*GP*V*V=832 adds so 4928 flops
                  LHS[iCV_E+A] = LHS[iCV_E+A]  
                                + FEW[iEGP]*CF[E-1]              *FES[iGP+V]  *FES[iGP+VV]  *RHS[iCVV_E+A]    // F
                                + FEW[iEGP]*CU[E-1]*OM[A-1]      *FED[iEGP1+V]*FES[iGP+VV]  *RHS[iCVV_E+A]    // U^T[1]
                                + FEW[iEGP]*CU[E-1]*OM[iNA+A-1]  *FED[iEGP2+V]*FES[iGP+VV]  *RHS[iCVV_E+A]    // U^T[2]    
                                + FEW[iEGP]*CU[E-1]*OM[2*iNA+A-1]*FED[iEGP3+V]*FES[iGP+VV]  *RHS[iCVV_E+A]    // U^T[3]
                                + FEW[iEGP]*CP[E-1]*OO[A-1]      *FED[iEGP1+V]*FED[iEGP1+VV]*RHS[iCVV_E+A]    // P[1,1]
                                + FEW[iEGP]*CP[E-1]*OO[iNA+A-1]  *FED[iEGP2+V]*FED[iEGP2+VV]*RHS[iCVV_E+A]    // P[2,2]
                                + FEW[iEGP]*CP[E-1]*OO[2*iNA+A-1]*FED[iEGP3+V]*FED[iEGP3+VV]*RHS[iCVV_E+A]    // P[3,3]
                                + FEW[iEGP]*CP[E-1]*OO[3*iNA+A-1]*FED[iEGP1+V]*FED[iEGP2+VV]*RHS[iCVV_E+A]    // P[1,2]
                                + FEW[iEGP]*CP[E-1]*OO[4*iNA+A-1]*FED[iEGP1+V]*FED[iEGP3+VV]*RHS[iCVV_E+A]    // P[1,3]
                                + FEW[iEGP]*CP[E-1]*OO[5*iNA+A-1]*FED[iEGP2+V]*FED[iEGP3+VV]*RHS[iCVV_E+A]    // P[2,3]
                                + FEW[iEGP]*CP[E-1]*OO[3*iNA+A-1]*FED[iEGP2+V]*FED[iEGP1+VV]*RHS[iCVV_E+A]    // P[2,1]
                                + FEW[iEGP]*CP[E-1]*OO[4*iNA+A-1]*FED[iEGP3+V]*FED[iEGP1+VV]*RHS[iCVV_E+A]    // P[3,1]
                                + FEW[iEGP]*CP[E-1]*OO[5*iNA+A-1]*FED[iEGP3+V]*FED[iEGP2+VV]*RHS[iCVV_E+A] ;  // P[3,2]
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CF[E-1]             *FES[iGP+V]  *FES[iGP+VV]  *RHS[iCVV_E+A]   // F         //  4  1
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CU[E-1]*OM[A-1]     *FED[iEGP1+V]*FES[iGP+VV]  *RHS[iCVV_E+A]   // U^T[1]    // [5  1]*3
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CU[E-1]*OM[iNA+A-1]  *FED[iEGP2+V]*FES[iGP+VV]  *RHS[iCVV_E+A]   // U^T[2]    
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CU[E-1]*OM[2*iNA+A-1]*FED[iEGP3+V]*FES[iGP+VV]  *RHS[iCVV_E+A]   // U^T[3]
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CP[E-1]*OO[A-1]     *FED[iEGP1+V]*FED[iEGP1+VV]*RHS[iCVV_E+A]   // P[1,1]    // [5  1]*9
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CP[E-1]*OO[iNA+A-1]  *FED[iEGP2+V]*FED[iEGP2+VV]*RHS[iCVV_E+A]   // P[2,2]
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CP[E-1]*OO[2*iNA+A-1]*FED[iEGP3+V]*FED[iEGP3+VV]*RHS[iCVV_E+A]   // P[3,3]
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CP[E-1]*OO[3*iNA+A-1]*FED[iEGP1+V]*FED[iEGP2+VV]*RHS[iCVV_E+A]   // P[1,2]
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CP[E-1]*OO[4*iNA+A-1]*FED[iEGP1+V]*FED[iEGP3+VV]*RHS[iCVV_E+A]   // P[1,3]
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CP[E-1]*OO[5*iNA+A-1]*FED[iEGP2+V]*FED[iEGP3+VV]*RHS[iCVV_E+A]   // P[2,3]
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CP[E-1]*OO[3*iNA+A-1]*FED[iEGP2+V]*FED[iEGP1+VV]*RHS[iCVV_E+A]   // P[2,1]
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CP[E-1]*OO[4*iNA+A-1]*FED[iEGP3+V]*FED[iEGP1+VV]*RHS[iCVV_E+A]   // P[3,1]
                  //LHS[iCV_E+A] = LHS[iCV_E+A] + FEW[iEGP]*CP[E-1]*OO[5*iNA+A-1]*FED[iEGP3+V]*FED[iEGP2+VV]*RHS[iCVV_E+A]   // P[3,2]
               } // A
            } // VV
         } // GP
      } // V
   } // E
#ifdef WITHOMP
 } // iTask
#pragma omp barrier
#else
   break;
#endif
} // iColor

} // END SUBROUTINE ApplyA_AVE1_Tet_SUPG
