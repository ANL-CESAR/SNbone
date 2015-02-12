! This subroutine performs the combined F,UT,U matrix application for a bunch of elements and angles for a set of tets
! The matrix operation can be viewed as: LHS(A,E) = C(A,E)*RHS(A,E)
! ------------
! FLOP SUMMARY
! There is significant pointer math which I ignore
!  832 additions/subtractions /E/A
! 4096 multiplications /E/A
! 4928 flops per element per angle
! ----
SUBROUTINE ApplyA_AVE1_Tet_SUPG(NumElements,NumAngles,NumVertices, C,                                 &
                                AS_NumColors,AS_NumThreads,TasksPerThread,MyThreadID,ThreadWiseWork,  &
                                    CF, CU, CP, FES,FED,FEW, OM,OO,                                   &
                                    LHS,RHS)
IMPLICIT NONE 
#include "PROTEUS_Preprocess.h"
! Problem size arguments and mesh data
PROTEUS_Int,     INTENT(IN) :: NumElements            ! Number of Es
PROTEUS_Int,     INTENT(IN) :: NumAngles              ! Number of angles
PROTEUS_Int,     INTENT(IN) :: NumVertices                ! Number of vertices in the mesh
PROTEUS_Int,     INTENT(IN) :: C(4,NumElements)       ! The connectivity matrix
! Thread size arguments and work seperation
PROTEUS_Int,     INTENT(IN) :: AS_NumColors,AS_NumThreads,TasksPerThread,MyThreadID
PROTEUS_Int,     INTENT(IN) :: ThreadWiseWork(2,AS_NumColors,AS_NumThreads) ! The element start/stop for each thread in each color
! Cross section vectors
PROTEUS_Real,    INTENT(IN) :: CF(NumElements)
PROTEUS_Real,    INTENT(IN) :: CU(NumElements)
PROTEUS_Real,    INTENT(IN) :: CP(NumElements)
! FE shape functions
PROTEUS_FE_Real, INTENT(IN) :: FES(4,4)               ! Vertices,GP     ! Shape
PROTEUS_FE_Real, INTENT(IN) :: FED(4,3,4,NumElements) ! Vertices,Dim,GP ! Derivatives in real space
PROTEUS_FE_Real, INTENT(IN) :: FEW(4,NumElements)     ! GP              ! Weight * |Jac|
! Angular cubature
PROTEUS_FE_Real, INTENT(IN) :: OM(NumAngles,3)                   ! Angular cubature evaluation
PROTEUS_FE_Real, INTENT(IN) :: OO(NumAngles,6)                   ! Basically it is OM(A,*) * OM(A,*)
! Vectors
PROTEUS_Real, INTENT(INOUT) :: LHS(NumAngles,NumVertices) ! LHS (to be modified)
PROTEUS_Real,    INTENT(IN) :: RHS(NumAngles,NumVertices) ! RHS (to be multiplied)

! Local
PROTEUS_Int :: E,A,V,VV,GP,iColor,iStart,iEnd
PROTEUS_Int :: AS_Thread,iTask

DO iColor = 1,AS_NumColors
#ifdef WITHOMP
 DO iTask = 1,TasksPerThread
   AS_Thread = (MyThreadID-1)*TasksPerThread + iTask
   iStart   = ThreadWiseWork(1,iColor,AS_Thread)
   iEnd     = ThreadWiseWork(2,iColor,AS_Thread)
#else
   iStart = 1
   iEnd   = NumElements
#endif
   DO E = iStart,iEnd
      DO V = 1,4
         DO GP = 1,4
            DO VV = 1,4
               DO A = 1,NumAngles
!               LHS(A,C(V,E)) = LHS(A,C(V,E)) &
!+ FEW(GP,E)*RHS(A,C(VV,E))*(CF(E)*FES(V,GP)*FES(VV,GP) &  ! F
!                        +CU(E)*(OM(A,1)*FED(V,1,GP,E)*FES(VV,GP)+OM(A,2)*FED(V,2,GP,E)*FES(VV,GP)+OM(A,3)*FED(V,3,GP,E)*FES(VV,GP))
!+CP(E)*(OO(A,1)*FED(V,1,GP,E)*FED(VV,1,GP,E)+OO(A,1)*FED(V,2,GP,E)*FED(VV,2,GP,E)+OO(A,3)*FED(V,3,GP,E)*FED(VV,3,GP,E)
!       +OO(A,4)*FED(V,1,GP,E)*FED(VV,2,GP,E)+OO(A,4)*FED(V,2,GP,E)*FED(VV,1,GP,E)+OO(A,5)*FED(V,1,GP,E)*FED(VV,3,GP,E) &
!       +OO(A,5)*FED(V,3,GP,E)*FED(VV,1,GP,E)+OO(A,6)*FED(V,2,GP,E)*FED(VV,3,GP,E)+OO(A,6)*FED(V,3,GP,E)*FED(VV,2,GP,E))  )
               ! There are 4+15+45*GP*V*V=64*4*4*4=4096 mults per element-angle and 13*GP*V*V=832 adds so 4928 flops
                  LHS(A,C(V,E)) = LHS(A,C(V,E)) &
                                + FEW(GP,E)*CF(E)        *FES(V,GP)    *FES(VV,GP)    *RHS(A,C(VV,E)) &  ! F
                                + FEW(GP,E)*CU(E)*OM(A,1)*FED(V,1,GP,E)*FES(VV,GP)    *RHS(A,C(VV,E)) &  ! U^T(1)
                                + FEW(GP,E)*CU(E)*OM(A,2)*FED(V,2,GP,E)*FES(VV,GP)    *RHS(A,C(VV,E)) &  ! U^T(2)    
                                + FEW(GP,E)*CU(E)*OM(A,3)*FED(V,3,GP,E)*FES(VV,GP)    *RHS(A,C(VV,E)) &  ! U^T(3)
                                + FEW(GP,E)*CP(E)*OO(A,1)*FED(V,1,GP,E)*FED(VV,1,GP,E)*RHS(A,C(VV,E)) &  ! P(1,1)
                                + FEW(GP,E)*CP(E)*OO(A,2)*FED(V,2,GP,E)*FED(VV,2,GP,E)*RHS(A,C(VV,E)) &  ! P(2,2)
                                + FEW(GP,E)*CP(E)*OO(A,3)*FED(V,3,GP,E)*FED(VV,3,GP,E)*RHS(A,C(VV,E)) &  ! P(3,3)
                                + FEW(GP,E)*CP(E)*OO(A,4)*FED(V,1,GP,E)*FED(VV,2,GP,E)*RHS(A,C(VV,E)) &  ! P(1,2)
                                + FEW(GP,E)*CP(E)*OO(A,5)*FED(V,1,GP,E)*FED(VV,3,GP,E)*RHS(A,C(VV,E)) &  ! P(1,3)
                                + FEW(GP,E)*CP(E)*OO(A,6)*FED(V,2,GP,E)*FED(VV,3,GP,E)*RHS(A,C(VV,E)) &  ! P(2,3)
                                + FEW(GP,E)*CP(E)*OO(A,4)*FED(V,2,GP,E)*FED(VV,1,GP,E)*RHS(A,C(VV,E)) &  ! P(2,1)
                                + FEW(GP,E)*CP(E)*OO(A,5)*FED(V,3,GP,E)*FED(VV,1,GP,E)*RHS(A,C(VV,E)) &  ! P(3,1)
                                + FEW(GP,E)*CP(E)*OO(A,6)*FED(V,3,GP,E)*FED(VV,2,GP,E)*RHS(A,C(VV,E))    ! P(3,2)
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CF(E)        *FES(V,GP)    *FES(VV,GP)    *RHS(A,C(VV,E))   ! F         !  4  1
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CU(E)*OM(A,1)*FED(V,1,GP,E)*FES(VV,GP)    *RHS(A,C(VV,E))   ! U^T(1)    ! (5  1)*3
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CU(E)*OM(A,2)*FED(V,2,GP,E)*FES(VV,GP)    *RHS(A,C(VV,E))   ! U^T(2)    
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CU(E)*OM(A,3)*FED(V,3,GP,E)*FES(VV,GP)    *RHS(A,C(VV,E))   ! U^T(3)
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CP(E)*OO(A,1)*FED(V,1,GP,E)*FED(VV,1,GP,E)*RHS(A,C(VV,E))   ! P(1,1)    ! (5  1)*9
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CP(E)*OO(A,2)*FED(V,2,GP,E)*FED(VV,2,GP,E)*RHS(A,C(VV,E))   ! P(2,2)
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CP(E)*OO(A,3)*FED(V,3,GP,E)*FED(VV,3,GP,E)*RHS(A,C(VV,E))   ! P(3,3)
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CP(E)*OO(A,4)*FED(V,1,GP,E)*FED(VV,2,GP,E)*RHS(A,C(VV,E))   ! P(1,2)
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CP(E)*OO(A,5)*FED(V,1,GP,E)*FED(VV,3,GP,E)*RHS(A,C(VV,E))   ! P(1,3)
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CP(E)*OO(A,6)*FED(V,2,GP,E)*FED(VV,3,GP,E)*RHS(A,C(VV,E))   ! P(2,3)
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CP(E)*OO(A,4)*FED(V,2,GP,E)*FED(VV,1,GP,E)*RHS(A,C(VV,E))   ! P(2,1)
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CP(E)*OO(A,5)*FED(V,3,GP,E)*FED(VV,1,GP,E)*RHS(A,C(VV,E))   ! P(3,1)
                  !LHS(A,C(V,E)) = LHS(A,C(V,E)) + FEW(GP,E)*CP(E)*OO(A,6)*FED(V,3,GP,E)*FED(VV,2,GP,E)*RHS(A,C(VV,E))   ! P(3,2)
               END DO ! A
            END DO ! VV
         END DO ! GP
      END DO ! V
   END DO ! E
#ifdef WITHOMP
 END DO ! iTask
!$OMP BARRIER
#else
   EXIT
#endif
END DO ! iColor

END SUBROUTINE ApplyA_AVE1_Tet_SUPG
