! This subroutine performs the combined F,UT,U matrix application for a bunch of elements and angles for a set of tets
! The matrix operation can be viewed as: LHS(A,E) = C(A,E)*RHS(A,E)
! ------------
! FLOP SUMMARY
! There is significant pointer math which I ignore
!  832 additions/subtractions /E/A
! 4096 multiplications /E/A
! 4928 flops per element per angle
! ----
SUBROUTINE ApplyA_AVE2_Tet_SUPG(NumElements,NumAngles,NumVertices, C,  &
                                AS_NumColors,AS_NumThreads,TasksPerThread,MyThreadID,ThreadWiseWork,  &
                                    CF, CU, CP, FES,FED,FEW, OM,OO,        &
                                    LHS,RHS,SRHS,SLHS)
IMPLICIT NONE 
#include "PROTEUS_Preprocess.h"
!#define Local_TryLHStoo
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
! Scratch vectors
PROTEUS_Real, INTENT(INOUT) :: SRHS(NumAngles,4),SLHS(NumAngles,4)

! Local
PROTEUS_Int :: E,A,V,VV,GP,iColor,iStart,iEnd,II
PROTEUS_Real :: LocalConst1,LocalConst2,LocalConst3,Products(10)
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
      ! Pull a copy of the vectors
      DO V = 1,4
         VV = C(V,E) ! The assembled vector position
         DO A = 1,NumAngles
            SRHS(A,V) = RHS(A,VV)
#ifdef Local_TryLHStoo
            SLHS(A,V) = LHS(A,VV)
#endif
         END DO
      END DO
      DO V = 1,4
#ifndef Local_TryLHStoo
         II = C(V,E)       
#endif
         DO GP = 1,4
            LocalConst1 = FEW(GP,E)*CF(E)
            LocalConst2 = FEW(GP,E)*CU(E)
            LocalConst3 = FEW(GP,E)*CP(E)
            DO VV = 1,4
               Products( 1) = LocalConst1*FES(V,GP)    *FES(VV,GP)    
               Products( 2) = LocalConst2*FED(V,1,GP,E)*FES(VV,GP) 
               Products( 3) = LocalConst2*FED(V,2,GP,E)*FES(VV,GP)    
               Products( 4) = LocalConst2*FED(V,3,GP,E)*FES(VV,GP)    
               Products( 5) = LocalConst3*FED(V,1,GP,E)*FED(VV,1,GP,E)
               Products( 6) = LocalConst3*FED(V,2,GP,E)*FED(VV,2,GP,E)
               Products( 7) = LocalConst3*FED(V,3,GP,E)*FED(VV,3,GP,E)
               Products( 8) = LocalConst3*(FED(V,1,GP,E)*FED(VV,2,GP,E)+FED(V,2,GP,E)*FED(VV,1,GP,E))
               Products( 9) = LocalConst3*(FED(V,1,GP,E)*FED(VV,3,GP,E)+FED(V,3,GP,E)*FED(VV,1,GP,E))
               Products(10) = LocalConst3*(FED(V,2,GP,E)*FED(VV,3,GP,E)+FED(V,3,GP,E)*FED(VV,2,GP,E))
               DO A = 1,NumAngles
#ifdef Local_TryLHStoo
                  SLHS(A,V) = SLHS(A,V) + Products( 1)        *SRHS(A,VV) &  ! F
                                        + Products( 2)*OM(A,1)*SRHS(A,VV) &  ! U^T(1)
                                        + Products( 3)*OM(A,2)*SRHS(A,VV) &  ! U^T(2)    
                                        + Products( 4)*OM(A,3)*SRHS(A,VV) &  ! U^T(3)
                                        + Products( 5)*OO(A,1)*SRHS(A,VV) &  ! P(1,1)
                                        + Products( 6)*OO(A,2)*SRHS(A,VV) &  ! P(2,2)
                                        + Products( 7)*OO(A,3)*SRHS(A,VV) &  ! P(3,3)
                                        + Products( 8)*OO(A,4)*SRHS(A,VV) &  ! P(1,2)
                                        + Products( 9)*OO(A,5)*SRHS(A,VV) &  ! P(1,3)
                                        + Products(10)*OO(A,6)*SRHS(A,VV)   ! P(2,3)
#else
                  LHS(A,II) = LHS(A,II) + Products( 1)        *SRHS(A,VV) &  ! F
                                        + Products( 2)*OM(A,1)*SRHS(A,VV) &  ! U^T(1)
                                        + Products( 3)*OM(A,2)*SRHS(A,VV) &  ! U^T(2)    
                                        + Products( 4)*OM(A,3)*SRHS(A,VV) &  ! U^T(3)
                                        + Products( 5)*OO(A,1)*SRHS(A,VV) &  ! P(1,1)
                                        + Products( 6)*OO(A,2)*SRHS(A,VV) &  ! P(2,2)
                                        + Products( 7)*OO(A,3)*SRHS(A,VV) &  ! P(3,3)
                                        + Products( 8)*OO(A,4)*SRHS(A,VV) &  ! P(1,2)
                                        + Products( 9)*OO(A,5)*SRHS(A,VV) &  ! P(1,3)
                                        + Products(10)*OO(A,6)*SRHS(A,VV)   ! P(2,3)
#endif
                  !SLHS(A,V) = SLHS(A,V) + Products( 1)        *SRHS(A,VV)   ! F
                  !SLHS(A,V) = SLHS(A,V) + Products( 2)*OM(A,1)*SRHS(A,VV)   ! U^T(1)
                  !SLHS(A,V) = SLHS(A,V) + Products( 3)*OM(A,2)*SRHS(A,VV)   ! U^T(2)    
                  !SLHS(A,V) = SLHS(A,V) + Products( 4)*OM(A,3)*SRHS(A,VV)   ! U^T(3)
                  !SLHS(A,V) = SLHS(A,V) + Products( 5)*OO(A,1)*SRHS(A,VV)   ! P(1,1)
                  !SLHS(A,V) = SLHS(A,V) + Products( 6)*OO(A,2)*SRHS(A,VV)   ! P(2,2)
                  !SLHS(A,V) = SLHS(A,V) + Products( 7)*OO(A,3)*SRHS(A,VV)   ! P(3,3)
                  !SLHS(A,V) = SLHS(A,V) + Products( 8)*OO(A,4)*SRHS(A,VV)   ! P(1,2)
                  !SLHS(A,V) = SLHS(A,V) + Products( 9)*OO(A,5)*SRHS(A,VV)   ! P(1,3)
                  !SLHS(A,V) = SLHS(A,V) + Products(10)*OO(A,6)*SRHS(A,VV)   ! P(2,3)
               END DO ! A
            END DO ! VV
         END DO ! GP
      END DO ! V
      ! Store the LHS result
#ifdef Local_TryLHStoo
      DO V = 1,4
         VV = C(V,E) ! The assembled vector position
         DO A = 1,NumAngles
            LHS(A,VV) = SLHS(A,V)
         END DO
      END DO
#endif
   END DO ! E
#ifdef WITHOMP
 END DO ! iTask
#else
   EXIT
#endif
!$OMP BARRIER
END DO ! iColor

END SUBROUTINE ApplyA_AVE2_Tet_SUPG
