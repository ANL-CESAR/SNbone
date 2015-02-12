!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine sets up the spatial and angular matrices
! It also forms the assembled matrix-vector system
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE AssembleNZmatrix(Input_Scheme)
USE CommonBlock
IMPLICIT NONE
!#define Local_Debug
!#define Debug_DumpAssembledMatrix
!#define Local_DumpDebugNZS
!#define Local_DebugNZS

PROTEUS_Int :: Input_Scheme
! Local variables
PROTEUS_Int I,J,K
PROTEUS_Int I1,I2,I3,I4
PROTEUS_Real :: ConstCoef

IF ((Input_Scheme .EQ. 0) .OR. (Input_Scheme .EQ. 3)) THEN
! Threadable section
!   DO iColor = 1,AS_NumColors
!      DO iTask = 1,TasksPerThread
!         iThread = ThreadID*TasksPerThread + iTask
!         DO K = AS_ThreadWiseWork(1,iThread,iColor),AS_ThreadWiseWork(2,iThread,iColor)

   ! Form and store the matrix
   DO K = 1,NumElements
      ! Row 1 ---------------------------------------------------------------
      I = Conn(1,K) ! The row vertex
      DO J = NZS_RowLoc(I),NZS_RowLoc(I+1)-1
         IF (Conn(1,K) .EQ. NZS_ColNum(J)) I1 = J
         IF (Conn(2,K) .EQ. NZS_ColNum(J)) I2 = J
         IF (Conn(3,K) .EQ. NZS_ColNum(J)) I3 = J
         IF (Conn(4,K) .EQ. NZS_ColNum(J)) I4 = J
      END DO
      ! 1,1
      DO J = 1,NumAngles
         ConstCoef = 2.0D0*fcoef(K)*ConstF(K)                                                           &
         - ConstU(K)*((Omega(J,1)*ucoef(1,K) + Omega(J,2)*ucoef(2,K) + Omega(J,3)*ucoef(3,K))  &
                 +(Omega(J,1)*ucoef(4,K) + Omega(J,2)*ucoef(5,K) + Omega(J,3)*ucoef(6,K))  &
                 +(Omega(J,1)*ucoef(7,K) + Omega(J,2)*ucoef(8,K) + Omega(J,3)*ucoef(9,K))) &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(01,K)+OmegaOmega(J,2)*pcoef(02,K)+OmegaOmega(J,3)*pcoef(03,K)     &
                 +OmegaOmega(J,4)*pcoef(04,K)+OmegaOmega(J,5)*pcoef(05,K)+OmegaOmega(J,6)*pcoef(06,K))    &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(07,K)+OmegaOmega(J,2)*pcoef(08,K)+OmegaOmega(J,3)*pcoef(09,K)     &
                 +OmegaOmega(J,4)*pcoef(10,K)+OmegaOmega(J,5)*pcoef(11,K)+OmegaOmega(J,6)*pcoef(12,K))    &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(13,K)+OmegaOmega(J,2)*pcoef(14,K)+OmegaOmega(J,3)*pcoef(15,K)     &
                 +OmegaOmega(J,4)*pcoef(16,K)+OmegaOmega(J,5)*pcoef(17,K)+OmegaOmega(J,6)*pcoef(18,K))    &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(19,K)+OmegaOmega(J,2)*pcoef(20,K)+OmegaOmega(J,3)*pcoef(21,K)     &
                 +OmegaOmega(J,4)*pcoef(22,K)+OmegaOmega(J,5)*pcoef(23,K)+OmegaOmega(J,6)*pcoef(24,K))    &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(25,K)+OmegaOmega(J,2)*pcoef(26,K)+OmegaOmega(J,3)*pcoef(27,K)     &
                 +OmegaOmega(J,4)*pcoef(28,K)+OmegaOmega(J,5)*pcoef(29,K)+OmegaOmega(J,6)*pcoef(30,K))    &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(31,K)+OmegaOmega(J,2)*pcoef(32,K)+OmegaOmega(J,3)*pcoef(33,K)     &
                 +OmegaOmega(J,4)*pcoef(34,K)+OmegaOmega(J,5)*pcoef(35,K)+OmegaOmega(J,6)*pcoef(36,K))    
         NZS_Data(J,I1) = NZS_Data(J,I1) + ConstCoef
      END DO
      ! 1,2
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                                                 &
         + ConstU(K)* (Omega(J,1)*ucoef(1,K) + Omega(J,2)*ucoef(2,K) + Omega(J,3)*ucoef(3,K))  &
         - ConstTau(K)*(OmegaOmega(J,1)*pcoef(01,K)+OmegaOmega(J,2)*pcoef(02,K)+OmegaOmega(J,3)*pcoef(03,K)     &
                 +OmegaOmega(J,4)*pcoef(04,K)+OmegaOmega(J,5)*pcoef(05,K)+OmegaOmega(J,6)*pcoef(06,K))
         NZS_Data(J,I2) = NZS_Data(J,I2) + ConstCoef
      END DO
      ! 1,3
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                                                 &
         + ConstU(K)* (Omega(J,1)*ucoef(4,K) + Omega(J,2)*ucoef(5,K) + Omega(J,3)*ucoef(6,K))  &
         - ConstTau(K)*(OmegaOmega(J,1)*pcoef(07,K)+OmegaOmega(J,2)*pcoef(08,K)+OmegaOmega(J,3)*pcoef(09,K)     &
                 +OmegaOmega(J,4)*pcoef(10,K)+OmegaOmega(J,5)*pcoef(11,K)+OmegaOmega(J,6)*pcoef(12,K))    &
         - ConstTau(K)*(OmegaOmega(J,1)*pcoef(19,K)+OmegaOmega(J,2)*pcoef(20,K)+OmegaOmega(J,3)*pcoef(21,K)     &
                 +OmegaOmega(J,4)*pcoef(22,K)+OmegaOmega(J,5)*pcoef(23,K)+OmegaOmega(J,6)*pcoef(24,K))
         NZS_Data(J,I3) = NZS_Data(J,I3) + ConstCoef
      END DO
      ! 1,4
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                                                 &
         + ConstU(K)* (Omega(J,1)*ucoef(7,K) + Omega(J,2)*ucoef(8,K) + Omega(J,3)*ucoef(9,K))  &
         - ConstTau(K)*(OmegaOmega(J,1)*pcoef(13,K)+OmegaOmega(J,2)*pcoef(14,K)+OmegaOmega(J,3)*pcoef(15,K)     &
                 +OmegaOmega(J,4)*pcoef(16,K)+OmegaOmega(J,5)*pcoef(17,K)+OmegaOmega(J,6)*pcoef(18,K))    &
         - ConstTau(K)*(OmegaOmega(J,1)*pcoef(19,K)+OmegaOmega(J,2)*pcoef(20,K)+OmegaOmega(J,3)*pcoef(21,K)     &
                 +OmegaOmega(J,4)*pcoef(22,K)+OmegaOmega(J,5)*pcoef(23,K)+OmegaOmega(J,6)*pcoef(24,K))    &
         - ConstTau(K)*(OmegaOmega(J,1)*pcoef(25,K)+OmegaOmega(J,2)*pcoef(26,K)+OmegaOmega(J,3)*pcoef(27,K)    &
                 +OmegaOmega(J,4)*pcoef(28,K)+OmegaOmega(J,5)*pcoef(29,K)+OmegaOmega(J,6)*pcoef(30,K))   &
         - ConstTau(K)*(OmegaOmega(J,1)*pcoef(31,K)+OmegaOmega(J,2)*pcoef(32,K)+OmegaOmega(J,3)*pcoef(33,K)    &
                 +OmegaOmega(J,4)*pcoef(34,K)+OmegaOmega(J,5)*pcoef(35,K)+OmegaOmega(J,6)*pcoef(36,K))    
         NZS_Data(J,I4) = NZS_Data(J,I4) + ConstCoef
      END DO
      ! Row 2 ---------------------------------------------------------------
      I = Conn(2,K) ! The row vertex
      DO J = NZS_RowLoc(I),NZS_RowLoc(I+1)-1
         IF (Conn(1,K) .EQ. NZS_ColNum(J)) I1 = J
         IF (Conn(2,K) .EQ. NZS_ColNum(J)) I2 = J
         IF (Conn(3,K) .EQ. NZS_ColNum(J)) I3 = J
         IF (Conn(4,K) .EQ. NZS_ColNum(J)) I4 = J
      END DO
      ! 2,1
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                          &
         - ConstU(K)*((Omega(J,1)*ucoef(1,K) + Omega(J,2)*ucoef(2,K) + Omega(J,3)*ucoef(3,K))   &
                 +(Omega(J,1)*ucoef(4,K) + Omega(J,2)*ucoef(5,K) + Omega(J,3)*ucoef(6,K))   &
                 +(Omega(J,1)*ucoef(7,K) + Omega(J,2)*ucoef(8,K) + Omega(J,3)*ucoef(9,K)))  &
         - ConstTau(K)*(OmegaOmega(J,1)*pcoef(01,K)+OmegaOmega(J,2)*pcoef(02,K)+OmegaOmega(J,3)*pcoef(03,K)      &
                 +OmegaOmega(J,4)*pcoef(04,K)+OmegaOmega(J,5)*pcoef(05,K)+OmegaOmega(J,6)*pcoef(06,K))     &
         - ConstTau(K)*(OmegaOmega(J,1)*pcoef(19,K)+OmegaOmega(J,2)*pcoef(20,K)+OmegaOmega(J,3)*pcoef(21,K)      &
                 +OmegaOmega(J,4)*pcoef(22,K)+OmegaOmega(J,5)*pcoef(23,K)+OmegaOmega(J,6)*pcoef(24,K))     &
         - ConstTau(K)*(OmegaOmega(J,1)*pcoef(25,K)+OmegaOmega(J,2)*pcoef(26,K)+OmegaOmega(J,3)*pcoef(27,K)      &
                 +OmegaOmega(J,4)*pcoef(28,K)+OmegaOmega(J,5)*pcoef(29,K)+OmegaOmega(J,6)*pcoef(30,K))
         NZS_Data(J,I1) = NZS_Data(J,I1) + ConstCoef
      END DO
      ! 2,2
      DO J = 1,NumAngles
         ConstCoef = 2.0D0*fcoef(K)*ConstF(K)                                                                 &
         + ConstU(K)*(Omega(J,1)*ucoef(1,K) + Omega(J,2)*ucoef(2,K) + Omega(J,3)*ucoef(3,K))  &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(01,K)+OmegaOmega(J,2)*pcoef(02,K)+OmegaOmega(J,3)*pcoef(03,K)    &
                 +OmegaOmega(J,4)*pcoef(04,K)+OmegaOmega(J,5)*pcoef(05,K)+OmegaOmega(J,6)*pcoef(06,K))
         NZS_Data(J,I2) = NZS_Data(J,I2) + ConstCoef
      END DO
      ! 2,3
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                                                 &
         + ConstU(K)* (Omega(J,1)*ucoef(4,K) + Omega(J,2)*ucoef(5,K) + Omega(J,3)*ucoef(6,K))  &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(19,K)+OmegaOmega(J,2)*pcoef(20,K)+OmegaOmega(J,3)*pcoef(21,K)     &
                 +OmegaOmega(J,4)*pcoef(22,K)+OmegaOmega(J,5)*pcoef(23,K)+OmegaOmega(J,6)*pcoef(24,K))
         NZS_Data(J,I3) = NZS_Data(J,I3) + ConstCoef
      END DO
      ! 2,4
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                                                 &
         + ConstU(K)* (Omega(J,1)*ucoef(7,K) + Omega(J,2)*ucoef(8,K) + Omega(J,3)*ucoef(9,K))  &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(25,K)+OmegaOmega(J,2)*pcoef(26,K)+OmegaOmega(J,3)*pcoef(27,K)     &
                 +OmegaOmega(J,4)*pcoef(28,K)+OmegaOmega(J,5)*pcoef(29,K)+OmegaOmega(J,6)*pcoef(30,K))
         NZS_Data(J,I4) = NZS_Data(J,I4) + ConstCoef
      END DO
      ! Row 3 ---------------------------------------------------------------
      I = Conn(3,K) ! The row vertex
      DO J = NZS_RowLoc(I),NZS_RowLoc(I+1)-1
         IF (Conn(1,K) .EQ. NZS_ColNum(J)) I1 = J
         IF (Conn(2,K) .EQ. NZS_ColNum(J)) I2 = J
         IF (Conn(3,K) .EQ. NZS_ColNum(J)) I3 = J
         IF (Conn(4,K) .EQ. NZS_ColNum(J)) I4 = J
      END DO
      ! 3,1
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                          &
         - (ConstU(K)*(Omega(J,1)*ucoef(1,K) + Omega(J,2)*ucoef(2,K) + Omega(J,3)*ucoef(3,K))   &
           +ConstU(K)*(Omega(J,1)*ucoef(4,K) + Omega(J,2)*ucoef(5,K) + Omega(J,3)*ucoef(6,K))   &
           +ConstU(K)*(Omega(J,1)*ucoef(7,K) + Omega(J,2)*ucoef(8,K) + Omega(J,3)*ucoef(9,K)))  &
         -  ConstTau(K)*(OmegaOmega(J,1)*pcoef(07,K)+OmegaOmega(J,2)*pcoef(08,K)+OmegaOmega(J,3)*pcoef(09,K)     &
                  +OmegaOmega(J,4)*pcoef(10,K)+OmegaOmega(J,5)*pcoef(11,K)+OmegaOmega(J,6)*pcoef(12,K))    &
         -  ConstTau(K)*(OmegaOmega(J,1)*pcoef(31,K)+OmegaOmega(J,2)*pcoef(32,K)+OmegaOmega(J,3)*pcoef(33,K)     &
                  +OmegaOmega(J,4)*pcoef(34,K)+OmegaOmega(J,5)*pcoef(35,K)+OmegaOmega(J,6)*pcoef(36,K))
         NZS_Data(J,I1) = NZS_Data(J,I1) + ConstCoef
      END DO
      ! 3,2
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                         &
         + ConstU(K)*(Omega(J,1)*ucoef(1,K) + Omega(J,2)*ucoef(2,K) + Omega(J,3)*ucoef(3,K))
         NZS_Data(J,I2) = NZS_Data(J,I2) + ConstCoef
      END DO
      ! 3,3
      DO J = 1,NumAngles
         ConstCoef = 2.0D0*fcoef(K)*ConstF(K)                                   &
         + ConstU(K)*(Omega(J,1)*ucoef(4,K) + Omega(J,2)*ucoef(5,K) + Omega(J,3)*ucoef(6,K))   &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(07,K)+OmegaOmega(J,2)*pcoef(08,K)+OmegaOmega(J,3)*pcoef(09,K)     &
                 +OmegaOmega(J,4)*pcoef(10,K)+OmegaOmega(J,5)*pcoef(11,K)+OmegaOmega(J,6)*pcoef(12,K)) 
         NZS_Data(J,I3) = NZS_Data(J,I3) + ConstCoef
      END DO
      ! 3,4
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                         &
         + ConstU(K)*(Omega(J,1)*ucoef(7,K) + Omega(J,2)*ucoef(8,K) + Omega(J,3)*ucoef(9,K))   &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(31,K)+OmegaOmega(J,2)*pcoef(32,K)+OmegaOmega(J,3)*pcoef(33,K)     &
                 +OmegaOmega(J,4)*pcoef(34,K)+OmegaOmega(J,5)*pcoef(35,K)+OmegaOmega(J,6)*pcoef(36,K))
         NZS_Data(J,I4) = NZS_Data(J,I4) + ConstCoef
      END DO
      ! Row 4 ---------------------------------------------------------------
      I = Conn(4,K) ! The row vertex
      DO J = NZS_RowLoc(I),NZS_RowLoc(I+1)-1
         IF (Conn(1,K) .EQ. NZS_ColNum(J)) I1 = J
         IF (Conn(2,K) .EQ. NZS_ColNum(J)) I2 = J
         IF (Conn(3,K) .EQ. NZS_ColNum(J)) I3 = J
         IF (Conn(4,K) .EQ. NZS_ColNum(J)) I4 = J
      END DO
      !WRITE(6,'("I1-I4",10I6)') I,I1,I2,I3,I4
      ! 4,1
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                          &
         - (ConstU(K)*(Omega(J,1)*ucoef(1,K) + Omega(J,2)*ucoef(2,K) + Omega(J,3)*ucoef(3,K))   &
           +ConstU(K)*(Omega(J,1)*ucoef(4,K) + Omega(J,2)*ucoef(5,K) + Omega(J,3)*ucoef(6,K))   &
           +ConstU(K)*(Omega(J,1)*ucoef(7,K) + Omega(J,2)*ucoef(8,K) + Omega(J,3)*ucoef(9,K)))  &
         -  ConstTau(K)*(OmegaOmega(J,1)*pcoef(13,K)+OmegaOmega(J,2)*pcoef(14,K)+OmegaOmega(J,3)*pcoef(15,K)  &
                  +OmegaOmega(J,4)*pcoef(16,K)+OmegaOmega(J,5)*pcoef(17,K)+OmegaOmega(J,6)*pcoef(18,K))
         NZS_Data(J,I1) = NZS_Data(J,I1) + ConstCoef
      END DO
      ! 4,2
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                         &
         + ConstU(K)*(Omega(J,1)*ucoef(1,K) + Omega(J,2)*ucoef(2,K) + Omega(J,3)*ucoef(3,K))
         NZS_Data(J,I2) = NZS_Data(J,I2) + ConstCoef
      END DO
      ! 4,3
      DO J = 1,NumAngles
         ConstCoef = fcoef(K)*ConstF(K)                                                      &
         + ConstU(K)*(Omega(J,1)*ucoef(4,K) + Omega(J,2)*ucoef(5,K) + Omega(J,3)*ucoef(6,K))
         NZS_Data(J,I3) = NZS_Data(J,I3) + ConstCoef
      END DO
      ! 4,4
      DO J = 1,NumAngles
         ConstCoef = 2.0D0*fcoef(K)*ConstF(K)                                                                     &
         + ConstU(K)*(Omega(J,1)*ucoef(7,K) + Omega(J,2)*ucoef(8,K) + Omega(J,3)*ucoef(9,K))                &
         + ConstTau(K)*(OmegaOmega(J,1)*pcoef(13,K)+OmegaOmega(J,2)*pcoef(14,K)+OmegaOmega(J,3)*pcoef(15,K) &
                 +OmegaOmega(J,4)*pcoef(16,K)+OmegaOmega(J,5)*pcoef(17,K)+OmegaOmega(J,6)*pcoef(18,K))    
         NZS_Data(J,I4) = NZS_Data(J,I4) + ConstCoef
      END DO
   END DO ! Element loop

#ifdef Local_DumpDebugNZS
   DO I = 1,NumVertices
      DO K = NZS_RowLoc(I),NZS_RowLoc(I+1)-1
         DO J = 1,NumAngles
            WRITE(23,777) 0,I,NZS_ColNum(K),J,NZS_Data(J,K)
         END DO
      END DO
   END DO
#endif

#ifdef Local_DebugNZS
   DO I = 1,NumVertices
      WRITE(Output_Unit,'("[SN-KERNEL] Row ",I8," stores data between ",I6,":",I6)') I,NZS_RowLoc(I),NZS_RowLoc(I+1)-1
      WRITE(Output_Unit,7721) (K,K=1,NZS_RowLoc(I+1)-NZS_RowLoc(I))
      WRITE(Output_Unit,7721) (NZS_ColNum(K),K=NZS_RowLoc(I),NZS_RowLoc(I+1)-1)
      DO J = 1,NumAngles
         WRITE(Output_Unit,7722) (NZS_Data(J,K),K=NZS_RowLoc(I),NZS_RowLoc(I+1)-1)
      END DO
      7721 FORMAT(9999(4X,I3,4X))
      7722 FORMAT(9999(1X,F9.3,1X))
   END DO
#endif
   !DO I = 1,NumVertices
   !   WRITE(Output_Unit,'("[SN-KERNEL] Row ",I8," stores data between ",I6,":",I6)') I,NZS_RowLoc(I),NZS_RowLoc(I+1)-1
   !END DO
END IF ! Scheme 0 and 3 

END SUBROUTINE AssembleNZmatrix
