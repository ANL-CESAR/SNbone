! This program generates a solution for the solver to obtain
! A * x = b   where x = LHS and b = RHS
! ----------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE GenerateXb(Input_Scheme,LHS_C,LHS_Answer,RHS_C)
USE CommonBlock
#ifdef WITHOMP
  USE OMP_LIB
#endif
!#define Local_Debug
IMPLICIT NONE
PROTEUS_Int,  INTENT(IN)    :: Input_Scheme
! Vectors
PROTEUS_Real, INTENT(INOUT) :: LHS_C(NumAngles,NumVertices),LHS_Answer(NumAngles,NumVertices),RHS_C(NumAngles,NumVertices)

! Local
INTEGER I,J,K,Angle,Element,Space,MyThreadID
INTEGER Lfail
INTEGER StartStop(2)

! We need to define the desired solution
   DO I=1,NumVertices
      DO J = 1,NumAngles
         !RHS_C(J,I)=0.2d0
         LHS_C(J,I)=0.0d0
         !LHS_Answer(J,I) = 0.0d0

         RHS_C(J,I)=0.0d0
         LHS_Answer(J,I) = 0.2d0
      END DO
   END DO
! Construct the source that goes with this answer Ax = b storing it in RHS_C
!$OMP PARALLEL shared(Input_Scheme,Conn,AS_ThreadWiseWork,LHS_Answer,RHS_C) private(MyThreadID)
#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1
#else
   MyThreadID = 1
#endif
IF ((Input_Scheme .EQ. 0) .OR. (Input_Scheme .EQ. 3)) THEN
   CALL SolveWGS_PassThrough_AVE3_NoHPM(LHS_Answer,RHS_C)
ELSE IF (Input_Scheme .EQ. 2) THEN
   CALL SolveWGS_PassThrough_AVE2_NoHPM(LHS_Answer,RHS_C)
ELSE 
   CALL SolveWGS_PassThrough_AVE1_NoHPM(LHS_Answer,RHS_C)
END IF
!$OMP END PARALLEL

#ifdef Local_Debug
   WRITE(Output_Unit,*)' Answer and Source '
   DO J = 1,NumAngles
      DO I = 1,NumVertices
         WRITE(6,300) I,J,LHS_Answer(J,I),RHS_C(J,I)
      END DO
   END DO
   300 FORMAT(I6,2X,I6,3X, 1PE13.6,1X,1PE13.6)
#endif

END SUBROUTINE GenerateXb

