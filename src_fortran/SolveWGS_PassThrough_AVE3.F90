!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine is a pass through point for GMRES
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE SolveWGS_PassThrough_AVE3(RHS_C,LHS_C)
USE CommonBlock
#ifdef WITHOMP
  USE OMP_LIB
#endif
IMPLICIT NONE 
! Passed in
PROTEUS_Real LHS_C(NumAngles,NumVertices),RHS_C(NumAngles,NumVertices)
PROTEUS_Int MyThreadID
PROTEUS_Int I,J,K,iRowStart,iRowEnd

#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1
   I = NumVertices/NumThreads
   iRowStart = (MyThreadID-1)*I + 1
   IRowEnd   = MyThreadID*I
   IF (MyThreadID .EQ. NumThreads) IRowEnd = NumVertices
   !WRITE(6,*)'iRowStart=',iRowStart,' iRowEnd=',iRowEnd
#else
   MyThreadID = 1
   iRowStart = 1
   iRowEnd   = NumVertices
#endif
#ifdef WITHBGQHPM
   IF (MyThreadID .EQ. 1) call hpm_start('AVE3_ApplyA')
#endif
! This barrier ensures that the incoming vector is fully defined by all threads
!$OMP BARRIER
DO I = iRowStart,iRowEnd
   DO J = NZS_RowLoc(I),NZS_RowLoc(I+1)-1
      DO K = 1,NumAngles
         LHS_C(K,I) = LHS_C(K,I) + NZS_Data(K,J)*RHS_C(K,NZS_ColNum(J))
      END DO
   END DO
END DO
! This second barrier is needed because the FGMRES threadwise splitting is currently different than the above. Likely need to change that
!$OMP BARRIER
#ifdef WITHBGQHPM
   IF (MyThreadID .EQ. 1) call hpm_stop('AVE3_ApplyA') ! Stops the hardware counters
#endif
END SUBROUTINE SolveWGS_PassThrough_AVE3


SUBROUTINE SolveWGS_PassThrough_AVE3_NoHPM(RHS_C,LHS_C)
USE CommonBlock
#ifdef WITHOMP
  USE OMP_LIB
#endif
IMPLICIT NONE 
! Passed in
PROTEUS_Real LHS_C(NumAngles,NumVertices),RHS_C(NumAngles,NumVertices)
PROTEUS_Int MyThreadID
PROTEUS_Int I,J,K,iRowStart,iRowEnd

#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1
   I = NumVertices/NumThreads
   iRowStart = (MyThreadID-1)*I + 1
   IRowEnd   = MyThreadID*I
   IF (MyThreadID .EQ. NumThreads) IRowEnd = NumVertices
   !WRITE(6,*)'iRowStart=',iRowStart,' iRowEnd=',iRowEnd
#else
   MyThreadID = 1
   iRowStart = 1
   iRowEnd   = NumVertices
#endif
! This barrier ensures that the incoming vector is fully defined by all threads
!$OMP BARRIER
DO I = iRowStart,iRowEnd
   DO J = NZS_RowLoc(I),NZS_RowLoc(I+1)-1
      DO K = 1,NumAngles
         LHS_C(K,I) = LHS_C(K,I) + NZS_Data(K,J)*RHS_C(K,NZS_ColNum(J))
      END DO
   END DO
END DO
! This second barrier is needed because the FGMRES threadwise splitting is currently different than the above. Likely need to change that
!$OMP BARRIER
END SUBROUTINE SolveWGS_PassThrough_AVE3_NoHPM

