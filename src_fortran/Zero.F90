! This program checks the solution
! ----------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE Zero(iStart,iEnd,iSize,VectorToInitialize)
IMPLICIT NONE
PROTEUS_Int,  INTENT(IN)    :: iStart,iEnd,iSize
PROTEUS_Real, INTENT(INOUT) :: VectorToInitialize(iSize)
! Local
INTEGER I

DO I = iStart,iEnd
   VectorToInitialize(I) = 0.0d0
END DO

END SUBROUTINE Zero

SUBROUTINE Zero_Threaded(iSize,VectorToInitialize)
#ifdef WITHOMP
  USE OMP_LIB
#endif
IMPLICIT NONE
PROTEUS_Int,  INTENT(IN)    :: iSize
PROTEUS_Real, INTENT(INOUT) :: VectorToInitialize(iSize)
! Local
INTEGER I

!$OMP PARALLEL DO SCHEDULE(STATIC)
   DO I = 1,iSize
      VectorToInitialize(I) = 0.0d0
   END DO
!$OMP END PARALLEL DO

END SUBROUTINE Zero_Threaded
