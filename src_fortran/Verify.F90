! This program checks the solution
! ----------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE Verify(Output_Unit,iSize,Guess,Answer,SomeString)
IMPLICIT NONE
PROTEUS_Int,  INTENT(IN) :: Output_Unit,iSize
PROTEUS_Real, INTENT(IN) :: Guess(iSize),Answer(iSize)
CHARACTER*(*) :: SomeString
! Local
INTEGER I
INTEGER Lfail

Lfail = 0
DO I = 1,iSize
   IF (DABS(Guess(I) - Answer(I)) .GT. 1.0d-5) THEN
      Lfail = I
      EXIT
   END IF
END DO

IF (Lfail .NE. 0) THEN
   WRITE(Output_Unit,2000) SomeString,Lfail
   2000 FORMAT('[SN-KERNEL] ****ERROR**** check failed for ',A12,' at ',I8)
ELSE
   WRITE(Output_Unit,2001) SomeString
   2001 FORMAT('[SN-KERNEL] ***SUCCESS*** check passed for ',A12)
END IF

   !WRITE(6,300) (I,Answer(I),Guess(I),I=1,iSize)
   !300 FORMAT(I6,2X,1PE13.6,1X,1PE13.6)

END SUBROUTINE Verify

