!  Copyright(c) 2005 Argonne National Laboratory
!---------------------------------------------------------------------------------------------------------------------------------
!  Controls the process of error reporting from PN2ND
!---------------------------------------------------------------------------------------------------------------------------------
!  Copyright(c) 2005 Argonne National Laboratory
SUBROUTINE Basic_CheckError(Output_Unit,ErrorCondition,StringToPrint)
IMPLICIT NONE
#include "PROTEUS_Preprocess.h"

! Passed in
PROTEUS_Int     Output_Unit
PROTEUS_Int     ErrorCondition  ! Some error code returned from HDF5
CHARACTER(*) StringToPrint ! The string to print

IF (ErrorCondition .NE. 0) THEN
   WRITE(Output_Unit,200) ErrorCondition
   WRITE(Output_Unit,205) StringToPrint
   WRITE(Output_Unit,300) 
   CALL Basic_Abort
   200 FORMAT('[ERROR]...Returned the error code.....<',I5,'>')
   205 FORMAT('[ERROR]...In the code you told me to report.<',A40,'>')
   300 FORMAT('[ERROR]...Sorry, but I must stop')
END IF

END SUBROUTINE Basic_CheckError
