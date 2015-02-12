!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine sets up the spatial and angular matrices
! It also forms the assembled matrix-vector system
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE BuildAngleCubature()
USE CommonBlock
IMPLICIT NONE
!#define Local_Debug
!#define Debug_DumpAssembledMatrix
!#define Local_DumpDebugNZS
!#define Local_DebugNZS

PROTEUS_Int :: Input_Scheme
! Local variables
PROTEUS_Int I

PROTEUS_Real :: Om1,Om2,Om3,SizeOm
PROTEUS_Real :: ConstCoef

! Initialize Omega with made-up numbers
ConstCoef = NumAngles
ConstCoef = 1.0d0 / ConstCoef
DO I = 1,NumAngles
   Om1 = (I-1)**2
   Om2 = (I+1)**2
   Om3 = -2*(I**2)-1
   SizeOm = sqrt( Om1*Om1 + Om2*Om2 + Om3*Om3) 
   Omega(I,1) = Om1 / SizeOm 
   Omega(I,2) = Om2 / SizeOm
   Omega(I,3) = Om3 / SizeOm
   AngleWeights(I) = ConstCoef
   !WRITE(*,'(I,I4,3F10.3)') 'Direction',I,Omega(I,1), Omega(I,2), Omega(I,3)  
   OmegaOmega(I,1) = Omega(I,1) * Omega(I,1)  ! 1 is 1*1
   OmegaOmega(I,2) = Omega(I,2) * Omega(I,2)  ! 2 is 2*2
   OmegaOmega(I,3) = Omega(I,3) * Omega(I,3)  ! 3 is 3*3
   OmegaOmega(I,4) = Omega(I,1) * Omega(I,2)  ! 4 is 1*2
   OmegaOmega(I,5) = Omega(I,1) * Omega(I,3)  ! 5 is 1*3
   OmegaOmega(I,6) = Omega(I,2) * Omega(I,3)  ! 6 is 2*3
END DO

END SUBROUTINE BuildAngleCubature
