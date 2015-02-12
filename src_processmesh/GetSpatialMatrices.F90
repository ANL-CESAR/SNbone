!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine evaluates the spatial matrices
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE GetSpatialMatrices()
USE CommonBlock
IMPLICIT NONE

PROTEUS_Int :: Input_Scheme
! Local variables
PROTEUS_Int I,J,K,L,M
PROTEUS_Int I1,I2,I3,I4
PROTEUS_Int Istart,IBatches,II,FixedBatchSize
PROTEUS_Real :: Om1,Om2,Om3,SizeOm
PROTEUS_Real :: ConstCoef
! Arrays needed to evaluate the derivatives
PROTEUS_FE_Real :: Local_FEDerivatives(4,FENumDim,FEGaussPoints)
PROTEUS_FE_Real :: Local_FEWeights(FEGaussPoints)

! Evaluate the element shape functions
CALL Element_NatTet(FEShapeFunctions,Local_FEDerivatives,Local_FEWeights)

! Compute the surface information
DO I = 1,NumElements
   BCInfo(1,I) = 0
   BCInfo(2,I) = 0
END DO
IF (NumVacuum .GT. 0) THEN
   J = ElementsWithVacuumBCs(1)
   BCInfo(1,J) = 1
   BCInfo(2,J) = 1
   DO I = 2,NumVacuum
      IF (ElementsWithVacuumBCs(I) .EQ. J) THEN ! Same element
         BCInfo(1,J) = BCInfo(1,J) + 1
      ELSE ! Not the same
         J = ElementsWithVacuumBCs(I) ! The element
         BCInfo(1,J) = 1   ! The number of surfaces
         BCInfo(2,J) = I-1 ! The starting surface ID - 1
      END IF
   END DO
END IF

! Made up finite element matrix data to emulate an identity like matrix, not really an identity, but a simple scaling
DO I1=1,NumElements                          !  SUPG                    GLS
   I = ElementLocalToGlobal(I1) ! The original element ordering in the serial case
   DO J = 1,FEGaussPoints
      FEDetJacandWgt(J,I1) = Local_FEWeights(J)*120.0d0
      DO K = 1,FENumDim
         DO L = 1,FEVertices
            FEDerivatives(L,K,J,I1) = 1.1d-24
         END DO
      END DO
   END DO
   ConstTau(I1) = 1.0d0 
   fcoef(I1)    = 1.0d0
   ConstF(I1)   = 1.0d0                      ! sigma_t                  sigma_t + tau*sigma_t*sigma_t
   ConstU(I1)   = 1.12d0  *1.1d-24*(I**0.5)  ! 0.0                      tau*sigma_t
   ConstUT(I1)  = 1.123d0 *1.1d-24*(I**0.5)  ! 1 + tau * sigma_t        1 + tau * sigma_t
   DO J = 1,9
      ucoef(J,I1) = 1.1d-24*(I**0.5+J**0.5)
   END DO
   DO J = 1,36
      pcoef(J,I1) = 1.1d-24*(I**0.5+J**0.5)
   END DO
END DO
! Make up the surface boundary condition info 
NumUnitNormals = NumVacuum/3
IF (NumUnitNormals .LE. 0) NumUnitNormals = 1
J = 0
DO I = 1,NumVacuum
   J = J + 1
   IF (J .GT. NumUnitNormals) J = J - NumUnitNormals
   IF (J .LE. 0) J = 1
   IndexNormal(I) = J ! I am trying to simulate a non-ordered pattern of access for the outward normals.
   Vac_Normals(1,I) = 1.0D0 - (I*0.99d0)/(NumVacuum*1.0d0)
   Vac_Normals(2,I) = 1.0D0 + (I*0.33d0)/(NumVacuum*1.0d0)
   Vac_Normals(3,I) = 1.0D0 - (I*0.88d0)/(NumVacuum*1.0d0)
   SurfaceAxB(I) = 1.0d-24*(I**0.5) ! This eliminates the contribution from the surface term
END DO

END SUBROUTINE GetSpatialMatrices
