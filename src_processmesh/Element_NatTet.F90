!---------------------------------------------------------------------------------------------------------------------------------
!  Element_NatTet_ : Finite element shape functions for the "Natural" TETRAHEDRON element
!  Copyright(c) 2005 Argonne National Laboratory
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE Element_NatTet(FEShapeFunctions,FEDerivatives,FEWeights)
IMPLICIT NONE
! Passed in
PROTEUS_FE_Real, INTENT(INOUT) :: FEShapeFunctions(4,4),FEDerivatives(4,3,4),FEWeights(4)
! Local stuff
PROTEUS_Int, PARAMETER :: TET3D_04PNTS = 4
PROTEUS_Real   :: TET3D_04GAUSSABSX(TET3D_04PNTS) = (/0.585410196624969d0,0.138196601125011d0,0.138196601125011d0, &
                                                        0.138196601125011d0                                         /)
PROTEUS_Real   :: TET3D_04GAUSSABSY(TET3D_04PNTS) = (/0.138196601125011d0,0.585410196624969d0,0.138196601125011d0, &
                                                        0.138196601125011d0                                         /)
PROTEUS_Real   :: TET3D_04GAUSSABSZ(TET3D_04PNTS) = (/0.138196601125011d0,0.138196601125011d0,0.585410196624969d0, &
                                                        0.138196601125011d0                                         /)
PROTEUS_Real   :: TET3D_04GAUSSWGT(TET3D_04PNTS)  = (/0.04166666666666666667d0,0.04166666666666666667d0,           &
                                                        0.04166666666666666667d0,0.04166666666666666667d0           /)

PROTEUS_Int  FIXEDNPE,FIXEDJAC
PROTEUS_Real XXABS,YYABS,ZZABS
PROTEUS_Int  GP

DO GP = 1,TET3D_04PNTS
   XXABS = TET3D_04GAUSSABSX(GP)
   YYABS = TET3D_04GAUSSABSY(GP)
   ZZABS = TET3D_04GAUSSABSZ(GP)

   FEShapeFunctions(1,GP) = 1.0D0 - XXABS - YYABS - ZZABS
   FEDerivatives(1,1,GP)  = -1.0D0
   FEDerivatives(1,2,GP)  = -1.0D0
   FEDerivatives(1,3,GP)  = -1.0D0

   FEShapeFunctions(2,GP) = XXABS
   FEDerivatives(2,1,GP)  = 1.0D0
   FEDerivatives(2,2,GP)  = 0.0D0
   FEDerivatives(2,3,GP)  = 0.0D0

   FEShapeFunctions(3,GP) = YYABS
   FEDerivatives(3,1,GP)  = 0.0D0
   FEDerivatives(3,2,GP)  = 1.0D0
   FEDerivatives(3,3,GP)  = 0.0D0

   FEShapeFunctions(4,GP) = ZZABS
   FEDerivatives(4,1,GP)  = 0.0D0
   FEDerivatives(4,2,GP)  = 0.0D0
   FEDerivatives(4,3,GP)  = 1.0D0

   FEWeights(GP)          = TET3D_04GAUSSWGT(GP)
END DO

END SUBROUTINE Element_NatTet
