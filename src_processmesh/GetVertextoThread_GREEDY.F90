!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine is unnecessary as the GREEDY approach assigns ownership by an element adjacency rule
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE GetVertextoThread_GREEDY(Output_Unit,NumThreads,NumVertices,GlobalXYZ,NumElements,FEVertices,Conn,  iVertexAssignment)
IMPLICIT NONE 
!#define Local_Debug
! Input
PROTEUS_Int, INTENT(IN)    :: Output_Unit
PROTEUS_Int, INTENT(IN)    :: NumThreads,NumVertices
PROTEUS_Real, INTENT(IN)   :: GlobalXYZ(NumVertices,3)
PROTEUS_Int, INTENT(IN)    :: NumElements,FEVertices
PROTEUS_Int, INTENT(IN)    :: Conn(FEVertices,NumElements)
! Arrays to fill
PROTEUS_Int, INTENT(INOUT) :: iVertexAssignment(NumVertices)      ! The thread assignment by vertex

PROTEUS_Int  I,J,K,iThread
PROTEUS_Int  iNumX,iNumY,iNumZ
PROTEUS_Real MinXYZ(3),MaxXYZ(3),deltaX,deltaY,dminX,dmaxX,dminY,dmaxY

! Initialization
iVertexAssignment(1:NumVertices) = 1

END SUBROUTINE GetVertextoThread_GREEDY
