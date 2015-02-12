!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine defines the new vertex ordering with respect to thread wise assigned work
! The idea is that the global vertex stuff should be bunched together
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE GetTVertexOrder(Output_Unit,NumThreads,NumVertices,NumElements,FEVertices,AS_NumColors,   &
                                Conn,iNewElementIndex,AS_ThreadWiseWork,                                &
                                iNewVertexIndex,iReverseVertexIndex)
IMPLICIT NONE 
!#define Local_Debug
! Input
PROTEUS_Int, INTENT(IN)    :: Output_Unit
PROTEUS_Int, INTENT(IN)    :: NumThreads,NumVertices,NumElements,FEVertices,AS_NumColors
PROTEUS_Int, INTENT(IN)    :: Conn(FEVertices,NumElements)
PROTEUS_Int, INTENT(IN)    :: iNewElementIndex(NumElements)                   ! The new element index to use
PROTEUS_Int, INTENT(IN)    :: AS_ThreadWiseWork(2,AS_NumColors,NumThreads) ! Gives the starting (1) and stoping (2) element ids for each thread for the assembly operation
! Arrays to fill
PROTEUS_Int, INTENT(INOUT) :: iNewVertexIndex(NumVertices)                 ! Vertex I was originally at position Index(I)
PROTEUS_Int, INTENT(INOUT) :: iReverseVertexIndex(NumVertices)             ! Vertex I is taken from Index(I)
! Scratch arrays

! Local
PROTEUS_Int I,J,K,L
PROTEUS_Int Element,iElement,Thread
PROTEUS_Int Vertex

! Initialization
DO I = 1,NumVertices
   iNewVertexIndex(I) = 0
   iReverseVertexIndex(I) = 0
END DO

! Loop over each thread adding in the vertices for each thread first to try and make its data contiguous
Vertex = 0
DO Thread = 1,NumThreads
   DO I = 1,AS_NumColors
      DO Element = AS_ThreadWiseWork(1,I,Thread),AS_ThreadWiseWork(2,I,Thread)
         iElement = iNewElementIndex(Element) ! element iElement moves to position Element
         DO K = 1,FEVertices
            L = Conn(K,iElement) ! The global vertex of interest
            IF (iReverseVertexIndex(L) .EQ. 0) THEN
               Vertex = Vertex + 1
               iNewVertexIndex(Vertex) = L      ! vertex L is moved to position Vertex
               iReverseVertexIndex(L) = Vertex  ! vertex (L) was originally at Vertex
            END IF
         END DO
      END DO
   END DO
END DO

#ifdef Local_Debug
   DO I = 1,NumVertices
      WRITE(Output_Unit,'("Vertex ",I5," is replaced by ",I5," or reverse went to ",I5)') I,iNewVertexIndex(I),iReverseVertexIndex(I)
   END DO
#endif

IF (Vertex .NE. NumVertices) THEN 
   WRITE(Output_Unit,*)'Sorry, but not all vertices were assigned in GetTVertexOrder'
   STOP
END IF

END SUBROUTINE GetTVertexOrder
