!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine determines the threads connected to each element
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE GetAdjacencyTElement(Output_Unit,NumThreads,NumVertices,NumElements,iVperE,Conn,iVertexAssignment, &
                             iElementNumThreads,iElementThreads)
IMPLICIT NONE 
!#define Local_Debug
! Input
PROTEUS_Int, INTENT(IN)    :: Output_Unit
PROTEUS_Int, INTENT(IN)    :: NumThreads,NumVertices,NumElements,iVperE
PROTEUS_Int, INTENT(IN)    :: Conn(iVperE,NumElements)
PROTEUS_Int, INTENT(IN)    :: iVertexAssignment(NumVertices)                ! The thread assignment by vertex
! Arrays to fill
PROTEUS_Int, INTENT(INOUT) :: iElementNumThreads(NumElements)               ! The number of threads adjacent to each element
PROTEUS_Int, INTENT(INOUT) :: iElementThreads(iVperE,NumElements)           ! The threads connected to each element IN SEQUENTIAL ORDERING (1,2,4,7,...)
! Local
PROTEUS_Int  Element,I,J,K,L
PROTEUS_Int  iThread
PROTEUS_Log  Matched

! Process each element
DO Element = 1,NumElements
   iElementNumThreads(Element) = 0
   DO I = 1,iVperE
      iElementThreads(I,Element) = 0
   END DO
   DO I = 1,iVperE
      iThread = iVertexAssignment(Conn(I,Element))
      Matched = .FALSE.
      DO K = 1,iElementNumThreads(Element)
         IF (iThread .EQ. iElementThreads(K,Element)) THEN
            Matched = .TRUE. ! Already present
            EXIT ! The K loop
         END IF
      END DO
      IF (.NOT. Matched) THEN
         ! Put iThread in the correct position
         DO L = 1,iElementNumThreads(Element)
            IF (iThread .LT. iElementThreads(L,Element)) THEN ! Replace position L with iThread and update iThread value
               K=iElementThreads(L,Element);iElementThreads(L,Element)=iThread;iThread=K;
            END IF
         END DO
         iElementNumThreads(Element) = iElementNumThreads(Element) + 1
         iElementThreads(iElementNumThreads(Element),Element) = iThread  ! This is the complete list of threads connected to this element
      END IF
   END DO
END DO

#ifdef Local_Debug
   DO I = 1,NumElements
      WRITE(Output_Unit,'("Element ",I5," is connected to ",I5," threads:",20(I5,2X))') &
         I,iElementNumThreads(I),(iElementThreads(J,I),J=1,iElementNumThreads(I))
   END DO
#endif

END SUBROUTINE GetAdjacencyTElement
