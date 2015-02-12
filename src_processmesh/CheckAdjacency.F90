!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine checks the adjacency of seperable work to ensure the connectivity does not overlap with other seperable work
! iWorkOption = Work options
! 0) Initialize the checking array
! 1) Setup the checking array for seperable work between iTarget1 and iTarget2
! 2) Merge the seperable work iTarget2 with iTarget1
! 3) Check whether the seperable work between iTarget1 and iTarget2 intersects with the current seperable work ID
!    LReturnStatus = TRUE is they are seperate, = FALSE means they intersect
! 4) Merge the seperable work flagging info but do not merge the seperable work
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE CheckAdjacency(Output_Unit,NumVertices,NumElements,FEVertices,   &
                              iTarget1,iTarget2,iWorkOption,LReturnStatus, &
                              iSeperableID,Conn,LVertexFlags)
IMPLICIT NONE 
!#define Local_Debug
! Input
PROTEUS_Int, INTENT(IN)    :: Output_Unit
PROTEUS_Int, INTENT(IN)    :: NumVertices,NumElements,FEVertices  
PROTEUS_Int, INTENT(IN)    :: iTarget1,iTarget2                ! The targeted seperable pieces of work
PROTEUS_Int, INTENT(IN)    :: iWorkOption                      ! 
PROTEUS_Log, INTENT(INOUT) :: LReturnStatus                    ! Return logical variable identifying 
PROTEUS_Int, INTENT(INOUT) :: iSeperableID(NumElements)        ! The seperable work ID of each element
PROTEUS_Int, INTENT(IN)    :: Conn(FEVertices,NumElements)         ! The connectivity of each element
PROTEUS_Log, INTENT(INOUT) :: LVertexFlags(NumVertices)        ! Scratch flag storage

! Local
PROTEUS_Int I,J

#ifdef Local_Debug
WRITE(Output_Unit,'("iWorkOption=",I3," iTarget1 =",I3," iTarget2=",I3)') iWorkOption,iTarget1,iTarget2
WRITE(Output_Unit,'(10(10L1,1X))') (LVertexFlags(I),I=1,NumVertices)
#endif

IF (iWorkOption .EQ. 0) THEN ! Inialize the array
   LReturnStatus = .TRUE.
   DO I = 1,NumVertices
      LVertexFlags(I) = .FALSE.
   END DO
END IF
IF (iWorkOption .EQ. 1) THEN ! Highlight the vertices of the seperable piece of work we are targeting
   LReturnStatus = .TRUE.
   DO I = 1,NumElements
      IF ((iTarget1 .LE. iSeperableID(I)) .AND. (iSeperableID(I) .LE. iTarget2)) THEN ! This element is part of the seperable work to check
         DO J = 1,FEVertices
            LVertexFlags(Conn(J,I)) = .TRUE.
         END DO
      END IF
   END DO
END IF
IF (iWorkOption .EQ. 2) THEN ! Merge the targeted work in
   LReturnStatus = .TRUE.
   DO I = 1,NumElements
      IF (iSeperableID(I) .EQ. iTarget2) THEN
         iSeperableID(I) = iTarget1 ! We merge iSeperableID as part of iTarget1
         DO J = 1,FEVertices
            LVertexFlags(Conn(J,I)) = .TRUE.
         END DO
      END IF
   END DO
END IF
IF (iWorkOption .EQ. 3) THEN
   LReturnStatus = .TRUE. ! This will indicate that the targeted work is truly isolated from the targeted work
   ElementLoop: DO I = 1,NumElements
      IF ((iTarget1 .LE. iSeperableID(I)) .AND. (iSeperableID(I) .LE. iTarget2)) THEN ! This element is part of the targeted seperable work
         DO J = 1,FEVertices
            IF (LVertexFlags(Conn(J,I))) THEN
               LReturnStatus = .FALSE. ! This seperable work overlaps the targeted work
               EXIT ElementLoop ! We are done with the check
            END IF
         END DO
      END IF
   END DO ElementLoop
END IF
IF (iWorkOption .EQ. 4) THEN ! Merge the targeted work in
   LReturnStatus = .TRUE.
   DO I = 1,NumElements
      IF ((iTarget1 .LE. iSeperableID(I)) .AND. (iSeperableID(I) .LE. iTarget2)) THEN ! Update the flags with this piece of work
         DO J = 1,FEVertices
            LVertexFlags(Conn(J,I)) = .TRUE.
         END DO
      END IF
   END DO
END IF
IF (iWorkOption .LT. 0) THEN ! Check/Merge a given element in with the current set of work
   ! iTarget = 0 means check the adjacency and add update the flags if it is not adjacent
   ! iTarget = 1 means only check the adjacency
   ! iTarget = 2 means add the element to the current list
   LReturnStatus = .TRUE. ! Element can be added as it does not intersect with the current set
   I = -iWorkOption
#ifdef Local_Debug
   WRITE(Output_Unit,*)'Incoming Element =',I
   WRITE(Output_Unit,'(10(10L1,1X))') (LVertexFlags(Conn(J,I)),J=1,FEVertices)
#endif
   IF (iTarget1 .LE. 1) THEN ! We only do the check
      DO J = 1,FEVertices
         IF (LVertexFlags(Conn(J,I))) THEN
            LReturnStatus = .FALSE. ! This seperable work overlaps the targeted work
            EXIT ! We are done with the check
         END IF
      END DO
   END IF
   IF ((iTarget1 .NE. 1) .AND. (LReturnStatus)) THEN ! Update the flags
      DO J = 1,FEVertices
         LVertexFlags(Conn(J,I)) = .TRUE.
      END DO
   END IF
#ifdef Local_Debug
   WRITE(Output_Unit,*)'Outgoing',I
   WRITE(Output_Unit,'(10(10L1,1X))') (LVertexFlags(Conn(J,I)),J=1,FEVertices)
#endif
END IF

#ifdef Local_Debug
WRITE(Output_Unit,'(10(10L1,1X))') (LVertexFlags(I),I=1,NumVertices)
WRITE(Output_Unit,'("Matched=",L1)') LReturnStatus
#endif

END SUBROUTINE CheckAdjacency
