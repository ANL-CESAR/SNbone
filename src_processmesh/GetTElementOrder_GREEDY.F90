!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine seperates the elements into thread wise seperable pieces of work with regard to element assembly
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE GetTElementOrder_GREEDY(Output_Unit,NumThreads,NumVertices,NumElements,FEVertices,iMaxConnections,iMaxSeperable, &
                                iElementNumThreads,iElementThreads,Conn,                                             &
                                iNewElementIndex,iReverseElementIndex,                                               &
                                iNumSeperable,iSeperableWork,iThreadsPerSeperable,iReOrderSeperable,                 &
                                iNumColors,iColorGroups,iScratchV)
IMPLICIT NONE 
!#define Local_Debug
! Input
PROTEUS_Int, INTENT(IN)    :: Output_Unit
PROTEUS_Int, INTENT(IN)    :: NumThreads,NumVertices,NumElements,FEVertices
PROTEUS_Int, INTENT(IN)    :: iMaxConnections                              ! The maximum number of connections for a given thread (should not be more than 27)
PROTEUS_Int, INTENT(IN)    :: iMaxSeperable                                ! NumThreads volumes, ~NumThreads*2 surfaces, ~NumThreads corners; more edges -> less corners
PROTEUS_Int, INTENT(IN)    :: iElementNumThreads(NumElements)              ! The number of threads adjacent to each
PROTEUS_Int, INTENT(IN)    :: iElementThreads(FEVertices,NumElements)           ! The maximum number of threads connected to any element
PROTEUS_Int, INTENT(IN)    :: Conn(FEVertices,NumElements)
! Arrays to fill
PROTEUS_Int, INTENT(INOUT) :: iNewElementIndex(NumElements)                ! J = iNewElementIndex(I)     ! Element J is moved to position I
PROTEUS_Int, INTENT(INOUT) :: iReverseElementIndex(NumElements)            ! iReverseElementIndex(J)=I   ! Element J goes to position I
PROTEUS_Int, INTENT(INOUT) :: iNumSeperable                                ! The number of seperable pieces of work (i.e. number of threads worth of work)
PROTEUS_Int, INTENT(INOUT) :: iSeperableWork(iMaxSeperable+1)              ! The starting element position for each piece of seperable work
PROTEUS_Int, INTENT(INOUT) :: iThreadsPerSeperable(iMaxSeperable)          ! Scratch
PROTEUS_Int, INTENT(INOUT) :: iReOrderSeperable(iMaxSeperable)             ! Scratch
PROTEUS_Int, INTENT(INOUT) :: iNumColors                                   ! The number of iColors of seperable work (i.e. number of outer loops)
PROTEUS_Int, INTENT(INOUT) :: iColorGroups(iMaxConnections)                ! The number of seperable pieces by connection (i.e. # ones, # twos, # threes, ...)
PROTEUS_Log, INTENT(INOUT) :: iScratchV(NumVertices)                       ! A scratch array
PROTEUS_Int I,J,K,II,JJ,KK
PROTEUS_Int iElement
PROTEUS_Log Matched

! Initialization
DO I = 1,NumElements
   iNewElementIndex(I) = 0
   iReverseElementIndex(I) = 0
END DO
iColorGroups(1:iMaxConnections) = 0

! In this scheme we just add as many elements as we can to each color where each element must be completely isolated
iNumColors = 0
iNumSeperable = 0
iSeperableWork(1) = 1
iElement = 0
DO 
   JJ = 0
   iScratchV(1:NumVertices) = .FALSE. ! Initialize
   DO I = 1,NumElements
      IF (iReverseElementIndex(I) .EQ. 0) THEN ! Unassigned at the moment
         J = -I
         K = 0
         CALL CheckAdjacency(Output_Unit,NumVertices,NumElements,FEVertices,K,K,J,Matched,iReverseElementIndex,Conn,iScratchV)
         IF (Matched) THEN ! This element did not touch the current group
            JJ = JJ + 1 ! Track the number of elements in the current "color"
            iReverseElementIndex(I) = -1 ! Indicate that this element needs to be assigned a seperable status
         END IF
      END IF
   END DO
   IF (JJ .EQ. 0) EXIT ! We are clearly done as there were no more elements
   ! Break the work into seperable pieces=NumThreads if possible
   iNumColors = iNumColors + 1 ! Increment the number of colors
   KK = JJ/NumThreads          ! The number of elements to assign per thread
   IF (KK .EQ. 0) KK = 1       ! JJ<NumThreads and we will run out of elements before threads
#ifdef Local_Debug
   WRITE(Output_Unit,'("iNumColors = ",I6," iNumSeperable=",I6," iElement=",I6," elements/thread=",I6)') iNumColors,iNumSeperable,iElement,KK
   WRITE(Output_Unit,*)'Element-wise seperable work identification'
   WRITE(Output_Unit,'(2(10I6,1X))') (iReverseElementIndex(I),I=1,NumElements)
#endif
   DO II = 1,NumThreads
      K = 0
      IF (II .EQ. NumThreads) KK = JJ ! Grab everything left for the last thread
      DO I = 1,NumElements
         IF (iReverseElementIndex(I) .EQ. -1) THEN
            iElement = iElement + 1
            iReverseElementIndex(I) = iElement
            K = K + 1
            IF (K .EQ. KK) EXIT
         END IF
      END DO
      IF (K .GT. 0) THEN ! Elements were modified
         iNumSeperable = iNumSeperable + 1
         iSeperableWork(iNumSeperable+1) = iElement+1
         iColorGroups(iNumColors) = iColorGroups(iNumColors) + 1 ! Number of seperable pieces in this color
      END IF
#ifdef Local_Debug
   WRITE(Output_Unit,'("Thread = ",I6," NumElements=",I6," iElement=",I6," iNumSeperable=",I3)') II,K,iElement,iNumSeperable
#endif
   END DO
END DO
iSeperableWork(iNumSeperable+1) = NumElements + 1
! Reverse the element index array
DO I = 1,NumElements
   J = iReverseElementIndex(I)
   iNewElementIndex(J) = I
END DO

END SUBROUTINE GetTElementOrder_GREEDY
