!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine seperates the elements into thread wise seperable pieces of work with regard to element assembly
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE GetTElementOrder_METIS(Output_Unit,NumThreads,NumVertices,NumElements,FEVertices,iMaxConnections,iMaxSeperable, &
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
PROTEUS_Int, INTENT(IN)    :: iElementNumThreads(NumElements)              ! The number of threads adjacent to each element
PROTEUS_Int, INTENT(IN)    :: iElementThreads(FEVertices,NumElements)      ! The maximum number of threads connected to any element
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
PROTEUS_Log, INTENT(INOUT) :: iScratchV(NumVertices,2)                     ! Scratch arrays
! Local
PROTEUS_Int I,J,K,L,M,N,II,JJ,KK,LL,FakeNumThreads
PROTEUS_Int iElement,iListSize,MaxEtoFind,iElementsInPass
PROTEUS_Log Matched

! Initialization
DO I = 1,NumElements
   iNewElementIndex(I) = 0
   iReverseElementIndex(I) = 0
END DO
iNumSeperable = 0
DO I = 1,iMaxSeperable
   iSeperableWork(I)       = 0
   iThreadsPerSeperable(I) = 0
   iReOrderSeperable(I)    = 0
END DO

#ifdef Local_Debug
   WRITE(Output_Unit,*)
   WRITE(Output_Unit,*)'Entered GetTElementOrder'
#endif

! iReverseElementIndex is used to identify which seperable piece of work each element is assigned to
iElement = 0 
DO I = 1,NumElements
   IF ((iReverseElementIndex(I) .EQ. 0) .AND. (iElementNumThreads(I) .EQ. 1)) THEN 
#ifdef Local_Debug
      WRITE(Output_Unit,'("Element=",I6," NumThreads =",I6," iElement=",I6)') I,iElementNumThreads(I),iElement
#endif
      ! This becomes the next batch of seperable work that we can assign to a thread
      iElement = iElement + 1
      iNumSeperable = iNumSeperable + 1
      iSeperableWork(iNumSeperable) = iSeperableWork(iNumSeperable) + 1 ! Tracks the number of elements in each seperable piece of work
      iReverseElementIndex(I)=iElement
      DO J = I+1,NumElements ! Look for any thread adjacency matches in the remaining elements
         IF ((iReverseElementIndex(J) .EQ. 0) .AND. (iElementNumThreads(J) .EQ. 1) &
             .AND. (iElementThreads(1,I) .EQ. iElementThreads(1,J)) ) THEN ! They have the same thread id
            iElement = iElement + 1
            iReverseElementIndex(J)=iElement
            iSeperableWork(iNumSeperable) = iSeperableWork(iNumSeperable) + 1
#ifdef Local_Debug
            WRITE(Output_Unit,'(" Element=",I6," threads matched threads in Element=",I6," -> ",30(I3,1X))') &
                              J,I,(iElementThreads(K,I),K=1,iElementNumThreads(I))
#endif
         END IF
      END DO
   END IF
END DO

iColorGroups(1:iMaxConnections) = 0
iNumColors = 1
iColorGroups(iNumColors) = iNumSeperable
! We need to set up iSeperable correctly
J = iSeperableWork(1)
iSeperableWork(1) = 1
DO I = 2,iNumSeperable+1
   K = iSeperableWork(I)
   iSeperableWork(I) =  iSeperableWork(I-1) + J
   J = K
END DO

#ifdef Local_Debug
   WRITE(Output_Unit,'("Found ",I6," seperable pieces of work for the FE assembly operation")') iNumSeperable
   DO I = 1,iNumSeperable
      WRITE(Output_Unit,'("Seperable piece ",I8," has ",I6," elements in it")') I,iSeperableWork(I+1)-iSeperableWork(I)
   END DO
   WRITE(Output_Unit,*)'Element-wise new element ordering'
   WRITE(Output_Unit,'(5I6,3X,5I6,3X,5I6,3X,5I6)') (iReverseElementIndex(I),I=1,NumElements)
   WRITE(Output_Unit,*)
   WRITE(Output_Unit,*)'Starting phase 2 of the iColor work in GetTElementOrder'

   DO I = 1,NumElements
      IF (iReverseElementIndex(I) .EQ. 0) WRITE(Output_Unit,77) I,iReverseElementIndex(I),(Conn(J,I),J=1,FEVertices)
      77 FORMAT('Element ',I5,' Status ',I3,' Conn-> ',100(I5,1X))
   END DO
#endif

DO 
   iElementsInPass = 0
   KK = NumElements-iElement   ! The number of elements left to assign
   MaxEtoFind = KK/NumThreads  ! The number of elements we can potentially merge together
#ifdef Local_Debug
   WRITE(Output_Unit,*)'Next pass with ',KK,' left to assign and ',MaxEtoFind,' per thread'
#endif
   FakeNumThreads = 0
   ! Construct a list of unassigned elements
   iListSize = 0
   DO I = 1,NumElements
      IF (iReverseElementIndex(I) .EQ. 0) THEN
         iListSize = iListSize + 1
         iNewElementIndex(iListSize) = I
#ifdef Local_Debug
   !WRITE(Output_Unit,*)'Element ',I,' was unassigned'
#endif
      END IF
   END DO
   IF (iListSize .EQ. 0) EXIT ! We are clearly done as there were no more elements
   ! Try to find clumps of elements > 2 to assign to each thread
   IF ((MaxEtoFind .GT. 1) .AND. (KK .GT. NumThreads*10+10)) THEN
      IF (MaxEtoFind .GT. 100) THEN
         MaxEtoFind = MaxEtoFind/3
      ELSE 
         MaxEtoFind = MaxEtoFind/2
      END IF
      iScratchV(1:NumVertices,1) = .FALSE. ! Initialize
      DO II = 1,iListSize
         I = iNewElementIndex(II)
         IF (iReverseElementIndex(I) .EQ. 0) THEN ! This element is not assigned yet
M=-I; ! This will check to see if element I touches the current work but not add
L=1;CALL CheckAdjacency(Output_Unit,NumVertices,NumElements,FEVertices,L,L,M,Matched,iReverseElementIndex,Conn,iScratchV(1,1))
   !WRITE(Output_Unit,*)'Matched 0=',I,Matched
            IF (Matched) THEN ! Element I is not connected to whatever is already in the list
#ifdef Local_Debug
   WRITE(Output_Unit,'("Building a clump from element ",I5," with max size of ",I5)') I,MaxEtoFind
#endif
! This will add the element to the alternate check list
               iScratchV(1:NumVertices,2) = .FALSE. ! Initialize
L=2;CALL CheckAdjacency(Output_Unit,NumVertices,NumElements,FEVertices,L,L,M,Matched,iReverseElementIndex,Conn,iScratchV(1,2))
               FakeNumThreads = FakeNumThreads + 1
               iReverseElementIndex(I) = -FakeNumThreads
               iElementsInPass = 0;
               ClumpSearch: DO
                  K = iElementsInPass
                  DO JJ = 1,iListSize ! Check for adjacencies with element I
                     J = iNewElementIndex(JJ)
                     IF (iReverseElementIndex(J) .EQ. 0) THEN ! It might not be zero anymore
M=-J
L=1;CALL CheckAdjacency(Output_Unit,NumVertices,NumElements,FEVertices,L,L,M,Matched,iReverseElementIndex,Conn,iScratchV(1,2)) ! Check only
   !WRITE(Output_Unit,*)'Matched 1=',J,Matched
                        IF (.NOT. Matched) THEN ! Element J is connected to other elements in this clump
L=1;CALL CheckAdjacency(Output_Unit,NumVertices,NumElements,FEVertices,L,L,M,Matched,iReverseElementIndex,Conn,iScratchV(1,1)) ! Check only
   !WRITE(Output_Unit,*)'Matched 2=',J,Matched
                           IF (Matched) THEN ! Element J does not touch any other elements in this color
#ifdef Local_Debug
   !WRITE(Output_Unit,'("Element ",I5," is connected to ",I5," but not other elements of this color")') J,I
#endif
                              iElementsInPass = iElementsInPass + 1 ! The total number of elements that are adjacent to "I"
                              iReverseElementIndex(J) = -FakeNumThreads
L=2;CALL CheckAdjacency(Output_Unit,NumVertices,NumElements,FEVertices,L,L,M,Matched,iReverseElementIndex,Conn,iScratchV(1,2)) ! Add to list
                           END IF
                        END IF
                     END IF
                     ! Exit the search
                     IF (iElementsInPass .EQ. MaxEtoFind) EXIT ClumpSearch
                  END DO
#ifdef Local_Debug
   !WRITE(Output_Unit,*)'Element-wise new element ordering after pass'
   !WRITE(Output_Unit,'(5I6,3X,5I6,3X,5I6,3X,5I6)') (iReverseElementIndex(I),I=1,NumElements)
   IF (K .NE. iElementsInPass) THEN
      WRITE(Output_Unit,*)'Unassigned element status ',iElementsInPass,K
      WRITE(Output_Unit,'(5I6,3X,5I6,3X,5I6,3X,5I6)') (iReverseElementIndex(iNewElementIndex(JJ)),JJ=1,iListSize)
   END IF
#endif
                  IF (K .EQ. iElementsInPass) EXIT ClumpSearch ! We did nothing in the last pass so don't bother doing another one
               END DO ClumpSearch
               IF ((MaxEtoFind .GT. 4) .AND. (iElementsInPass .LE. MaxEtoFind/2)) THEN ! This element can likely be left until the last pass as they are orphaned
#ifdef Local_Debug
   WRITE(Output_Unit,'("There were an insufficient number of elements ",I5," associated with ",I5)') iElementsInPass,I
#endif
                  iElementsInPass = 0
                  DO JJ = 1,iListSize ! Check for adjacencies with element I
                     J = iNewElementIndex(JJ)
                     IF (iReverseElementIndex(J) .EQ. -FakeNumThreads) iReverseElementIndex(J) = 0 ! Set it back to zero
                  END DO
                  FakeNumThreads = FakeNumThreads - 1
               ELSE
                  DO JJ = 1,iListSize
                     J = iNewElementIndex(JJ)
                     IF (iReverseElementIndex(J) .EQ. -FakeNumThreads) THEN
#ifdef Local_Debug
   !WRITE(Output_Unit,*)'Adding elements in this clump'
#endif
M=-J ! Add the element to the set 1 list to indicate they are now clumped together for a given processor
L=2;CALL CheckAdjacency(Output_Unit,NumVertices,NumElements,FEVertices,L,L,M,Matched,iReverseElementIndex,Conn,iScratchV(1,1)) ! Add to list
                     END IF
                  END DO
               END IF
            END IF ! .NOT. Matched
         END IF ! iReverseElementIndex(I) .EQ. 0
         IF (FakeNumThreads .EQ. NumThreads) EXIT ! We have found enough work to do for each thread
      END DO ! II
   END IF ! (MaxEtoFind .GT. 2)
   IF (FakeNumThreads .EQ. 0) THEN ! No work was found to assign to the current color so just assign everything that is left to thread 0
#ifdef Local_Debug
      WRITE(Output_Unit,*)'Dumping all elements into a single thread'
#endif
      DO II = 1,iListSize
         I = iNewElementIndex(II)
         IF (iReverseElementIndex(I) .EQ. 0) iReverseElementIndex(I) = -1 ! Indicate that this element needs to be assigned a seperable status
      END DO
   END IF
   iNumColors = iNumColors + 1 ! Increment the number of colors
#ifdef Local_Debug
   WRITE(Output_Unit,'("iNumColors = ",I6," iNumSeperable=",I8," iElement=",I6)') iNumColors,iNumSeperable,iElement
   WRITE(Output_Unit,*)'Element ordering before assigning to threads'
   WRITE(Output_Unit,'(5I6,3X,5I6,3X,5I6,3X,5I6)') (iReverseElementIndex(I),I=1,NumElements)
#endif
   DO KK = 1,NumThreads
      K = 0
      DO II = 1,iListSize
         I = iNewElementIndex(II)
         IF (iReverseElementIndex(I) .EQ. -KK) THEN
            K = K + 1
            iElement = iElement + 1
            iReverseElementIndex(I) = iElement
         END IF
      END DO
      IF (K .GT. 0) THEN ! Elements were modified
         iNumSeperable = iNumSeperable + 1
         iSeperableWork(iNumSeperable+1) = iElement+1
         iColorGroups(iNumColors) = iColorGroups(iNumColors) + 1 ! Number of seperable pieces in this color
      END IF
#ifdef Local_Debug
   WRITE(Output_Unit,'("Thread = ",I6," NumElements=",I6," iElement=",I6," iNumSeperable=",I8)') II,K,iElement,iNumSeperable
#endif
   END DO
END DO
iSeperableWork(iNumSeperable+1) = NumElements + 1

#ifdef Local_Debug
   WRITE(Output_Unit,'(5I6,3X,5I6,3X,5I6,3X,5I6)') (iReverseElementIndex(I),I=1,NumElements)
#endif


! Reverse the element index array
DO I = 1,NumElements
   J = iReverseElementIndex(I)
   iNewElementIndex(J) = I
END DO

#ifdef Local_Debug
   WRITE(Output_Unit,'("Found ",I6," iColors of work for the FE assembly operation ")') iNumColors
   DO I = 1,iNumColors
      WRITE(Output_Unit,'("iColor ",I6," has ",I8," seperable pieces of work")') I,iColorGroups(I)
   END DO
   WRITE(Output_Unit,*)'Current seperable work'
   DO I = 1,iNumSeperable
      WRITE(Output_Unit,'("Seperable piece ",I6," has elements ",I6,":",I6," in it")') I,iSeperableWork(I),iSeperableWork(I+1)-1
   END DO
   WRITE(Output_Unit,*)'New element index'
   WRITE(Output_Unit,'(2(10I6,1X))') (iNewElementIndex(I),I=1,NumElements)
   WRITE(Output_Unit,*)'Reverse element index'
   WRITE(Output_Unit,'(2(10I6,1X))') (iReverseElementIndex(I),I=1,NumElements)
#endif

END SUBROUTINE GetTElementOrder_METIS
