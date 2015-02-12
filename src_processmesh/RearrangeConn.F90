!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine takes a grid of tetrahedral elements and reorders the global vertices and connectivity for optimal threading
! In 2D, we can partition the mesh into 4 pieces (threads) where there are 4 lines of overlap between the sets of elements
!  ---- ----  
! | 1  | 2  | 
!  ---- ----
! | 3  | 4  | 
!  ---- ----
! In this approach, we can conceptually think about three levels of threaded operations
! 1) the individually assigned volumes of elements that exclude elements with connections to the interior "edges"
! 2) the interior edges between each volume block which excludes any elements that interact at corners
! 3) the center vertices which connect 4 adjacent edges
! 
! In summary, this routine bunches the data into iColors of work which there are 1,2,3,4, threads connected.
! The idea being that in each situation, there is enough pieces of work that it can be distributed over the N threads
!---------------------------------------------------------------------------------------------------------------------------------
! These arrays are reordered in this subroutine
! GlobalXYZ, Conn, ElementsWithVacuumBCs
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE RearrangeConn(Input_Part)
USE CommonBlock
IMPLICIT NONE 
!#define Local_Debug
#define Local_TestVertexReorder
! Passed in
PROTEUS_Int Input_Part

! Local scratch
PROTEUS_Int :: iMaxConnections = 0, iMaxSeperable = 0, NumSeperable = 0
PROTEUS_Int,  ALLOCATABLE :: iVertexAssignment(:)
PROTEUS_Int,  ALLOCATABLE :: iNewElementIndex(:),iNewVertexIndex(:)         ! The new element and vertex ordering to use on the local process
PROTEUS_Int,  ALLOCATABLE :: iReverseElementIndex(:),iReverseVertexIndex(:) ! The reverse map of new element and vertex ordering to use on the local process
PROTEUS_Int,  ALLOCATABLE :: iElementLocalToGlobal(:),iVertexLocalToGlobal(:)
PROTEUS_Int,  ALLOCATABLE :: iElementNumThreads(:)    ! The number of threads per element
PROTEUS_Int,  ALLOCATABLE :: iElementThreads(:,:)   ! The threads ids that connect to each element
PROTEUS_Int,  ALLOCATABLE :: iSeperableWork(:,:)    ! (NumSeperable+1,3) position (I) gives the starting element for each seperable piece of work
PROTEUS_Int,  ALLOCATABLE :: iColorGroups(:)          ! (AS_NumColors)   position (I) gives the number of seperable pieces of work in each iColor
! Memory copy stuff
PROTEUS_Int,  ALLOCATABLE :: iConn(:,:)             ! (NumElements,4)  temp copy of Conn matrix
PROTEUS_Real, ALLOCATABLE :: dGlobalXYZ(:,:)        ! (NumVertices,3)  temp copy of global xyz
! Local
PROTEUS_Int Element,Surface
PROTEUS_Int I,J,K,L
PROTEUS_Real dXX,dYY,dZZ,rXX,rYY,rZZ

IF ((NumThreads .GT. 1) .AND. (NumThreads .GT. NumElements/4)) THEN
   WRITE(Output_Unit,*)'Sorry, but NumThreads must be small relative to the number of elements'
   STOP
END IF

! We make reasonable guesses as to the work
iMaxConnections = NumElements
iMaxSeperable = NumElements 
! Local scratch
ALLOCATE(iSeperableWork(iMaxSeperable+1,3),iColorGroups(iMaxConnections),                &
         iVertexAssignment(NumVertices),                                               &
         iNewElementIndex(NumElements),iNewVertexIndex(NumVertices),                   &
         iReverseElementIndex(NumElements),iReverseVertexIndex(NumVertices),           &
         iElementLocalToGlobal(NumElements),iVertexLocalToGlobal(NumVertices*2),       &
         iElementNumThreads(NumElements),iElementThreads(FEVertices,NumElements),          &
         iConn(FEVertices,NumElements),dGlobalXYZ(NumVertices,3))

! Partition the vertices to each thread: iVertexAssignment
IF (Input_Part .EQ. 1) THEN
   CALL GetVertextoThread_METIS(Output_Unit,NumThreads,NumVertices,NumElements,FEVertices,Conn,  iVertexAssignment)
   ! Identify the adjacency of the threads and threads touched per element: iNumAdjacent,iThreadAdjacency,iElementThreads
   CALL GetAdjacencyTElement(Output_Unit,NumThreads,NumVertices,NumElements,FEVertices,Conn,iVertexAssignment, &
                              iElementNumThreads,iElementThreads)
   ! Identify the new element indexing along with the NumThreads constrained seperable work iColors
   CALL GetTElementOrder_METIS(Output_Unit,NumThreads,NumVertices,NumElements,FEVertices,iMaxConnections,iMaxSeperable, &
                                iElementNumThreads,iElementThreads,Conn,                                       &
                                iNewElementIndex,iReverseElementIndex,                                         &
                                NumSeperable,iSeperableWork(1,1),iSeperableWork(1,2),iSeperableWork(1,3),      &
                                AS_NumColors,iColorGroups,iVertexLocalToGlobal) ! iVertexLocalToGlobal is scratch
ELSE
   CALL GetVertextoThread_GREEDY(Output_Unit,NumThreads,NumVertices,GlobalXYZ,NumElements,FEVertices,Conn,  iVertexAssignment)
   ! Identify the adjacency of the threads and threads touched per element: iNumAdjacent,iThreadAdjacency,iElementThreads
   CALL GetAdjacencyTElement(Output_Unit,NumThreads,NumVertices,NumElements,FEVertices,Conn,iVertexAssignment, &
                              iElementNumThreads,iElementThreads)
   ! Identify the new element indexing along with the NumThreads constrained seperable work iColors
   CALL GetTElementOrder_GREEDY(Output_Unit,NumThreads,NumVertices,NumElements,FEVertices,iMaxConnections,iMaxSeperable, &
                                iElementNumThreads,iElementThreads,Conn,                                       &
                                iNewElementIndex,iReverseElementIndex,                                         &
                                NumSeperable,iSeperableWork(1,1),iSeperableWork(1,2),iSeperableWork(1,3),      &
                                AS_NumColors,iColorGroups,iVertexLocalToGlobal) ! iVertexLocalToGlobal is scratch
END IF

! We can now define the thread wise work on the FE assembly, matrix-matrix, and disassembly operations
ALLOCATE(DA_ThreadWiseWork(2,NumThreads),MM_ThreadWiseWork(2,NumThreads),AS_ThreadWiseWork(2,AS_NumColors,NumThreads),STAT=I)
IF (I .NE. 0) THEN
   WRITE(6,*)'Allocation failure in rearrangeConn ',AS_NumColors
   STOP
ENDIF
K = NumElements/NumThreads
L = 1
DO I = 1,NumThreads
   DA_ThreadWiseWork(1,I) = L
   DA_ThreadWiseWork(2,I) = L+K-1
   MM_ThreadWiseWork(1,I) = L
   MM_ThreadWiseWork(2,I) = L+K-1
   L = L + K
   DO J = 1,AS_NumColors
      AS_ThreadWiseWork(1,J,I) = 1  ! Rather than add an additional if test on this, we just rely upon the do loop to close immediately
      AS_ThreadWiseWork(2,J,I) = 0
   END DO
END DO
DA_ThreadWiseWork(2,NumThreads) = NumElements
MM_ThreadWiseWork(2,NumThreads) = NumElements

! Fill the assembly array
L = 1
DO I = 1,AS_NumColors
   K = 0
   DO J = L,L - 1 + iColorGroups(I)  ! iColors(I) is the number of seperable pieces in group I
      K = K + 1
      AS_ThreadWiseWork(1,I,K) = iSeperableWork(J,1)     ! Thread K group work starts at this element
      AS_ThreadWiseWork(2,I,K) = iSeperableWork(J+1,1)-1 ! Thread K group work ends   at this element
   END DO
   L = L + iColorGroups(I)
END DO

! Compute the reverse element index
DO I = 1,NumElements
   J = iNewElementIndex(I)
   iReverseElementIndex(J) = I
END DO

#ifdef Local_Debug
   WRITE(Output_Unit,'("The assembly operation is broken into the following steps of element assigned work")')
   WRITE(Output_Unit,'("Threads:",2X, 500(6X,I5,6X))') (I,I=1,NumThreads)
   DO I = 1,AS_NumColors
     WRITE(Output_Unit,'("Step ",I5,500("|",I7,":",I7,"|"))') I,(AS_ThreadWiseWork(1,I,J),AS_ThreadWiseWork(2,I,J),J=1,NumThreads)
   END DO
#endif
WRITE(Output_Unit,'("The assembly operation is broken into the following steps of elements")')
WRITE(Output_Unit,'("Threads:",2X, 500(2X,I4,2X))') (I,I=1,NumThreads)
DO I = 1,AS_NumColors
  WRITE(Output_Unit,'("Step ",I5,500("|",I6,"|"))') I,(AS_ThreadWiseWork(2,I,J)-AS_ThreadWiseWork(1,I,J)+1,J=1,NumThreads)
END DO

#ifdef Local_TestVertexReorder
   ! Identify the new global vertex ordering
   CALL GetTVertexOrder(Output_Unit,NumThreads,NumVertices,NumElements,FEVertices,AS_NumColors, &
                             Conn,iNewElementIndex,AS_ThreadWiseWork,                          &
                             iNewVertexIndex,iReverseVertexIndex) 
   DO I = 1,NumVertices
      J = iNewVertexIndex(I)     ! Vertex J moves to position I
      iVertexLocalToGlobal(I) = VertexLocalToGlobal(J)
      dGlobalXYZ(I,1) = GlobalXYZ(J,1)
      dGlobalXYZ(I,2) = GlobalXYZ(J,2)
      dGlobalXYZ(I,3) = GlobalXYZ(J,3)
   END DO
#else
   ! In real applications, the vertex data will not be reordered as we need to perform scatter-gather operations
   !  and would prefer some amount of data adjacency
   DO I = 1,NumVertices
      iVertexLocalToGlobal(I) = I
      iReverseVertexIndex(I)  = I
   END DO
#endif

! Reorder the global connectivity
DO I = 1,NumElements
   J = iNewElementIndex(I) ! Element J is moved to position I
   iElementLocalToGlobal(I) = ElementLocalToGlobal(J)
   DO K = 1,FEVertices
      L = Conn(K,J) ! The vertex assignment for the original element
      iConn(K,I) = iReverseVertexIndex(L) ! The new vertex ID for the reordered element
   END DO
END DO

! Reorder the element surface info
DO I = 1,NumVacuum
   J = ElementsWithVacuumBCs(I) ! The original element ID
   ElementsWithVacuumBCs(I) = iReverseElementIndex(J) ! The new element ID for that element
END DO

! Copy back the vertex and connectivity data
DO I = 1,NumVertices
   VertexLocalToGlobal(I) = iVertexLocalToGlobal(I)
   GlobalXYZ(I,1) = dGlobalXYZ(I,1)
   GlobalXYZ(I,2) = dGlobalXYZ(I,2)
   GlobalXYZ(I,3) = dGlobalXYZ(I,3)
END DO
DO I = 1,NumElements
   ElementLocalToGlobal(I) = iElementLocalToGlobal(I)
END DO
DO I = 1,NumElements
   DO K = 1,FEVertices
      Conn(K,I) = iConn(K,I)
   END DO
END DO

CALL ExportVTK(iVertexAssignment,iNewElementIndex,iNewVertexIndex)

! Free the local scratch
DEALLOCATE(iSeperableWork,iColorGroups,                      &
         iVertexAssignment,                                  &
         iNewElementIndex,iNewVertexIndex,                   &
         iReverseElementIndex,iReverseVertexIndex,           &
         iElementLocalToGlobal,iVertexLocalToGlobal,         &
         iElementNumThreads,iElementThreads,                 &
         iConn,dGlobalXYZ,STAT=I)
   !WRITE(Output_Unit,*)'position 10'
IF (I .NE. 0) THEN
   WRITE(6,*)'Deallocation failure in rearrangeConn ',AS_NumColors
   STOP
ENDIF

END SUBROUTINE RearrangeConn
