!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine assigns each global vertex to a thread by using METIS
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE GetVertextoThread_METIS(Output_Unit,NumThreads,NumVertices,NumElements,FEVertices,Conn,iVertexAssignment)
IMPLICIT NONE 
!#define Local_Debug
! Input
PROTEUS_Int, INTENT(IN)    :: Output_Unit
PROTEUS_Int, INTENT(IN)    :: NumThreads,NumVertices
PROTEUS_Int, INTENT(IN)    :: NumElements,FEVertices
PROTEUS_Int, INTENT(IN)    :: Conn(FEVertices,NumElements)
! Arrays to fill
PROTEUS_Int, INTENT(INOUT) :: iVertexAssignment(NumVertices)      ! The thread assignment by vertex

! Local
PROTEUS_Int  I,J,K,IOS
PROTEUS_Int  Element,rVertex,cVertex,RowDOF,ColDOF
PROTEUS_Int  ImaxLength,IcurLength
PROTEUS_Int  edge_weights,wgtflag,numflag,options(5),edgecut
PROTEUS_Int  IIsize,IImax,IIcVer
PROTEUS_Int, POINTER :: RowLoc(:),ColNum(:),tempColNum(:),Vertex_Weights(:)
PROTEUS_Int, POINTER :: iNumCols(:),ColNumbers(:,:),tempColNumbers(:,:)

! Memory allocation
IIsize = 10 !5000 ! The number of global vertices we work with at a time (if we have 300,000 vertices, then we will do 60 passes)
IF (IIsize .GT. NumVertices) IIsize = NumVertices
IImax = 20 ! Initial assumption about connections per vertex
IcurLength = 0
ImaxLength = NumVertices * IImax ! A guess at the current maximum storage
ALLOCATE(ColNumbers(IImax,IIsize),iNumCols(IIsize),                              &
         RowLoc(NumVertices+1),ColNum(ImaxLength),Vertex_Weights(NumVertices),   &
         STAT=IOS)

! Initialization
DO I = 1,NumVertices
   iVertexAssignment(I) = 1
   Vertex_Weights(I)    = 1
END DO

! Fill the vertex adjacency
IIcVer = 1
DO 
   DO I = 1,IIsize
      iNumCols(I) = 0
   END DO
   DO Element = 1,NumElements
      DO rVertex = 1,FEVertices
         RowDOF = Conn(rVertex,Element)
         IF ((IIcVer .LE. RowDOF) .AND. (RowDOF .LE. IIcVer+IIsize-1)) THEN ! This global vertex is within the current storage range
            DO cVertex = 1,FEVertices
               IF (cVertex .EQ. rVertex) CYCLE
               ColDOF = Conn(cVertex,Element)
               K = RowDOF - IIcVer + 1
               ! See if this column is already in the list
               J = 0
               DO I = 1,iNumCols(K)
                  IF (ColNumbers(I,K) .EQ. ColDOF) THEN
                     J = 1
                     EXIT
                  END IF
               END DO
               IF (J .EQ. 0) THEN ! Column was not present
                  IF (iNumCols(K) .EQ. IImax) THEN ! Resize the storage array
                     ALLOCATE(tempColNumbers(IImax*2,IIsize),STAT=IOS)
                     tempColNumbers(1:IImax,1:IIsize) = ColNumbers(1:IImax,1:IIsize)
                     DEALLOCATE(ColNumbers,STAT=IOS)
                     ALLOCATE(ColNumbers(IImax*2,IIsize),STAT=IOS)
                     ColNumbers(1:IImax,1:IIsize) = tempColNumbers(1:IImax,1:IIsize)
                     DEALLOCATE(tempColNumbers,STAT=IOS)
                     IImax = IImax * 2
                  END IF ! Resize
                  iNumCols(K) = iNumCols(K) + 1
                  ColNumbers(iNumCols(K),K) = ColDOF
               END IF
            END DO ! cVertex
         END IF ! valid global vertex for current search
      END DO ! row vertex
   END DO ! Element

#ifdef Local_Debug
   WRITE(Output_Unit,*)'Partial for vertices ',IIcVer,' : ',IIcVer+IIsize-1
   DO I = 1,IIsize
      WRITE(Output_Unit,7722) IIcVer+I-1,iNumCols(I),(ColNumbers(J,I),J=1,iNumCols(I))
      7722 FORMAT('Vertex ',I5,' is connected to vertices -> ',100(I4,1X))
   END DO
#endif

   ! Store the data
   J = 0
   DO I = 1,IIsize
      J = J + iNumCols(I)
   END DO
   IF (IcurLength+J .GT. ImaxLength) THEN ! Resize the storage
      ALLOCATE(tempColNum(ImaxLength*2),STAT=IOS)
      tempColNum(1:ImaxLength) = ColNum(1:ImaxLength)
      DEALLOCATE(ColNum,STAT=IOS)
      ALLOCATE(ColNum(ImaxLength*2),STAT=IOS)
      ColNum(1:ImaxLength) = tempColNum(1:ImaxLength)
      DEALLOCATE(tempColNum,STAT=IOS)
      ImaxLength = ImaxLength * 2
   END IF
   ! Store the data
   DO I = 1,IIsize
      K = IIcVer+I-1
      IF (K .GT. NumVertices) EXIT
      RowLoc(K) = IcurLength+1
      DO J = 1,iNumCols(I)
        ColNum(IcurLength+J) = ColNumbers(J,I)
      END DO
      IcurLength = IcurLength + iNumCols(I)
   END DO
   ! Increment the search range
   IIcVer = IIcVer + IIsize
   IF (IIcVer .GT. NumVertices) EXIT
END DO
! Finalize the storage
RowLoc(NumVertices+1) = IcurLength+1
#ifdef Local_Debug
   WRITE(Output_Unit,*)'Targeting ',NumThreads
   DO Element = 1,NumElements
      WRITE(Output_Unit,7744) Element,(Conn(J,Element),J=1,FEVertices)
      7744 FORMAT('Element ',I5,' Conn-> ',100(I4,1X))
   END DO
   DO I = 1,NumVertices
      WRITE(Output_Unit,7733) I,RowLoc(I+1)-RowLoc(I),(ColNum(J),J=RowLoc(I),RowLoc(I+1)-1)
      7733 FORMAT('Row ',I5,' nonzeros=',I5,' columns=',1000(I4,1X))
   END DO
#endif
! Put in c++ ordering
!DO I = 1,NumVertices+1
!   RowLoc(I) = RowLoc(I) - 1
!END DO
!DO I = 1,IcurLength
!   ColNum(I) = ColNum(I) - 1
!END DO

#ifdef USEMETIS
IF (NumThreads .GT. 1) THEN
   ! Use METIS to define the mesh partitioning
   edge_weights   = 0
   wgtflag        = 0
   numflag        = 1
   options(1)     = 0
   options(2)     = 0
   options(3)     = 0
   options(4)     = 0
   options(5)     = 0
   edgecut        = 0
   CALL METIS_PartGraphRecursive(NumVertices,RowLoc,ColNum,Vertex_Weights,edge_weights,&
                                 wgtflag,numflag,NumThreads,options,edgecut,iVertexAssignment)
END IF
#else
STOP
#endif

DO K = 1,NumVertices
#ifdef Local_Debug
   WRITE(Output_Unit,'("Vertex ",I5," is assigned to thread ",I5)') K,iVertexAssignment(K)
#endif
END DO

DEALLOCATE(ColNumbers,iNumCols,  &
           RowLoc,ColNum,Vertex_Weights,STAT=IOS)

END SUBROUTINE GetVertextoThread_METIS
