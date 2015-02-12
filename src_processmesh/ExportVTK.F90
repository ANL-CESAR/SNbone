!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine exports the rearranged mesh to a VTK file for visual confirmation
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE ExportVTK(iVertexAssignment,iNewElementIndex,iNewVertexIndex)
USE CommonBlock
IMPLICIT NONE 
! Passed in
PROTEUS_Int :: iNewElementIndex(NumElements)
PROTEUS_Int :: iVertexAssignment(NumVertices),iNewVertexIndex(NumVertices)
! Local
PROTEUS_Int :: File_Unit=10
PROTEUS_Real, ALLOCATABLE :: dPlotValues(:,:)    ! (NumElements,10)  1=DA assignment 2=MM assignment 3=AS assignment 4=E ordering
PROTEUS_Int I,J,K,L,Element

! Export to a VTK file
ALLOCATE(dPlotValues(NumElements,10))
DO I = 1,NumThreads
   DO K = DA_ThreadWiseWork(1,I),DA_ThreadWiseWork(2,I)
      dPlotValues(K,1) = I
   END DO
   DO K = MM_ThreadWiseWork(1,I),MM_ThreadWiseWork(2,I)
      dPlotValues(K,2) = I
      dPlotValues(iNewElementIndex(K),9) = I
   END DO
END DO
L = 0
DO J = 1,AS_NumColors
   DO I = 1,NumThreads
      L = L + 1
      DO K = AS_ThreadWiseWork(1,J,I),AS_ThreadWiseWork(2,J,I)
         dPlotValues(K,3) = I ! Thread
         dPlotValues(K,6) = J ! iColor
         dPlotValues(K,7) = L ! seperable
      END DO
   END DO
END DO
DO I = 1,NumElements
   dPlotValues(I,4) = iNewElementIndex(I) ! Because we have now reordered the list, this should be the old ordering id
   dPlotValues(I,5) = I
   IF (iNewElementIndex(I) .LE. 5) THEN
      dPlotValues(I,8) = iNewElementIndex(I) - ((iNewElementIndex(I)/6))*5.0d0
   ELSE
      dPlotValues(I,8) = iNewElementIndex(I) - ((iNewElementIndex(I)/5)+1.0d0)*5.0d0 + 5.0d0
      IF (dPlotValues(I,8) .EQ. 0.0d0) dPlotValues(I,8) = 5.0d0
   END IF
END DO

! Open the input file for writing
OPEN(UNIT=File_Unit,IOSTAT=I,FILE='pmesh.vtk',ACCESS='SEQUENTIAL',FORM='FORMATTED')
WRITE(File_Unit,'("# vtk DataFile Version 3.0")')                          ! Generic junk about the file type
WRITE(File_Unit,'("ANL NE devision neutronics export routine 08-2006")')   ! A title description of the type of output
WRITE(File_Unit,'("ASCII")')                                               ! This will be an ascii export option
WRITE(File_Unit,'("DATASET UNSTRUCTURED_GRID")')                           ! Indicates that we are defining a unstructured grid
WRITE(File_Unit,'("POINTS ",I12," float")')  NumVertices                   ! The number of mesh points which are of type float
DO I = 1,NumVertices
   !RadiusXY  = DSQRT(0.3d0+GlobalXYZ(I,1)*GlobalXYZ(I,1)+GlobalXYZ(I,2)*GlobalXYZ(I,2))+1.0D-15
   !RadiusXYZ = DSQRT(0.4d0+GlobalXYZ(I,1)*GlobalXYZ(I,1)+GlobalXYZ(I,2)*GlobalXYZ(I,2)+GlobalXYZ(I,3)*GlobalXYZ(I,3))+1.0D-15
   !Phi = ACOS((GlobalXYZ(I,1)+0.1d0)/RadiusXY)
   !Theta = ACOS((GlobalXYZ(I,3)*0.8d0+0.2d0)/RadiusXYZ)
   !IF (GlobalXYZ(I,2) .LT. 0.1d0) RadiusXYZ = RadiusXYZ - 0.1d0
   !IF (GlobalXYZ(I,3) .GT. 0.5d0) RadiusXYZ = RadiusXYZ + 0.2d0
   !WRITE(File_Unit,'(3(1PE12.5,1X))') RadiusXYZ*((COS(Phi))**0.4d0)*SIN(Theta),  &
   !                                   RadiusXYZ*((SIN(Phi*1.5d0))**0.6d0)*((COS(Phi*0.5d0))**0.1d0)*SIN(Theta),  &
   !                                   RadiusXYZ*(COS(Theta))**0.4d0
   ! We use a FE transformation
   !dXX = (GlobalXYZ(I,1) - 0.5d0)*2.0d0;dYY = (GlobalXYZ(I,2) - 0.5d0)*2.0d0;dZZ = (GlobalXYZ(I,3) - 0.5d0)*2.0d0;
   !CALL f_1_Element_Quadratic_Tet(rXX,rYY,rZZ,dXX,dYY,dZZ)
   !WRITE(File_Unit,'(3(1PE12.5,1X))') rXX,rYY,rZZ
   WRITE(File_Unit,'(3(1PE12.5,1X))') GlobalXYZ(I,1),GlobalXYZ(I,2),GlobalXYZ(I,3)
END DO

WRITE(File_Unit,*) ! Blank line

WRITE(File_Unit,'("CELLS ",I12,1X,I12)') NumElements,NumElements*5

DO Element = 1,NumElements
   WRITE(File_Unit,'(30I12)') 4,Conn(1,Element)-1,Conn(2,Element)-1,Conn(3,Element)-1,Conn(4,Element)-1
END DO

WRITE(File_Unit,*) ! Blank line

WRITE(File_Unit,'("CELL_TYPES ",I12)') NumElements

WRITE(File_Unit,'(40I3)') (10, Element=1,NumElements)

WRITE(File_Unit,*) ! Blank line

WRITE(File_Unit,'("POINT_DATA ",I12)') NumVertices
WRITE(File_Unit,'("FIELD Vertex_data 3 ")')
WRITE(File_Unit,'("ThreadAssignment 1 ",I12," float")') NumVertices
WRITE(File_Unit,889) (iVertexAssignment(iNewVertexIndex(I))+0.0d0,I=1,NumVertices)
WRITE(File_Unit,'("VertexID 1 ",I12," float")') NumVertices
WRITE(File_Unit,889) (I+0.0d0,I=1,NumVertices)
WRITE(File_Unit,'("OldVertexID 1 ",I12," float")') NumVertices
WRITE(File_Unit,889) (iNewVertexIndex(I)+0.0d0,I=1,NumVertices)
889 FORMAT(10(1PE12.5,1X))

WRITE(File_Unit,*) ! Blank line

WRITE(File_Unit,'("CELL_DATA ",I12)') NumElements
WRITE(File_Unit,'("FIELD MeshInformation 9")')
WRITE(File_Unit,'("DA_ThreadID 1 ",I12," float")') NumElements
WRITE(File_Unit,888) dPlotValues(1:NumElements,1)
WRITE(File_Unit,'("MM_ThreadID 1 ",I12," float")') NumElements
WRITE(File_Unit,888) dPlotValues(1:NumElements,2)
WRITE(File_Unit,'("AS_ThreadID 1 ",I12," float")') NumElements
WRITE(File_Unit,888) dPlotValues(1:NumElements,3)
WRITE(File_Unit,'("OldElementID 1 ",I12," float")') NumElements
WRITE(File_Unit,888) dPlotValues(1:NumElements,4)
WRITE(File_Unit,'("NewElementID 1 ",I12," float")') NumElements
WRITE(File_Unit,888) dPlotValues(1:NumElements,5)
WRITE(File_Unit,'("ASiColor 1 ",I12," float")') NumElements
WRITE(File_Unit,888) dPlotValues(1:NumElements,6)
WRITE(File_Unit,'("ASseperable 1 ",I12," float")') NumElements
WRITE(File_Unit,888) dPlotValues(1:NumElements,7)
WRITE(File_Unit,'("BaseTetID 1 ",I12," float")') NumElements
WRITE(File_Unit,888) dPlotValues(1:NumElements,8)
WRITE(File_Unit,'("AltEOrderThreadID 1 ",I12," float")') NumElements
WRITE(File_Unit,888) dPlotValues(1:NumElements,9)
888 FORMAT(10(1PE12.5,1X))

CLOSE(File_Unit)

! Free the local scratch
DEALLOCATE(dPlotValues,STAT=I)

END SUBROUTINE ExportVTK
