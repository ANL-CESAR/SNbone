!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine imports a tetrahedral mesh
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE Import_Mesh()
USE CommonBlock
IMPLICIT NONE 
! Local
PROTEUS_Int Element,Vertex,File_Unit,IOS,I

File_Unit = 777
OPEN(UNIT=File_Unit,IOSTAT=IOS,STATUS='OLD',FILE='inputmesh.ascii')
IF (IOS .NE. 0) THEN
   WRITE(Output_Unit,*)'Could not open inputmesh.ascii'
   CALL Abort
END IF
READ(File_Unit,*,IOSTAT=IOS)        ! Generic junk about the file type
IF (IOS .NE. 0) THEN
   WRITE(Output_Unit,*)'Read failure of ascii string'
   CALL Abort
END IF
READ(File_Unit,*,IOSTAT=IOS) NumVertices,NumElements,NumVacuum
IF (IOS .NE. 0) THEN
   WRITE(Output_Unit,*)'Read failure of control data card'
   CALL Abort
END IF
CALL CommonBlock_Allocate()

DO Vertex = 1,NumVertices
   READ(File_Unit,*,IOSTAT=IOS) GlobalXYZ(Vertex,1),GlobalXYZ(Vertex,2),GlobalXYZ(Vertex,3)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of vertex data line # ',Vertex
      CALL Abort
   END IF
END DO
DO Element = 1,NumElements
   READ(File_Unit,*,IOSTAT=IOS) Conn(1,Element),Conn(2,Element),Conn(3,Element),Conn(4,Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of element data line # ',Element
      CALL Abort
   END IF
END DO
DO Element = 1,NumVacuum
   READ(File_Unit,*,IOSTAT=IOS) ElementsWithVacuumBCs(Element),LocalSurfaceIndex(Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of surface data line # ',Element
      CALL Abort
   END IF
END DO

CLOSE(File_Unit)

OPEN(UNIT=File_Unit,IOSTAT=IOS,FILE='inputmesh.vtk',ACCESS='SEQUENTIAL',FORM='FORMATTED')
WRITE(File_Unit,'("# vtk DataFile Version 3.0")')                          ! Generic junk about the file type
WRITE(File_Unit,'("ANL NE devision neutronics export routine 08-2006")')   ! A title description of the type of output
WRITE(File_Unit,'("ASCII")')                                               ! This will be an ascii export option
WRITE(File_Unit,'("DATASET UNSTRUCTURED_GRID")')                           ! Indicates that we are defining a unstructured grid
WRITE(File_Unit,'("POINTS ",I12," float")')  NumVertices                   ! The number of mesh points which are of type float
DO I = 1,NumVertices
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

WRITE(File_Unit,'("CELL_DATA ",I12)') NumElements
WRITE(File_Unit,'("FIELD MeshInformation 1")')
WRITE(File_Unit,'("Element_ID  1 ",I12," float")') NumElements
WRITE(File_Unit,888) (I*1.0d0,I=1,NumElements)
888 FORMAT(10(1PE12.5,1X))
CLOSE(File_Unit)

END SUBROUTINE Import_Mesh
