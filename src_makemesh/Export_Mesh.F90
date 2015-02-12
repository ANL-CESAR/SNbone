!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine exports the tetrahedral mesh
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE Export_Mesh(NumVertices,NumElements,NumVacuum,     &
                       GlobalXYZ,Conn,ElementsWithVacuumBCs,LocalSurfaceIndex)
IMPLICIT NONE 
! Parameters
PROTEUS_Int :: NumVertices
PROTEUS_Int :: NumElements
PROTEUS_Int :: NumVacuum
! Arrays to fill
PROTEUS_Real    GlobalXYZ(NumVertices,3)
PROTEUS_Int     Conn(4,NumElements)
PROTEUS_Int     ElementsWithVacuumBCs(NumVacuum)
PROTEUS_Int     LocalSurfaceIndex(NumVacuum)
! Local
PROTEUS_Int Element,Vertex,File_Unit,IOS,I
PROTEUS_Real dXX,dYY,dZZ

! We do a FE transformation here to make the geometry look like something other than a brick :-)
! This is done to highlight the fact that CFE does not care about element-angle ordering
DO I = 1,NumVertices
   dXX = (GlobalXYZ(I,1) - 0.5d0)*2.0d0
   dYY = (GlobalXYZ(I,2) - 0.5d0)*2.0d0
   dZZ = (GlobalXYZ(I,3) - 0.5d0)*2.0d0
   CALL Transform_Hex20(GlobalXYZ(I,1),GlobalXYZ(I,2),GlobalXYZ(I,3),dXX,dYY,dZZ)
END DO

200 FORMAT(10I8)
210 FORMAT(1P,6(E13.6,1X))
OPEN(UNIT=File_Unit,IOSTAT=IOS,FILE='grid_tet_mesh.ascii')
WRITE(File_Unit,'("# Version 1.0 file format for CFE SN miniapp")')        ! Generic junk about the file type
WRITE(File_Unit,200) NumVertices,NumElements,NumVacuum
DO Vertex = 1,NumVertices
   WRITE(File_Unit,210) GlobalXYZ(Vertex,1),GlobalXYZ(Vertex,2),GlobalXYZ(Vertex,3)
END DO
DO Element = 1,NumElements
   WRITE(File_Unit,200) Conn(1,Element),Conn(2,Element),Conn(3,Element),Conn(4,Element)
END DO
DO Element = 1,NumVacuum
   WRITE(File_Unit,200) ElementsWithVacuumBCs(Element),LocalSurfaceIndex(Element)
END DO
CLOSE(File_Unit)


OPEN(UNIT=File_Unit,IOSTAT=IOS,FILE='grid_tet_mesh.vtk',ACCESS='SEQUENTIAL',FORM='FORMATTED')
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

END SUBROUTINE Export_Mesh
