! This module is used to emulate a common block behavior in the mini-app
! We would normally do this as a data structure save rather than individual arrays, 
!  but I chose to do it this way in the mini-app to keep it as simple as possible
! If you recieve an allocate error in the mini-app, ummm...try a smaller problem?
MODULE CommonBlock
IMPLICIT NONE
#include "PROTEUS_Preprocess.h"

! Output controls
PROTEUS_Int,PARAMETER :: Output_Unit = 6   ! Output unit

! Problem size
PROTEUS_Int :: NumElements          ! Total number of tetrahedral elements (5 per structured grid)
PROTEUS_Int :: NumVertices          ! Total number of vertices in structured grid
PROTEUS_Int :: NumVacuum            ! Number of tetrahedral surfaces (triangles) with vacuum boundary (2 tri/surf and 2 global surfs)

! Mesh related data structure
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: GlobalXYZ             ! The global element vertex positions
PROTEUS_Int, DIMENSION(:,:),ALLOCATABLE,SAVE :: Conn                  ! Connectivity matrix (lists 4 global vertices for 1st element, 4 global vertices for 2nd element, etc.)
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: ElementsWithVacuumBCs ! Elements which have a vacuum boundary condition (listed once for each vacuum surface)
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: LocalSurfaceIndex     ! Identifies the reference surface with the vacuum bc (1,2,3 or 4)

CONTAINS

SUBROUTINE CommonBlock_Allocate()

! Mesh related info
ALLOCATE( GlobalXYZ(NumVertices,3) )
ALLOCATE( Conn(4,NumElements)   )  ! Connectivity matrix (lists 4 global vertices for 1st element, 4 global vertices for 2nd element, etc.)
ALLOCATE( ElementsWithVacuumBCs(NumVacuum))  ! Elements which have a vacuum boundary condition (listed once for each vacuum surface)
ALLOCATE( LocalSurfaceIndex(NumVacuum)    )  ! Identifies the reference surface with the vacuum bc (1,2,3 or 4)

END SUBROUTINE CommonBlock_Allocate

SUBROUTINE CommonBlock_Deallocate()

! Deallocate the local memory
DEALLOCATE(GlobalXYZ,Conn,ElementsWithVacuumBCs,LocalSurfaceIndex)

END SUBROUTINE CommonBlock_Deallocate

END MODULE CommonBlock
