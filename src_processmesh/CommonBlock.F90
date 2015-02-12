MODULE CommonBlock
IMPLICIT NONE
#include "PROTEUS_Preprocess.h"

! Output controls
PROTEUS_Int,PARAMETER :: Output_Unit = 6   ! Output unit

! Problem size
PROTEUS_Int :: NumElements          ! Total number of tetrahedral elements (5 per structured grid)
PROTEUS_Int :: NumVertices          ! Total number of vertices in structured grid
PROTEUS_Int :: NumVacuum            ! Number of tetrahedral surfaces (triangles) with vacuum boundary (2 tri/surf and 2 global surfs)
PROTEUS_Int :: NumUnitNormals       ! Number of vacuum boundary unit normals
! Threaded related variables
PROTEUS_Int :: NumThreads           ! Total number of threads
PROTEUS_Int :: AS_NumColors = 0     ! The number of non-overlapping (in memory) pieces of work that the FE assembly operation can be broken into

! Mesh related data 
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: GlobalXYZ             ! The global element vertex positions
PROTEUS_Int, DIMENSION(:,:),ALLOCATABLE,SAVE :: Conn                  ! Connectivity matrix (lists 4 global vertices for 1st element, 4 global vertices for 2nd element, etc.)
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: ElementsWithVacuumBCs ! Elements which have a vacuum boundary condition (listed once for each vacuum surface)
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: LocalSurfaceIndex     ! Identifies the reference surface with the vacuum bc (1,2,3 or 4)
! Element-wise "cross section" values which change with respect to energy
PROTEUS_Real,DIMENSION(:), ALLOCATABLE,SAVE :: ConstTau   ! The element-wise stabilization factor
PROTEUS_Real,DIMENSION(:), ALLOCATABLE,SAVE :: ConstF     ! The coefficient for the F_element matrix
PROTEUS_Real,DIMENSION(:), ALLOCATABLE,SAVE :: ConstU     ! The coefficient for the U_element matrix
PROTEUS_Real,DIMENSION(:), ALLOCATABLE,SAVE :: ConstUT    ! The coefficient for the UT_element matrix
! These are the spatial "matrices" that result when using the tetrahedral elements
PROTEUS_Real,DIMENSION(:),  ALLOCATABLE,SAVE :: fcoef       ! 1  coeff/ele
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: pcoef       ! 36 coeff/ele: 1st 6 -> p11c(1:6), 2nd 6 -> p22(1:6), ..
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: ucoef       ! 9  coeff/ele: (/ iJ11,iJ12,iJ13,iJ21,iJ22,iJ23,iJ31,iJ32,iJ33 /)
PROTEUS_Real,DIMENSION(:),  ALLOCATABLE,SAVE :: SurfaceAxB ! The surface coefficient for the W matrix
! These are required to identify the boundary surfaces
PROTEUS_Int, DIMENSION(:,:),ALLOCATABLE,SAVE :: BCInfo      ! (2,NumElements) (1=)the number of boundary surfaces (2=) the first surface
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: IndexNormal ! Specifies the index of the unit normal for this surface (unnecessary since we're just storing them all uniquely, but would occur in a more complex code)
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: Vac_Normals ! X,Y,Z components of the unit normal for each vacuum bc surface
! Spatial matrix integration approach
PROTEUS_Int, PARAMETER                              :: FEVertices    = 4
PROTEUS_Int, PARAMETER                              :: FEGaussPoints = 4
PROTEUS_Int, PARAMETER                              :: FENumDim      = 3
! These are needed for the conventional finite element implementation where the spatial matrices are not stored by computed during each iteration
PROTEUS_FE_Real,DIMENSION(:,:),    ALLOCATABLE,SAVE :: FEShapeFunctions
PROTEUS_FE_Real,DIMENSION(:,:,:,:),ALLOCATABLE,SAVE :: FEDerivatives  ! 4,FENumDim,FEGaussPoints,NumElements
PROTEUS_FE_Real,DIMENSION(:,:)    ,ALLOCATABLE,SAVE :: FEDetJacandWgt ! FEGaussPoints,NumElements
! These allow us to translate the solution back to the serial space
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: VertexLocalToGlobal   ! The global id of each vertex
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: ElementLocalToGlobal  ! The global id of each element
! Threading setup information
PROTEUS_Int, ALLOCATABLE,SAVE :: DA_ThreadWiseWork(:,:)     ! (2,NumThreads) Gives the starting (1) and stoping (2) element ids for each thread for the disassembly operation
PROTEUS_Int, ALLOCATABLE,SAVE :: MM_ThreadWiseWork(:,:)     ! (2,NumThreads) Gives the starting (1) and stoping (2) element ids for each thread for the Matrix-matrix product
PROTEUS_Int, ALLOCATABLE,SAVE :: AS_ThreadWiseWork(:,:,:)   ! (2,AS_NumColors,NumThreads) Gives the starting (1) and stoping (2) element ids for each thread for the assembly operation

CONTAINS

SUBROUTINE CommonBlock_Allocate()

! Mesh related info
ALLOCATE( GlobalXYZ(NumVertices,3) )
ALLOCATE( Conn(FEVertices,NumElements)   )  ! Connectivity matrix (lists 4 global vertices for 1st element, 4 global vertices for 2nd element, etc.)
ALLOCATE( ElementsWithVacuumBCs(NumVacuum))  ! Elements which have a vacuum boundary condition (listed once for each vacuum surface)
ALLOCATE( LocalSurfaceIndex(NumVacuum)    )  ! Identifies the reference surface with the vacuum bc (1,2,3 or 4)
! Cross sections
ALLOCATE( ConstTau(NumElements)) ! space dependent element coefficient
ALLOCATE( ConstF(NumElements) )  ! method-dependent coefficient for the F_element matrix (=sigt+tau*a1*sigt**2)
ALLOCATE( ConstU(NumElements) )  ! method-dependent coefficient for the U_element matrix (=tau*a3*sigt)
ALLOCATE( ConstUT(NumElements))  ! method-dependent coefficient for the UT_element matrix (=1+tau*a2*sigt)
! Spatial "matrices"
ALLOCATE( fcoef(NumElements)   ) ! 1  coeff/ele
ALLOCATE( pcoef(36,NumElements)) ! 36 coeff/ele: 1st 6 -> p11c(1:6), 2nd 6 -> p22(1:6), ..
ALLOCATE( ucoef(9,NumElements) ) ! 9  coeff/ele: (/ iJ11,iJ12,iJ13,iJ21,iJ22,iJ23,iJ31,iJ32,iJ33 /)
ALLOCATE( SurfaceAxB(NumVacuum)) ! method-dependent coefficient for the W_element matrix (=surface area)
! Numerical evaluation of the element
ALLOCATE( FEShapeFunctions(FEVertices,FEGaussPoints))
ALLOCATE( FEDerivatives(FEVertices,FENumDim,FEGaussPoints,NumElements))
ALLOCATE( FEDetJacandWgt(FEGaussPoints,NumElements))
! Boundary surface info
ALLOCATE( BCInfo(2,NumElements)   )   ! The boundary condition setup
ALLOCATE( Vac_Normals(3,NumVacuum))   ! X,Y,Z components of the unit normal for each vacuum bc surface
ALLOCATE( IndexNormal(NumVacuum)  )   ! Specifies the index of the unit normal for this surface (unnecessary since we're just storing them all uniquely, but would occur in a more complex code)
! Remapping arrays
ALLOCATE( VertexLocalToGlobal(NumVertices))
ALLOCATE( ElementLocalToGlobal(NumElements))

END SUBROUTINE CommonBlock_Allocate

SUBROUTINE CommonBlock_Deallocate()
! Deallocate the local memory
DEALLOCATE( ConstTau,ConstF,ConstU,ConstUT,&
            GlobalXYZ,Conn,ElementsWithVacuumBCs,LocalSurfaceIndex,  &
            fcoef,pcoef,ucoef,SurfaceAxB, &
            BCInfo,Vac_Normals,IndexNormal, &
            FEShapeFunctions,FEDerivatives,FEDetJacandWgt, &
            VertexLocalToGlobal,ElementLocalToGlobal,       &
            DA_ThreadWiseWork,MM_ThreadWiseWork,AS_ThreadWiseWork)
END SUBROUTINE CommonBlock_Deallocate

END MODULE CommonBlock
