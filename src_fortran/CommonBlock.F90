! This module is used to emulate a common block behavior in the mini-app
! We would normally do this as a data structure save rather than individual arrays, 
!  but I chose to do it this way in the mini-app to keep it as simple as possible
! If you recieve an allocate error in the mini-app, try a smaller problem
MODULE CommonBlock
IMPLICIT NONE
#include "PROTEUS_Preprocess.h"

! Output controls
PROTEUS_Int,PARAMETER :: Output_Unit = 6   ! Output unit

! Communicator setup
PROTEUS_Int :: ParallelRank  = 0
PROTEUS_Int :: ParallelSize  = 0
PROTEUS_Int :: ParallelComm  = 0
PROTEUS_Int :: ParallelGroup = 0

! Control variables
PROTEUS_Int :: NumElements          ! Total number of tetrahedral elements (5 per structured grid)
PROTEUS_Int :: NumVertices          ! Total number of vertices in structured grid
PROTEUS_Int :: NumAngles	         ! Number of angles visible to the local process
PROTEUS_Int :: NumVacuum            ! Number of tetrahedral surfaces (triangles) with vacuum boundary (2 tri/surf and 2 global surfs)
PROTEUS_Int :: NumUnitNormals       ! Number of vacuum boundary unit normals
PROTEUS_Int :: NumThreads           ! Total number of threads
PROTEUS_Int :: AS_NumColors         ! The number of non-overlapping (in memory) pieces of work that the FE assembly operation can be broken into
PROTEUS_Int :: AS_NumThreads        ! The number of threads worth of data the mesh was broken into (Ideally this matches NumThreads)
PROTEUS_Int :: FEGaussPoints        ! The number of GaussPoints in the cubature
PROTEUS_Int :: FEVertices           ! The number of vertices per element (4)
PROTEUS_Int :: FENumDim             ! The number of dimensions in the problem (3)

! Threading work breakdown for starting (1) and stopping (2) element id.
PROTEUS_Int :: TasksPerThread = 0                           ! This variable accounds for the AS_NumThreads <> NumThreads
PROTEUS_Int, ALLOCATABLE,SAVE :: DA_ThreadWiseWork(:,:)     ! (2,AS_NumThreads) The disassembly operation
PROTEUS_Int, ALLOCATABLE,SAVE :: MM_ThreadWiseWork(:,:)     ! (2,AS_NumThreads) The matrix-matrix product
PROTEUS_Int, ALLOCATABLE,SAVE :: AS_ThreadWiseWork(:,:,:)   ! (2,AS_NumColors,AS_NumThreads) The assembly operation
!    AS_ThreadWiseWork(2,AS_NumColors,AS_NumThreads) 
! ThreadID = omp_get_thread_num() + 1
! AS_ThreadID = 
! Mesh related data structure
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: GlobalXYZ             ! The global element vertex positions
PROTEUS_Int, DIMENSION(:,:),ALLOCATABLE,SAVE :: Conn                  ! Connectivity matrix (lists 4 global vertices for 1st element, 4 global vertices for 2nd element, etc.)
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: ElementsWithVacuumBCs ! Elements which have a vacuum boundary condition (listed once for each vacuum surface)
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: LocalSurfaceIndex     ! Identifies the reference surface with the vacuum bc (1,2,3 or 4)
! Element-wise "cross section" values which change with respect to energy
PROTEUS_Real,DIMENSION(:), ALLOCATABLE,SAVE :: ConstTau ! The element-wise stabilization factor
PROTEUS_Real,DIMENSION(:), ALLOCATABLE,SAVE :: ConstF   ! The coefficient for the F_element matrix
PROTEUS_Real,DIMENSION(:), ALLOCATABLE,SAVE :: ConstU   ! The coefficient for the U_element matrix
PROTEUS_Real,DIMENSION(:), ALLOCATABLE,SAVE :: ConstUT  ! The coefficient for the UT_element matrix
! These are the spatial "matrices" that result when using the tetrahedral elements
PROTEUS_Real,DIMENSION(:),  ALLOCATABLE,SAVE :: fcoef       ! 1  coeff/ele
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: pcoef       ! 36 coeff/ele: 1st 6 -> p11c(1:6), 2nd 6 -> p22(1:6), ..
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: ucoef       ! 9  coeff/ele: (/ iJ11,iJ12,iJ13,iJ21,iJ22,iJ23,iJ31,iJ32,iJ33 /)
PROTEUS_Real,DIMENSION(:),  ALLOCATABLE,SAVE :: SurfaceAxB ! The surface coefficient for the W matrix
! These are needed for the conventional finite element implementation where the spatial matrices are not stored by computed during each iteration
PROTEUS_FE_Real,DIMENSION(:,:),    ALLOCATABLE,SAVE :: FEShapeFunctions
PROTEUS_FE_Real,DIMENSION(:,:,:,:),ALLOCATABLE,SAVE :: FEDerivatives  ! FEVertices,FENumDim,FEGaussPoints,NumElements
PROTEUS_FE_Real,DIMENSION(:,:)    ,ALLOCATABLE,SAVE :: FEDetJacandWgt ! FEGaussPoints,NumElements
! Angular related data structure which makes up the first dimension of the matrix-matrix product
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: Omega          ! Omega stores O1, O2, O3 components of each angle
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: OmegaOmega     ! OmegaOmega stores O1*O1,O2*O2,O3*O3,O1*O2,O1*O3,O2*O3
PROTEUS_Real,DIMENSION(:),  ALLOCATABLE,SAVE :: AngleWeights   ! This stores the angular weights used to integrate spherical harmonics
! These are required to identify the boundary surfaces
PROTEUS_Int, DIMENSION(:,:),ALLOCATABLE,SAVE :: BCInfo      ! (2,NumElements) (1=)the number of boundary surfaces (2=) the first surface
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: IndexNormal ! Specifies the index of the unit normal for this surface (unnecessary since we're just storing them all uniquely, but would occur in a more complex code)
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE,SAVE :: Vac_Normals ! X,Y,Z components of the unit normal for each vacuum bc surface
! These allow us to translate the solution back to the serial space
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: VertexLocalToGlobal   ! The global id of each vertex
PROTEUS_Int, DIMENSION(:),  ALLOCATABLE,SAVE :: ElementLocalToGlobal  ! The global id of each element
! Scratch vector storage
PROTEUS_Real,DIMENSION(:,:,:),  ALLOCATABLE,SAVE :: Scratch_V1,Scratch_V2 ! (NumAngles,FEVertices,Threads)
! Sparse matrix storage 
PROTEUS_Int                                    :: NZS_NonZeros = 0
PROTEUS_Int, DIMENSION(:),    ALLOCATABLE,SAVE :: NZS_RowLoc ! (NumVertices+1)  The NZS_Data storage position
PROTEUS_Int, DIMENSION(:),    ALLOCATABLE,SAVE :: NZS_ColNum ! (NZS_NonZeros)   The column numbers associated with NZS_Data
PROTEUS_Real,DIMENSION(:,:),  ALLOCATABLE,SAVE :: NZS_Data   ! (NumAngles,NZS_NonZeros)

CONTAINS

SUBROUTINE CommonBlock_Allocate()

! Threaded coloring info
ALLOCATE(DA_ThreadWiseWork(2,AS_NumThreads))
ALLOCATE(MM_ThreadWiseWork(2,AS_NumThreads))
ALLOCATE(AS_ThreadWiseWork(2,AS_NumColors,AS_NumThreads))
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
! Angular Omega matrices
ALLOCATE( Omega(NumAngles,3)     )  ! Omega stores O1, O2, O3 components of each angle
ALLOCATE( OmegaOmega(NumAngles,6))  ! OmegaOmega stores O1*O1,O2*O2,O3*O3,O1*O2,O1*O3,O2*O3
ALLOCATE( AngleWeights(NumAngles))  ! Store the angular weights
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
ALLOCATE( BCInfo(2,NumElements)   )       ! The boundary condition setup
ALLOCATE( Vac_Normals(3,NumUnitNormals))  ! X,Y,Z components of the unit normal for each vacuum bc surface
ALLOCATE( IndexNormal(NumVacuum)  )       ! Specifies the index of the unit normal for this surface (unnecessary since we're just storing them all uniquely, but would occur in a more complex code)
! Remapping arrays
ALLOCATE( VertexLocalToGlobal(NumVertices))
ALLOCATE( ElementLocalToGlobal(NumElements))
! Scratch arrays needed for the threaded GMRES solve process
ALLOCATE( Scratch_V1(NumAngles,FEVertices,NumThreads),&
          Scratch_V2(NumAngles,FEVertices,NumThreads) )

END SUBROUTINE CommonBlock_Allocate

SUBROUTINE CommonBlock_Deallocate()

! Deallocate the local memory
DEALLOCATE( ConstTau,ConstF,ConstU,ConstUT,&
            Omega,OmegaOmega, &
            Conn,ElementsWithVacuumBCs,LocalSurfaceIndex,Vac_Normals,IndexNormal, &
            SurfaceAxB, &
            DA_ThreadWiseWork,MM_ThreadWiseWork,AS_ThreadWiseWork, &
            NZS_RowLoc,NZS_ColNum,NZS_Data)

DEALLOCATE( Scratch_V1,Scratch_V2)

END SUBROUTINE CommonBlock_Deallocate

END MODULE CommonBlock
