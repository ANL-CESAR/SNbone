! This program is meant to emulate the on-node work of A*x which forms the bulk of the flop work in CFE methods
! It demonstrates the performance of simple linear tetrahedral matrix-vector products in SUPG
! ----------------------------------------------------------------------------
! Element Type: Linear Tetrahedron
! Grid:         Structured grid where each block is composed of 12 tetrahedrons
PROGRAM Driver
! Module inclusions
USE CommonBlock
USE Method_Krylov
#ifdef WITHOMP
  USE OMP_LIB
#endif
#ifdef WITHBGQHPM
#include "mpif.h"
#endif
IMPLICIT NONE
! preprocessing inclusions
#include "PROTEUS_Preprocess.h"

! Key Problem Size Variables - Specified via user command line input
PROTEUS_Int :: Input_Scheme          ! AVE-1, AVE-2, AVE-3, ... 0=all of them
PROTEUS_Int :: Input_Angles          ! Total Angles to simulate
PROTEUS_Int :: Input_Iterations      ! User-defined # of iterations (statistics)
PROTEUS_Int :: Input_Nthreads        ! number of threads to use in openmp
PROTEUS_Int :: Input_BackVectors     ! Number of back vectors to use in the solver portion

! Local
PROTEUS_Int  :: I,J,K,L,N, ReturnedError,NumArguments,ii_start,ii_end,I_SizeVec
PROTEUS_Int  :: Element,Angle,Space,MyThreadID
PROTEUS_Int,PARAMETER  :: iNumMethods = 3
CHARACTER*8  :: MyHPMname(iNumMethods)=(/ 'AVE-1','AVE-2','AVE-3'/)
PROTEUS_Real :: MyFlopCnt(iNumMethods)=(/ 4929.0d0,4929.0d0,0.0d0/) ! Flops / element/angle
PROTEUS_Int  :: MyIterCnt(iNumMethods)
PROTEUS_Real :: MyTime(iNumMethods) ! 1 is the cumulative time, 2 is the current timer setting
CHARACTER*32 :: Some_String

! Entire DOF vectors with duplications of vertices
TYPE (Method_Krylov_Type)  SN_GMRESdata   ! (Grids) at some point

PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE :: LHS_C,LHS_Answer
PROTEUS_Real,DIMENSION(:,:),ALLOCATABLE :: RHS_C

INTEGER,parameter            :: PP_AllCodes = 15
INTEGER,parameter            :: PP_NumSegments = 6 ! We will repeat the entire set of calculations this many times
CHARACTER*8 :: PP_AllEventNames(PP_AllCodes)
INTEGER*8   :: PP_values(PP_AllCodes,iNumMethods)
PROTEUS_Real :: SomeReal1,SomeReal2,AssemblyTime

#ifdef WITHBGQHPM
   CALL MPI_INIT(ReturnedError)
#endif

! Read command line input
NumArguments = COMMAND_ARGUMENT_COUNT()   ! # User defined command line arguments
IF ((NumArguments .LT. 4) .OR. (NumArguments .GT. 5)) THEN
200 FORMAT('[SN-KERNEL]',109('.'))
499 FORMAT('[SN-KERNEL] The list of arguments was incomplete....................',52('.'))
500 FORMAT('[SN-KERNEL] Version 1.0 SNaCFE mini-app to study on node performance',52('.'))
501 FORMAT('[SN-KERNEL] This mini-app is a test of the within-group FGMRES solver for a CFE SUPG based SN methodology',15('.'))
502 FORMAT('[SN-KERNEL] Usage:   snacfe.x  Scheme Iter BackV Angles Threads',52('.'))
503 FORMAT('[SN-KERNEL] Example: snacfe.x  1      100  30    32     1      ',52('.'))
504 FORMAT('[SN-KERNEL] Scheme          specifies which scheme to use for the study (0=all)...........',15('.'))
505 FORMAT('[SN-KERNEL] Iter(ation)     specifies the maximum FGMRES iterations to allow..............',30('.'))
506 FORMAT('[SN-KERNEL] Back V(ectors)  specifies the maximum back vectors to use in FGMRES...........',30('.'))
507 FORMAT('[SN-KERNEL] Angles          specifies the number of angles assigned to the local process..',30('.'))
510 FORMAT('[SN-KERNEL] T(hreads)       specifies the number of threads to use during the execution...',30('.'))
   WRITE(Output_Unit,499)
   WRITE(Output_Unit,500)
   WRITE(Output_Unit,501)
   WRITE(Output_Unit,200)
   WRITE(Output_Unit,502)
   WRITE(Output_Unit,503)
   WRITE(Output_Unit,200)
   WRITE(Output_Unit,504)
   WRITE(Output_Unit,505)
   WRITE(Output_Unit,506)
   WRITE(Output_Unit,507)
   WRITE(Output_Unit,510)
   WRITE(Output_Unit,200)
   STOP
ELSE
   WRITE(Output_Unit,500)
   WRITE(Output_Unit,501)
   WRITE(Output_Unit,200)
   WRITE(Output_Unit,502)
   WRITE(Output_Unit,503)
   WRITE(Output_Unit,200)
   ! Assign command line input to variables
   CALL GET_COMMAND_ARGUMENT(1,Some_String)
   READ(Some_String,*,IOSTAT=ReturnedError) Input_Scheme
   CALL GET_COMMAND_ARGUMENT(2,Some_String)
   READ(Some_String,*,IOSTAT=ReturnedError) Input_Iterations
   CALL GET_COMMAND_ARGUMENT(3,Some_String)
   READ(Some_String,*,IOSTAT=ReturnedError) Input_BackVectors
   CALL GET_COMMAND_ARGUMENT(4,Some_String)
   READ(Some_String,*,IOSTAT=ReturnedError) Input_Angles
   IF (NumArguments .EQ. 5) THEN
      CALL GET_COMMAND_ARGUMENT(5,Some_String)
      READ(Some_String,*,IOSTAT=ReturnedError) Input_Nthreads
   ELSE
      Input_Nthreads = 0
   END IF
   WRITE(Output_Unit,200)
END IF

#ifdef WITHOMP
   IF (Input_Nthreads .NE. 0) THEN
      call omp_set_num_threads(Input_Nthreads)
   ELSE
      Input_Nthreads = 1 !omp_get_num_threads() ! This should be unnecessary, but...
      call omp_set_num_threads(Input_Nthreads) ! We always impose serial when no option is given
   END IF
!$OMP PARALLEL
  WRITE(Output_Unit,'("[SN-KERNEL] Thread id ",I5," of ",I5)') omp_get_thread_num(),omp_get_num_threads()
!$OMP END PARALLEL
#else
   IF (Input_Nthreads .EQ. 0) Input_Nthreads = 1
#endif

! Correct stupid inputs
ii_start = Input_Scheme
ii_end   = Input_Scheme
IF ((Input_Scheme .LE. 0) .OR. (Input_Scheme .GT. iNumMethods)) THEN
   Input_Scheme = 0
   ii_start = 1
   ii_end   = iNumMethods
END IF

IF (Input_Iterations .LT. 1)  Input_Iterations = 1
IF (Input_BackVectors .LT. 3) Input_BackVectors = 3 ! This would be CG equivalent
IF (Input_BackVectors .GT. 50) Input_BackVectors = 50
IF (Input_Angles  .LT. 1) Input_Angles  = 1
NumThreads=Input_Nthreads
NumAngles=Input_Angles

! Inform the user as to the problem size we will be running
WRITE(Output_Unit,'("[SN-KERNEL] Running schemes ",I2,":",I2)') ii_start,ii_end
WRITE(Output_Unit,'("[SN-KERNEL] Number of Iterations    ",I8)') Input_Iterations
WRITE(Output_Unit,'("[SN-KERNEL] Number of Back Vectors  ",I8)') Input_BackVectors
WRITE(Output_Unit,'("[SN-KERNEL] Number of Angles        ",I8)') NumAngles
WRITE(Output_Unit,'("[SN-KERNEL] Number of Threads       ",I8)') NumThreads

! Import the mesh
WRITE(Output_Unit,'("[SN-KERNEL] Importing the processed mesh")')
CALL Import_PMesh()

! Setup the angular cubature
WRITE(Output_Unit,'("[SN-KERNEL] Setting up the angle cubature")')
CALL BuildAngleCubature()

WRITE(Output_Unit,'("[SN-KERNEL] Stenciling the spatial NZ matrix")')
CALL StencilNZmatrix(Input_Scheme)

WRITE(Output_Unit,'("[SN-KERNEL] Number of Elements         ",I8)') NumElements
WRITE(Output_Unit,'("[SN-KERNEL] Number of Vertices         ",I8)') NumVertices
WRITE(Output_Unit,'("[SN-KERNEL] Vector Size Assembled      ",I8)') NumVertices*NumAngles
WRITE(Output_Unit,'("[SN-KERNEL] Vector Size Dis-assem      ",I8)') NumElements*NumAngles*FEVertices
WRITE(Output_Unit,'("[SN-KERNEL] Number of Non-zeros        ",I8)') NZS_NonZeros
WRITE(Output_Unit,'("[SN-KERNEL] Average Connections/Vertex ",I8)') NZS_NonZeros/NumVertices

! GMRES and solution vectors
WRITE(Output_Unit,'("[SN-KERNEL] Allocating FGMRES memory and solution vectors")')
ALLOCATE( LHS_C(NumAngles,NumVertices))
ALLOCATE( LHS_Answer(NumAngles,NumVertices))
ALLOCATE( RHS_C(NumAngles,NumVertices))
I_SizeVec = NumAngles*NumVertices
CALL Method_Krylov_Define(SN_GMRESdata,I_SizeVec,Input_BackVectors,Method_Krylov_FGMRES)

! Initialize the timers for all methods
MyTime(1:iNumMethods)  = 0.0d0

! Assemble the matrix noting that it is part of the method 3 timing
WRITE(Output_Unit,'("[SN-KERNEL] Building the non-zero space-angle matrices")')
CALL GETTHETIME(SomeReal1)
CALL AssembleNZmatrix(Input_Scheme)
CALL GETTHETIME(SomeReal2)
AssemblyTime = SomeReal2 - SomeReal1
!IF ((Input_Scheme .EQ. 0) .OR. (Input_Scheme .EQ. 3)) &
!   WRITE(Output_Unit,'("[SN-KERNEL] Took ",F13.6," seconds to assemble")') AssemblyTime

! Get a source RHS_C, the solution LHS_Answer, and the guess for the LHS
I_SizeVec = NumVertices*NumAngles
WRITE(Output_Unit,'("[SN-KERNEL] Generating an answer and its associated source with size ",I9)') I_SizeVec
CALL GenerateXb(Input_Scheme,LHS_C,LHS_Answer,RHS_C)

MyIterCnt = 0

! -----------------------------------------------------
! This is the part of the code which we wish to measure
! -----------------------------------------------------
#ifdef WITHBGQHPM
   call summary_start() ! This is the MPI tracking stuff
#endif

DO I = ii_start,ii_end,1
   CALL GETTHETIME(SomeReal1)
#ifdef WITHBGQHPM
   call hpm_start(MyHPMname(I)) ! Initializes and starts the hardware counters
#endif
   CALL SolveWGS(Input_Iterations,MyIterCnt(I),I,LHS_C,RHS_C,SN_GMRESdata)
#ifdef WITHBGQHPM
      call hpm_stop(MyHPMname(I)) ! Stops the hardware counters
#endif
   CALL GETTHETIME(SomeReal2)
   MyTime(I) = SomeReal2 - SomeReal1
   ! Verify that the outgoing vector is accurate
   J = NumVertices*NumAngles
   CALL Verify(Output_Unit,J,LHS_C,LHS_Answer,MyHPMname(I))
END DO ! methods

MyFlopCnt(3) = NZS_NonZeros*NumAngles ! This is the number of mults per application of A
I = SN_GMRESdata%BackVectors
CALL PrintSummary(Output_Unit,NumElements,NumAngles,NumVertices,Input_Iterations,I, &
                  iNumMethods,ii_start,ii_end,AssemblyTime,   &
                  MyTime,MyFlopCnt,MyIterCnt,MyHPMname,PP_AllCodes,PP_NumSegments,PP_AllEventNames,PP_values)
WRITE(Output_Unit,200)

#ifdef WITHBGQHPM
   CALL summary_stop() ! This finishes the MPI measurements
   CALL MPI_Finalize()
#endif

CALL CommonBlock_Deallocate()

END PROGRAM Driver


SUBROUTINE GETTHETIME(SomeReal)
#ifdef WITHOMP
  USE OMP_LIB
#endif
PROTEUS_Real SomeReal

#ifdef WITHOMP
   SomeReal = omp_get_wtime()
#else
   CALL CPU_Time(SomeReal)
#endif

END SUBROUTINE GETTHETIME