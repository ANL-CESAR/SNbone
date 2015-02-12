! This program takes a given tetrahedral mesh and processes it for use in the CFE SN mini-app
! ----------------------------------------------------------------------------
PROGRAM Driver
! Module inclusions
USE CommonBlock
IMPLICIT NONE
! preprocessing inclusions
#include "PROTEUS_Preprocess.h"

! Key Problem Size Variables - Specified via user command line input
PROTEUS_Int :: Input_Part          ! The type of thread decomposition scheme to use
PROTEUS_Int :: Input_nthreads        ! number of threads to use in openmp
! Local
CHARACTER*32 :: Some_String
PROTEUS_Int  :: NumArguments
PROTEUS_Int ReturnedError

! Read command line input
NumArguments = COMMAND_ARGUMENT_COUNT()   ! # User defined command line arguments
IF ((NumArguments .LT. 2) .OR. (NumArguments .GT. 2)) THEN
200 FORMAT('[SN-KERNEL]',109('.'))
499 FORMAT('[SN-KERNEL] The list of arguments was incomplete....................',52('.'))
500 FORMAT('[SN-KERNEL] Version 1.0 mesh processing utility.....................',52('.'))
502 FORMAT('[SN-KERNEL] Usage:   processmesh.x  Partition Threads...............',52('.'))
503 FORMAT('[SN-KERNEL] Example: processmesh.x  1         32     ...............',52('.'))
504 FORMAT('[SN-KERNEL] Scheme   Use (1) METIS, (other) GREEDY..................',52('.'))
510 FORMAT('[SN-KERNEL] Threads  Maximum threads to setup for in the mini-app...',52('.'))
   WRITE(Output_Unit,200)
   WRITE(Output_Unit,499)
   WRITE(Output_Unit,500)
   WRITE(Output_Unit,200)
   WRITE(Output_Unit,502)
   WRITE(Output_Unit,503)
   WRITE(Output_Unit,504)
   WRITE(Output_Unit,510)
   WRITE(Output_Unit,200)
   STOP
ELSE
   WRITE(Output_Unit,200)
   WRITE(Output_Unit,500)
   WRITE(Output_Unit,502)
   WRITE(Output_Unit,503)
   WRITE(Output_Unit,200)
   ! Assign command line input to variables
   CALL GET_COMMAND_ARGUMENT(1,Some_String)
   READ(Some_String,*,IOSTAT=ReturnedError) Input_Part
   CALL GET_COMMAND_ARGUMENT(2,Some_String)
   READ(Some_String,*,IOSTAT=ReturnedError) Input_Nthreads
   WRITE(Output_Unit,200)
END IF

! Correct stupid inputs
IF (Input_Nthreads .EQ. 0) Input_Nthreads = 1
NumThreads = Input_Nthreads

! Inform the user as to the problem size we will be running
IF (Input_Part .EQ. 1) THEN
   WRITE(Output_Unit,'("[SN-KERNEL] Preparing METIS  decomposition")') 
ELSE
   WRITE(Output_Unit,'("[SN-KERNEL] Preparing GREEDY decomposition")')
END IF
WRITE(Output_Unit,'("[SN-KERNEL] Number of Threads       ",I8)') Input_Nthreads

WRITE(Output_Unit,'("[SN-KERNEL] Importing the mesh")')
CALL Import_Mesh

WRITE(Output_Unit,'("[SN-KERNEL] Number of Elements    ",I8)') NumElements
WRITE(Output_Unit,'("[SN-KERNEL] Number of Vertices    ",I8)') NumVertices
WRITE(Output_Unit,'("[SN-KERNEL] Number of Threads     ",I8)') NumThreads

WRITE(Output_Unit,'("[SN-KERNEL] Re-arranging the connectivity by coloring the domain")')
CALL RearrangeConn(Input_Part)

WRITE(Output_Unit,'("[SN-KERNEL] Obtaining the spatial matrices")')
CALL GetSpatialMatrices()

WRITE(Output_Unit,'("[SN-KERNEL] Exporting the processed mesh")')
CALL Export_PMesh()

CALL CommonBlock_Deallocate()

END PROGRAM Driver
