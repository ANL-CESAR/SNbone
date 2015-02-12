!---------------------------------------------------------------------------------------------------------------------------------
! Assembles the table of timing data from all of the processors and prints it to the screen
!  Copyright(c) 2005 Argonne National Laboratory
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Basic_Abort()
USE Basic_CommonBlock
IMPLICIT NONE
! MPI preprocessing information handled via PETSc interface
#ifdef USEMPI
#include "mpif.h"
#endif
#include "PROTEUS_Preprocess.h"

PROTEUS_Int ReturnedError,I(10)

! Make certain everyone is ready to quit at the same time

#ifdef USEMPI
CALL  MPI_ABORT(Basic_Parallel_Comm,I,ReturnedError)
#else
CALL ABORT
#endif

END SUBROUTINE Basic_Abort
