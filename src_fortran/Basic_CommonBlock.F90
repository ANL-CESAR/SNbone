!---------------------------------------------------------------------------------------------------------------------------------
! This is a common block used to work around mpi abort and barrier issues
!  Copyright(c) 2005 Argonne National Laboratory
!---------------------------------------------------------------------------------------------------------------------------------
MODULE Basic_CommonBlock
IMPLICIT NONE
#include "PROTEUS_Preprocess.h"

PROTEUS_Int, SAVE  :: Basic_Output_Unit = 6
PROTEUS_Int, SAVE  :: Basic_Parallel_Comm = 6

CONTAINS

SUBROUTINE Basic_CommonBlock_SetValues(Output_Unit,Parallel_Comm)
PROTEUS_Int Output_Unit,Parallel_Comm
Basic_Output_Unit   = Output_Unit
Basic_Parallel_Comm = Parallel_Comm
END SUBROUTINE Basic_CommonBlock_SetValues

END MODULE Basic_CommonBlock
