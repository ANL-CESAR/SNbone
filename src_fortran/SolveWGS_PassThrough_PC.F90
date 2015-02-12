!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine is a pass through point for GMRES
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE SolveWGS_PassThrough_PC(RHS_C,LHS_C)
USE CommonBlock
#ifdef WITHOMP
  USE OMP_LIB
#endif
IMPLICIT NONE 
! Passed in
PROTEUS_Real LHS_C(NumAngles*NumVertices),RHS_C(NumAngles*NumVertices)
PROTEUS_Int I,J
PROTEUS_Int Istart,Iend

#ifdef WITHOMP
   I = omp_get_thread_num() + 1
   J = NumAngles*NumVertices/NumThreads
   Istart = (I-1)*J + 1
   Iend   =  I*J
   IF (I .EQ. NumThreads) Iend = NumAngles*NumVertices
#else
   Istart = 1
   Iend   = NumAngles*NumVertices
#endif
! No incoming barrier is needed as this preconditioner does not have any overlap, but we are by rule required to impose a barrier
DO I = Istart,Iend
   LHS_C(I) = RHS_C(I) ! The identity operation
END DO
!$OMP BARRIER

END SUBROUTINE SolveWGS_PassThrough_PC
