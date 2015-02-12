!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine is a pass through point for GMRES
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE SolveWGS_PassThrough_AVE2(RHS_C,LHS_C)
USE CommonBlock
#ifdef WITHOMP
  USE OMP_LIB
#endif
IMPLICIT NONE 
! Passed in
PROTEUS_Real LHS_C(NumAngles,NumVertices),RHS_C(NumAngles,NumVertices)
PROTEUS_Int MyThreadID

#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1
#else
   MyThreadID = 1
#endif
#ifdef WITHBGQHPM
   IF (MyThreadID .EQ. 1) call hpm_start('AVE2_ApplyA')
#endif
! This barrier ensures that the incoming vector is fully defined by all threads
!$OMP BARRIER
   CALL ApplyA_AVE2_Tet_SUPG(NumElements,NumAngles,NumVertices,Conn,   &
                        AS_NumColors,AS_NumThreads,TasksPerThread,MyThreadID,AS_ThreadWiseWork, &
                        ConstF, ConstU, ConstUT,                        &
                        FEShapeFunctions,FEDerivatives,FEDetJacandWgt,  &
                        Omega,OmegaOmega,                               &
                        LHS_C,RHS_C,   &
                        Scratch_V1(1,1,MyThreadID),Scratch_V2(1,1,MyThreadID))
#ifdef WITHBGQHPM
   IF (MyThreadID .EQ. 1) call hpm_stop('AVE2_ApplyA') ! Stops the hardware counters
#endif
END SUBROUTINE SolveWGS_PassThrough_AVE2



#include "PROTEUS_Preprocess.h"
SUBROUTINE SolveWGS_PassThrough_AVE2_NoHPM(RHS_C,LHS_C)
USE CommonBlock
#ifdef WITHOMP
  USE OMP_LIB
#endif
IMPLICIT NONE 
! Passed in
PROTEUS_Real LHS_C(NumAngles,NumVertices),RHS_C(NumAngles,NumVertices)
PROTEUS_Int MyThreadID

#ifdef WITHOMP
   MyThreadID = omp_get_thread_num() + 1
#else
   MyThreadID = 1
#endif
! This barrier ensures that the incoming vector is fully defined by all threads
!$OMP BARRIER
   CALL ApplyA_AVE2_Tet_SUPG(NumElements,NumAngles,NumVertices,Conn,   &
                        AS_NumColors,AS_ThreadWiseWork(1,1,MyThreadID), &
                        ConstF, ConstU, ConstUT,                        &
                        FEShapeFunctions,FEDerivatives,FEDetJacandWgt,  &
                        Omega,OmegaOmega,                               &
                        LHS_C,RHS_C,   &
                        Scratch_V1(1,1,MyThreadID),Scratch_V2(1,1,MyThreadID))
END SUBROUTINE SolveWGS_PassThrough_AVE2_NoHPM
