!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine drives the solution of the SN equations via GMRES
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE SolveWGS(Input_Iterations,IterationCount,iMethod,LHS_C,RHS_C,SN_GMRESdata)
USE CommonBlock
USE Method_Krylov
#ifdef WITHOMP
  USE OMP_LIB
#endif
IMPLICIT NONE 
! Passed in
PROTEUS_Int Input_Iterations,IterationCount,iMethod
PROTEUS_Real LHS_C(NumAngles,NumVertices),RHS_C(NumAngles,NumVertices)

TYPE (Method_Krylov_Type)  SN_GMRESdata

! Local
PROTEUS_Int  I,J,I_SizeVec
PROTEUS_Int  MyThreadID,iStart,iEnd
PROTEUS_Log  GuessIsNonZero
PROTEUS_Int  ReasonForConvergence       ! Divergence (<0), MaxIterationCount (=0), Convergence (>0)
PROTEUS_Real ResidualNorm               ! Norm of the residual
! Additional Threaded stuff
PROTEUS_Real VectorNorm
PROTEUS_Real, ALLOCATABLE :: VectorNorm_Local(:,:)
PROTEUS_Int, ALLOCATABLE :: MyStart(:), MyEnd(:) ! Buffers for starting and ending coordinates

! Function declarations
EXTERNAL SolveWGS_PassThrough_AVE1,SolveWGS_PassThrough_AVE2,SolveWGS_PassThrough_AVE3,&
         SolveWGS_PassThrough_AVE4,SolveWGS_PassThrough_AVE5,SolveWGS_PassThrough_PC

I_SizeVec = NumAngles*NumVertices

GuessIsNonZero = .FALSE.
SN_GMRESdata%Maximum_Iterations = Input_Iterations
SN_GMRESdata%Absolute_Tolerance = 1.0d-12
SN_GMRESdata%Relative_Tolerance = 1.0d-8
IF      (iMethod .EQ. 1) THEN
   WRITE(Output_Unit,'("[SN-KERNEL]...Calling FGMRES solver for AVE1")')
ELSE IF (iMethod .EQ. 2) THEN
   WRITE(Output_Unit,'("[SN-KERNEL]...Calling FGMRES solver for AVE2")')
ELSE IF (iMethod .EQ. 3) THEN
   WRITE(Output_Unit,'("[SN-KERNEL]...Calling FGMRES solver for AVE3")')
ELSE IF (iMethod .EQ. 4) THEN
   WRITE(Output_Unit,'("[SN-KERNEL]...Calling FGMRES solver for AVE4")')
ELSE 
   WRITE(Output_Unit,'("[SN-KERNEL]...Calling FGMRES solver for AVE5")')
END IF

ALLOCATE(VectorNorm_Local(NumThreads,SN_GMRESdata%BackVectors))
DO I = 1,SN_GMRESdata%BackVectors
   DO J = 1,NumThreads
      VectorNorm_Local(J,I) = 0.0d0
   END DO
END DO

ALLOCATE(MyStart(NumThreads))
ALLOCATE(MyEnd(NumThreads))
I = (NumVertices*NumAngles)/NumThreads
DO J=1,NumThreads
  MyStart(J) = (J-1)*I + 1
  IF (J .EQ. NumThreads) THEN
    MyEnd(J) = NumVertices*NumAngles
  ELSE
    MyEnd(J) = J*I
  END IF
END DO

IF      (iMethod .EQ. 1) THEN
   CALL FGMRES_Threaded(SN_GMRESdata,LHS_C,RHS_C,                                           &
                      GuessIsNonZero,ReasonForConvergence,IterationCount, &
                      ResidualNorm,VectorNorm,VectorNorm_Local,VectorNorm_Local,                               &
                      MyStart,MyEnd, &
                      SolveWGS_PassThrough_AVE1,SolveWGS_PassThrough_PC)
ELSE IF (iMethod .EQ. 2) THEN
   CALL FGMRES_Threaded(SN_GMRESdata,LHS_C,RHS_C,                                           &
                      GuessIsNonZero,ReasonForConvergence,IterationCount, &
                      ResidualNorm,VectorNorm,VectorNorm_Local,VectorNorm_Local,                               &
                      MyStart,MyEnd, &
                      SolveWGS_PassThrough_AVE2,SolveWGS_PassThrough_PC)
ELSE IF (iMethod .EQ. 3) THEN
   CALL FGMRES_Threaded(SN_GMRESdata,LHS_C,RHS_C,                                           &
                      GuessIsNonZero,ReasonForConvergence,IterationCount, &
                      ResidualNorm,VectorNorm,VectorNorm_Local,VectorNorm_Local,                               &
                      MyStart,MyEnd, &
                      SolveWGS_PassThrough_AVE3,SolveWGS_PassThrough_PC)
END IF

DEALLOCATE(VectorNorm_Local, MyStart, MyEnd)

WRITE(Output_Unit,'("[SN-KERNEL]...FGMRES returned with an error of ",1PE13.6," after ",I6," iterations")') &
   ResidualNorm,IterationCount

END SUBROUTINE SolveWGS
