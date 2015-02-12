!--------------------------------------------------------------------------------------------------------------------------------- 
!  Method_Krylov : data module for setting up the memory for Krylov Solver
!  Copyright(c) 2005 Argonne National Laboratory
!---------------------------------------------------------------------------------------------------------------------------------
! Method_Krylov_SetUnit      ! This routine sets the output unit and debug print settings
! Method_Krylov_Define       ! This routine allocates the data structure
! Method_Krylov_Void         ! This routine voids the data structure
! Method_Krylov_Print        !  Prints the data type
! ----------------------------
! Method_Krylov_SetUnit(Output_Unit)
! Method_Krylov_Define(User_Krylov,Local_Owned,BackVectors,Krylov_Method)
! Method_Krylov_Void(User_Krylov)
! Method_Krylov_CheckError(ErrorCondition,StringToPrint)
! Method_Krylov_Print(User_Krylov)
!---------------------------------------------------------------------------------------------------------------------------------

#include "PROTEUS_Preprocess.h"

MODULE Method_Krylov
IMPLICIT NONE

TYPE Method_Krylov_Type
   ! Global informative variables
   PROTEUS_Log           :: Defined       = .FALSE. ! Logical flag indicating that the data structure is setup
   PROTEUS_Real          :: MemoryUsage   = 0.0D0   ! The total memory usage (MB)
   PROTEUS_Int           :: Local_Owned   = 0       ! The number of dof assigned to this processor
   PROTEUS_Int_64bit     :: Local_Visible = 0       ! The number of dof visible to this processor
   PROTEUS_Int           :: Method        = 0       ! A flag to identify the method used (None,CG,GMRES,FGMRES,...)
   ! Solution vector
   PROTEUS_Int           :: AllocationSize= 0       ! Duplicate Local_Owned to store the allocated size
   PROTEUS_Real, POINTER :: Storage_Local(:)        ! (Local_Owned)      Locally assigned parallel flux vector
   ! Krylov solver variable
   PROTEUS_Int           :: BackVectors          = 0       ! Number of BackVectors to be used
   PROTEUS_Int           :: Maximum_Iterations   = 0       ! Maximum number of iterations
   PROTEUS_Int           :: Iterations           = 0       ! Current number of iteration performed to solve the problem
   PROTEUS_Real          :: Absolute_Tolerance   = 1.0d-4   ! Absolute tolerance for convergence
   PROTEUS_Real          :: Relative_Tolerance   = 1.0d-3   ! Relative tolerance for convergence (relative to norm of the first residual)
   PROTEUS_Real          :: Divergence_Tolerance = 1.0d-3   ! Relative tolerance to determine divergence (relative to norm of the first residual)
   PROTEUS_Real, POINTER :: Basis(:,:)                     ! (Local_Owned,BackVectors) Krylov subspace orthonormal basis vector
   PROTEUS_Real, POINTER :: Hessenberg(:,:)                ! (BackVectors,BackVectors+1) Hessenberg matrix coefficients (dot products)
   PROTEUS_Real, POINTER :: Givens(:,:)                    ! (BackVectors,2) Givens rotation (cos,sin) coefficients
   PROTEUS_Real, POINTER :: PC_Basis(:,:)                  ! (Local_Owned,BackVectors) Used by FGMRES to store the preconditionned basis vector
   PROTEUS_Real, POINTER :: Modified_RHS(:)                ! (BackVectors) Modified right hand side (used by GMRES and FGMRES) 
END TYPE Method_Krylov_Type

PROTEUS_Int :: MODULE_OUT = 6

! Define some constant to identify the different Krylov methodology
PROTEUS_Int, PARAMETER :: Method_Krylov_None   = 0
PROTEUS_Int, PARAMETER :: Method_Krylov_CG     = 1
PROTEUS_Int, PARAMETER :: Method_Krylov_GMRES  = 2
PROTEUS_Int, PARAMETER :: Method_Krylov_FGMRES = 3

PRIVATE MODULE_OUT

CONTAINS

! --------------------------------------------------------------------------------------------------------------------------------
! Method_Krylov_SetUnit(Output_Unit)
! This routine sets the output unit and debug print settings
! --------------------------------------------------------------------------------------------------------------------------------
!  Copyright(c) 2005 Argonne National Laboratory
#undef __FUNCT__
#define __FUNCT__ "Method_Krylov_SetUnit"
SUBROUTINE Method_Krylov_SetUnit(Output_Unit)
IMPLICIT NONE
PROTEUS_Int Output_Unit                             ! Reset the object output unit
MODULE_OUT = Output_Unit
END SUBROUTINE Method_Krylov_SetUnit

! --------------------------------------------------------------------------------------------------------------------------------
! Method_Krylov_Define(User_Krylov,Local_Owned,BackVectors,Krylov_Method)
! This routine allocates the data structure
! --------------------------------------------------------------------------------------------------------------------------------
!  Copyright(c) 2005 Argonne National Laboratory
#undef __FUNCT__
#define __FUNCT__ "Method_Krylov_Define"
SUBROUTINE Method_Krylov_Define(User_Krylov,Local_Owned,BackVectors,Krylov_Method)
IMPLICIT NONE
! Passed in
TYPE (Method_Krylov_Type)  User_Krylov
PROTEUS_Int                   Local_Owned
PROTEUS_Int                   BackVectors
PROTEUS_Int                   Krylov_Method

! Local
PROTEUS_Int IOS
PROTEUS_Int I,J

100 FORMAT('[Method]...Sorry, but I must stop')
105 FORMAT('[Method]...There was a fatal error that occured in (Define)',38('.'))

! Void the data structure if it is already defined
IF (User_Krylov%Defined) CALL Method_Krylov_Void(User_Krylov)

IF (Local_Owned .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[Method]...Local specification was bogus",75("."))')
   WRITE(MODULE_OUT,'("[Method]...Number of local     dofs....",75("."),I16)') Local_Owned
   WRITE(MODULE_OUT,100)
   CALL Basic_Abort
END IF

IF (BackVectors .LE. 0) THEN
   WRITE(MODULE_OUT,105)
   WRITE(MODULE_OUT,'("[Method]...Local specification was bogus",75("."))')
   WRITE(MODULE_OUT,'("[Method]...Number of local backvectors....",75("."),I16)') BackVectors
   WRITE(MODULE_OUT,100)
   CALL Basic_Abort
END IF

SELECT CASE (Krylov_Method)
   CASE(Method_Krylov_None)
      ALLOCATE(User_Krylov%Storage_Local(Local_Owned),User_Krylov%Basis(1,1),User_Krylov%Hessenberg(1,1), &
               User_Krylov%Givens(1,1),User_Krylov%PC_Basis(1,1),User_Krylov%Modified_RHS(1),STAT=IOS)
      IF (IOS .NE. 0) THEN
         WRITE(MODULE_OUT,105)
         WRITE(MODULE_OUT,'("[Method]...Failed to allocate memory",75("."))')
         WRITE(MODULE_OUT,100)
         CALL Basic_Abort
      END IF
      User_Krylov%MemoryUsage = (Local_Owned + 1 + 1 + 1 + 1 + 1)*PROTEUS_Real_Size
      DO I = 1,Local_Owned
         User_Krylov%Storage_Local(I) = 0.0d0
      END DO
   CASE(Method_Krylov_CG)
      ALLOCATE(User_Krylov%Storage_Local(Local_Owned),User_Krylov%Basis(Local_Owned,4),User_Krylov%Hessenberg(1,1), &
               User_Krylov%Givens(1,1),User_Krylov%PC_Basis(1,1),User_Krylov%Modified_RHS(1),STAT=IOS)
      IF (IOS .NE. 0) THEN
         WRITE(MODULE_OUT,105)
         WRITE(MODULE_OUT,'("[Method]...Failed to allocate memory",75("."))')
         WRITE(MODULE_OUT,100)
         CALL Basic_Abort
      END IF
      User_Krylov%MemoryUsage = (Local_Owned + Local_Owned*BackVectors + 1 + 1 + 1 + 1)*PROTEUS_Real_Size
      DO I = 1,Local_Owned
         User_Krylov%Storage_Local(I) = 0.0d0
         DO J = 1,4
            User_Krylov%Basis(I,J) = 0.0d0
         END DO 
      END DO
    CASE(Method_Krylov_GMRES)
      ! In GMRES we use PC_Basis as a scratch storage 
      ! PC_Basis(:,1) will receive the preconditioned RHS
      ! PC_Basis(:,2) will be a scratch vector 
      ALLOCATE(User_Krylov%Storage_Local(Local_Owned),User_Krylov%Basis(Local_Owned,BackVectors+1),&
               User_Krylov%Hessenberg(BackVectors+1,BackVectors+1), &
               User_Krylov%Givens(BackVectors+1,2),User_Krylov%PC_Basis(Local_Owned,1),&
               User_Krylov%Modified_RHS(BackVectors+1),STAT=IOS)
      IF (IOS .NE. 0) THEN
         WRITE(MODULE_OUT,105)
         WRITE(MODULE_OUT,'("[Method]...Failed to allocate memory",75("."))')
         WRITE(MODULE_OUT,100)
         CALL Basic_Abort
      END IF
      User_Krylov%MemoryUsage = (Local_Owned + Local_Owned*(BackVectors+1) + (BackVectors+1)*(BackVectors+1) &
            + 2*(BackVectors+1) + 1 + BackVectors+1)*PROTEUS_Real_Size
      DO I = 1,Local_Owned
         User_Krylov%Storage_Local(I) = 0.0d0
         User_Krylov%PC_Basis(I,1) = 0.0d0
         DO J = 1,BackVectors+1
            User_Krylov%Basis(I,J) = 0.0d0
         END DO
      END DO
      DO I = 1,BackVectors+1
         User_Krylov%Modified_RHS(I) = 0.0d0
         DO J = 1,BackVectors+1
            User_Krylov%Hessenberg(I,J) = 0.0d0
         END DO
         DO J = 1,2
            User_Krylov%Givens(I,J) = 0.0d0
         END DO
      END DO
    CASE(Method_Krylov_FGMRES)
      ALLOCATE(User_Krylov%Storage_Local(Local_Owned),User_Krylov%Basis(Local_Owned,BackVectors+1),&
               User_Krylov%Hessenberg(BackVectors+1,BackVectors+1), &
               User_Krylov%Givens(BackVectors+1,2),User_Krylov%PC_Basis(Local_Owned,BackVectors+1),&
               User_Krylov%Modified_RHS(BackVectors+1),STAT=IOS)
      IF (IOS .NE. 0) THEN
         WRITE(MODULE_OUT,105)
         WRITE(MODULE_OUT,'("[Method]...Failed to allocate memory",75("."))')
         WRITE(MODULE_OUT,100)
         CALL Basic_Abort
      END IF
      User_Krylov%MemoryUsage = (Local_Owned + Local_Owned*(BackVectors+1.0d0) &
                              + (BackVectors+1)*(BackVectors+1.0d0) + 2.0d0*(BackVectors+1.0d0) &
                              + Local_Owned*(BackVectors+1.0d0) + BackVectors+1.0d0)*PROTEUS_Real_Size
      DO I = 1,Local_Owned
         User_Krylov%Storage_Local(I) = 0.0d0
         DO J = 1,BackVectors+1
            User_Krylov%Basis(I,J) = 0.0d0
            User_Krylov%PC_Basis(I,J) = 0.0d0
         END DO
      END DO
      DO I = 1,BackVectors+1
         User_Krylov%Modified_RHS(I) = 0.0d0
         DO J = 1,BackVectors+1
            User_Krylov%Hessenberg(I,J) = 0.0d0
         END DO
         DO J = 1,2
            User_Krylov%Givens(I,J) = 0.0d0
         END DO
      END DO

    CASE DEFAULT
       WRITE(MODULE_OUT,105)
       WRITE(MODULE_OUT,'("[Method]...Local specification was bogus",75("."))')
       WRITE(MODULE_OUT,'("[Method]...Unknown krylov methodology",67("."),I16)') Krylov_Method 
       WRITE(MODULE_OUT,100)
       CALL Basic_Abort
END SELECT
User_Krylov%Defined       = .TRUE.
User_Krylov%Local_Owned   = Local_Owned
User_Krylov%AllocationSize= Local_Owned
User_Krylov%Local_Visible = 0
User_Krylov%BackVectors   = BackVectors
User_Krylov%Method        = Krylov_Method 

END SUBROUTINE Method_Krylov_Define

! --------------------------------------------------------------------------------------------------------------------------------
! Method_Krylov_Void(User_Krylov)  
! This routine voids the data structure
! --------------------------------------------------------------------------------------------------------------------------------
!  Copyright(c) 2005 Argonne National Laboratory
#undef __FUNCT__
#define __FUNCT__ "Method_Krylov_Void"
SUBROUTINE Method_Krylov_Void(User_Krylov)
IMPLICIT NONE
! Passed in
TYPE (Method_Krylov_Type) :: User_Krylov
! Local
PROTEUS_Int ReturnedError

100 FORMAT('[Method]...Sorry, but I must stop')
105 FORMAT('[Method]...There was a fatal error that occured in (Void)',42('.'))

ReturnedError = 0
IF (User_Krylov%Defined) THEN
   DEALLOCATE(User_Krylov%Storage_Local,User_Krylov%Basis,User_Krylov%Hessenberg,User_Krylov%Givens,User_Krylov%PC_Basis, &
              User_Krylov%Modified_RHS,STAT=ReturnedError)
   CALL Basic_CheckError(MODULE_OUT,ReturnedError,'Method_Krylov deallocation in void')
END IF

User_Krylov%Defined       = .FALSE.
User_Krylov%MemoryUsage   = 0.0D0
User_Krylov%Local_Owned   = 0
User_Krylov%AllocationSize= 0
User_Krylov%Local_Visible = 0
User_Krylov%BackVectors   = 0
User_Krylov%Method        = 0
User_Krylov%Maximum_Iterations   = 0
User_Krylov%Iterations           = 0
User_Krylov%Absolute_Tolerance   = 1.0d-4
User_Krylov%Relative_Tolerance   = 1.0d-3
User_Krylov%Divergence_Tolerance = 1.0d-3
END SUBROUTINE Method_Krylov_Void

!---------------------------------------------------------------------------------------------------------------------------------
! Method_Krylov_Print(User_Krylov)
! Prints the data type
!---------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Method_Krylov_Print(User_Krylov)
IMPLICIT NONE
! Passed in
TYPE (Method_Krylov_Type) User_Krylov
! Local
PROTEUS_Int I
CHARACTER(16) Krylov_Method(0:10)

100 FORMAT('[Method]...Sorry, but I must stop')
105 FORMAT('[Method]...There was a non-fatal error that occured in (Krylov_Print)',42('.'))

!PROTEUS_Int :: Method_Krylov_None   = 0
!PROTEUS_Int :: Method_Krylov_CG     = 1
!PROTEUS_Int :: Method_Krylov_GMRES  = 2
!PROTEUS_Int :: Method_Krylov_FGMRES = 3
Krylov_Method(Method_Krylov_None)   = "NONE"
Krylov_Method(Method_Krylov_CG)     = "CG"
Krylov_Method(Method_Krylov_GMRES)  = "GMRES"
Krylov_Method(Method_Krylov_FGMRES) = "FGMRES"

IF (User_Krylov%Defined) THEN
   WRITE(MODULE_OUT,200) User_Krylov%Local_Owned
   WRITE(MODULE_OUT,201) User_Krylov%AllocationSize
   WRITE(MODULE_OUT,205) User_Krylov%Local_Visible
   WRITE(MODULE_OUT,210) User_Krylov%MemoryUsage
   WRITE(MODULE_OUT,215) User_Krylov%BackVectors
   WRITE(MODULE_OUT,220) Krylov_Method(User_Krylov%Method)
   200 FORMAT('[Method]...The number of local moments........',34('.'),I16)
   201 FORMAT('[Method]...The allocation size = local moments',34('.'),I16)
   205 FORMAT('[Method]...The number of visible moments......',34('.'),I16)
   210 FORMAT('[Method]...The total memory usage.............',34('.'),F16.9)
   215 FORMAT('[Method]...The number of backvectors..........',34('.'),I16)
   220 FORMAT('[Method]...The Krylov methodology used........',34('.'),A16)
   DO I = 1,User_Krylov%Local_Owned
      WRITE(MODULE_OUT,310) I,User_Krylov%Storage_Local(I)
   END DO
   310 FORMAT('[Method]...',I8,1X,1PE16.9)
END IF

END SUBROUTINE Method_Krylov_Print

END MODULE Method_Krylov
