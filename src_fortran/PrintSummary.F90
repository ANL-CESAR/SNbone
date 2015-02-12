! This program checks the solution
! ----------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE PrintSummary(Output_Unit,NumElements,NumAngles,NumVertices,Input_Iterations,BackVectors, &
                        NumMethods,ii_start,ii_end,AssemblyTime,   &
                        MyTime,MyFlopCnt,MyIterCnt,MyHPMname,PP_AllCodes,PP_NumSegments,PP_AllEventNames,PP_values)
IMPLICIT NONE
PROTEUS_Int,  INTENT(IN) :: Output_Unit
PROTEUS_Int,  INTENT(IN) :: NumElements,NumAngles,NumVertices,Input_Iterations,BackVectors
PROTEUS_Int,  INTENT(IN) :: NumMethods,ii_start,ii_end
PROTEUS_Real, INTENT(IN) :: AssemblyTime
PROTEUS_Real, INTENT(IN) :: MyTime(NumMethods)
PROTEUS_Real, INTENT(INOUT) :: MyFlopCnt(NumMethods)
PROTEUS_Int,  INTENT(IN) :: MyIterCnt(NumMethods)
CHARACTER*8,  INTENT(IN) :: MyHPMname(NumMethods)
PROTEUS_Int,  INTENT(IN) :: PP_AllCodes,PP_NumSegments
CHARACTER*(*),INTENT(IN) :: PP_AllEventNames(PP_AllCodes)
INTEGER*8,    INTENT(IN) :: PP_values(PP_AllCodes,NumMethods)
! Local
PROTEUS_Int  :: I,J,K,L
PROTEUS_Real :: Temp1,Temp2,Temp3,Temp4,Temp5,Temp6,ddtemp
PROTEUS_Real :: TempSc1(100),TempSc2(100),TempSc3(100),TempSc4(100),TempSc5(100)
PROTEUS_Real :: MyMemory(100) ! This will just be the memory required to apply the A matrix
PROTEUS_Real :: EstimatedFlops(100) ! This is the estimated flop work that was done

! -------------------------------------------------------------------------------------------------
! Convert the flops estimate per element/angle into flops for this particular problem
Temp1 = NumElements*NumAngles
MyFlopCnt(1) = MyFlopCnt(1) * Temp1
MyFlopCnt(2) = MyFlopCnt(2) * Temp1
!MyFlopCnt(3) = MyFlopCnt(3) * Temp1 No change on this one as it uses

! -------------------------------------------------------------------------------------------------
! Compute the amount of memory and flops used in each methodology
! Method 1 uses gaussian integration
Temp1 = NumElements
Temp2 = Temp1*NumAngles
Temp3 = NumVertices*NumAngles
Temp4 = 4.d0*Temp1
MyMemory = 0.0d0
MyMemory(1) = Temp4*0.5d0                                & ! Connectivity
            + 3.0d0 * Temp1                              & ! ConstF UT and Tau
            + 4.d0*4.d0 + (4.d0*3.d0*4.d0)*Temp1 + 4.0d0 & ! Basis and derivative function evaluations
            + NumAngles * 9.0d0                            ! Omega and OmegaOmega
! Method 2 is a copy
MyMemory(2) = MyMemory(1)
! Method 5 uses assembled stored matrices
Temp5 = NumAngles
MyMemory(3) = MyFlopCnt(3) * 1.5d0 + (NumVertices+1.0d0)  ! Total storage of the space-angle assembled matrix

! -------------------------------------------------------------------------------------------------
! Intermediate memory print out
!Temp4 = BackVectors
Temp4 = MIN(MyIterCnt(ii_start),BackVectors) ! The number of iterations used to solve the problem
ddtemp =   Temp3                 & ! Solution vector
           + Temp3*Temp4*2.d0    & ! Basis and PC basis
           + Temp4*(Temp4+1.0d0) & ! Hessenberg
           + (Temp4+1.0d0)*3.d0    ! Givens and Modified RHS
!WRITE(Output_Unit,1770) (4.d0*NumElements)*0.5d0/131072.d0
!WRITE(Output_Unit,1771) (NumElements*1.0d0)*(1.d0 + 9.d0 + 36.d0)/131072.d0
WRITE(Output_Unit,1773) ddtemp/131072.d0
!1770 FORMAT('[SN-KERNEL] Connectivity memory   (MB) ',F13.3)
!1771 FORMAT('[SN-KERNEL] Spatial matrix data   (MB) ',F13.3)
1773 FORMAT('[SN-KERNEL] Average FGMRES memory (MB) ',F13.3)
!WRITE(Output_Unit,*)'MyMemory=',MyMemory(1:6)/131072.d0

! -------------------------------------------------------------------------------------------------
! Include the FGMRES memory used in each case and the FGMRES flop work
Temp2 = NumVertices*NumAngles
Temp3 = NumVertices*NumAngles
Temp4 = 0.0d0
DO I = 1,BackVectors
   Temp4 = Temp4 + I*2.0d0 + 2.0d0 ! The number of dot products and scaling operations needed for each inner in FGMRES
END DO
EstimatedFlops = 0.0d0 
DO I = 1,NumMethods
   !Temp1 = BackVectors ! The number of iterations used to solve the problem
   Temp1 = MIN(MyIterCnt(I),BackVectors) ! The number of iterations used to solve the problem
   !WRITE(6,*) Temp1,MyMemory(I)
   MyMemory(I) = MyMemory(I) + Temp2               & ! Solution vector
                             + Temp2*Temp1*2.d0    & ! Basis and PC basis
                             + Temp1*(Temp1+1.0d0) & ! Hessenberg
                             + (Temp1+1.0d0)*3.d0    ! Givens and Modified RHS
   J = MyIterCnt(I)/BackVectors ! The number of times we restarted in FGMRES
   K = MyIterCnt(I) - J*BackVectors ! The number of iterations past the last restart used
   Temp5 = 0.0d0
   DO L = 1,K
      Temp5 = Temp5 + L*2.0d0 + 2.0d0 ! The number of dot products and scaling operations needed for each inner in FGMRES
   END DO
   Temp6 = J  ! The number of outers
   EstimatedFlops(I) = 2.d0*Temp3 + MyFlopCnt(I)                        &    ! The work performed upon entry
                     + 1.d0*Temp3                                       &    ! The vector scaling before the outer loop
                     + (Temp6+1.0d0)*(BackVectors*(BackVectors+1.0d0)   &    ! The hessian matrix after the outer loop
                                     + Temp3*BackVectors                &    ! The Modified RHS scaling of each vector
                                     + MyFlopCnt(I))                    &    ! The Apply A outside of the outer
                     +     Temp6*MyFlopCnt(I) + Temp4                   &    ! The Apply A inside the inner and cost of the inner work
                     + (K*1.0d0)*MyFlopCnt(I) + Temp5                        ! The Apply A inside the last pass through the inner
END DO

! -------------------------------------------------------------------------------------------------
ddtemp = 1.0d0
Temp1 = 1.0d0 / (1024.d0*1024.d0*1024.d0) ! Convert to Giga flops
#ifdef WITHPAPI
   ddtemp = PP_NumSegments
   ddtemp = 1.0d0 / ddtemp
#endif
DO I = 1,NumMethods
   TempSc1(I)= MyMemory(I)/131072.d0 ! * 8.0d0 / 1024d0/1024d0
   TempSc2(I)= EstimatedFlops(I)*Temp1  ! In Gflops
   TempSc3(I)= MyTime(I)*ddtemp
   TempSc4(I)= TempSc2(I)/(TempSc3(I)+1.d-24)
   TempSc5(I)= TempSc2(I)/TempSc1(I)*1024.d0 ! Gflops/GByte
END DO
! Say something about the matrix assembly
IF (AssemblyTime .GT. 0.0d0) THEN
   !TempSc3(3) = TempSc3(3) + AssemblyTime
   WRITE(Output_Unit,300) AssemblyTime,AssemblyTime / MyTime(3) * MyIterCnt(3)
300 FORMAT('[SN-KERNEL] Assembly Time =',F13.6,' seconds or ',F13.1,' FGMRES iterations ')
END IF

#ifdef WITHPAPI
                        !123456789*  123456789*123  123456789*123  123456789*123  123456789*123  123456789*123
2000 FORMAT('[SN-KERNEL] Method           Time       Est. Mem(MB)   Est. GFlops   Est. GFlops/s   GFlops/GByte  ',&
      'PAPI->',30(2X,A12,2X))')
2001 FORMAT('=SN-KERNEL] ',A10,2X,F13.5,2X,F13.3,2X,F13.4,2X,F13.3,2X,F13.2,2X,30I16)
WRITE(Output_Unit,2000) (PP_AllEventNames(I),I=1,PP_AllCodes)
DO J = ii_start,ii_end
   WRITE(Output_Unit,2001) MyHPMname(J),TempSc3(J),TempSc1(J),TempSc2(J),TempSc4(J),TempSc5(J),(PP_values(I,J),I=1,PP_AllCodes)
END DO
#else
                        !123456789*  123456789*123  123456789*123  123456789*123  123456789*123  123456789*123
2000 FORMAT('[SN-KERNEL] Method           Time       Est. Mem(MB)   Est. GFlops   Est. GFlops/s   GFlops/GByte')
2001 FORMAT('=SN-KERNEL] ',A10,2X,F13.5,2X,F13.3,2X,F13.4,2X,F13.3,2X,F13.2)
WRITE(Output_Unit,2000)
DO I = ii_start,ii_end
   WRITE(Output_Unit,2001) MyHPMname(I),TempSc3(I),TempSc1(I),TempSc2(I),TempSc4(I),TempSc5(I)
END DO

#endif

END SUBROUTINE PrintSummary

