!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine stencils the spatially connected FE system
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE StencilNZmatrix(Input_Scheme)
USE CommonBlock
IMPLICIT NONE
!#define Local_Debug
!#define Debug_DumpAssembledMatrix
!#define Local_DumpDebugNZS
!#define Local_DebugNZS

PROTEUS_Int :: Input_Scheme
! Local variables
PROTEUS_Int I,J,K,L,M
PROTEUS_Int I1,I2,I3,I4
PROTEUS_Int Istart,IBatches,II,FixedBatchSize
! Arrays needed to evaluate the derivatives
PROTEUS_Log,ALLOCATABLE     :: Local_VertexFlag(:,:)
PROTEUS_Int,ALLOCATABLE     :: Temp_NZS_ColNum(:)

IF ((Input_Scheme .EQ. 0) .OR. (Input_Scheme .EQ. 3)) THEN
   FixedBatchSize = 30
   ALLOCATE(Local_VertexFlag(NumVertices,FixedBatchSize))
   ALLOCATE(NZS_RowLoc(NumVertices+1))
   NZS_RowLoc(1) = 1

   ! NON-THREADABLE SECTION (or of questionable value)
   ! We need to stencil the spatial matrix and do so by sweeping the connectivity matrix Batch times
   IBatches = NumVertices/FixedBatchSize
   IF (IBatches * FixedBatchSize .NE. NumVertices) IBatches = IBatches + 1
   Istart    = 0
   DO I = 1,IBatches
      IF (Istart+FixedBatchSize .GT. NumVertices) FixedBatchSize = NumVertices - Istart
      DO II = 1,FixedBatchSize
         DO J = 1,NumVertices
            Local_VertexFlag(J,II) = .FALSE.
         END DO
      END DO
      ! Flag all of the vertices that are connected
      DO K = 1,NumElements
         DO L = 1,FEVertices
            DO II = 1,FixedBatchSize
               IF (Conn(L,K) .EQ. Istart+II) THEN
                  DO M = 1,FEVertices
                     Local_VertexFlag(Conn(M,K),II) = .TRUE.
                  END DO
               END IF
            END DO
         END DO
      END DO
      ! Update NZS_NonZeros and store the rowloc information
      K = NZS_NonZeros
      DO II = 1,FixedBatchSize
         DO J = 1,NumVertices
            IF (Local_VertexFlag(J,II)) K = K + 1
         END DO
         !WRITE(Output_Unit,*)'Row ',Istart+II,' K=',K
         NZS_RowLoc(Istart+II+1) = K+1 ! Start of the next row
      END DO
      ! Resize ColNum
      IF (NZS_NonZeros .EQ. 0) THEN
         ALLOCATE(NZS_ColNum(K))
      ELSE
         ALLOCATE(Temp_NZS_ColNum(NZS_NonZeros))
         DO J = 1,NZS_NonZeros
            Temp_NZS_ColNum(J) = NZS_ColNum(J)
         END DO
         DEALLOCATE(NZS_ColNum)
         ALLOCATE(NZS_ColNum(K))
         DO J = 1,NZS_NonZeros
            NZS_ColNum(J) = Temp_NZS_ColNum(J)
         END DO
         DEALLOCATE(Temp_NZS_ColNum)
      END IF
      ! Store ColNum for this set of data
      K = NZS_NonZeros
      DO II = 1,FixedBatchSize
         DO J = 1,NumVertices
            IF (Local_VertexFlag(J,II)) THEN
               K = K + 1
               NZS_ColNum(K) = J
            END IF
         END DO
      END DO
      NZS_NonZeros = K
      Istart = Istart + FixedBatchSize
   END DO ! End of Batch of vertices
   ALLOCATE(NZS_Data(NumAngles,NZS_NonZeros))
   DO I = 1,NZS_NonZeros
      DO J = 1,NumAngles
         NZS_Data(J,I) = 0.0d0
      END DO
   END DO

   DEALLOCATE(Local_VertexFlag)
ELSE
   ALLOCATE(NZS_RowLoc(1),NZS_ColNum(1),NZS_Data(1,1))
END IF ! Scheme 0 and 3 

END SUBROUTINE StencilNZmatrix
