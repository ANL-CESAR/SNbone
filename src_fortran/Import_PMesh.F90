!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine imports the processed tetrahedral mesh
! In includes the original mesh and the spatial matrices
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
!#define Local_Debug
SUBROUTINE Import_PMesh()
USE CommonBlock
IMPLICIT NONE 
! Local
PROTEUS_Int Element,Vertex,File_Unit,IOS,I,J,K

File_Unit = 777

200 FORMAT(10I8)
210 FORMAT(1P,10(E13.6,1X))
OPEN(UNIT=File_Unit,IOSTAT=IOS,STATUS='OLD',FILE='pmesh.ascii')
IF (IOS .NE. 0) THEN
   WRITE(Output_Unit,*,IOSTAT=IOS)'Could not open inputmesh.ascii'
   CALL Abort
END IF
!READ(File_Unit,*,IOSTAT=IOS)
!IF (IOS .NE. 0) THEN
!   WRITE(Output_Unit,*)'Read failure of ascii string'
!   CALL Abort
!END IF
READ(File_Unit,*,IOSTAT=IOS) NumVertices,NumElements,NumVacuum,NumUnitNormals,   &
                  AS_NumThreads,AS_NumColors,FEVertices,FEGaussPoints,FENumDim
IF (IOS .NE. 0) THEN
   WRITE(Output_Unit,*)'Read failure of control data'
   CALL Abort
END IF

IF (AS_NumThreads .LT. NumThreads) THEN
   WRITE(Output_Unit,'("Fatal error as processed mesh ",I8," is not prepared for ",I8," threads")') AS_NumThreads,NumThreads
   CALL Abort
END IF

TasksPerThread = AS_NumThreads/NumThreads
IF (TasksPerThread*NumThreads .NE. AS_NumThreads) THEN
   WRITE(Output_Unit,'("Fatal error as processed mesh ",I8," is not divisible by ",I8," threads")') AS_NumThreads,NumThreads
   CALL Abort
END IF

CALL CommonBlock_Allocate()

! Threading setup information
DO I = 1,AS_NumThreads
   READ(File_Unit,*,IOSTAT=IOS) DA_ThreadWiseWork(1,I),DA_ThreadWiseWork(2,I),MM_ThreadWiseWork(1,I),MM_ThreadWiseWork(2,I)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of DA and MM thread controls'
      CALL Abort
   END IF
   DO J = 1,AS_NumColors
      READ(File_Unit,*,IOSTAT=IOS) AS_ThreadWiseWork(1,J,I),AS_ThreadWiseWork(2,J,I)
      IF (IOS .NE. 0) THEN
         WRITE(Output_Unit,*)'Read failure of AS thread controls'
         CALL Abort
      END IF
   END DO
END DO
DO Vertex = 1,NumVertices
   READ(File_Unit,*,IOSTAT=IOS) GlobalXYZ(Vertex,1),GlobalXYZ(Vertex,2),GlobalXYZ(Vertex,3)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of XYZ data'
      CALL Abort
   END IF
END DO
DO I = 1,NumVertices
   READ(File_Unit,*,IOSTAT=IOS) VertexLocalToGlobal(I)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of vertex local to global'
      CALL Abort
   END IF
END DO
DO Element = 1,NumElements
   READ(File_Unit,*,IOSTAT=IOS) Conn(1,Element),Conn(2,Element),Conn(3,Element),Conn(4,Element),&
                     BCInfo(1,Element),BCInfo(2,Element),ElementLocalToGlobal(Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of conn and boundary conditions'
      CALL Abort
   END IF
END DO
DO Element = 1,NumElements
   READ(File_Unit,*,IOSTAT=IOS) ConstTau(Element),ConstF(Element),ConstU(Element),ConstUT(Element),fcoef(Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of cross section data'
      CALL Abort
   END IF
   READ(File_Unit,*,IOSTAT=IOS) pcoef( 1:6,Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of pcoef1'
      CALL Abort
   END IF
   READ(File_Unit,*,IOSTAT=IOS) pcoef( 7:12,Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of pcoef2'
      CALL Abort
   END IF
   READ(File_Unit,*,IOSTAT=IOS) pcoef(13:18,Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of pcoef3'
      CALL Abort
   END IF
   READ(File_Unit,*,IOSTAT=IOS) pcoef(19:24,Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of pcoef4'
      CALL Abort
   END IF
   READ(File_Unit,*,IOSTAT=IOS) pcoef(25:30,Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of pcoef5'
      CALL Abort
   END IF
   READ(File_Unit,*,IOSTAT=IOS) pcoef(31:36,Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of pcoef6'
      CALL Abort
   END IF
   READ(File_Unit,*,IOSTAT=IOS) ucoef(1:6,Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of ucoef1'
      CALL Abort
   END IF
   READ(File_Unit,*,IOSTAT=IOS) ucoef(7:9,Element)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of ucoef2'
      CALL Abort
   END IF
END DO
DO I = 1,NumVacuum
   READ(File_Unit,*,IOSTAT=IOS) LocalSurfaceIndex(I),IndexNormal(I)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of LocalSurfaceIndex'
      CALL Abort
   END IF
END DO
DO I = 1,NumVacuum
   READ(File_Unit,*,IOSTAT=IOS) Vac_Normals(1,I),Vac_Normals(2,I),Vac_Normals(3,I),SurfaceAxB(I)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of normals'
      CALL Abort
   END IF
END DO
DO I = 1,FEGaussPoints
   READ(File_Unit,*,IOSTAT=IOS) (FEShapeFunctions(J,I),J=1,FEVertices)
   IF (IOS .NE. 0) THEN
      WRITE(Output_Unit,*)'Read failure of FE shape functions'
      CALL Abort
   END IF
END DO
DO Element = 1,NumElements
   DO I = 1,FEGaussPoints
      READ(File_Unit,*,IOSTAT=IOS) FEDetJacandWgt(I,Element)
      IF (IOS .NE. 0) THEN
         WRITE(Output_Unit,*)'Read failure of jacobian and weights'
         CALL Abort
      END IF
      DO K = 1,FENumDim
         READ(File_Unit,*,IOSTAT=IOS) (FEDerivatives(J,K,I,Element),J=1,FEVertices)
         IF (IOS .NE. 0) THEN
            WRITE(Output_Unit,*)'Read failure of derivatives'
            CALL Abort
         END IF
      END DO
   END DO
END DO

CLOSE(File_Unit)

#ifdef Local_Debug
   ! Print the connectivity matrix out
   DO Element = 1,NumElements
      WRITE(6,'("Element ",I5," vertex numbers -> ",4(I6,1X))') &
         Element, Conn(1,Element),Conn(2,Element),Conn(3,Element),Conn(4,Element)
   END DO
   DO I = 1,NumVertices
      WRITE(6,'("Vertex ",I5," positions -> ",3(F13.6,1X))') &
         I, GlobalXYZ(I,1),GlobalXYZ(I,2),GlobalXYZ(I,3)
   END DO
#endif

END SUBROUTINE Import_PMesh
