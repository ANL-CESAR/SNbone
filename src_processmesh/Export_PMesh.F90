!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine exports the processed tetrahedral mesh
! In includes the original mesh and the spatial matrices
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE Export_PMesh()
USE CommonBlock
IMPLICIT NONE 
! Local
PROTEUS_Int Element,Vertex,File_Unit,IOS,I,J,K


File_Unit = 777

200 FORMAT(10I8)
210 FORMAT(1P,10(E13.6,1X))
OPEN(UNIT=File_Unit,IOSTAT=IOS,FILE='pmesh.ascii')
!WRITE(File_Unit,'("# Version 1.0 file format for CFE SN miniapp")')
WRITE(File_Unit,200) NumVertices,NumElements,NumVacuum,NumUnitNormals,  &
                     NumThreads,AS_NumColors,FEVertices,FEGaussPoints,FENumDim

! Threading setup information
DO I = 1,NumThreads
   WRITE(File_Unit,200) DA_ThreadWiseWork(1,I),DA_ThreadWiseWork(2,I),MM_ThreadWiseWork(1,I),MM_ThreadWiseWork(2,I)
   DO J = 1,AS_NumColors
      WRITE(File_Unit,200) AS_ThreadWiseWork(1,J,I),AS_ThreadWiseWork(2,J,I)
   END DO
END DO
DO Vertex = 1,NumVertices
   WRITE(File_Unit,210) GlobalXYZ(Vertex,1),GlobalXYZ(Vertex,2),GlobalXYZ(Vertex,3)
END DO
DO I = 1,NumVertices
   WRITE(File_Unit,200) VertexLocalToGlobal(I)
END DO
DO Element = 1,NumElements
   WRITE(File_Unit,200) Conn(1,Element),Conn(2,Element),Conn(3,Element),Conn(4,Element),&
                        BCInfo(1,Element),BCInfo(2,Element),ElementLocalToGlobal(Element)
END DO
DO Element = 1,NumElements
   WRITE(File_Unit,210) ConstTau(Element),ConstF(Element),ConstU(Element),ConstUT(Element),fcoef(Element)
   WRITE(File_Unit,210) pcoef( 1:6,Element)
   WRITE(File_Unit,210) pcoef( 7:12,Element)
   WRITE(File_Unit,210) pcoef(13:18,Element)
   WRITE(File_Unit,210) pcoef(19:24,Element)
   WRITE(File_Unit,210) pcoef(25:30,Element)
   WRITE(File_Unit,210) pcoef(31:36,Element)
   WRITE(File_Unit,210) ucoef(1:6,Element)
   WRITE(File_Unit,210) ucoef(7:9,Element)
END DO
DO I = 1,NumVacuum
   WRITE(File_Unit,200) LocalSurfaceIndex(I),IndexNormal(I)
END DO
DO I = 1,NumVacuum
   WRITE(File_Unit,210) Vac_Normals(1,I),Vac_Normals(2,I),Vac_Normals(3,I),SurfaceAxB(I)
END DO
DO I = 1,FEGaussPoints
   WRITE(File_Unit,210) (FEShapeFunctions(J,I),J=1,FEVertices)
END DO
DO Element = 1,NumElements
   DO I = 1,FEGaussPoints
      WRITE(File_Unit,210) FEDetJacandWgt(I,Element)
      DO K = 1,FENumDim
         WRITE(File_Unit,210) (FEDerivatives(J,K,I,Element),J=1,FEVertices)
      END DO
   END DO
END DO

CLOSE(File_Unit)

END SUBROUTINE Export_PMesh
