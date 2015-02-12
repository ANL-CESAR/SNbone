!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine creates a structured grid of tetrahedral elements 
! The connectivity pattern that is created should be comparable to an unstructured FE based access pattern
!
! List of 9 vertices in a given box:     List of 12 tetrahedrons inside box, each described by 4 vertices:
!       0  0  0 -> Vertex 1 coords                1 2 4 9     -> The first tetrahedron has vertices 1,2,4, and 5.
!       1  0  0 -> Vertex 2 coords                1 4 3 9     -> The second tetrahedron has vertices 1 4 3 9
!       0  1  0 ...                               5 7 8 9     ...
! XYZ=  1  1  0                            Conn=  5 8 6 9          +Y            +X            +Z
!       0  0  1                                   1 6 2 9     5 - 6 7 - 8   7 - 5 8 - 6   1 - 2  5 - 6
!       1  0  1                                   1 5 6 9     |   | |   |   |   | |   |   |   |  |   |
!       0  1  1                                   3 4 8 9     1 - 2 3 - 4   3 - 1 4 - 2   3 - 4  7 - 8
!       1  1  1                                   3 8 7 9     
!      .5 .5 .5                                   3 5 1 9     Elements 1- 2 are on -Z   3- 4 are on +Z
!                                                 3 7 5 9              5- 6 are on -Y   7- 8 are on +Y
!                                                 4 2 6 9              9-10 are on -X  11-12 are on +X
!                                                 4 6 8 9
! 
! Each structured X-Y-Z position corresponds to 12 tetrahedrons.
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE Build_Grid_Tet(Input_BCmodel,Input_GridX,Input_GridY,Input_GridZ,           &
                          NumVertices,NumElements,NumVacuum,                           &
                          GlobalXYZ,Conn,ElementsWithVacuumBCs,LocalSurfaceIndex)
IMPLICIT NONE 
!#define Local_Debug
!#define Local_ExportVTK
! Parameters
PROTEUS_Int :: Input_BCmodel
PROTEUS_Int :: Input_GridX
PROTEUS_Int :: Input_GridY
PROTEUS_Int :: Input_GridZ
PROTEUS_Int :: NumVertices
PROTEUS_Int :: NumElements
PROTEUS_Int :: NumVacuum
! Arrays to fill
PROTEUS_Real    GlobalXYZ(NumVertices,3)
PROTEUS_Int     Conn(4,NumElements)
PROTEUS_Int     ElementsWithVacuumBCs(NumVacuum)
PROTEUS_Int     LocalSurfaceIndex(NumVacuum)
! Local
PROTEUS_Int Element,Surface,I,J,K,L,III_1,III_3,III_5,III_7,III_9,III,JJJ
PROTEUS_Int :: File_Unit=10
PROTEUS_Real gScaleX,gScaleY,gScaleZ
PROTEUS_Real dXX,dYY,dZZ,rXX,rYY,rZZ

! Connectivity is correct as far as I can tell
gScaleX = 1.0d0/(Input_GridX+0.0d0)
gScaleY = 1.0d0/(Input_GridY+0.0d0)
gScaleZ = 1.0d0/(Input_GridZ+0.0d0)
III_1 = 1                                       ! The index position of the 1 vertices
III_3 = III_1 + Input_GridX + 1                 ! The index position of the 3 vertices
III_5 = III_1 + (Input_GridX+1)*(Input_GridY+1) ! The index position of the 5 vertices
III_7 = III_5 + Input_GridX + 1                 ! The index posiiton of the 7 vertices
III_9 = (Input_GridX+1)*(Input_GridY+1)*(Input_GridZ+1) + 1 ! The index position of the 9 vertices
Element = 0
DO K = 1,Input_GridZ
   DO J = 1,Input_GridY
      DO I = 1,Input_GridX
         Conn(1,Element+1) = III_1+0;Conn(2,Element+1) = III_1+1;Conn(3,Element+1) = III_3+1;Conn(4,Element+1) = III_9;  ! 1 2 4 9
         Conn(1,Element+2) = III_1+0;Conn(2,Element+2) = III_3+1;Conn(3,Element+2) = III_3+0;Conn(4,Element+2) = III_9;  ! 1 4 3 9
         Conn(1,Element+3) = III_5+0;Conn(2,Element+3) = III_7+0;Conn(3,Element+3) = III_7+1;Conn(4,Element+3) = III_9;  ! 5 7 8 9
         Conn(1,Element+4) = III_5+0;Conn(2,Element+4) = III_7+1;Conn(3,Element+4) = III_5+1;Conn(4,Element+4) = III_9;  ! 5 8 6 9
         Conn(1,Element+5) = III_1+0;Conn(2,Element+5) = III_5+1;Conn(3,Element+5) = III_1+1;Conn(4,Element+5) = III_9;  ! 1 6 2 9
         Conn(1,Element+6) = III_1+0;Conn(2,Element+6) = III_5+0;Conn(3,Element+6) = III_5+1;Conn(4,Element+6) = III_9;  ! 1 5 6 9
         Conn(1,Element+7) = III_3+0;Conn(2,Element+7) = III_3+1;Conn(3,Element+7) = III_7+1;Conn(4,Element+7) = III_9;  ! 3 4 8 9
         Conn(1,Element+8) = III_3+0;Conn(2,Element+8) = III_7+1;Conn(3,Element+8) = III_7+0;Conn(4,Element+8) = III_9;  ! 3 8 7 9
         Conn(1,Element+9) = III_3+0;Conn(2,Element+9) = III_5+0;Conn(3,Element+9) = III_1+0;Conn(4,Element+9) = III_9;  ! 3 5 1 9
         Conn(1,Element+10)= III_3+0;Conn(2,Element+10)= III_7+0;Conn(3,Element+10)= III_5+0;Conn(4,Element+10)= III_9;  ! 3 7 5 9
         Conn(1,Element+11)= III_3+1;Conn(2,Element+11)= III_1+1;Conn(3,Element+11)= III_5+1;Conn(4,Element+11)= III_9;  ! 4 2 6 9
         Conn(1,Element+12)= III_3+1;Conn(2,Element+12)= III_5+1;Conn(3,Element+12)= III_7+1;Conn(4,Element+12)= III_9;  ! 4 6 8 9
         ! There are only 8 unique vertices per pass
         GlobalXYZ(III_1+0,1) = (I-1)*gScaleX;GlobalXYZ(III_1+0,2) = (J-1)*gScaleY;GlobalXYZ(III_1+0,3) = (K-1)*gScaleZ;
         GlobalXYZ(III_1+1,1) =  I   *gScaleX;GlobalXYZ(III_1+1,2) = (J-1)*gScaleY;GlobalXYZ(III_1+1,3) = (K-1)*gScaleZ;
         GlobalXYZ(III_3+0,1) = (I-1)*gScaleX;GlobalXYZ(III_3+0,2) =  J   *gScaleY;GlobalXYZ(III_3+0,3) = (K-1)*gScaleZ;
         GlobalXYZ(III_3+1,1) =  I   *gScaleX;GlobalXYZ(III_3+1,2) =  J   *gScaleY;GlobalXYZ(III_3+1,3) = (K-1)*gScaleZ;
         GlobalXYZ(III_5+0,1) = (I-1)*gScaleX;GlobalXYZ(III_5+0,2) = (J-1)*gScaleY;GlobalXYZ(III_5+0,3) =  K   *gScaleZ;
         GlobalXYZ(III_5+1,1) =  I   *gScaleX;GlobalXYZ(III_5+1,2) = (J-1)*gScaleY;GlobalXYZ(III_5+1,3) =  K   *gScaleZ;
         GlobalXYZ(III_7+0,1) = (I-1)*gScaleX;GlobalXYZ(III_7+0,2) =  J   *gScaleY;GlobalXYZ(III_7+0,3) =  K   *gScaleZ;
         GlobalXYZ(III_7+1,1) =  I   *gScaleX;GlobalXYZ(III_7+1,2) =  J   *gScaleY;GlobalXYZ(III_7+1,3) =  K   *gScaleZ;
         GlobalXYZ(III_9  ,1) =(I-.5)*gScaleX;GlobalXYZ(III_9  ,2) =(J-.5)*gScaleY;GlobalXYZ(III_9  ,3) =(K-.5)*gScaleZ;
         ! We increment the counters
         Element = Element + 12
         III_1 = III_1 + 1
         III_3 = III_3 + 1
         III_5 = III_5 + 1
         III_7 = III_7 + 1
         III_9 = III_9 + 1
      END DO
      III_1 = III_1 + 1
      III_3 = III_3 + 1
      III_5 = III_5 + 1
      III_7 = III_7 + 1
      !III_9 = III_9 + 1 This is unnecessary as they will go sequentially
   END DO
   III_1 = III_1 + 1 + Input_GridX
   III_3 = III_3 + 1 + Input_GridX
   III_5 = III_5 + 1 + Input_GridX
   III_7 = III_7 + 1 + Input_GridX
   !III_9 = III_9 + 1 + Input_GridX This is unnecessary
END DO

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
   ! Create a VTK file
#endif

#ifdef Local_ExportVTK
   ! Open the input file for writing
   OPEN(UNIT=File_Unit,IOSTAT=I,FILE='f_ApplyA.vtk',ACCESS='SEQUENTIAL',FORM='FORMATTED')
   WRITE(File_Unit,'("# vtk DataFile Version 3.0")')                          ! Generic junk about the file type
   WRITE(File_Unit,'("ANL NE devision neutronics export routine 08-2006")')   ! A title description of the type of output
   WRITE(File_Unit,'("ASCII")')                                               ! This will be an ascii export option
   WRITE(File_Unit,'("DATASET UNSTRUCTURED_GRID")')                           ! Indicates that we are defining a unstructured grid
   WRITE(File_Unit,'("POINTS ",I12," float")')  NumVertices                   ! The number of mesh points which are of type float
   DO I = 1,NumVertices
      !RadiusXY  = DSQRT(0.3d0+GlobalXYZ(I,1)*GlobalXYZ(I,1)+GlobalXYZ(I,2)*GlobalXYZ(I,2))+1.0D-15
      !RadiusXYZ = DSQRT(0.4d0+GlobalXYZ(I,1)*GlobalXYZ(I,1)+GlobalXYZ(I,2)*GlobalXYZ(I,2)+GlobalXYZ(I,3)*GlobalXYZ(I,3))+1.0D-15
      !Phi = ACOS((GlobalXYZ(I,1)+0.1d0)/RadiusXY)
      !Theta = ACOS((GlobalXYZ(I,3)*0.8d0+0.2d0)/RadiusXYZ)
      !IF (GlobalXYZ(I,2) .LT. 0.1d0) RadiusXYZ = RadiusXYZ - 0.1d0
      !IF (GlobalXYZ(I,3) .GT. 0.5d0) RadiusXYZ = RadiusXYZ + 0.2d0
      !WRITE(File_Unit,'(3(1PE12.5,1X))') RadiusXYZ*((COS(Phi))**0.4d0)*SIN(Theta),  &
      !                                   RadiusXYZ*((SIN(Phi*1.5d0))**0.6d0)*((COS(Phi*0.5d0))**0.1d0)*SIN(Theta),  &
      !                                   RadiusXYZ*(COS(Theta))**0.4d0
      ! We use a FE transformation
      !dXX = (GlobalXYZ(I,1) - 0.5d0)*2.0d0;dYY = (GlobalXYZ(I,2) - 0.5d0)*2.0d0;dZZ = (GlobalXYZ(I,3) - 0.5d0)*2.0d0;
      !CALL f_1_Transform_Hex20(rXX,rYY,rZZ,dXX,dYY,dZZ)
      !WRITE(File_Unit,'(3(1PE12.5,1X))') rXX,rYY,rZZ
      WRITE(File_Unit,'(3(1PE12.5,1X))') GlobalXYZ(I,1),GlobalXYZ(I,2),GlobalXYZ(I,3)
   END DO

   WRITE(File_Unit,*) ! Blank line

   WRITE(File_Unit,'("CELLS ",I12,1X,I12)') NumElements,NumElements*5

   DO Element = 1,NumElements
      WRITE(File_Unit,'(30I12)') 4,Conn(1,Element)-1,Conn(2,Element)-1,Conn(3,Element)-1,Conn(4,Element)-1
   END DO

   WRITE(File_Unit,*) ! Blank line

   WRITE(File_Unit,'("CELL_TYPES ",I12)') NumElements

   WRITE(File_Unit,'(40I3)') (10, Element=1,NumElements)

   WRITE(File_Unit,*) ! Blank line

   WRITE(File_Unit,'("CELL_DATA ",I12)') NumElements
   WRITE(File_Unit,'("FIELD MeshInformation 1")')
   WRITE(File_Unit,'("Base_Tet_ID 1 ",I12," float")') NumElements
   WRITE(File_Unit,'(10(F8.3,1X))') ( (I,I=1,12),J=1,Input_GridX*Input_GridY*Input_GridZ)

   CLOSE(File_Unit)
#endif

! 03-25-2013 ---Verified these to assign correct ElementsWithVacuumBCS
! 03-25-2013 ---Not sure about correct LocalSurfaceIndex but Mike said it doesn't matter;
! 03-25-2013 --- Also, the SurfaceAxB data are made-up
Surface = 0
! Boundary conditions on Z-plane surfaces (top and bottom)
DO J = 1,Input_GridY
   DO I = 1,Input_GridX
      III = ((1          -1)*Input_GridY*Input_GridX + (J-1)*Input_GridX + I)*12 - 3    ! Lower Z Surface
      JJJ = ((Input_GridZ-1)*Input_GridY*Input_GridX + (J-1)*Input_GridX + I)*12 - 1    ! Upper Z Surface
      ! Lower surface
      IF (Input_BCmodel .GE. 3) THEN
         Surface = Surface + 1 ! Vacuum b.c.      
         ElementsWithVacuumBCs(Surface) = III     ! First triangle on face
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0  ! made-up
         !Vac_Normals(1,Surface)=-0.999d0;Vac_Normals(2,Surface)=0.001d0;Vac_Normals(3,Surface)=0.002d0;
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = III+1   ! Second triangle on face
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=-0.999d0;Vac_Normals(2,Surface)=0.001d0;Vac_Normals(3,Surface)=0.002d0;
      END IF
      ! Upper surface
      IF (Input_BCmodel .GE. 6) THEN
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = JJJ
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=0.999d0;Vac_Normals(2,Surface)=0.001d0;Vac_Normals(3,Surface)=0.002d0;
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = JJJ+1
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=0.999d0;Vac_Normals(2,Surface)=0.001d0;Vac_Normals(3,Surface)=0.002d0;
      END IF
   END DO
END DO

! Boundary conditions on Y-plane surfaces (Left and Right)
DO K = 1,Input_GridZ
   DO I = 1,Input_GridX
      III = ((K-1)*Input_GridY*Input_GridX + (1          -1)*Input_GridX + I)*12 - 7
      JJJ = ((K-1)*Input_GridY*Input_GridX + (Input_GridY-1)*Input_GridX + I)*12 - 5
      ! Front Surface
      IF (Input_BCmodel .GE. 2) THEN
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = III
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=-0.999d0;Vac_Normals(2,Surface)=-0.001d0;Vac_Normals(3,Surface)=0.002d0;
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = III+1
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=-0.999d0;Vac_Normals(2,Surface)=-0.001d0;Vac_Normals(3,Surface)=0.002d0;
      END IF
      ! Back surface
      IF (Input_BCmodel .GE. 5) THEN
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = JJJ
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=0.999d0;Vac_Normals(2,Surface)=0.001d0;Vac_Normals(3,Surface)=0.002d0;
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = JJJ+1
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=0.999d0;Vac_Normals(2,Surface)=0.001d0;Vac_Normals(3,Surface)=0.002d0;
      END IF
   END DO
END DO

! Boundary conditions on X-plane surfaces (Left and right)
DO K = 1,Input_GridZ
   DO J = 1,Input_GridY
      III = ((K-1)*Input_GridY*Input_GridX + (J-1)*Input_GridX +           1)*12 - 11
      JJJ = ((K-1)*Input_GridY*Input_GridX + (J-1)*Input_GridX + Input_GridX)*12 -  9    
      ! Left surface
      IF (Input_BCmodel .GE. 1) THEN
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = III
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=-0.999d0;Vac_Normals(2,Surface)=-0.001d0;Vac_Normals(3,Surface)=-0.002d0;
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = III+1
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=-0.999d0;Vac_Normals(2,Surface)=-0.001d0;Vac_Normals(3,Surface)=-0.002d0;
      END IF
      ! Upper surface
      IF (Input_BCmodel .GE. 4) THEN
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = JJJ
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=0.999d0;Vac_Normals(2,Surface)=0.001d0;Vac_Normals(3,Surface)=0.002d0;
         Surface = Surface + 1 ! Vacuum b.c.
         ElementsWithVacuumBCs(Surface) = JJJ+1
         LocalSurfaceIndex(Surface)     = 4
         !IndexNormal(Surface)           = Surface
         !SurfaceAxB(Surface)            = 1.0d0 + I*0.001d0 + J*0.002d0
         !Vac_Normals(1,Surface)=0.999d0;Vac_Normals(2,Surface)=0.001d0;Vac_Normals(3,Surface)=0.002d0;
      END IF
   END DO
END DO

NumVacuum = Surface

END SUBROUTINE Build_Grid_Tet
