! This program creates a mesh of tetradrons by building a structured grid
! ----------------------------------------------------------------------------
! Element Type: Linear Tetrahedron
! Grid:         Structured grid where each block is composed of 12 tetrahedrons
PROGRAM Driver
! Module inclusions
USE CommonBlock
IMPLICIT NONE
! preprocessing inclusions
#include "PROTEUS_Preprocess.h"

! Key Problem Size Variables - Specified via user command line input
PROTEUS_Int :: Input_GridX           ! Structured grid cells in X  
PROTEUS_Int :: Input_GridY           ! Structured grid cells in Y
PROTEUS_Int :: Input_GridZ           ! Structured grid cells in Z
PROTEUS_Int :: Input_BCmodel         ! The boundary condition model to assume (0-3, 0=no b.c.s  1= X surface bc   2= X&Y surface bc  3=X&Y&Z surface bc)
! Local
CHARACTER*32 :: Some_String
PROTEUS_Int ReturnedError,NumArguments,I

! Read command line input
NumArguments = COMMAND_ARGUMENT_COUNT()   ! # User defined command line arguments
IF ((NumArguments .LT. 4) .OR. (NumArguments .GT. 4)) THEN
200 FORMAT('[SN-KERNEL]',109('.'))
499 FORMAT('[SN-KERNEL] The list of arguments was incomplete.....................',51('.'))
500 FORMAT('[SN-KERNEL] Version 1.0 Linear tetrahedral mesh generation...........',51('.'))
502 FORMAT('[SN-KERNEL] Usage:   makemesh.x   X    Y    Z  BCmodel...............',51('.'))
503 FORMAT('[SN-KERNEL] Example: makemesh.x  10   20   10  0      ...............',51('.'))
509 FORMAT('[SN-KERNEL] X, Y, Z specifies the number of X-Y-Z meshes (E=X*Y*Z*12)',51('.'))
508 FORMAT('[SN-KERNEL] BC 0-6  specifies how many box surfaces have vacuum b.c.s',51('.'))
   WRITE(Output_Unit,200)
   WRITE(Output_Unit,499)
   WRITE(Output_Unit,500)
   WRITE(Output_Unit,200)
   WRITE(Output_Unit,502)
   WRITE(Output_Unit,503)
   WRITE(Output_Unit,200)
   WRITE(Output_Unit,509)
   WRITE(Output_Unit,508)
   WRITE(Output_Unit,200)
   STOP
ELSE
   WRITE(Output_Unit,200)
   WRITE(Output_Unit,500)
   WRITE(Output_Unit,502)
   WRITE(Output_Unit,503)
   WRITE(Output_Unit,200)
   ! Assign command line input to variables
   CALL GET_COMMAND_ARGUMENT(1,Some_String)
   READ(Some_String,*,IOSTAT=ReturnedError) Input_GridX
   CALL GET_COMMAND_ARGUMENT(2,Some_String)
   READ(Some_String,*,IOSTAT=ReturnedError) Input_GridY
   CALL GET_COMMAND_ARGUMENT(3,Some_String)
   READ(Some_String,*,IOSTAT=ReturnedError) Input_GridZ
END IF

IF (Input_GridX   .LT. 1) Input_GridX   = 1
IF (Input_GridY   .LT. 1) Input_GridY   = 1
IF (Input_GridZ   .LT. 1) Input_GridZ   = 1

! Inform the user as to the mesh we will be building
WRITE(Output_Unit,'("[SN-KERNEL] Number of X meshes      ",I8)') Input_GridX
WRITE(Output_Unit,'("[SN-KERNEL] Number of Y meshes      ",I8)') Input_GridY
WRITE(Output_Unit,'("[SN-KERNEL] Number of Z meshes      ",I8)') Input_GridZ
WRITE(Output_Unit,'("[SN-KERNEL] Number of BC surfaces   ",I8)') Input_BCmodel

! These are geometrical formulas related to linear tetrahedrons inside a structured mesh
NumElements = Input_GridX*Input_GridY*Input_GridZ*12            ! Total number of tetrahedral elements (12 per structured grid)
NumVertices = (Input_GridX+1)*(Input_GridY+1)*(Input_GridZ+1) &
            + Input_GridX*Input_GridY*Input_GridZ ! Total number of vertices in structured grid
I =     Input_GridX*Input_GridY   ! Number of structured grid surfaces on an XY external (vacuum) boundary
I = I + Input_GridX*Input_GridZ   ! Number of structured grid surfaces on an XZ external (vacuum) boundary
I = I + Input_GridY*Input_GridZ   ! Number of structured grid surfaces on an YZ external (vacuum) boundary
NumVacuum      = I*4              ! Number of tetrahedral surfaces (triangles) with vacuum boundary (2 tri/surf and 2 global surfs)

WRITE(Output_Unit,'("[SN-KERNEL] Number of Elements    ",I8)') NumElements
WRITE(Output_Unit,'("[SN-KERNEL] Number of Vertices    ",I8)') NumVertices

! Allocate matrices without checking as we should not be working with ridiculous amounts of memory per node
CALL CommonBlock_Allocate()

! Build a mesh with boundary conditions
WRITE(Output_Unit,'("[SN-KERNEL] Building the Tetrahedral mesh ")')
CALL Build_Grid_Tet(Input_BCmodel,Input_GridX,Input_GridY,Input_GridZ,  &
                    NumVertices,NumElements,NumVacuum,                  &
                    GlobalXYZ,Conn,ElementsWithVacuumBCs,LocalSurfaceIndex)

WRITE(Output_Unit,'("[SN-KERNEL] Exporting the Tetrahedral mesh ")')
CALL Export_Mesh(NumVertices,NumElements,NumVacuum,     &
                 GlobalXYZ,Conn,ElementsWithVacuumBCs,LocalSurfaceIndex)

CALL CommonBlock_Deallocate()

END PROGRAM Driver
