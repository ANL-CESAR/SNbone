//---------------------------------------------------------------------------------------------------------------------------------
// This subroutine imports the processed tetrahedral mesh
// In includes the original mesh and the spatial matrices
//---------------------------------------------------------------------------------------------------------------------------------
#include "ApplyA_functions.h"
#include "ApplyA_common.h"
#include <stdio.h>
#include <stdlib.h>
#ifdef WITHOMP
#include <omp.h>
#endif
// Fortran interface routines
void import_pmesh() {Import_PMesh();}
void import_pmesh_() {Import_PMesh();}
void Import_PMesh() {
//PROTEUS_Int Input_Iterations,IterationCount,iMethod
//PROTEUS_Real LHS_C(NumAngles,NumVertices),RHS_C(NumAngles,NumVertices)

// Local
FILE *fp;
int  Element,Vertex,I,J,K,itemp1,istupid;
double dstupid;
char line[45];

// Open the input file
if ((fp = fopen("pmesh.ascii", "r")) == NULL) {
     printf("[SN-KERNEL] Error: pmesh.ascii cannot be opened \n");
     abort();
   }

// The header info for this file is:
// WRITE(File_Unit,'("# Version 1.0 file format for CFE SN miniapp")')
// WRITE(File_Unit,200) NumVertices,NumElements,NumVacuum,NumUnitNormals,  &
//                     NumThreads,AS_NumColors,FEVertices,FEGaussPoints,FENumDim

// Read the character string on the first line
// fscanf(fp, "%s", line);
// Read the control characters
fscanf(fp, "%d", &NumVertices);
fscanf(fp, "%d", &NumElements);
fscanf(fp, "%d", &NumVacuum);
fscanf(fp, "%d", &NumUnitNormals);
fscanf(fp, "%d", &AS_NumThreads);
fscanf(fp, "%d", &AS_NumColors);
fscanf(fp, "%d", &FEVertices);
fscanf(fp, "%d", &FEGaussPoints);
fscanf(fp, "%d", &FENumDim);

if (AS_NumThreads < NumThreads) {
   printf("[SN-KERNEL]...Fatal error as processed mesh %d is not prepared for %d threads \n",AS_NumThreads,NumThreads);
   abort();
   }

TasksPerThread = AS_NumThreads/NumThreads;
if (TasksPerThread*NumThreads < AS_NumThreads) {
   printf("[SN-KERNEL]...Fatal error as processed mesh %d is not divisible by %d threads \n",AS_NumThreads,NumThreads);
   abort();
   }

CommonBlock(); // allocates all of the arrays we are about to read in

// Threading setup information
for (I = 1;I < AS_NumThreads+1; I++) { // DO I = 1,AS_NumThreads
   itemp1 = (I-1)*2 + 1 - 1;   
   //READ(File_Unit,*,IOSTAT=IOS) DA_ThreadWiseWork[1,I],DA_ThreadWiseWork[2,I],MM_ThreadWiseWork[1,I],MM_ThreadWiseWork[2,I];
   fscanf(fp, "%d", &istupid); DA_ThreadWiseWork[itemp1]=istupid;
   fscanf(fp, "%d", &istupid); DA_ThreadWiseWork[itemp1+1]=istupid;
   fscanf(fp, "%d", &istupid); MM_ThreadWiseWork[itemp1]=istupid;
   fscanf(fp, "%d", &istupid); MM_ThreadWiseWork[itemp1+1]=istupid;

   for (J = 1;J < AS_NumColors+1; J++) { // DO J = 1,AS_NumColors
      //READ(File_Unit,*,IOSTAT=IOS) AS_ThreadWiseWork(1,J,I),AS_ThreadWiseWork(2,J,I)
      itemp1 = (I-1)*AS_NumColors*2 + (J-1)*2 + 1 - 1;
      fscanf(fp, "%d", &istupid); AS_ThreadWiseWork[itemp1]=istupid;
      fscanf(fp, "%d", &istupid); AS_ThreadWiseWork[itemp1+1]=istupid;
      }
   }
// mesh vertices
for (Vertex = 1;Vertex < NumVertices+1; Vertex++) { // DO Vertex = 1,NumVertices
   //READ(File_Unit,*,IOSTAT=IOS) GlobalXYZ(Vertex,1),GlobalXYZ(Vertex,2),GlobalXYZ(Vertex,3)
   fscanf(fp, "%lf", &dstupid); GlobalXYZ[              Vertex-1]=dstupid;
   fscanf(fp, "%lf", &dstupid); GlobalXYZ[NumVertices  +Vertex-1]=dstupid;
   fscanf(fp, "%lf", &dstupid); GlobalXYZ[NumVertices*2+Vertex-1]=dstupid;
   }

for (I = 1;I < NumVertices+1; I++) { // DO Vertex = 1,NumVertices
   //READ(File_Unit,*,IOSTAT=IOS) VertexLocalToGlobal(I)
   fscanf(fp, "%d", &istupid); VertexLocalToGlobal[I-1]=istupid;
   }

for (Element = 1;Element < NumElements+1; Element++) { // DO Element = 1,NumElements
   //READ(File_Unit,*,IOSTAT=IOS) Conn(1,Element),Conn(2,Element),Conn(3,Element),Conn(4,Element),&
   //                             BCInfo(1,Element),BCInfo(2,Element),ElementLocalToGlobal(Element)
   itemp1 = (Element-1)*4 + 1 - 1; 
   fscanf(fp, "%d", &istupid); Conn[itemp1  ]=istupid;
   fscanf(fp, "%d", &istupid); Conn[itemp1+1]=istupid;
   fscanf(fp, "%d", &istupid); Conn[itemp1+2]=istupid;
   fscanf(fp, "%d", &istupid); Conn[itemp1+3]=istupid;
   itemp1 = (Element-1)*2 + 1 - 1; 
   fscanf(fp, "%d", &istupid); BCInfo[itemp1  ]=istupid;
   fscanf(fp, "%d", &istupid); BCInfo[itemp1+1]=istupid;
   fscanf(fp, "%d", &istupid); ElementLocalToGlobal[Element-1]=istupid;
   }

for (Element = 1;Element < NumElements+1; Element++) { // DO Element = 1,NumElements
   //READ(File_Unit,*,IOSTAT=IOS) ConstTau(Element),ConstF(Element),ConstU(Element),ConstUT(Element),fcoef(Element)
   fscanf(fp, "%lf", &dstupid); ConstTau[Element-1]=dstupid;
   fscanf(fp, "%lf", &dstupid); ConstF[Element-1]=dstupid;
   fscanf(fp, "%lf", &dstupid); ConstU[Element-1]=dstupid;
   fscanf(fp, "%lf", &dstupid); ConstUT[Element-1]=dstupid;
   fscanf(fp, "%lf", &dstupid); fcoef[Element-1]=dstupid;
   
   //READ(File_Unit,*,IOSTAT=IOS) pcoef( 1:6,Element)
   //READ(File_Unit,*,IOSTAT=IOS) pcoef( 7:12,Element)
   //READ(File_Unit,*,IOSTAT=IOS) pcoef(13:18,Element)
   //READ(File_Unit,*,IOSTAT=IOS) pcoef(19:24,Element)
   //READ(File_Unit,*,IOSTAT=IOS) pcoef(25:30,Element)
   //READ(File_Unit,*,IOSTAT=IOS) pcoef(31:36,Element)
   itemp1 = (Element-1)*36 + 1 - 1;
   for (I = 1; I<37; I++) {
      fscanf(fp, "%lf", &dstupid); pcoef[itemp1+I-1]=dstupid;
      }

   //READ(File_Unit,*,IOSTAT=IOS) ucoef(1:6,Element)
   //READ(File_Unit,*,IOSTAT=IOS) ucoef(7:9,Element)
   itemp1 = (Element-1)*9 + 1 - 1;
   for (I = 1; I<10; I++) {
      fscanf(fp, "%lf", &dstupid); ucoef[itemp1+I-1]=dstupid;
      }
   }

for (I = 1;I < NumVacuum+1; I++) { // DO I = 1,NumVacuum
   //READ(File_Unit,*,IOSTAT=IOS) LocalSurfaceIndex(I),IndexNormal(I)
   fscanf(fp, "%d", &istupid); LocalSurfaceIndex[I-1]=istupid;
   fscanf(fp, "%d", &istupid); IndexNormal[I-1]=istupid;
   }

for (I = 1;I < NumVacuum+1; I++) { // DO I = 1,NumVacuum
   //READ(File_Unit,*,IOSTAT=IOS) Vac_Normals(1,I),Vac_Normals(2,I),Vac_Normals(3,I),SurfaceAxB(I)
   itemp1 = (I-1)*3 + 1 - 1;
   fscanf(fp, "%lf", &dstupid); Vac_Normals[itemp1  ]=dstupid;
   fscanf(fp, "%lf", &dstupid); Vac_Normals[itemp1+1]=dstupid;
   fscanf(fp, "%lf", &dstupid); Vac_Normals[itemp1+2]=dstupid;
   fscanf(fp, "%lf", &dstupid); SurfaceAxB[I-1]=dstupid;
   }

for (I = 1;I < FEGaussPoints*FEVertices+1; I++) { // DO I = 1,FEGaussPoints
   //READ(File_Unit,*,IOSTAT=IOS) (FEShapeFunctions(J,I),J=1,FEVertices)
   fscanf(fp, "%lf", &dstupid); FEShapeFunctions[I-1]=dstupid;
   }

for (Element = 1;Element < NumElements+1; Element++) { // DO Element = 1,NumElements
   for (I = 1;I < FEGaussPoints+1; I++) { // DO I = 1,FEGaussPoints
      //READ(File_Unit,*,IOSTAT=IOS) FEDetJacandWgt(I,Element)
      itemp1 = (Element-1)*FEGaussPoints + 1 - 1;
      fscanf(fp, "%lf", &dstupid); FEDetJacandWgt[itemp1+I-1]=dstupid;

      itemp1 = (Element-1)*FEGaussPoints*FENumDim*FEVertices + (I-1)*FENumDim*FEVertices + 1 - 1;

      for (K = 1;K < FENumDim*FEVertices+1; K++) { // DO K = 1,FENumDim
         //READ(File_Unit,*,IOSTAT=IOS) (FEDerivatives(J,K,I,Element),J=1,FEVertices)
         fscanf(fp, "%lf", &dstupid); FEDerivatives[itemp1+K-1]=dstupid;
         }
      }
   }

fclose(fp);

#ifdef Local_Debug
   // Print the connectivity matrix out
   for (Element = 1;Element < NumElements+1; Element++) { // DO Element = 1,NumElements
   //   WRITE(6,'("Element ",I5," vertex numbers -> ",4(I6,1X))') &
   //      Element, Conn(1,Element),Conn(2,Element),Conn(3,Element),Conn(4,Element)
      itemp1 = (Element-1)*4 + 1 - 1;
      printf("[SN-KERNEL]...Element %d vertex numbers -> %d %d %d %d \n",
             Element,Conn[itemp1],Conn[itemp1+1],Conn[itemp1+2],Conn[itemp1+3]);
      }
   for (I = 1;I < NumVertices+1; I++) { // DO I = 1,NumVertices
      //WRITE(6,'("Vertex ",I5," positions -> ",3(F13.6,1X))') I, GlobalXYZ(I,1),GlobalXYZ(I,2),GlobalXYZ(I,3)
      printf("[SN-KERNEL]...Vertex %d positions -> %13.6f %13.6f %13.6f \n",
             I,GlobalXYZ[I-1],GlobalXYZ[NumVertices+I-1],GlobalXYZ[2*NumVertices+I-1]);
      }
#endif

} // END SUBROUTINE Import_PMesh
