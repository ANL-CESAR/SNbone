// This program is meant to emulate the on-node work of A*x which forms the bulk of the flop work in CFE methods
// It demonstrates the performance of simple linear tetrahedral matrix-vector products in SUPG
// ----------------------------------------------------------------------------
// Element Type: Linear Tetrahedron
// Grid:         Structured grid where each block is composed of 12 tetrahedrons
#include "ApplyA_common.h"
#include "ApplyA_functions.h"
#include <stdlib.h>
#include <stdio.h>
#ifdef WITHOMP
#include <omp.h>
#endif
int main(int argc, char *argv[]) {
// Module inclusions
#ifdef WITHBGQHPM
#include "mpif.h"
#endif
// Key Problem Size Variables - Specified via user command line input
int Output_Unit;
int Input_Scheme;          // The specific scheme to test
int Input_Angles;          // Total Angles to simulate
int Input_Iterations;      // User-defined # of iterations (statistics)
int Input_Nthreads;        // number of threads to use in openmp
int Input_BackVectors;     // Number of back vectors to use in the solver portion

// Local
int I,J,K,L,N,ReturnedError,NumArguments,ii_start,ii_end,I_SizeVec;
int Element,Angle,Space,MyThreadID;
int iNumMethods = 3;
char MyHPMname[] = "AVE-1AVE-2AVE-3";
unsigned int Len_MyHPMname=5;
double  MyFlopCnt[]={4929.0,4929.0,0.0}; // Flops / element/angle
int MyIterCnt[]={0,0,0};
double  MyTime[]={0.0,0.0,0.0};

double *LHS_C, *LHS_Answer;
double *RHS_C;

int    PP_AllCodes = 15;
int    PP_NumSegments = 6; // We will repeat the entire set of calculations this many times
char   PP_AllEventNames[]={"umm"};
signed long long  PP_values[PP_AllCodes*iNumMethods];
unsigned int Len_AllEventNames=1;
double SomeReal2,SomeReal1,AssemblyTime;

#ifdef WITHBGQHPM
   CALL MPI_INIT(ReturnedError)
#endif

// Read command line input
NumArguments = argc - 1; // # User defined command line arguments
//printf("NumArguments= %8d",NumArguments);
if ((NumArguments < 4) || (NumArguments > 5)) {
printf("[SN-KERNEL].............................................................................................................\n");
printf("[SN-KERNEL] The list of arguments was incomplete........................................................................\n");
printf("[SN-KERNEL] Version 1.0 SNaCFE mini-app to study on node performance....................................................\n");
printf("[SN-KERNEL] This mini-app is a test of the within-group FGMRES solver for a CFE SUPG based SN methodology...............\n");
printf("[SN-KERNEL] Usage:   snacfe.x  Scheme Iter BackV Angles Threads.........................................................\n");
printf("[SN-KERNEL] Example: snacfe.x  1      100  30    32     1      .........................................................\n");
printf("[SN-KERNEL] Scheme          specifies which scheme to use for the study (0=all).........................................\n");
printf("[SN-KERNEL] Iter(ation)     specifies the maximum FGMRES iterations to allow............................................\n");
printf("[SN-KERNEL] Back V(ectors)  specifies the maximum back vectors to use in FGMRES.........................................\n");
printf("[SN-KERNEL] Angles          specifies the number of angles assigned to the local process................................\n");
printf("[SN-KERNEL] T(hreads)       specifies the number of threads to use during the execution.................................\n");
printf("[SN-KERNEL].............................................................................................................\n");
abort();
}
else {
printf("[SN-KERNEL].............................................................................................................\n");
printf("[SN-KERNEL] Version 1.0 SNaCFE mini-app to study on node performance....................................................\n");
printf("[SN-KERNEL] This mini-app is a test of the within-group FGMRES solver for a CFE SUPG based SN methodology...............\n");
printf("[SN-KERNEL] Usage:   snacfe.x  Scheme Iter BackV Angles Threads.........................................................\n");
printf("[SN-KERNEL] Example: snacfe.x  1      100  30    32     1      .........................................................\n");
printf("[SN-KERNEL].............................................................................................................\n");
   // Assign command line input to variables
   Input_Scheme      = atoi(argv[1]);
   Input_Iterations  = atoi(argv[2]);
   Input_BackVectors = atoi(argv[3]);
   Input_Angles      = atoi(argv[4]);
   Input_Nthreads = 0;
   if (NumArguments == 5) {Input_Nthreads = atoi(argv[5]); }
}

#ifdef WITHOMP
   if (Input_Nthreads != 0) {omp_set_num_threads(Input_Nthreads);}
                       else {Input_Nthreads = omp_get_num_threads();} // This should be unnecessary, but...
#pragma omp parallel
{
  printf("[SN-KERNEL] Thread id %5d of %5d \n",omp_get_thread_num(),omp_get_num_threads());
}
#else
   if (Input_Nthreads == 0) {Input_Nthreads = 1;}
#endif

// Correct stupid inputs
ii_start = Input_Scheme;
ii_end   = Input_Scheme;
if ((Input_Scheme <= 0) || (Input_Scheme > iNumMethods)) {
   Input_Scheme = 0;
   ii_start = 1;
   ii_end   = iNumMethods;
}

if (Input_Iterations < 1) {Input_Iterations = 1;}
if (Input_BackVectors < 3) {Input_BackVectors  = 3;}
if (Input_BackVectors > 50) {Input_BackVectors = 50;}
if (Input_Angles  < 1) {Input_Angles  = 1;}

NumThreads=Input_Nthreads;
NumAngles=Input_Angles;

// Inform the user as to the problem size we will be running
printf("[SN-KERNEL] Running schemes %2d :  %2d \n",ii_start,ii_end);
printf("[SN-KERNEL] Number of Iterations   %8d \n", Input_Iterations);
printf("[SN-KERNEL] Number of Back Vectors %8d \n", Input_BackVectors);
printf("[SN-KERNEL] Number of Angles       %8d \n", Input_Angles);
printf("[SN-KERNEL] Number of Threads      %8d \n", Input_Nthreads);

// Import the mesh
printf("[SN-KERNEL] Importing the processed mesh \n");
Import_PMesh();

// Setup the angular cubature
printf("[SN-KERNEL] Setting up the angle cubature \n");
BuildAngleCubature(&NumAngles,Omega,OmegaOmega,AngleWeights);

printf("[SN-KERNEL] Stenciling the spatial NZ matrix \n");
StencilNZmatrix(&Input_Scheme);

printf("[SN-KERNEL] Number of Elements         %8d \n", NumElements);
printf("[SN-KERNEL] Number of Vertices         %8d \n", NumVertices);
printf("[SN-KERNEL] Vector Size Assembled      %8d \n", NumVertices*NumAngles);
printf("[SN-KERNEL] Vector Size Dis-assem      %8d \n", NumElements*NumAngles*FEVertices);
printf("[SN-KERNEL] Number of Non-zeros        %8d \n", NZS_NonZeros);
printf("[SN-KERNEL] Average Connections/Vertex %8d \n", NZS_NonZeros/NumVertices);

// GMRES and solution vectors
printf("[SN-KERNEL] Allocating FGMRES memory and solution vectors \n");
LHS_C      = (double*) malloc(NumAngles*NumVertices*sizeof(double));
LHS_Answer = (double*) malloc(NumAngles*NumVertices*sizeof(double));
RHS_C      = (double*) malloc(NumAngles*NumVertices*sizeof(double));
I_SizeVec = NumAngles*NumVertices;

// Assemble the matrix noting that it is part of the method 3 timing
printf("[SN-KERNEL] Building the non-zero space-angle matrices \n");
// -----------------------------------------------------
// This is the part of the code which we wish to measure
// -----------------------------------------------------
SomeReal1 = 0.0;
SomeReal2 = 0.0;
SomeReal1 = GETTHETIME();
AssembleNZmatrix(&Input_Scheme);
SomeReal2 = GETTHETIME();
AssemblyTime = SomeReal2 - SomeReal1;
//IF ((Input_Scheme .EQ. 0) .OR. (Input_Scheme .EQ. 3)) &
//   WRITE(Output_Unit,'("[SN-KERNEL] Took ",F13.6," seconds to assemble")') AssemblyTime
//printf("[SN-KERNEL] Assembly time %f %f %f \n",SomeReal1,SomeReal2,AssemblyTime);

// Get a source RHS_C, the solution LHS_Answer, and the guess for the LHS
printf("[SN-KERNEL] Generating an answer and its associated source with size %8d \n",I_SizeVec);
GenerateXb(&Input_Scheme,LHS_C,LHS_Answer,RHS_C);

// -----------------------------------------------------
// This is the part of the code which we wish to measure
// -----------------------------------------------------
#ifdef WITHBGQHPM
   summary_start() // This is the MPI tracking stuff
#endif

for (I = ii_start;I < ii_end+1; I++) { // DO I = ii_start,ii_end,1
   SomeReal1 = GETTHETIME();
#ifdef WITHBGQHPM
   call hpm_start(MyHPMname[(I-1)*Len_MyHPMname],Len_MyHPMname); // Initializes and starts the hardware counters
#endif
   SolveWGS(&Input_Iterations,&MyIterCnt[I-1],&I,LHS_C,RHS_C);
#ifdef WITHBGQHPM
   call hpm_stop(MyHPMname[(I-1)*Len_MyHPMname],Len_MyHPMname); // Stops the hardware counters
#endif
   SomeReal2 = GETTHETIME();
   MyTime[I-1] = SomeReal2 - SomeReal1;
   // Verify that the outgoing vector is accurate
   J = NumVertices*NumAngles;
   Verify(&Output_Unit,&J,LHS_C,LHS_Answer,&MyHPMname[(I-1)*Len_MyHPMname],Len_MyHPMname);
   } // methods

MyFlopCnt[2] = NZS_NonZeros*NumAngles; // This is the number of mults per application of A
PrintSummary(&Output_Unit,&NumElements,&NumAngles,&NumVertices,&Input_Iterations,&Krylov_BackVectors,
                  &iNumMethods,&ii_start,&ii_end,&AssemblyTime,
                  MyTime,MyFlopCnt,MyIterCnt,MyHPMname,&PP_AllCodes,&PP_NumSegments,PP_AllEventNames,PP_values,
                  Len_MyHPMname,Len_AllEventNames);
printf("[SN-KERNEL].............................................................................................................\n");

#ifdef WITHBGQHPM
   summary_stop(); // This finishes the MPI measurements
   MPI_Finalize();
#endif
}

double GETTHETIME() {
double SomeReal;
#ifdef WITHOMP
   SomeReal = omp_get_wtime();
   //printf("[SN-KERNEL] omp time %f \n",SomeReal);
#else
   SomeReal = clock();
#endif
return SomeReal;
} // END SUBROUTINE GETTHETIME