#include <stdlib.h>
// This is the main c common block
#define APPLYA_COMMON
#include "ApplyA_common.h"

// The end target setup..or close to it
void CommonBlock() {

// Threading work breakdown for starting (1) and stopping (2) element id.
DA_ThreadWiseWork = (int*) malloc(2*AS_NumThreads*sizeof(int));  // (2,AS_NumThreads) The disassembly operation
MM_ThreadWiseWork = (int*) malloc(2*AS_NumThreads*sizeof(int));  // (2,AS_NumThreads) The matrix-matrix product
AS_ThreadWiseWork = (int*) malloc(2*AS_NumColors*AS_NumThreads*sizeof(int));  // (2,AS_NumColors,AS_NumThreads) The assembly operation

// Mesh related info
GlobalXYZ             = (double*) malloc(NumVertices*3*sizeof(double)); // The global element vertex positions
Conn                  = (int*)    malloc(FEVertices*NumElements*sizeof(int));  // Connectivity matrix (lists 4 global vertices for 1st element, 4 global vertices for 2nd element, etc.)
ElementsWithVacuumBCs = (int*)    malloc(NumVacuum*sizeof(int));  // Elements which have a vacuum boundary condition (listed once for each vacuum surface)
LocalSurfaceIndex     = (int*)    malloc(NumVacuum*sizeof(int));  // Identifies the reference surface with the vacuum bc (1,2,3 or 4)

// Method-dependent arrays
ConstTau = (double*) malloc(NumElements*sizeof(double));  // space dependent element coefficient
ConstF   = (double*) malloc(NumElements*sizeof(double));  // method-dependent coefficient for the F_element matrix (=sigt+tau*a1*sigt**2)
ConstU   = (double*) malloc(NumElements*sizeof(double));  // method-dependent coefficient for the U_element matrix (=tau*a3*sigt)
ConstUT  = (double*) malloc(NumElements*sizeof(double));  // method-dependent coefficient for the UT_element matrix (=1+tau*a2*sigt)

// Angular Omega matrices
Omega        = (double*) malloc(NumAngles*3*sizeof(double));  // Omega stores O1, O2, O3 components of each angle
OmegaOmega   = (double*) malloc(NumAngles*6*sizeof(double));  // OmegaOmega stores O1*O1,O2*O2,O3*O3,O1*O2,O1*O3,O2*O3
AngleWeights = (double*) malloc(NumAngles*sizeof(double));  // Store the angular weights

// These are the spatial "matrices" that result when using the tetrahedral elements
fcoef        = (double*) malloc(NumElements*sizeof(double)); // 1  coeff/ele
pcoef        = (double*) malloc(36*NumElements*sizeof(double)); // 36 coeff/ele: 1st 6 -> p11c(1:6), 2nd 6 -> p22(1:6), ..
ucoef        = (double*) malloc(9*NumElements*sizeof(double)); // 9  coeff/ele: (/ iJ11,iJ12,iJ13,iJ21,iJ22,iJ23,iJ31,iJ32,iJ33 /)
SurfaceAxB   = (double*) malloc(NumVacuum*sizeof(double)); // Surfaces

// These are needed for the conventional finite element implementation where the spatial matrices are not stored by computed during each iteration
FEShapeFunctions  = (double*) malloc(FEVertices*FEGaussPoints*sizeof(double));
FEDerivatives     = (double*) malloc(FEVertices*FENumDim*FEGaussPoints*NumElements*sizeof(double));
FEDetJacandWgt    = (double*) malloc(FEGaussPoints*NumElements*sizeof(double));

// These are required to identify the boundary surfaces
BCInfo       = (int*)    malloc(2*NumElements*sizeof(int));
Vac_Normals  = (double*) malloc(3*NumUnitNormals*sizeof(double));  // X,Y,Z components of the unit normal for each vacuum bc surface
IndexNormal  = (int*)    malloc(NumVacuum*sizeof(int));  // Specifies the index of the unit normal for this surface (unnecessary since we're just storing them all uniquely, but would occur in a more complex code)

// These allow us to translate the solution back to the serial space
VertexLocalToGlobal   = (int*) malloc(NumVertices*sizeof(int));
ElementLocalToGlobal  = (int*) malloc(NumElements*sizeof(int));

// Scratch vector storage
Scratch_V1            = (double*) malloc(NumElements*FEVertices*NumThreads*sizeof(double));
Scratch_V2            = (double*) malloc(NumElements*FEVertices*NumThreads*sizeof(double));

// Krylov storage
Krylov_Local_Owned         = NumVertices*NumAngles;
Krylov_BackVectors         = 30;
Krylov_Maximum_Iterations  = 100;
Krylov_Iterations          = 0;
Krylov_Absolute_Tolerance  = 1.0e-12;
Krylov_Relative_Tolerance  = 1.0e-8;
Krylov_Divergence_Tolerance = 1.0e-3;
Krylov_Storage_Local = (double*) malloc(Krylov_Local_Owned*sizeof(double));                              // (Local_Owned)      Locally assigned parallel flux vector
Krylov_Basis         = (double*) malloc(Krylov_Local_Owned*(Krylov_BackVectors+1)*sizeof(double));       // (Local_Owned,BackVectors) Krylov subspace orthonormal basis vector
Krylov_Hessenberg    = (double*) malloc((Krylov_BackVectors+1)*(Krylov_BackVectors+1)*sizeof(double));   // (BackVectors,BackVectors+1) Hessenberg matrix coefficients (dot products)
Krylov_Givens        = (double*) malloc((Krylov_BackVectors+1)*2*sizeof(double));                        // (BackVectors,2) Givens rotation (cos,sin) coefficients
Krylov_PC_Basis      = (double*) malloc(Krylov_Local_Owned*(Krylov_BackVectors+1)*sizeof(double));       // (Local_Owned,BackVectors) Used by FGMRES to store the preconditionned basis vector
Krylov_Modified_RHS  = (double*) malloc((Krylov_BackVectors+1)*sizeof(double));                          // (BackVectors) Modified right hand side (used by GMRES and FGMRES) 

}

