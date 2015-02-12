#include <stdlib.h>
#include "ApplyA_common.h"
// This function is temporary as I transition from all fortran to all c
void setcommonfromfortran1_(
 int *fParallelRank,   int *fParallelSize,int *fParallelComm,int *fParallelGroup,
 int *fNumElements,    int *fNumVertices, int *fNumAngles,   int *fNumVacuum,
 int *fNumUnitNormals, int *fNumThreads, int *fAS_NumColors, int *fAS_NumThreads,
 int *fFEGaussPoints,  int *fFEVertices, int *fFENumDim, int *fTasksPerThread,
 int *fDA_ThreadWiseWork, int *fMM_ThreadWiseWork,int *fAS_ThreadWiseWork,
 double *fGlobalXYZ,int *fConn,int *fElementsWithVacuumBCs,int *fLocalSurfaceIndex,
 double *fConstTau,double *fConstF,double *fConstU,double *fConstUT,
 double *ffcoef,double *fpcoef, double *fucoef, double *fSurfaceAxB,
 double *fFEShapeFunctions, double *fFEDerivatives, double *fFEDetJacandWgt,
 double *fOmega,double *fOmegaOmega,double *fAngleWeights,
 int *fBCInfo,int *fIndexNormal,double *fVac_Normals,
 int *fVertexLocalToGlobal, int *fElementLocalToGlobal,
 double *fScratch_V1,double *fScratch_V2) {
// set the c integers and pointers
ParallelRank          =*fParallelRank;
ParallelSize          =*fParallelSize;
ParallelComm          =*fParallelComm;
ParallelGroup         =*fParallelGroup;
NumElements           =*fNumElements;
NumVertices           =*fNumVertices;
NumAngles             =*fNumAngles;
NumVacuum             =*fNumVacuum,
NumUnitNormals        =*fNumUnitNormals;
NumThreads            =*fNumThreads;
AS_NumColors          =*fAS_NumColors;
AS_NumThreads         =*fAS_NumThreads;
FEGaussPoints         =*fFEGaussPoints;
FEVertices            =*fFEVertices;
FENumDim              =*fFENumDim;
TasksPerThread        =*fTasksPerThread;

DA_ThreadWiseWork     = fDA_ThreadWiseWork;
MM_ThreadWiseWork     = fMM_ThreadWiseWork;
AS_ThreadWiseWork     = fAS_ThreadWiseWork;

GlobalXYZ             = fGlobalXYZ;
Conn                  = fConn;
ElementsWithVacuumBCs = fElementsWithVacuumBCs;
LocalSurfaceIndex     = fLocalSurfaceIndex;

ConstTau              = fConstTau;
ConstF                = fConstF;
ConstU                = fConstU;
ConstUT               = fConstUT;

fcoef                 = ffcoef;
pcoef                 = fpcoef;
ucoef                 = fucoef;
SurfaceAxB            = fSurfaceAxB;

FEShapeFunctions      = fFEShapeFunctions;
FEDerivatives         = fFEDerivatives;
FEDetJacandWgt        = fFEDetJacandWgt;

Omega                 = fOmega;
OmegaOmega            = fOmegaOmega;
AngleWeights          = fAngleWeights;

BCInfo                = fBCInfo;
IndexNormal           = fIndexNormal;
Vac_Normals           = fVac_Normals;

VertexLocalToGlobal   = fVertexLocalToGlobal;
ElementLocalToGlobal  = fElementLocalToGlobal;

Scratch_V1            = fScratch_V1;
Scratch_V2            = fScratch_V2;
}

void setcommonfromfortran2_(int *fNZS_NonZeros, int *fNZS_RowLoc, int *fNZS_ColNum, double *fNZS_Data) {
NZS_NonZeros          =*fNZS_NonZeros;
NZS_RowLoc            = fNZS_RowLoc;
NZS_ColNum            = fNZS_ColNum;
NZS_Data              = fNZS_Data;
}

void setcommonfromfortran3_(
int    *fKrylov_Local_Owned,int    *fKrylov_BackVectors,int    *fKrylov_Maximum_Iterations,int    *fKrylov_Iterations,
double *fKrylov_Absolute_Tolerance,double *fKrylov_Relative_Tolerance,double *fKrylov_Divergence_Tolerance,
double *fKrylov_Storage_Local,double *fKrylov_Basis,double *fKrylov_Hessenberg,double *fKrylov_Givens,
double *fKrylov_PC_Basis,double *fKrylov_Modified_RHS) {
Krylov_Local_Owned            =*fKrylov_Local_Owned           ;
Krylov_BackVectors            =*fKrylov_BackVectors           ;
Krylov_Maximum_Iterations     =*fKrylov_Maximum_Iterations   ;
Krylov_Iterations             =*fKrylov_Iterations           ;
Krylov_Absolute_Tolerance     =*fKrylov_Absolute_Tolerance   ;
Krylov_Relative_Tolerance     =*fKrylov_Relative_Tolerance   ;
Krylov_Divergence_Tolerance   =*fKrylov_Divergence_Tolerance ;
Krylov_Storage_Local          =fKrylov_Storage_Local       ;
Krylov_Basis                  =fKrylov_Basis               ;
Krylov_Hessenberg             =fKrylov_Hessenberg          ;
Krylov_Givens                 =fKrylov_Givens              ;
Krylov_PC_Basis               =fKrylov_PC_Basis            ;
Krylov_Modified_RHS           =fKrylov_Modified_RHS        ;
}

