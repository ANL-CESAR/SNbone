#ifdef APPLYA_COMMON
  #define EXTERN
#else
  #define EXTERN extern
#endif
// So this is how you do common block variables in c
// Communicator setup
EXTERN int ParallelRank  ;
EXTERN int ParallelSize  ;
EXTERN int ParallelComm  ;
EXTERN int ParallelGroup ;

// Problem size
EXTERN int  NumElements    ;     // Total number of tetrahedral elements (5 per structured grid)
EXTERN int  NumVertices    ;     // Total number of vertices in structured grid
EXTERN int  NumAngles	   ;     // Number of angles visible to the local process
EXTERN int  NumVacuum      ;     // Number of tetrahedral surfaces (triangles) with vacuum boundary (2 tri/surf and 2 global surfs)
EXTERN int  NumUnitNormals ;     // Number of vacuum boundary unit normals
EXTERN int  NumThreads     ;     // Total number of threads
EXTERN int  AS_NumColors   ;     // The number of non-overlapping (in memory) pieces of work that the FE assembly operation can be broken into
EXTERN int  AS_NumThreads  ;     // The number of threads worth of data the mesh was broken into (Ideally this matches NumThreads)
EXTERN int  FEGaussPoints  ;     // The number of GaussPoints in the cubature
EXTERN int  FEVertices     ;     // The number of vertices per element (4)
EXTERN int  FENumDim       ;     // The number of dimensions in the problem (3)
// Threading work breakdown for starting (1) and stopping (2) element id.
EXTERN int  TasksPerThread     ;    // This variable accounds for the AS_NumThreads <> NumThreads
EXTERN int  *DA_ThreadWiseWork ;    // (2,AS_NumThreads) The disassembly operation
EXTERN int  *MM_ThreadWiseWork ;    // (2,AS_NumThreads) The matrix-matrix product
EXTERN int  *AS_ThreadWiseWork ;    // (2,AS_NumColors,AS_NumThreads) The assembly operation
// Mesh related data structure
EXTERN double  *GlobalXYZ            ; // The global element vertex positions
EXTERN int     *Conn                 ; // Connectivity matrix (lists 4 global vertices for 1st element 4 global vertices for 2nd element etc.)
EXTERN int     *ElementsWithVacuumBCs; // Elements which have a vacuum boundary condition (listed once for each vacuum surface)
EXTERN int     *LocalSurfaceIndex    ; // Identifies the reference surface with the vacuum bc (123 or 4)
// Element-wise "cross section" values which change with respect to energy
EXTERN double   *ConstTau; // The element-wise stabilization factor
EXTERN double   *ConstF  ; // The coefficient for the F_element matrix
EXTERN double   *ConstU  ; // The coefficient for the U_element matrix
EXTERN double   *ConstUT ; // The coefficient for the UT_element matrix
// These are the spatial "matrices" that result when using the tetrahedral elements
EXTERN double  *fcoef      ; // 1  coeff/ele
EXTERN double  *pcoef      ; // 36 coeff/ele: 1st 6 -> p11c(1:6) 2nd 6 -> p22(1:6) ..
EXTERN double  *ucoef      ; // 9  coeff/ele: (/ iJ11iJ12iJ13iJ21iJ22iJ23iJ31iJ32iJ33 /)
EXTERN double  *SurfaceAxB ; // Surfaces
// These are needed for the conventional finite element implementation where the spatial matrices are not stored by computed during each iteration
EXTERN double *FEShapeFunctions;
EXTERN double *FEDerivatives   ; // FEVerticesFENumDimFEGaussPointsNumElements
EXTERN double *FEDetJacandWgt  ; // FEGaussPointsNumElements
// Angular related data structure
EXTERN double    *Omega       ;  // Omega stores O1 O2 O3 components of each angle
EXTERN double    *OmegaOmega  ;  // OmegaOmega stores O1*O1O2*O2O3*O3O1*O2O1*O3O2*O3
EXTERN double    *AngleWeights;  // This stores the angular weights used to integrate spherical harmonics
// These are required to identify the boundary surfaces
EXTERN int       *BCInfo      ;  // (2,NumElements) (1=)the number of boundary surfaces (2=) the first surface
EXTERN int       *IndexNormal ;  // Specifies the index of the unit normal for this surface (unnecessary since we're just storing them all uniquely, but would occur in a more complex code)
EXTERN double    *Vac_Normals ;  // X,Y,Z components of the unit normal for each vacuum bc surface
// These allow us to translate the solution back to the serial space
EXTERN int     *VertexLocalToGlobal  ; // The global id of each vertex
EXTERN int     *ElementLocalToGlobal ; // The global id of each element
// Scratch vector storage
EXTERN double    *Scratch_V1,*Scratch_V2; // (NumAngles,FEVertices,Threads)
// Sparse matrix storage 
EXTERN int       NZS_NonZeros;
EXTERN int      *NZS_RowLoc  ; // (NumVertices+1)  The NZS_Data storage position
EXTERN int      *NZS_ColNum  ; // (NZS_NonZeros)   The column numbers associated with NZS_Data
EXTERN double   *NZS_Data    ; // (NumAngles*NZS_NonZeros)

// Krylov storage
EXTERN int    Krylov_Local_Owned;           // The number of dof assigned to this processor
EXTERN int    Krylov_BackVectors;           // Number of BackVectors to be used
EXTERN int    Krylov_Maximum_Iterations;    // Maximum number of iterations
EXTERN int    Krylov_Iterations;            // Current number of iteration performed to solve the problem
EXTERN double Krylov_Absolute_Tolerance;           // Absolute tolerance for convergence
EXTERN double Krylov_Relative_Tolerance;           // Relative tolerance for convergence (relative to norm of the first residual)
EXTERN double Krylov_Divergence_Tolerance;         // Relative tolerance to determine divergence (relative to norm of the first residual)
EXTERN double *Krylov_Storage_Local;               // (Local_Owned)      Locally assigned parallel flux vector
EXTERN double *Krylov_Basis;                       // (Local_Owned,BackVectors) Krylov subspace orthonormal basis vector
EXTERN double *Krylov_Hessenberg;                  // (BackVectors,BackVectors+1) Hessenberg matrix coefficients (dot products)
EXTERN double *Krylov_Givens;                      // (BackVectors,2) Givens rotation (cos,sin) coefficients
EXTERN double *Krylov_PC_Basis;                    // (Local_Owned,BackVectors) Used by FGMRES to store the preconditionned basis vector
EXTERN double *Krylov_Modified_RHS;                // (BackVectors) Modified right hand side (used by GMRES and FGMRES) 
