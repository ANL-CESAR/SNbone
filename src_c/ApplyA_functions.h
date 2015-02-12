// The various c function headers

void Verify(int *Output_Unit,int *iSize,double *Guess,double *Answer,char *SomeString,unsigned int Len_MyHPMname);
void Zero(int *iStart,int *iEnd,int *iSize,double *VectorToInitialize);
void Zero_Threaded(int *iSize,double *VectorToInitialize);
void BuildAngleCubature(int *NumAngles,double *Omega,double *OmegaOmega,double *AngleWeights);
void GenerateXb(int *Input_Scheme,double *LHS_C, double *LHS_Answer, double *RHS_C);
void PrintSummary(int *Output_Unit,int *NumElements,int *NumAngles,int *NumVertices,int *Input_Iterations,int *BackVectors, 
                      int *NumMethods,int *ii_start,int *ii_end,double *AssemblyTime, 
                      double *MyTime,double *MyFlopCnt,int *MyIterCnt,char *MyHPMname,
                      int *PP_AllCodes,int *PP_NumSegments,char *PP_AllEventNames,signed long long *PP_values,
                      unsigned int iLC1, unsigned int iLC2);
void FGMRES_Threaded(int *Output_Unit,
  int *Krylov_Local_Owned,int *Krylov_BackVectors,int *Krylov_Maximum_Iterations,int *Krylov_Iterations,
  double *Krylov_Absolute_Tolerance,double *Krylov_Relative_Tolerance,double *Krylov_Divergence_Tolerance,
  double *Krylov_Basis,double *Krylov_Hessenberg,double *Krylov_Givens,double *Krylov_PC_Basis, double *Krylov_Modified_RHS,
  double *Solution,double *RightHandSide,
  int *MyThreadID,int *MyStart,int *MyEnd,
  int *NumThreads,int *GuessIsNonZero,int *ReasonForConvergence,int *IterationCount,int *ParallelComm,int *ParallelRank,
  double *ResidualNorm,double *VectorNorm,double *VectorNorm_Local,double *HessenNorm_Local,
  int *iMethod);
double GETTHETIME();
void SolveWGS(int *Input_Iterations,int *IterationCount,int *iMethod,double *LHS_C,double *RHS_C);
void LookDownBelow(int *iMethod,double *RHS_C,double *LHS_C);
void StencilNZmatrix(int *Input_Scheme);
void Import_PMesh();

void ApplyA_AVE1_Tet_SUPG(int *NE, int *NA, int *NV, int *C,
                              int *AS_NumColors,int *AS_NumThreads,int *TasksPerThread,int *MyThreadID,int *ThreadWiseWork,
                              double *CF, double *CU, double *CP, double *FES, double *FED, double *FEW, double *OM, double *OO,
                              double *LHS, double *RHS);
void ApplyA_AVE2_Tet_SUPG(int *NE, int *NA, int *NV, int *C,
                           int *AS_NumColors,int *AS_NumThreads,int *TasksPerThread,int *MyThreadID,int *ThreadWiseWork,
                              double *CF, double *CU, double *CP, double *FES, double *FED, double *FEW, double *OM, double *OO,
                              double *LHS, double *RHS,double *SRHS,double *SLHS);

void SolveWGS_PassThrough_AVE1(double *RHS_C, double *LHS_C);
void SolveWGS_PassThrough_AVE1_NoHPM(double *RHS_C, double *LHS_C);
void SolveWGS_PassThrough_AVE2(double *RHS_C, double *LHS_C);
void SolveWGS_PassThrough_AVE2_NoHPM(double *RHS_C, double *LHS_C);
void SolveWGS_PassThrough_AVE3(double *RHS_C, double *LHS_C);
void SolveWGS_PassThrough_AVE3_NoHPM(double *RHS_C, double *LHS_C);
void SolveWGS_PassThrough_PC(double *RHS_C, double *LHS_C);

void AssembleNZmatrix(int *Input_Scheme);