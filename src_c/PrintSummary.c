// This program checks the solution
// ----------------------------------------------------------------------------
#include "ApplyA_functions.h"
#include <stdio.h>
// Fortran interface routines
void printsummary(int *Output_Unit,int *NumElements,int *NumAngles,int *NumVertices,int *Input_Iterations,int *BackVectors, 
                      int *NumMethods,int *ii_start,int *ii_end,double *AssemblyTime, 
                      double *MyTime,double *MyFlopCnt,int *MyIterCnt,char *MyHPMname,
                      int *PP_AllCodes,int *PP_NumSegments,char *PP_AllEventNames,signed long long *PP_values,
                      unsigned int Len_MyHPMname, unsigned int Len_AllEventNames) {
PrintSummary(Output_Unit,NumElements,NumAngles,NumVertices,Input_Iterations,BackVectors, 
                 NumMethods,ii_start,ii_end,AssemblyTime, 
                 MyTime,MyFlopCnt,MyIterCnt,MyHPMname,
                 PP_AllCodes,PP_NumSegments,PP_AllEventNames,PP_values,
                 Len_MyHPMname,Len_AllEventNames);
}
void printsummary_(int *Output_Unit,int *NumElements,int *NumAngles,int *NumVertices,int *Input_Iterations,int *BackVectors, 
                      int *NumMethods,int *ii_start,int *ii_end,double *AssemblyTime, 
                      double *MyTime,double *MyFlopCnt,int *MyIterCnt,char *MyHPMname,
                      int *PP_AllCodes,int *PP_NumSegments,char *PP_AllEventNames,signed long long *PP_values,
                      unsigned int Len_MyHPMname, unsigned int Len_AllEventNames) {
PrintSummary(Output_Unit,NumElements,NumAngles,NumVertices,Input_Iterations,BackVectors, 
                 NumMethods,ii_start,ii_end,AssemblyTime, 
                 MyTime,MyFlopCnt,MyIterCnt,MyHPMname,
                 PP_AllCodes,PP_NumSegments,PP_AllEventNames,PP_values,
                 Len_MyHPMname,Len_AllEventNames);
}
// Main subroutine header
void PrintSummary(int *Output_Unit,int *NumElements,int *NumAngles,int *NumVertices,int *Input_Iterations,int *BackVectors, 
                      int *NumMethods,int *ii_start,int *ii_end,double *AssemblyTime, 
                      double *MyTime,double *MyFlopCnt,int *MyIterCnt,char *MyHPMname,
                      int *PP_AllCodes,int *PP_NumSegments,char *PP_AllEventNames,signed long long *PP_values,
                      unsigned int iLC1, unsigned int iLC2) {
//  MyTime(NumMethods+1),MyFlopCnt(NumMethods),MyIterCnt(NumMethods),CHAR*8 MyHPMname(NumMethods),PP_AllEventNames(PP_AllCodes)
//  int*8 PP_values(PP_AllCodes*NumMethods)
// Local
int     I,J,K,L;
double  Temp1,Temp2,Temp3,Temp4,Temp5,Temp6,ddtemp;
double  TempSc1[100],TempSc2[100],TempSc3[100],TempSc4[100],TempSc5[100];
double  MyMemory[100]; // This will just be the memory required to apply the A matrix
double  EstimatedFlops[100]; // This is the estimated flop work that was done

// -------------------------------------------------------------------------------------------------
// Convert the flops estimate per element/angle into flops for this particular problem
Temp1 = *NumElements * *NumAngles;
MyFlopCnt[0] = MyFlopCnt[0] * Temp1;
MyFlopCnt[1] = MyFlopCnt[1] * Temp1;

// -------------------------------------------------------------------------------------------------
// Compute the amount of memory and flops used in each methodology
// Method 1 uses gaussian integration
Temp1 = *NumElements;
Temp2 = *NumAngles * Temp1;
Temp3 = *NumAngles * *NumVertices;
Temp4 = 4.0*Temp1;
for (I=1; I<100; I++) {MyMemory[I-1] = 0.0;}
MyMemory[0] = Temp4*0.5                                  // Connectivity
            + 3.0 * Temp1                                // ConstF UT and Tau
            + 4.0*4.0 + (4.0*3.0*4.0)*Temp1 + 4.0        // Basis and derivative function evaluations
            + *NumAngles * 9.0 ;                         // Omega and OmegaOmega
// Method 2 is a copy
MyMemory[1] = MyMemory[0];
// Method 5 uses assembled stored matrices
Temp5 = *NumAngles;
MyMemory[2] = MyFlopCnt[2] * 1.5 + (*NumVertices+1.0);  // Total storage of the space-angle assembled matrix

// -------------------------------------------------------------------------------------------------
// Intermediate memory print out
Temp4 = MyIterCnt[*ii_start-1];
if (*BackVectors < Temp4) {Temp4=*BackVectors; }
ddtemp =   Temp3                // Solution vector
           + Temp3*Temp4*2.0    // Basis and PC basis
           + Temp4*(Temp4+1.0)  // Hessenberg
           + (Temp4+1.0)*3.0;   // Givens and Modified RHS
printf("[SN-KERNEL] Average FGMRES memory (MB) %13.3f \n",ddtemp/131072.0);

// -------------------------------------------------------------------------------------------------
// Include the FGMRES memory used in each case and the FGMRES flop work
Temp2 = *NumVertices* *NumAngles;
Temp3 = *NumVertices* *NumAngles;
Temp4 = 0.0;
for (I = 1;I<*BackVectors+1;I++) {
   Temp4 = Temp4 + I*2.0 + 2.0; // The number of dot products and scaling operations needed for each inner in FGMRES
}
for (I=1; I<100; I++) {EstimatedFlops[I-1] = 0.0;}
for (I = 1;I < *NumMethods+1;I++) {
   Temp1 = MyIterCnt[I-1]; // The number of iterations used to solve the problem
   if (Temp1 < *BackVectors) {Temp1 = *BackVectors;}  // The number of iterations used to solve the problem
   MyMemory[I-1] = MyMemory[I-1] + Temp2               // Solution vector
                                 + Temp2*Temp1*2.0     // Basis and PC basis
                                 + Temp1*(Temp1+1.0)   // Hessenberg
                                 + (Temp1+1.0)*3.0;    // Givens and Modified RHS
   J = MyIterCnt[I-1] / *BackVectors;  // The number of times we restarted in FGMRES
   K = MyIterCnt[I-1] - *BackVectors*J; // The number of iterations past the last restart used
   Temp5 = 0.0;
   for (L = 1;L<K+1;L++) {
      Temp5 = Temp5 + L*2.0 + 2.0; // The number of dot products and scaling operations needed for each inner in FGMRES
   }
   Temp6 = J;  // The number of outers
   EstimatedFlops[I-1] = 2.0*Temp3 + MyFlopCnt[I-1]                           // The work performed upon entry
                       + 1.0*Temp3                                            // The vector scaling before the outer loop
                       + (Temp6+1.0)*(*BackVectors*(*BackVectors+1.0)         // The hessian matrix after the outer loop
                                     + *BackVectors*Temp3                     // The Modified RHS scaling of each vector
                                     + MyFlopCnt[I-1])                        // The Apply A outside of the outer
                     +     Temp6*MyFlopCnt[I-1] + Temp4                       // The Apply A inside the inner and cost of the inner work
                     + (K*1.0)*MyFlopCnt[I-1] + Temp5;                        // The Apply A inside the last pass through the inner
}

// -------------------------------------------------------------------------------------------------
ddtemp = 1.0;
Temp1 = 1.0 / (1024.0*1024.0*1024.0); // Convert to Giga flops
#ifdef WITHPAPI
   ddtemp = PP_NumSegments;
   ddtemp = 1.0 / ddtemp;
#endif
for (I = 1;I<*NumMethods+1;I++) {
   TempSc1[I-1]= MyMemory[I-1]/131072.0; // * 8.0 / 1024d0/1024d0
   TempSc2[I-1]= EstimatedFlops[I-1]*Temp1;  // In Gflops
   TempSc3[I-1]= MyTime[I-1]*ddtemp;
   TempSc4[I-1]= TempSc2[I-1]/(TempSc3[I-1]+1.e-24);
   TempSc5[I-1]= TempSc2[I-1]/TempSc1[I-1]*1024.0; // Gflops/GByte
}

// Say something about the matrix assembly
if (*AssemblyTime > 0.0) {
printf("[SN-KERNEL] Assembly Time =%13.6f seconds or %13.1f FGMRES iterations \n",*AssemblyTime,*AssemblyTime/MyTime[2] * MyIterCnt[2]);
}

#ifdef WITHPAPI
                        //123456789*  123456789*123  123456789*123  123456789*123  123456789*123  123456789*123
printf("[SN-KERNEL] Method           Time       Est. Mem(MB)   Est. GFlops   Est. GFlops/s   GFlops/GByte  ");
// probably should fix this
//for (I = 1;I < *PP_AllCodes+1; I++){
//    printf("PAPI-> %16d",PP_AllEventNames[I-1]);
//}
printf("\n");

for (J = *ii_start;J < *ii_end+1; J++) {
   printf("=SN-KERNEL] %10c  %13.5f  %13.3f  %13.4f  %13.3f  %13.2f  ",
          MyHPMname[J-1],TempSc3[J-1],TempSc1[J-1],TempSc2[J-1],TempSc4[J-1],TempSc5[J-1]);
   for (I = 1,I < *PP_AllCodes+1, I++) {printf("%16d",PP_values[(J-1)*PP_AllCodes+I-1]);}
   printf("\n");
}
#else
                        //123456789*  123456789*123  123456789*123  123456789*123  123456789*123  123456789*123
printf("[SN-KERNEL] Method           Time       Est. Mem(MB)   Est. GFlops   Est. GFlops/s   GFlops/GByte  \n");
for (J = *ii_start;J < *ii_end+1; J++) {
   if (TempSc4[J-1] > 1000000.0) {TempSc4[J-1]=-1.0;}
   printf("=SN-KERNEL] ");
   K = (J-1)*iLC1;
   for (I = 1;I <iLC1+1;I++) {printf("%c",MyHPMname[K+I-1]);}
   printf("  %13.5f  %13.3f %13.4f  %13.3f  %13.2f  \n",TempSc3[J-1],TempSc1[J-1],TempSc2[J-1],TempSc4[J-1],TempSc5[J-1]);
}

#endif
} // END SUBROUTINE PrintSummary

