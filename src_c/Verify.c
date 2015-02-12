// This program checks the solution
// ----------------------------------------------------------------------------
#include "ApplyA_functions.h"
#include <stdio.h>
// Fortran interface routines
void verify(int *Output_Unit,int *iSize,double *Guess,double *Answer,char *SomeString,unsigned int Len_MyHPMname) {
Verify(Output_Unit,iSize,Guess,Answer,SomeString,Len_MyHPMname);
}
void verify_(int *Output_Unit,int *iSize,double *Guess,double *Answer,char *SomeString,unsigned int Len_MyHPMname) {
Verify(Output_Unit,iSize,Guess,Answer,SomeString,Len_MyHPMname);
}
// Main subroutine header
void Verify(int *Output_Unit,int *iSize,double *Guess,double *Answer,char *SomeString,unsigned int Len_MyHPMname) {
// int     Output_Unit,iSize
// double   Guess(iSize),Answer(iSize)
// CHARACTER*8  SomeString
int I;
int Lfail;
double dummy;

Lfail = 0;
for (I = 1;I<*iSize+1;I++) {
   dummy = Guess[I-1] - Answer[I-1];
   if (dummy < 0.0) {dummy = -dummy;}
   if (dummy > 0.00001) {
      Lfail = I;
      break;
   }
}

if (Lfail != 0) {
   printf("[SN-KERNEL] ****ERROR**** check failed for ");
   for (I=1;I<Len_MyHPMname+1;I++) {printf("%c",SomeString[I-1]);}
   printf(" at %8d \n",Lfail); }
else {
   printf("[SN-KERNEL] ***SUCCESS*** check passed for ");
   for (I=1;I<Len_MyHPMname+1;I++) {printf("%c",SomeString[I-1]);}
   printf("\n");
}

} // END SUBROUTINE f_0_Verify

