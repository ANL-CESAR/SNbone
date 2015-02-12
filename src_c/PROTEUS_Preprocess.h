!  Copyright(c) 2005 Argonne National Laboratory
!  This header file is used to control the data set assignments in the UNIC code
!     8 bit INTEGER -128 to 128
!    16 bit INTEGER -32,768 to 32,767
!    32 bit INTEGER -2,147,483,648 to 2,147,483,647
!    64 bit INTEGER -9,223,372,036,854,775,808 to 9,223,372,036,854,775,807
!   128 bit INTEGER why?
!    32 bit REAL 1.17549435E-38 to 3.40282347E38
!    64 bit REAL 2.2250738585072013D-308 to 1.7976931348623158D308
!   128 bit REAL 6.4751751194380251109244389582276465524996Q-4966 to 1.189731495357231765085759326628007016196477Q4932
!
!
! UNIC wide
! --------------
!#define PROTEUS_Debug
!#define Debug_ImposeIdentityMatrix

#define PROTEUS_R32_Tolerance       1.0d-6
#define PROTEUS_R64_Tolerance       1.0d-13
#define PROTEUS_R64_LowerCutoff     8.232d-11
#define PROTEUS_R64_HigherCutoff    1.304d19
#define PROTEUS_PeakToAverage       1.0D-7
#define PROTEUS_Definition_PI       3.141592653589793d0
#define PROTEUS_Definition_Avogadro 0.602214179d0
#define PROTEUS_Definition_J2MeV    6.24150974D12
#define PROTEUS_Convert_To_MB       9.53674316D-7
#define PROTEUS_Tchebychev_DR_low   0.1d0
#define PROTEUS_Tchebychev_DR_high  1.0d0
#define PROTEUS_MaxFileUnits        1000
#define PROTEUS_Definition_Yes      '      Yes       '
#define PROTEUS_Definition_No       '      No        '
#define PROTEUS_BC_Reflected        'REFLECTIVE      '
#define PROTEUS_BC_Vacuum           'VOID            '
#define PROTEUS_BC_GenericInterface 'GENERICINTERFACE'
#define PROTEUS_ULO_Error_Eigenvalue 1.0D-3
#define PROTEUS_LLO_Error_Eigenvalue 1.0D-10
#define PROTEUS_DO_Error_Eigenvalue  1.0D-6
#define PROTEUS_ULO_Error_Fission    1.0D-3
#define PROTEUS_LLO_Error_Fission    1.0D-10
#define PROTEUS_DO_Error_Fission     5.0D-6
#define PROTEUS_ULO_Error_Flux       1.0D-3
#define PROTEUS_LLO_Error_Flux       1.0D-10
#define PROTEUS_DO_Error_Flux        5.0D-7
!  Variable definitions
#define PROTEUS_Int_08bit     INTEGER(KIND=1)
#define PROTEUS_Int_16bit     INTEGER(KIND=2)
#define PROTEUS_Int_32bit     INTEGER(KIND=4)
#define PROTEUS_Int_64bit     INTEGER(KIND=8)
#define PROTEUS_Int           INTEGER(KIND=4)
#define PROTEUS_Real_32bit    REAL(KIND=4)
#define PROTEUS_Real_64bit    REAL(KIND=8)
#define PROTEUS_Real          REAL(KIND=8)
#define PROTEUS_Log           LOGICAL(KIND=4)
#define PROTEUS_FileNames     CHARACTER*128
#define PROTEUS_Char          CHARACTER*16
#define PROTEUS_Char32        CHARACTER*32
! This determines the pointer size (integer kind)
#define PROTEUS_ADDRESS_SIZE 8
#define PROTEUS_ADDRESS integer(kind=PROTEUS_ADDRESS_SIZE)

! Processor Angle-Space decomposition style
! #define PROTEUS_UseOldCommSetup

!  The matching byte size definitions
#define PROTEUS_Int_08bit_Size  1
#define PROTEUS_Int_16bit_Size  2
#define PROTEUS_Int_32bit_Size  4
#define PROTEUS_Int_64bit_Size  8
#define PROTEUS_Int_Size        4
#define PROTEUS_Real_32bit_Size 4
#define PROTEUS_Real_64bit_Size 8
#define PROTEUS_Real_Size       8
#define PROTEUS_Log_Size        4
#define PROTEUS_FileNames_Size  128
#define PROTEUS_Char_Size       16
#define PROTEUS_Char32_Size     32
#define PROTEUS_MessageLength  120

! Thread-wise operation settings
#define PROTEUS_AssemblyElementSize 1000


! NonZeroStorage
! --------------
#define PROTEUS_NZS_CatchErrors
#define PROTEUS_NZS_Int       INTEGER(KIND=4)
#define PROTEUS_NZS_MaxSize   INTEGER(KIND=8)
!  The matching byte size definitions
#define PROTEUS_NZS_Int_Size     4
#define PROTEUS_NZS_MaxSize_Size 8
#define PROTEUS_NZS_Matlab_View  .FALSE.

!  YLM_Matrices
! --------------
!#define PROTEUS_Debug_YLM
#define PROTEUS_Ylm_Real        REAL(KIND=8)
#define PROTEUS_Ylm_NZRWS_local TYPE (NZRWStorage_D_RowWiseData)
#define PROTEUS_Ylm_NZRWS_type  TYPE (NZRWStorage_D_Data)
#define PROTEUS_Ylm_NZS_type    TYPE (NZS_D_NonZeroStorage)

!  FE_Matrices
! --------------
!  #define PROTEUS_Debug_FE
#define PROTEUS_FE_Real       REAL(KIND=8)
#define PROTEUS_FE_Real_Size  8
