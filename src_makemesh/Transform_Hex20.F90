!---------------------------------------------------------------------------------------------------------------------------------
! This subroutine allows a series of points to be transformed using a hardwired Quadratic FE interpolation scheme
!---------------------------------------------------------------------------------------------------------------------------------
#include "PROTEUS_Preprocess.h"
SUBROUTINE Transform_Hex20(ReturnX,ReturnY,ReturnZ,XXABS,YYABS,ZZABS)
IMPLICIT NONE
! Passed in/out
PROTEUS_Real ReturnX,ReturnY,ReturnZ
PROTEUS_Real XXABS,YYABS,ZZABS ! The reference position we want to evaluate the basis functions at
! Local junk
PROTEUS_Int, PARAMETER :: FIXEDNPE = 20, FIXEDJAC = 3
PROTEUS_Real XI(FIXEDNPE),ETA(FIXEDNPE),ZETA(FIXEDNPE)
DATA   XI/-1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0,-1.0D0,-1.0D0,-1.0D0, 1.0D0, 1.0D0,-1.0D0,-1.0D0, 0.0D0, 1.0D0, &
           1.0D0, 1.0D0, 0.0D0,-1.0D0,-1.0D0/
DATA  ETA/-1.0D0,-1.0D0,-1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0,-1.0D0,-1.0D0, 1.0D0, 1.0D0,-1.0D0,-1.0D0,-1.0D0, &
           0.0D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0/
DATA ZETA/-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0, &
           1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/
! FE transformation of the entire domain (Bottom 8 vertices, then middle 4, then top 8)
! No change
!PROTEUS_Real :: CoorX(FIXEDNPE)=(/ -1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0,-1.0D0,-1.0D0,  -1.0D0, 1.0D0, 1.0D0,-1.0D0, &
!                                   -1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0,-1.0D0,-1.0D0/)
!PROTEUS_Real :: CoorY(FIXEDNPE)=(/ -1.0D0,-1.0D0,-1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0,  -1.0D0,-1.0D0, 1.0D0, 1.0D0, &
!                                   -1.0D0,-1.0D0,-1.0D0, 0.0D0, 1.0D0, 1.0D0, 1.0D0, 0.0D0/)
!PROTEUS_Real :: CoorZ(FIXEDNPE)=(/-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,-1.0d0,    0.0d0, 0.0d0, 0.0d0, 0.0d0, &
!                                   1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0, 1.0d0/)
! 30 degree twisted annulus
!PROTEUS_Real :: CoorX(FIXEDNPE)=(/ 1.0D0, 1.5D0, 2.0D0, 1.41D0,0.0D0, 0.0D0, 0.0D0,0.71D0,   1.43D0, 2.37D0, 0.0D0, 0.25D0, &
!                                   1.87D0,2.3d0, 2.73D0,1.52D0,0.0D0,0.25D0, 0.5D0,1.26D0/)
!PROTEUS_Real :: CoorY(FIXEDNPE)=(/ 0.0D0, 0.0D0, 0.0D0, 1.41D0,2.0D0, 1.5D0, 1.0D0,0.71D0,   0.0D0,  0.25D0, 1.62D0, 0.68D0, &
!                                   0.0D0, 0.25D0,0.5D0, 1.43D0,1.23D0,0.8D0, 0.37D0,0.47D0/)
!PROTEUS_Real :: CoorZ(FIXEDNPE)=(/-1.17d0,-1.12d0,-1.07d0,-1.02d0,-0.98d0,-0.93d0,-0.88d0,-0.83d0,-0.1d0, 0.1d0,-0.1d0, 0.1d0, &
!                                   0.83d0, 0.88d0, 0.93d0, 0.98d0, 1.02d0, 1.07d0, 1.12d0, 1.17d0/)
! 45 degree twisted annulus
!PROTEUS_Real :: CoorX(FIXEDNPE)=(/ 1.5D0, 1.75d0,2.0D0, 1.00D0,0.0D0,0.25D0,0.50D0,1.00D0,   1.25D0, 2.00D0, 0.00D0,0.25D0, &
!                                   1.0D0, 1.5D0, 2.0D0, 1.41D0,0.0D0, 0.0D0, 0.0D0,0.71D0/)
!PROTEUS_Real :: CoorY(FIXEDNPE)=(/ 0.0D0, 0.55D0,1.1D0, 2.00D0,1.1D0,0.55D0,0.0D0,0.46D0,   0.0D0,  0.55D0, 1.55D0, 0.5D0, &
!                                   0.0D0, 0.0D0, 0.0D0, 1.41D0,2.0D0, 1.5D0, 1.0D0,0.71D0/)
!PROTEUS_Real :: CoorZ(FIXEDNPE)=(/-1.17d0,-1.12d0,-1.07d0,-1.02d0,-0.98d0,-0.93d0,-0.88d0,-0.83d0,-0.1d0, 0.1d0,-0.1d0, 0.1d0, &
!                                   0.83d0, 0.88d0, 0.93d0, 0.98d0, 1.02d0, 1.07d0, 1.12d0, 1.17d0/)
! 90 degree twisted annulus
PROTEUS_Real :: CoorX(FIXEDNPE)=(/ 1.0D0, 1.5D0, 2.0D0, 1.41D0,0.0D0, 0.0D0, 0.0D0,0.71D0,   1.5D0, 2.0D0, 0.0D0, 0.5D0, &
                                   2.0D0, 2.0d0, 2.0D0, 0.59D0,0.0D0, 0.5D0, 1.0D0,1.29D0/)
PROTEUS_Real :: CoorY(FIXEDNPE)=(/ 0.0D0, 0.0D0, 0.0D0, 1.41D0,2.0D0, 1.5D0, 1.0D0,0.71D0,   0.5D0, 1.0D0, 1.0D0, 0.5D0, &
                                   1.0D0, 1.5D0, 2.0D0, 1.41D0,0.0D0, 0.0D0, 0.0D0,0.71D0/)
PROTEUS_Real :: CoorZ(FIXEDNPE)=(/-1.17d0,-1.12d0,-1.07d0,-1.02d0,-0.98d0,-0.93d0,-0.88d0,-0.83d0,-0.1d0, 0.1d0,-0.1d0, 0.1d0, &
                                   0.83d0, 0.88d0, 0.93d0, 0.98d0, 1.02d0, 1.07d0, 1.12d0, 1.17d0/)
PROTEUS_Real AR,BR,CR,DR,ER,HHH
PROTEUS_Int I,GP

ReturnX = 0.0D0
ReturnY = 0.0D0
ReturnZ = 0.0D0

    ! 4 bottom corners
   DO I = 1,8,2
      AR = 1.0D0 + XI(I)   * XXABS
      BR = 1.0D0 + ETA(I)  * YYABS
      CR = 1.0D0 + ZETA(I) * ZZABS
      HHH = 0.125D0 * AR * BR * CR * (XI(I)*XXABS + ETA(I)*YYABS + ZETA(I)*ZZABS - 2.0D0)
      ReturnX = ReturnX + CoorX(I) * HHH
      ReturnY = ReturnY + CoorY(I) * HHH
      ReturnZ = ReturnZ + CoorZ(I) * HHH
   END DO
   ! 4 top corners
   DO I = 13,20,2
      AR = 1.0D0 + XI(I)   * XXABS
      BR = 1.0D0 + ETA(I)  * YYABS
      CR = 1.0D0 + ZETA(I) * ZZABS
      HHH = 0.125D0 * AR * BR * CR * (XI(I)*XXABS + ETA(I)*YYABS + ZETA(I)*ZZABS - 2.0D0)
      ReturnX = ReturnX + CoorX(I) * HHH
      ReturnY = ReturnY + CoorY(I) * HHH
      ReturnZ = ReturnZ + CoorZ(I) * HHH
   END DO
   ! 4 bottom edges
   DO I = 2,8,2
      AR = 1.0D0 + XI(I)   * XXABS
      BR = 1.0D0 + ETA(I)  * YYABS
      CR = 1.0D0 + ZETA(I) * ZZABS

      DR = (1.0D0 - XXABS * XXABS) * ETA(I) * ETA(I)
      ER = (1.0D0 - YYABS * YYABS) * XI(I)  * XI(I)
      HHH = 0.25D0 * ( DR * BR * CR  +  ER * AR * CR )
      ReturnX = ReturnX + CoorX(I) * HHH
      ReturnY = ReturnY + CoorY(I) * HHH
      ReturnZ = ReturnZ + CoorZ(I) * HHH
   END DO   
   ! 4 top edges
   DO I = 14,20,2
      AR = 1.0D0 + XI(I)   * XXABS
      BR = 1.0D0 + ETA(I)  * YYABS
      CR = 1.0D0 + ZETA(I) * ZZABS

      DR = (1.0D0 - XXABS * XXABS) * ETA(I) * ETA(I)
      ER = (1.0D0 - YYABS * YYABS) * XI(I)  * XI(I)
      HHH = 0.25D0 * ( DR * BR * CR  +  ER * AR * CR )
      ReturnX = ReturnX + CoorX(I) * HHH
      ReturnY = ReturnY + CoorY(I) * HHH
      ReturnZ = ReturnZ + CoorZ(I) * HHH
   END DO 
   ! 4 medians
   DO I = 9,12,1
      AR = 1.0D0 + XI(I)   * XXABS
      BR = 1.0D0 + ETA(I)  * YYABS
      CR = 1.0D0 + ZETA(I) * ZZABS
      DR = 1.0D0 - ZZABS * ZZABS
      HHH = 0.25D0 * DR * AR * BR
      ReturnX = ReturnX + CoorX(I) * HHH
      ReturnY = ReturnY + CoorY(I) * HHH
      ReturnZ = ReturnZ + CoorZ(I) * HHH
   END DO

END SUBROUTINE Transform_Hex20

