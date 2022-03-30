Copyright 1996 Massachusetts Institute of Technology and The Regents of the
C              University of California, San Diego.  All Rights Reserved.
C
      subroutine atm_grad_part(elevation,azimuth,partial)
C
C PURPOSE    : To compute N/S and E/W atmospheric gradient partials 
C              wrt phase observations.
C              Ref: G Chen and T.A. Herring Effects of Atmospheric 
C              Asymmetry in the Analysis of Space Geodetic Data 1996.
C
C PARAMETERS :
C         IN :  ELEVATION : Elevation angle of the Satellite     R*8
C               AZIMUTH   : Azimuth of the Satellite             R*8
C
C        OUT :  PARTIAL   : N/S and E/W atmospheric gradients 
C                           partial derivatives wrt phase        R*8(2) 
C
C SUBROUTINES CALLED :
C                                                       
C CREATED    :  96/09/27              LAST MODIFIED :  96/09/27
C
C AUTHOR     :  SIMON MCCLUSKY
C
      implicit none
C
      real*8 elevation, azimuth, partial(2)
C
C     N/S partial
      partial(1) = cos(azimuth)/(sin(elevation)*tan(elevation)+0.003d0)
C
C     E/W partial
      partial(2) = sin(azimuth)/(sin(elevation)*tan(elevation)+0.003d0)           
C
      return
      end
