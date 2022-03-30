Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
       SUBROUTINE HOFRAD(ANGLE,IDMS)
C
C      CONVERT FROM RADIANS TO DEGREES, MINUTES, SECONDS AND DECIMAL SECONDS
C
       implicit none 

       include '../includes/dimpar.h'   
       include '../includes/arc.h'


       real*8 angles,angle,scrat
       integer idms
       DIMENSION IDMS(4)
C
       ANGLES = ANGLE*360.D0/TWOPI
       IDMS(1) = ANGLES
       SCRAT = (ANGLES-IDMS(1))*60.D0
       IDMS(2) = SCRAT
       SCRAT = (SCRAT-IDMS(2))*60.D0
       IDMS(3) = SCRAT
C
C      NOTICE THE ROUNDUP IN THE FOLLOWING
       IDMS(4) = (SCRAT-IDMS(3))*10.D0 +0.5D0
       IF (IDMS(4).GE.10) IDMS(3)=IDMS(3)+1
       IF (IDMS(4).GE.10) IDMS(4)=IDMS(4)-10
C
       RETURN
       END
