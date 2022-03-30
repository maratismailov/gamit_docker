c     Common for GNSS type and models used in both ARC and MODEL
c     R. King 180321

      character*5 frame,precmod,gravmod,nutmod
      character*8 speopmod
      common/framemod/frame,precmod,gravmod,nutmod,speopmod 
                                             
      character*5 srpmod,eradmod,antradmod 
      common/orbmod/srpmod,eradmod,antradmod 

c   frame     c*5   inertial frame J2000 or B1950 from t-file
c   precmod   c*5   precession IAU68 or IAU76 from t-file
c   nutmod    c*5   nutation  IAU80 or IAU00 from t-file
c   srpmod    c*5   solar radiation pressure BERNE etc. from t-file
c   gravmod   c*5   gravity field, EGR08 etc from t-file
c   eradmod   c*5   Earth radiation pressure NCLE1 etc from t-file
c   antradmod c*5   antenna radiation  ANTBK etc from t-file  
c   speopmod  c*8   model for short-period UT1 and pole, IERS96, IERS10, etc

      character*1 gnss
* MOD TAH 200511: Added gnsslf  
      character*2 gnsslf  ! GNSS with lower frequency designation     
      common /systyp/ gnss, gnsslf
      
c    GNSS constellation
c      G GPS
c      R Glonass
c      C Beidou
c      E Galileo
c      I INRSS 






