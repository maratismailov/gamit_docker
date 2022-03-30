Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
      SUBROUTINE DAYJUL( JD,iyear,IDOY )
C     WRITTEN BY R.KING 1 AUG 83
C     GIVEN JD, RETURNS DAY OF YEAR AND YEAR    
c     Modified for Y2K R King 27 Jul 99
C
      implicit none
      integer*4 jd,jd1900,id,idoy,iyear
      data jd1900/2415020/
C
      ID= (JD - jd1900)*100
      iyear= int(ID/36525) 
      IDOY= ID/100 - iyear*365 - (iyear-1)/4   
c     earliest GPS launch date is 1978; earliest space-geodetic data about 1960 
      if( iyear.gt.60 ) then
         iyear = iyear + 1900
      else
         iyear = iyear + 2000
      endif  
      RETURN
      END
