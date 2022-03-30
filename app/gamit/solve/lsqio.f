Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.
C
      Subroutine LSQIO

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'     
      include 'parameters.h' 

      integer i3,ik,ictyp,jj,ix1,i
 
      logical debug/.false./

c     Open M-file, read its header and get station coordinates
C
      CHARACTER*4 BUF4,BUF42,UPPERC


      real*8 phi,along,radius,height,geodrad,x,y,z
           
      if(debug) print *,'LSQIO obfiln(1) cfiln(1) rlabel(1) '
     .  ,obfiln(1),cfiln(1),rlabel(1)
       
c     initialization to avoid compiler warning
      i3 = 0
      ik = 0

      ICTYP=1
      IF(PREVAL(3).LT.6000.D0) ICTYP=2
C
C     CHECK FOR STATION PARTIALS
      IF(ISLOT1(1).GT.300) GO TO 101
C
C     CONSTRUCT CARTESIAN COORDINATES (FOR WEIGHTING) (UNITS-KM)
      DO 100 JJ=1,nsite                
        BUF4=UPPERC(cfiln(JJ)(2:5)) 
        DO 120 I=1,nsite
        IX1=(I-1)*3
        I3=IX1+1
        BUF42=UPPERC(RLABEL(I3)(1:4))  
        if(debug) print *,'LSQIO jj i ix1 i3, cfiln rlabel buf4 buf42 '
     .     ,  jj,i,ix1,i3,cfiln(jj),rlabel(i3),buf4,buf42 
        IF (BUF4.EQ.BUF42) THEN
          IK=IX1
          I3=(JJ-1)*3   
          GOTO 130
        ENDIF
  120   CONTINUE        
 130    PHI=PREVAL(IK+1)
        ALONG=PREVAL(IK+2)
        GO TO (400,401), ICTYP
  400   radius=preval(ik+3)
        CALL SPHXYZ(PHI,ALONG,radius,X,Y,Z,1)
        if(debug) print *,'LSQIO nsite jj ik i3 phi along radius '
     .       ,nsite,jj,ik,i3,phi,along,radius 
        GO TO 402  
  401   height = preval(ik+3)
        CALL GEOXYZ(semi,finv,phi,along,height,geodrad,x,y,z,1)
  402   COORDS(I3+1)=X
        COORDS(I3+2)=Y
        COORDS(I3+3)=Z      
        if(debug) print *,'LSQIO i3 coords '
     .      ,i3,coords(i3+1),coords(i3+2),coords(i3+3)
  100 CONTINUE                 

  101 RETURN
      END
