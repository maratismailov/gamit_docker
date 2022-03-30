C
      Subroutine WSTAT
C
C     ADD STATION COORDINATE WEIGHTS TO NORMAL MATRIX

      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer ictyp,indxs,ind,jj,i,j,k

      real*8 covplr(6),curvm,curvn,f,e2,coslat,sinlat,term,radius

c     in parameters.h
c     real*8 preval(maxprm),sigma(maxprm)
c     integer*4 idms(maxprm),islot1(maxprm),free(maxprm)
c     character*20 rlabel(maxprm)


c     initialization to avoid compiler warning
      curvm = 0.d0
      curvn = 0.d0

c      WRITE(6,999)
c  999 FORMAT(' WEIGHTING STATION COORDINATE PARAMETERS')
C
C type of coordinates - default is geodetic (not geocentric)
C
      ICTYP=2
      F=1.D0/FINV
      E2=2.D0*F-F*F
C
C  assumption is that site coordinates are first in parameter list
C
C  loop over all stations
      DO 100 I=1,nsite
       INDXS=3*(I-1)
CD      WRITE(6,*) I,INDXS
C If station coordinates not free parameters - skip
C  assumption here that all coordinates of a site are fixed or free
       IF(FREE(INDXS+1).EQ.0) GO TO 100
C
C  form lat-long-radius(height) covariance matrix (diagonal)
C
       COSLAT=DCOS(PREVAL(INDXS+1))
       SINLAT=DSIN(PREVAL(INDXS+2))
       RADIUS=PREVAL(INDXS+3)
C set geocentric coordinates flag
       IF(RADIUS.GT.6000.D0) ICTYP=1
       IF(ICTYP.EQ.2) then
         TERM=DSQRT(1.D0-E2*SINLAT*SINLAT)
C units are meters
         CURVM=(SEMI*(1.D0-E2))/(TERM*TERM*TERM)
         CURVN=SEMI/TERM
       end if

      DO 200 J=1,3
         JJ=(J*J-J)/2
       DO 200 K=1,J
       IND=K+JJ
       COVPLR(IND)=0.D0
C  convert latitude and longitude to radians
       IF(J.NE.K) GO TO 200
       IF(J.EQ.1.AND.ICTYP.EQ.2) RADIUS=CURVM
       IF(J.EQ.2.AND.ICTYP.EQ.2) RADIUS=CURVN
       IF(J.EQ.1) COVPLR(IND)=stat_apr(I,J)/RADIUS
       IF(J.EQ.2)
     1    COVPLR(IND)=stat_apr(I,J)/(RADIUS*COSLAT)
       IF(J.EQ.3) COVPLR(IND)=stat_apr(I,J)
C  convert error to variance
       COVPLR(IND)=COVPLR(IND)*COVPLR(IND)
C  invert to get weight
       COVPLR(IND)=1.D0/COVPLR(IND)
  200 CONTINUE
C
CD      WRITE(6,500)
CD 500  FORMAT(/)

C  Add weight matrix to normal matrix

      call addwgt( indxs,6,3,covplr )
      
c      DO 300 J=1,3
c      JJ1=(J*J-J)/2
c      JJ2=J+INDXS
c      JJ3=(JJ2*JJ2-JJ2)/2
c      DO 300 K=1,J
c      IND1=K+JJ1
c      IND2=K+INDXS+JJ3
c      a(IND2)=a(IND2)+COVPLR(IND1)
c      ALC(IND2)=ALC(IND2)+COVPLR(IND1)
CD      WRITE(6,*) IND1,IND2,COVPLR(IND1)
c  300 CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
