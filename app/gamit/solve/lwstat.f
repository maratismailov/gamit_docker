C
      SUBROUTINE LWSTAT
C
C     Add station weight constraints to normal matrix for loose solution
C
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'
      include 'parameters.h'

      integer ictyp,indxs,ind,jj,i,j,k

      real*8 covplr(6),curvm,curvn,f,e2,coslat,sinlat,radius,term,temp      

c     debug
c      integer*4 jj1,jj2,jj3,ind1,ind2

c     initialization to avoid compiler warning
      curvm = 0.d0
      curvn = 0.d0           

C type of coordinates - default is geodetic (not geocentric)
C
      ICTYP=2
      F=1.D0/FINV
      E2=2.D0*F-F*F


c     print *,'LWSTAT: ',(stat_apr2(1,j),j=1,3)

c print and write to q-file loose-solution constraints

          if( logprt ) write(6,110) (i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .                 ,(stat_apr2(i,j),j=1,3),i=1,nsite)
          write(10,110) (i,rlabel((i-1)*3+1)(1:4),sitnam(i)
     .                  ,(stat_apr2(i,j),j=1,3),i=1,nsite)
  110     format(
     .    ' A priori coordinate errors in kilometers',/,
     .    '                          Latitude  Longitude  Radius',/
     .    , 100(i3,1x,a4,2x,a12,2x,3f10.5/))


C  assumption is that site coordinates are first in parameter list

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
C
       DO 200 J=1,3
          JJ=(J*J-J)/2
       DO 200 K=1,J
       IND=K+JJ
       COVPLR(IND)=0.D0
C  convert latitude and longitude to radians
       IF(J.NE.K) GO TO 200
       IF(J.EQ.1.AND.ICTYP.EQ.2) RADIUS=CURVM
       IF(J.EQ.2.AND.ICTYP.EQ.2) RADIUS=CURVN
       IF(J.EQ.1) COVPLR(IND)=stat_apr2(I,J)/RADIUS
       IF(J.EQ.2)
     1    COVPLR(IND)=stat_apr2(I,J)/(RADIUS*COSLAT)
       IF(J.EQ.3) COVPLR(IND)=stat_apr2(I,J)
C  convert error to variance
       COVPLR(IND)=COVPLR(IND)*COVPLR(IND)
       if(stat_apr(i,j).gt.0.d0) then
         temp = stat_apr2(i,j)/stat_apr(i,j)
       else
c        no station weights in constrained solutions
         temp=1.d0
       endif

c   weight increment = (new weight) - (old weight)
c                    = 1.0/(new variance) - 1.0/(old variance)
c                    = (1.0 - temp*temp)/(new variance)
c   new weight = old weight + weight increment

c      print *,'LWSTAT temp: ',temp
       COVPLR(IND)=(1.D0-TEMP*TEMP)/COVPLR(IND)
  200  CONTINUE
C
CD      WRITE(6,500)
CD 500  FORMAT(/)

C  Add weight increment to normal matrix
       
      call addwgt( indxs,6,3,covplr )
     
c      DO 300 J=1,3
c         JJ1=(J*J-J)/2
c         JJ2=J+INDXS
c         JJ3=(JJ2*JJ2-JJ2)/2
c         DO 300 K=1,J
c            IND1=K+JJ1
c            IND2=K+INDXS+JJ3  
c            print *,'LWSTAT: data  ind2 alc ',ind2,alc(ind2)
c            if (l2flag.eq.1.or.l2flag.eq.3.or.l2flag.eq.4 ) then
c               ALC(IND2)=ALC(IND2)+COVPLR(IND1)  
c           print*,'LWSTAT: indx ind1 ind2 alc ',indx,ind1,ind2,alc(ind2)
c            else
c               A(IND2)=A(IND2)+COVPLR(IND1)
c            endif
c      print*,'LWSTAT: j,k,ind1,covplr,ind2,alc = '
c     .                ,j,k,ind1,covplr(ind1),ind2,alc(ind2)
c  300 CONTINUE
C
  100 CONTINUE
C
      RETURN
      END
