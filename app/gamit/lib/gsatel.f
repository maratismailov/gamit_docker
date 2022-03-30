Copyright (c) Massachusetts Institute of Technology,1986. All rights reserved.

      Subroutine GSATEL ( ispeed,trun,iunit,ksat,rvec,ytrp,yy
     .                  , nsat,sdelt,nintrs,ji0,jil,iy1,iy2
     .                  , jlast,jnow,iendf,nepcht )

c     Based on routine SATEL in MODEL.  This version documents the
c     interpolation logic and corrects a few bugs.  It is intended
c     to replace SATEL eventually.   Note that the meaning of the control
c     variable ISPEED is different from the current (91.12.26) SATEL, and
c     that the code for partials of velocity has been removed.
c     The routine is still inefficient in the sense that if partials are
c     on the file, interpolation Y-vectors are set up for them even if they
c     are not needed.  This is to avoid have to know on the first call (e.g.
c     in MODEL) whether velocities and partials will be needed in a subsequent
c     call.    R. King 911226

c       ISPEED = 1   Position only  (first call in MODEL)
c              = 2   Position and velocity (TTOG, TTOICS, TTONGS, ORBFIT)
c              = 3   Velocity only after position
c              = 4   Position, velocity, and partials (ORBFIT)
c              = 5   Velocity and partials after position (second call in MODEL)

c       NIINTRS = 3  T-file has only position (no partials)
c               = 30 T-file has position (1-3) and 9 partials

      implicit none

      include '../includes/dimpar.h'

      character*80 prog_name

      integer*4 iunit,isat,ksat,ispeed,nintrs,nepcht
     .        , ji0,jil,jnow,jlast,iendf,iy1,iy2,kup,kdwn,nsat
     .        , ir,jyvc,jcrd,jstart,jend,j,k,jj,len,rcpar,i


      real*8 sdelt,sdelti,fact,yy,rvec,trun,p,q,p2,q2,sum1,sum2
     .        , evcf,ytrp

      logical debug/.false./

      DIMENSION RVEC(MAXYT2)
      DIMENSION YTRP(5,2,MAXYT2,MAXSAT),YY(10,MAXYTP,MAXSAT)
C
C
C-------------------------------------------------------------------------
      DIMENSION EVCF(5,5),FACT(9)
C       INTERPOLATION COEFFICIENTS
      DATA EVCF/
     1 1.7873015873015873D0,
     1-0.4960317460317460D0, 0.1206349206349206D0,
     2-0.1984126984126984D-01, 0.1587301587301587D-02,
     3-0.9359567901234568D0,
     4 0.6057098765432098D0,-0.1632716049382716D0,
     5 0.2779982363315696D-01,-0.2259700176366843D-02,
     6 0.1582175925925926D0,
     7-0.1171296296296296D0, 0.4606481481481481D-01,
     8-0.8796296296296296D-02, 0.7523148148148148D-03,
     9-0.9755291005291005D-02,
     1 0.7605820105820106D-02,-0.3505291005291005D-02,
     1 0.8597883597883598D-03,-0.8267195767195767D-04,
     2 0.1929012345679012D-03,
     3-0.1543209876543210D-03, 0.7716049382716048D-04,
     4-0.2204585537918871D-04, 0.2755731922398589D-05/
C
      DATA FACT
     1 /1.0D0,3.0D0,5.0D0,7.0D0,9.0D0,2.0D0,4.0D0,6.0D0,18.0D0/
C-------------------------------------------------------------------------

c     get calling program name and m-file name for report_stat
      len = rcpar(0,prog_name)
                      

c Cannot interpolate within 5 epochs of the end of the T-file

      if( trun .lt. 4.999d0*sdelt .or.
     .    trun .gt. (nepcht-4.999d0)*sdelt ) then
          print *,'trun sdelt nepcht ',trun,sdelt,nepcht 
         call report_stat('FATAL',prog_name,'lib/gsatel',' '
     .       ,'Requested time out of T-file range',0)    
      endif

CRIA Comment: TRUN is time since first epoch of ephemeris
CRWK Comment: No, it is time + one epoch
CRWK Changed in MODEL, TTOG, TTONGS, and TTOICS to be time since first epoch,
c    but offset retained in GSATEL and SATEL at least temporarily, for safety
c    6/13/88; 12/28/91
C
      SDELTI= 1.D0/SDELT
      P=(TRUN+SDELT)*SDELTI

C RWK Comment:  JNOW is one less if P is exactly an integer
C               but the results seem to be OK (Check this)
      JNOW=P


c  Come here after each read to see if additional records are needed in storage

990   CONTINUE
      IF(IENDF.NE.0) GO TO 1000
      IF((JNOW+5.LE.JLAST).AND.(JLAST.GE.10)) GO TO 1000


c  A read is necessary:  update the coordinates-vector (YY) indices

      JLAST=JLAST+1
      JI0=MOD(JI0,10)+1
      JIL=MOD(JI0+3,10)+1
      READ(IUNIT,END=995)
     $  ((YY(JIL,IR,ISAT),IR=1,NINTRS),ISAT=1,NSAT)
c     if( trun.gt.13000.d0 ) stop
c     print *,'jlast,jnow,jil ',jlast,jnow,jil
c      WRITE (6,8000) yy(jil,1,1)
c8000  FORMAT (1X,'yy(jil,1,1):',f10.3)
c     next statement added by rwk 911227 : necessary to avoid undefined
c     yy(1-4,6-10) at beginning of T-file
      if( jlast.lt.10 ) goto 990
      IF((JNOW+8.LE.JLAST).AND.(JLAST.GE.12)) GO TO 990



c  We now have sufficient values in storage:
c  setup the interpolation Y-vector

c     JYVC here is the component pointer within both the coordinate
c     (YY) vector and the interpolation Y-vector

c     The maximum Y-vector index to be set up
c         = 3               for position only
c     or  = NINTRS + 3 (33) for position and (velocity or partials)

c     The velocity Y-vectors are stored in slots NINTRS+1 to NINTRS+3 (nominally 31-33).
c     We could avoid setting up all the Y-vectors if we either 1) knew that velocities
c     and/or partials would not be needed later, or 2) assigned a variable to remember
c     whether the non-position Y-vectors have been set up.

      IY1=IY2
      IY2=MOD(IY2,2)+1

      if(debug.and.ksat.eq.1 ) then   
        write(*,'(a,f8.1,f10.5,i4)') 'GSATEL trun p ispeed '
     .      ,trun,p,ispeed
        write(*,'(a)') 'yy-vectors :'
        do i=1,10 
          write(*,'(3f16.8)') (yy(i,j,ksat),j=1,3)
        enddo 
        if( ispeed.gt.3 ) then      
           write(*,'(a)') 'yy vectors for partials '
           do j=4,nintrs
             write(*,'(i3,10d15.4)') j,(yy(i,j,ksat),i=1,10)  
           enddo
        endif
      endif 

c     print *,'jlast,jnow,jil,iy1,iy2 ',jlast,jnow,jil,iy1,iy2

      DO 520 ISAT=1,NSAT
      DO 520 jyvc=1,nintrs
        DO 500 J=1,5
        YTRP(J,IY2,jyvc,ISAT)=0.0D0
        YTRP(J,IY2,jyvc+NINTRS,ISAT)=0.0D0
        DO 400 K=5,2,-1
          KUP=MOD(JI0+K+8,10)+1
          KDWN=MOD(JI0-K+10,10)+1
          YTRP(J,IY2,jyvc,ISAT) = YTRP(J,IY2,jyvc,ISAT)
     .               + EVCF(K,J)*(YY(KUP,jyvc,ISAT)+YY(KDWN,jyvc,ISAT))
  400     continue
        YTRP(J,IY2,jyvc,ISAT)=
     .             YTRP(J,IY2,jyvc,ISAT) + EVCF(1,J)*YY(JI0,jyvc,ISAT)
      if (ytrp(j,iy2,jyvc,isat).gt.1.d7) then
         print *,j,iy2,jyvc,isat,ytrp(j,iy2,jyvc,isat)
         print *,kup,kdwn,ji0,evcf(1,j),yy(ji0,jyvc,isat)
         print *,trun,jnow
         print *,(yy(jj,jyvc,isat),jj=1,10)
         print *,(yy(jj,1,isat),jj=1,10)
         print *,'Possible bad value in Y-vector, stop in GSATEL'
         print *,' See Bob King'
         call report_stat('FATAL',prog_name,'lib/gsatel',' '
     .      ,'Possible bad value in Y-vector--see screen output',0)
       endif     
cd       if( debug.and.ksat.eq.1 ) 
cd     .   print *,'j,iy2,jyvc,nintrs,isat,fact(j),sdelt,ytrp: '
cd     .         , j,iy2,jyvc,isat,fact(j),sdelt,ytrp(j,iy2,jyvc,isat)
cd     .         , ytrp(j,iy2,jyvc+nintrs,isat)
        YTRP(J,IY2,jyvc+NINTRS,ISAT) =
     .             YTRP(J,IY2,jyvc,ISAT)*FACT(J)*SDELTI
  500   continue
  520 continue
      GO TO 990
c     Go see if more records need to be read in

c     An end-of-file was reached on the tape, don't read any more records
c     (This shouldn't happen because the interpolation can no longer use 10 points)
  995 IENDF=1
      JNOW=JNOW-1



c  Now the interpolation Y-vector is set up:  do the interpolation

 1000 P=P-JLAST+5
      Q=1.D0-P
      P2=P*P
      Q2=Q*Q

c     We're within the range of the setup Y-vectors
c     Determine what components are needed

c                   Control flag   RVEC index    Y-vector index
c                      ISPEED       JCRD            JYVC
c  Position only          1        1 - 3            1 - 3
c  Pos + Velocity         2        1 - 6            1 - 3, NINTRS+1 - NINTRS+3
c  Vel. only              3        4 - 6            NINTRS+1 - NINTRS+3
c  Pos, Vel., Partials    4        1 - NINTRS+3     1 - NINTRS+3
c  Vel. and Partials      5        4 - NINTRS+3     1 - NINTRS+3

      jstart = 1
      if( ispeed.eq.3 .or. ispeed.eq.5 ) jstart = 4
      jend   = 3
      if( ispeed.eq.2 .or. ispeed.eq.3 ) jend = 6
      if( ispeed.ge.4 ) jend = nintrs + 3

      do 2500 jcrd = jstart,jend

        if( jcrd.le.3 ) jyvc = jcrd
        if( jcrd.ge.4 .and. jcrd.le.6 ) jyvc = jcrd - 3 + nintrs
        if( jcrd.gt.6 ) jyvc = jcrd -3

        SUM1=0.0D0
        SUM2=0.0D0
        DO 2000 J=5,2,-1
         SUM1=P2*(SUM1+YTRP(J,IY2,jyvc,KSAT))
         SUM2=Q2*(SUM2+YTRP(J,IY1,jyvc,KSAT))
 2000    continue

      if (jcrd.ge.4 .and. jcrd.le.6 ) then
c       velocity
        RVEC(jcrd) = (YTRP(1,IY2,jyvc,KSAT)+SUM1)
     .             - (YTRP(1,IY1,jyvc,KSAT)+SUM2)
       else
c       position or partial derivative of position
        RVEC(jcrd) = P * (YTRP(1,IY2,jyvc,KSAT)+SUM1)
     .             + Q * (YTRP(1,IY1,jyvc,KSAT)+SUM2)
       endif

      if( debug.and.jnow.eq.6 .and. ksat.eq.1 )
     .   write(6,'(a,/,2i3,7f16.8)') 
     .    'jcrd,jyvc,p,q,sum1,sum2,ytrp(1,1-2,jyvc,1),rvec(jcrd)'
     .   ,jcrd,jyvc,p,q,sum1,sum2,(ytrp(1,j,jyvc,1),j=1,2),rvec(jcrd)

 2500 continue

      RETURN
      END
