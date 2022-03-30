      SUBROUTINE SHWPRT ( PJD, LAMBDA, ECLTYP, beta )
C
C     Print out the time of beginning and end of satellite eclipsing
C     R. W. King  30 March 1987
c     Mods B. Oral November 1991
c     Mods for svnav_read R King June 1998
C
C     PJD      PEP Julian date
C     LAMBDA   Multiplication factor for solar pressure (1.= no eclipse)
c     ecltyp   Type of eclipse: Earth (E) or Moon (M) 
c     beta     Angle between Sun and orbital plane: printed for information only

      implicit none

      include '../includes/dimpar.h'
      include '../includes/arc.h'
               
c     In common /radiation/ lamprt: Flag for printout (F = not eclipsing;  T = eclipsing )

      integer*4 iyr,idoy,imonth,iday,ihr,imin

      real*8 lambda,pjd,beta

      character*1 ecltyp
                         
      real*8 rad2deg
      data rad2deg/57.29577951308232d0/

c     Trap two cases:
c       1) Second interation (correction of predictor-corrector) at same time:
c           don't repeat printout
      if( pjdlast.eq.pjd) goto 50
c       2) Backward integration has finished and integrator has skipped back
c          to the initial epoch:  identify by seeing if the current epoch
c          differs from the last by more than a step-size and treat the case
c          like a new beginning, leaving end_eclipse at zero to be trapped later
      if( dabs(pjdlast-pjd).gt.(diint+.1)/86400.d0)  lamprt = .false.

      if( lamprt ) goto 20

c     Last value of LAMBDA was 1.0 (not eclipsing), see if now eclipsing

      if( lambda.eq.1.0d0)  goto 50
C     Now eclipsing, so print
      call pjdhms(PJD,iyr,idoy,imonth,iday,ihr,imin)
      if ( first ) then
        write(iarh,5)
    5   format(22X,'PJD',6x,'Lambda',3x,'Yr',1x,'DOY',2x,'Mo',2x,
     .      'Dy',2x,'Hr',1x,'Min',1x,'PRN',4x,'Earth/Moon')
      else
        write(iarh,'(a80)')
      endif   
      neclipse = neclipse + 1
      eclipse_start(neclipse) = pjd
      eclipse_end(neclipse) = 0. 
      eclipse_type(neclipse) = ecltyp    
      eclipse_beta(neclipse) = beta * rad2deg 
      write(iarh,10) pjd,lambda,iyr,idoy,imonth,iday,ihr,imin,iprn
     .           ,eclipse_type(neclipse),eclipse_beta(neclipse)
   10 FORMAT(1X,'Start eclipse: ',F14.5,1x,F6.3,i5,5i4,1x,'PRN ',i2
     .   ,2x,a1,' beta=',f6.2)
      lamprt= .true.
      first = .false.
      goto 50

C     Last value of LAMBDA was < 1.0 (eclipsing)

   20 if( lambda.eq.0.d0 ) goto 50
C       If in total eclipse, no printout
      if( lambda.eq.1.d0 ) goto 30
      call pjdhms( pjd,iyr,idoy,imonth,iday,ihr,imin )
      write(iarh,25) pjd,lambda,iyr,idoy,imonth,iday,ihr,imin,iprn
     .          , eclipse_type(neclipse)
   25 format(1x,'Eclipsing:     ',F14.5,1x,F6.3,i5,5i4,1x,'PRN ',i2
     .    ,2x,a1)
      goto 50
   30 call pjdhms( pjd,iyr,idoy,imonth,iday,ihr,imin ) 
      eclipse_end(neclipse) = pjd 
      if( neclipse.eq.0 ) neclipse = 1
      eclipse_type(neclipse) = ecltyp
      write(iarh,35) pjd,lambda,iyr,idoy,imonth,iday,ihr,imin,iprn
     .          , eclipse_type(neclipse)
   35 format(1x,'End eclipse:   ',F14.5,1x,F6.3,i5,5i4,1x,'PRN ',i2
     .      ,2x,a1)
      lamprt= .false.
C
   50 pjdlast = pjd
      return
      end


