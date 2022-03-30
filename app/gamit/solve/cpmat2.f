      Subroutine CPMAT2( icall, phi,elev, pmat )
c
c PURPOSE: subroutine to construct the carrier beat covariance
c          matrix, given as input several different weighting
c          functions/strategies.....
c
c Strategy (1) BASELINE DEPENDENT WEIGHTING:
C              CONSTRUCT CARRIER-BEAT PHASE COVARIANCE MATRIX
C              USING THE MODEL : A**2 + (B**2)(DISTANCE**2)
C              A IS GIVEN IN MM, B IN MM/KM, DISTANCE IN KM
c
c Strategy (2) UNIFORM WEIGHTING
c              Construct a carrier-beat phase covariance matrix
c              using the model: A**2
c
c Strategy (3) ELEVATION ANGLE DEPENDENT WEIGHTING:
c              Construct a carrier-beat phase covariance matrix
c              using the model: A**2 + C**2/(sinE)**2
c
c  PARAMETERS:
C          IN: PHI -    MATRIX OF CARRIER-BEAT PHASES USEFUL FOR
C                       AN EPOCH OF DOUBLE DIFFERENCES
C              COORDS - VECTOR OF CARTESIAN COORDINATES
C              nsite  - NUMBER OF STATIONS (ROWS IN PHI)
C              nsat   - NUMBER OF SATELLITES (COLUMNS IN PHI)
C              APHI   - CONSTANT TERM IN ERROR MODEL (L1 cycles)
C              BPHI   - DISTANCE PROPORTIONAL TERM IN ERROR MODEL 
c                       (bphi or bkappa) (unitless)
c              ELEV   - Elevation angles from stations to each satellite
c              ICALL  - data error or ionosphere VCV switch
C
C         OUT: PMAT   - COVARIANCE MATRIX OF CARRIER-BEAT-PHASES,
C                       PMAT IS TRI-DIAGONAL IF OBSERVATIONS ARE
C                       ORDERED BY STATION
C                       I.E., STATION 1, SATELLITES 1,2,3, ETC.
c                       Stored in calling program OPERA as CPHI
c                       or in calling program OPERA2 as CPHIK
C
c SUBROUTINES CALLED: distnc, matprt
c
c CREATED: Don't Know               LAST MODIFIED:  1st SEPT 1995
c
c AUTHOR: someone else and S McClusky.
c
      implicit none

      include '../includes/dimpar.h'
      include 'solve.h'        
      include 'parameters.h'
                                                  
      integer*4 istat1,isat1,ielem,istat2,isat2,iobs,ir,i,j

      real*8 lmax,alpha,beta2,cc,radius,smax,term,sechmx,aterm,bterm
     .     , dist,dummy

C     DIMENSIONED FOR MAXSIT x MAXSAT SIMULTANEOUS CARRIER-BEAT-PHASES
      real*8 PHI(MAXSIT,MAXSAT)
      integer IKEY(MAXOBS,2)
      real*8 PMAT(MAXWM1),elev(maxsit,maxsat)
                    
c     a priori data noise and model

      character*4 icall
      character*256 message

      logical debug/.false./

      if(debug) print *,'CPMAT2 nsite nsat,elev(2,4) ',elev(2,4) 

c these defaults (see Schraffrin and Bock), may be changed below
      alpha = 1.d0
      beta2 = 0.d0     

c set this to avoid a compiler warning
      cc = 0.d0

c earth radius in km
        radius=6378.137d0

      if (err_mod(1) .eq. 'baseline' .or. icall .eq. 'iono') then

c maximum distance in km for ion-constraint formulation 
        smax=10000.d0 
c       maximum arc length
        lmax=smax/radius
        cc=dlog(1.d0+dsqrt(2.d0))/lmax
        term=cc*lmax
        sechmx=2.d0/(dexp(term)+1.d0/dexp(term))

c constant and distance-proportional terms in L1 cycles  
        aterm=aphi
c       bphi is dimensionless; convert max distance to cycles    
c        print *,'bphi2 ',bphi**2
c        print *,'smax cyc ',smax/190.29367d-6
c        print *,'sechterm ',1.d0/(1.d0-sechmx)
        bterm = (bphi**2) * (smax/190.29367d-6)**2 / (1.d0-sechmx)
c** Old formulation with bkappa (=bphi here) converted to cycles in read_bfl_sess 
c**        bterm=(bphi*1.d4)/(1.d0-sechmx)/conv 
c**        beta2=aterm*aterm+bterm*bterm  
c** Tthere are two errors:  1) a gratuitous factor of 1/95.146 applied
c**                         2) (1. - sechmx) should not be squared in beta2
c**   With these changes, beta2 gets smaller by a factor of 50.
        beta2 = aterm*aterm + bterm  
c        print *,'CPMAT2 ION  cc term sechmx ',cc,term,sechmx
c        print *,' aphi bphi bterm beta2 ',aphi,bphi,bterm,beta2  

c Case of zero constraint on ionosphere
        if(aphi.gt.0.d0.or.bphi.gt.0.d0) then
          alpha=1.d0-(aphi*aphi)/(beta2)
        else
          alpha=1.d0
        end if

      endif
c
C DETERMINE IKEY ARRAY (KEY TO OBSERVATIONS WHICH ARE STORED IN PHI)
C   SORTED BY STATION (BY ROW OF PHI)
C   FIRST INDEX - STATION
C   SECOND INDEX - SATELLITE
      IOBS=0
      DO 5 I=1,nsite
        DO 6 J=1,nsat
          IF(PHI(I,J).EQ.0.D0) GO TO 6
          IOBS=IOBS+1
          IKEY(IOBS,1)=I
          IKEY(IOBS,2)=J
    6   CONTINUE
    5 CONTINUE
       if(debug) then
         print *,'CPMAT2 iobs ikey 1 ',iobs,(ikey(i,1),i=1,iobs)       
         print *,'CPMAT2 iobs ikey 2 ',iobs,(ikey(i,2),i=1,iobs)
       endif
C
C LOOP OVER NON-ZERO OBSERVATIONS
      DO 10 I=1,IOBS
        IR=(I*I-I)/2
        ISTAT1=IKEY(I,1)
        ISAT1=IKEY(I,2)
        DO 11 J=1,I
          IELEM=IR+J     
cd          print *,'i ir j ielem ',i,ir,j,ielem 
          PMAT(IELEM)=0.D0
          IF(I.EQ.J) GO TO 12
          if (i.ne.j .and. icall.eq.'data'.and.
     .        err_mod(istat1).ne.'baseline') goto 11
          ISTAT2=IKEY(J,1)
          ISAT2=IKEY(J,2)
          IF(ISTAT1.EQ.ISTAT2) GO TO 11
          IF(ISAT1.NE.ISAT2) GO TO 11
C OFF DIAGONAL, DISTANCE-DEPENDENT ELEMENTS
          CALL DISTNC(COORDS,ISTAT1,ISTAT2,DIST)
CD        WRITE(6,*) 'DIST',DIST
CD        WRITE(26,*) 'DIST',DIST
c         PMAT(IELEM)=-(BPHI*DIST)*(BPHI*DIST)
          PMAT(IELEM)=(alpha*beta2)/DCOSH(CC*(DIST/RADIUS)) 
c          If( ielem.eq.4) 
c     .      print *,' Off-diag ielem dist pmat ',ielem,dist,pmat(ielem)
          GO TO 11
C DIAGONAL ELEMENTS
   12     CONTINUE
          if (err_mod(istat1).eq.'baseline' .or. icall.eq. 'iono') then
            PMAT(IELEM)=beta2  
c            print *,' ION ielem pmat ',ielem,pmat(ielem)
          else
            pmat(ielem) = sit_err(istat1)**2 + sat_err(isat1)**2
c apply elevation angle dependant weighting
            if ( err_mod(istat1) .eq. 'uniform' ) then
              pmat(ielem) = pmat(ielem)   
c              print *,' DATA ielem pmat ',ielem,pmat(ielem)
            elseif ( err_mod(istat1) .eq. 'elevation' ) then
              pmat(ielem) = pmat(ielem) +
     .        sit_elv(istat1)**2/(dsin(elev(istat1,isat1)))**2   
              if(debug) print *
     .         ,' istat1 isat1 err_mod ielem elev sit_elv pmat '
     .           ,istat1,isat1,err_mod(istat1),ielem,elev(istat1,isat1)
     .            ,sit_elv(istat1),pmat(ielem)
            else                               
              write(message,'(a,a10,a)') 'Unknown data weighting model('
     .        ,err_mod(istat1),') uniform data weighting applied'
              call report_stat('WARNING','SOLVE','cpmat2',' ',message,0)
            endif
          endif
   11   CONTINUE      
cd       print *,'i ir istat1,isat1 ielem ',i,ir,istat1,isat1,ielem
   10 CONTINUE    
C                 
c**  The array 'work' is used only for debug; to avoid a compiler warning for
c    an unused calling argument, add this meaningless statement:
      dummy = work(1)
CD      CALL MATPRT(PMAT,WORK,IOBS,1)
c      if(debug) stop 1
      RETURN
      END
