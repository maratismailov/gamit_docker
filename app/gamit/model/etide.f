      Subroutine ETIDE( jd,t,sidtm,tdtgpst  
     .                , evec,sun,moon,tidevec,atmvec )

c     Compute tidal corrections to station coordinates.
c     Written by R. King 14 Apr 88 with substantial modifications
c     by P. Morgan (Oct 93), S. McClusky (Dec 93), P. Tregoning (May 95),
c     Last modified by R. King Dec 98.          

c     This version assumes that the reference system for the calculations in 
c     MODEL is inertial since this makes it simpler and Earth-fixed doesn't
c     yet work.  If Earth-fixed is allowed, see the older code in       
c     --rwk 981222
c
c     PT040107: added calculations for atmospheric pressure loading. While it
c               isn't strictly a "tidal" effect, the calculations are pretty
c               much the same and it seems appropriate to do it here!
c
c     CW050118: added IERS2003 compliant solid Earth tide using code from
c               Veronique Dehant, ftp://omaftp.oma.be/dist/astro/dehant/IERS/  
c
c     RWK050217: solid-Earth selection now made in FIXDRV:
c                IERS92  : IERS Conventions (1996), McCarthy et al., TN 21)
c                          except without latitude-dependent love numbers;
c                          GAMIT default pre 10.2
c                IERS03  : IERS Conventions (2003) McCarthy and Petit (2004)
c                          This produces a "conventional tide free" solution.
c                          The effect of the "permanent tide" is not
c                          "uncorrected" to achieve a "mean tide" solution.
c     
c     PT050509: add the atm tidal displacements (S1 and S2). These will be computed
c               separately but then added to the tidevec(3) output variable by way
c               of adding datotal to the summation of E-fixed tidal corrections.

c     Input Parameters
c     ****************

c     jd, t    PEP Julian day and seconds-of-day for evaluation (GPST)
c     sidtm    Siderial time (radians) (also passed back by calls to ROTSNP if called
c     tdtgpst  TDT - GPST (s)   

c     ietide   binary coded integer for which tides are to be computed     (model.h)
c              1 = solid earth tides
c              2 = K1 frequency dependent earth tide
c              4 = Pole tide
c              8 = Ocean Tide
c             16 = Remove mean pole position in computing pole tide
c             32 = Atmospheric tides
c                Hence, ietide = 63 means all tides are turned on, with the
c                mean pole removed from the pole tide.
c                There is a  LOGICAL function KBIT(ivariable,ibit) in the library 
c                which return the value .true. if ibit is set.  

  
c     frame    inertial frame (J2000 or B1950)
c     precmod  precession model (IAU68 or IAU76)  
c     isptide  binary coded control for diurnal pole (bit 1) and UT1 (bit 2)
c     simulation  T/F to indicate whether this is simulatation mode

c     Units from common /units/ in model.h     
c     iscrn    logical unit for debug to screen  
c     iut1     logical unit for ut1. file
c     ipole    logical unit for pole. file
c     inut     logical unit for nutation (nutabl.) file 
      
c     etidemod  c*8   'IERS92  ' for IERS 1992, 'IERS03  ' for IERS2003  (model.h)
c     otide_source c*3   'OSO    ' for Scherneck, 'NAO    ' for Matsumoto (model.h)
c     otides(54,6) ocean tide components for station (3 amplitudes and 3 phases) 
c                  11 constituents for Scherneck
c                  54 constituents for Matsumoto   
                 
c     evec0(3,2) Earth-fixed station coordinates (km) L1, L2 (model.h)
c     evecf(3)   Earth-fixed station coordinates in meters
c     evec(6)  inertial station coordinates (km)  
c     sun(6)   coordinates of Sun wrt Earth
c     moon(6)   coordinates of Moon wrt Earth

c     atmlod(3)  actual NEU displacements in mm (R*4, assigned to atmneu(3) R*8 here). (model.h)
c     atides(2,6)  (from common in model.h)  cosine and sine amplitudes of 
c                     atm tidal loading at S1 and S2 frequencies


c     Output
c     *******

c     tidevec(3)  tidal corrections in inertial system (km)
c     atmvec(3) atmospheric pressure loading corrections in inertial system (km)

      implicit none
                                    
      include '../includes/dimpar.h'   
      include '../includes/units.h'   
      include '../includes/global.h'
      include '../includes/model.h'

      integer*4 i,j,jd,idir   

c     debug
      integer*4 icounter

* MOD TAH 2002119: Renoved lmptid (logical) and replaced with
*     model name
      logical  kbit,debug

      character*8 ptd_model  ! Name of poletide model (NONE, IERS10, 
                             ! IERS20). Raplaces lmptid
             
      real*8 evec(6),evecf(3),tidevecf(3),tidevec(3),r(3),rmag
      real*8 atmtotal(3),atmvec(3),atmneu(3)
     .      ,atmtidneu(3),datotal(3),dtatm,freq(2)
      real*8 gme,gmm,gms,dot,pi,t,xjd,xmjd,oangle(54)
     .     , eloveh,elovel,gm,convd,told,tdtgpst,utc
      real*8 factor1,factor2,factor3, dptotal(3),dototal(3),docmc(3)
      real*8 ur(3), uevec(3), dftotal(3), dstotal(3)
      real*8 dttotal(3),dttotali(3)
      real*8 xrot(3,3), xdtrot(3,3),geoc(3),dh,etidefixed(3)
      real*8 sidtm,sidtm1,evecsph(3),xpole,ypole
      real*8 rotmat(3,3),evecfsph(3)
      real*8 eveclat, eveclon,etideneu(3),geocn(3),ftideneu(3)
      real*8 ptideneu(3),otideneu(3),totalneu(3)
      real*8 sun(6),moon(6)
      real*8 taiutc

c     CW050118:  Solid Earth tide variables in existing GAMIT code are:
c       dstotal       : Inertial Solid Earth tide
c       etidefixed(3) : ECEF Solid Earth tide components
c       etideneu(3)   : NEU Solid Earth tide components
c       ftideneu(3)   : K1 Frequency dependent NEU components
c
c     Additional solid Earth tide variables for IERS2003 are:
c       sunf          : Position of sun in ECEF
c       moonf         : Position of moon in ECEF
c       suni          : Sun position in intertial space
c       mooni         : Moon position in intertial space
c       dstotalf      : Solid Earth tide in ECEF
c       fhr           : Fractional hour of the day 
c       date          : stores year, month, day, hour, minute
c       seconds       : stores output from jd_to_ymdhms routine
c       utctdt        : difference b/w UTC and TDT, should be 
c                       ~64.184 seconds at 1999.0
c
      integer*4 date(5)
      real*8 sunf(3), moonf(3), suni(3), mooni(3), dstotalf(3)
      real*8 fhr, seconds
      real*8 utctdt,utctai

c       Gravitational constants for Sun, Moon, Earth
      data gms/1.32712440D11/,gmm/4.9027993d3/,gme/3.98600448d5/
c       Conversion from degrees to radians
      data convd/0.017453292519943296d0/
c       Love numbers h and l
      data eloveh/.6090d0/,elovel/.0852d0/   
c      initialize trigger for debug plots
      data told/0.d0/

c      determine pi
      pi = convd*180.0d0

      debug = .false.
        
      if(debug ) then
        print *,'Entering ETIDE ietide simulation isptide '
     .    ,ietide,simulation,isptide
      endif
                              
c      Dummy statements to avoid warning until simulation code finished
c       if( simulation ) then
c         write(*,'(a)') 'Simulation on in ETIDE'
c       endif

c---- Uncomment below to debug ----
cd      debug = .true.
cd      data icounter /0/  
cd      icounter = icounter + 1
cd      if (icounter.gt. 30) stop
cd      open(unit=79,file='tides.plt',status='unknown')
cd    CW 050118: Changed plot file (for NEU components in mm)
cd      open(unit=79,file='tides.plt',status='unknown',access='append')
c---- Uncomment above to debug ----

c      Zero out tidal correction term variables

      do j=1,3
         dstotal(j)=0.0d0
         dftotal(j)=0.0d0
         dptotal(j)=0.0d0
         dototal(j)=0.0d0 
         datotal(j)=0.0d0 
         dttotal(j)=0.0d0
         dttotali(j)=0.0d0
      enddo             
              
c      Compute quantities needed for the tides

c        - solid-Earth components need only inertial quantities but Earth-fixed
c          quantities used for debug 
c        - pole and ocean tides need only Earth-fixed quantities 
c        - compute most quantities in both systems since easy
      
c     Get e-fixed site vector in meters from vector in km stored in model.h
      do i=1,3
        evecf(i) = evec0(i,1)*1.d3
      enddo

c     Inertial unit vector
      rmag = dsqrt(evec(1)**2+evec(2)**2+evec(3)**2)
      do i=1,3
        uevec(i) = evec(i)/rmag   
      enddo  
      
c     Spherical equivalents in both the inertial and E-fixed systems
      evecsph(2) = datan2(evec(2),evec(1))
      evecsph(3) = dsqrt(evec(3)**2 +(evec(1)*dcos(evecsph(2)) +
     .             evec(2)*dsin(evecsph(2)))**2)
      evecsph(1) = dasin(evec(3)/evecsph(3))
      evecfsph(2) = datan2(evec0(2,1),evec0(1,1))
      evecfsph(3) = dsqrt(evec0(3,1)**2 +(evec0(1,1)*dcos(evecfsph(2)) +
     .              evec0(2,1)*dsin(evecfsph(2)))**2)
      evecfsph(1) = dasin(evec0(3,1)/evecfsph(3))
                  
c     Convert spherical E-fixed to geodetic for freq-dependent, pole, or ocean  tide
      if ( kbit(ietide,2). or. kbit(ietide,3) .or. kbit(ietide,4) ) then
        call xyz_to_geod(rotmat, evecf, geoc)
c       convert co-latitude to latitude; note geoc(3) in meters
        geoc(1) =  pi/2.d0 - geoc(1)
      endif
      

c     Compute the solid-Earth tide components (vector formulation) 
      if (kbit(ietide,1)) then
      
c       OK, Solid Earth tides are to be applied, now distinguish between models
        if( etidemod(1:6).eq.'IERS92' 
c**     add this temporily to cover misnamed model 
     .       .or. etidemod(1:6).eq.'IERS96' ) then
c       OK, IERS1992 has been selected, use existing GAMIT code
        do ibody =1,2
          if (ibody.eq.1 ) then
            gm = gmm                         
            do i=1,3
              r(i) = moon(i)
            enddo
          else
            gm = gms
            do i=1,3
              r(i) = sun(i)
            enddo
          endif
          rmag = dsqrt(r(1)**2+r(2)**2+r(3)**2)
          do i=1,3
            ur(i) = r(i)/rmag 
          enddo
          factor1 = (gm * (dot(evec,evec))**2.0d0 )
     .                 / (gme * (dot(r,r))**1.5d0)
          factor2 = 3.0d0*elovel*(dot(ur,uevec))
          factor3 = 3.0d0*(eloveh/2.d0 -elovel) 
     .                 * ((dot(ur,uevec))**2) - eloveh/2.0d0
          do i=1,3
           dstotal(i)=dstotal(i) +factor1*factor2*ur(i)
     .                           +factor1*factor3*uevec(i)
          enddo
        enddo
        
c-------Debug for IERS1992 primary solid-Earth component
        if( debug ) then

          write(iscrn,'(2a)') '**Debug for IERS1992 tides:'
     .           ,'Coords and totals in km; components in mm'
          write(iscrn,'(a,i8,f12.3)') '  JD, sod: ',jd,t
          write(iscrn,'(a,3f16.8)') 'Inertial stn vector: '
     .           ,(evec(i),i=1,3)
          write(iscrn,'(a,3f16.8)') '  in spherical units:'  
     .           ,evecsph(1)/convd,evecsph(2)/convd, evecsph(3)  
          write(iscrn,'(a,a)') 'Solid Earth Tide Model :     '
     .      ,etidemod
          write(iscrn,'(a,3x,3f8.1)') 'Inertial solid-E dXYZ :'
     .           ,(dstotal(i)*1.d6,i=1,3)
c         rotate to Earth-fixed
          idir = -1                
          call rotsnp( idir,jd,t,tdtgpst,isptide
     .               , frame,precmod,iut1,ipole,inut
     .               , xrot,xdtrot,sidtm1,xpole,ypole ) 
          call matmpy (xrot,dstotal,etidefixed,3,3,1)
          write(iscrn,'(a,4x,3f8.1)') 
     .           'E-fixed solid-E dXYZ :'
     .           ,(etidefixed(i)*1.d6,i=1,3)
c         convert earth fixed x,y,z corrections to NEU components
c         kf routine rotate_geod expects units of meters: convert from km
          call rotate_geod (etidefixed, etideneu,'XYZ','NEU',evecf
     .                   ,geocn,rotmat)    
          do i = 1,3
             etidefixed(i) = etidefixed(i)/1000.d0 
             etideneu(i) = etideneu(i)/1000.d0
          enddo
          write(iscrn,'(a,10x,3f8.1)') 'Solid-E dNEU :  '
     .                         ,(etideneu(i)*1.d6,i=1,3) 
        endif
c-------End debug for IERS1992 primary solid-Earth component


c       Now, check to see if K1 freq dependend term is required under IERS1992
        if (kbit(ietide,2)) then
        
c         Compute the frequency-dependent (K1) component of the solid tide
          eveclat = geoc(1)
          eveclon = geoc(2)
          dh = -0.0000253d0*dsin(eveclat)*dcos(eveclat)*
     .          dsin(sidtm+eveclon)
          dftotal(1) = dh*dcos(evecfsph(2))
          dftotal(2) = dh*dsin(evecfsph(2))
          dftotal(3) = dh*dtan(evecfsph(1))  

c-------Debug for IERS1992 frequency-dependent component
        if( debug ) then
          write(iscrn,'(a,3f16.8)') 'Geodetic lat long ht(m): '
     .                , geoc(1)/convd,geoc(2)/convd,geoc(3)   
          write(iscrn,'(a,f16.8)') 'Sidereal time: ',sidtm
          write(iscrn,'(a,20x,f8.1)') 'Freq-dep radial tide: ',dh*1.d6
          write(iscrn,'(a,2x,3f8.1)') 'E-fixed freq-dep dXYZ : '
     .               , (dftotal(i)*1.d6,i=1,3) 
c         convert to NEU for comparison with other terms
          do i = 1,3
            dftotal(i) = dftotal(i)*1000.0d0
          enddo 
          call rotate_geod (dftotal, ftideneu,'XYZ','NEU',evecf
     .                     ,geocn,rotmat)
          do i = 1,3
            dftotal(i) = dftotal(i)/1000.0d0 
            ftideneu(i) = ftideneu(i)/1000.d0
          enddo
          write(iscrn,'(a,10x,3f8.1)') 'Freq dep dNEU : '
     .                           ,(ftideneu(i)*1.d6,i=1,3)
        endif
c-------End debug for IERS1992 frequency dependent component

        endif
c       End IERS1992 code, now check for IERS2003

        elseif ( etidemod(1:6).eq.'IERS03' .or. 
     .           etidemod(1:6).eq.'IERS10') then
        
c         OK, IERS2003 has been selected, get the job done using Dehant code
c         Firslty ensure freq-dependent components are set to zero as they
c         are not required, however still used in the formation of dttotal
          do i = 1,3
            ftideneu(i) = 0.0d0
            dftotal(i) = 0.0d0
          enddo

c         Make a copy of sun and moon position variables (velocity not required)
c         IERS2003 routine expects meters so convert from km
          do i = 1,3
            suni(i) = sun(i)*1000.0d0
            mooni(i) = moon(i)*1000.0d0
          enddo
         
c         OK, rotate to Earth-fixed frame for iers2003_etide routine
          idir = -1
          call rotsnp( idir,jd,t,tdtgpst,isptide
     .               , frame,precmod,iut1,ipole,inut
     .               , xrot,xdtrot,sidtm1,xpole,ypole ) 
          call matmpy (xrot,suni,sunf,3,3,1)
          call matmpy (xrot,mooni,moonf,3,3,1)
 
c         Need year, month, day and fractional hour of the day for IERS2003
c         routine.  Use the existing routine /gamit/lib/jd_to_ymdhms.f.
c         jd_to_ymdhms expects true Julian date (= PEP_JD - 0.5)
          xjd = jd - 0.5d0
          call JD_to_YMDHMS( xjd, date, seconds )
          
c         Compute the fractional hour of the day, FHR
c         Use t (which is seconds of day, GPST).
c         FHR must be in UTC for iers2003_etide.f.  
          fhr = (t - taiutc(jd) + 19.d0)/3600.0d0

          if ( etidemod(1:6).eq.'IERS03' ) then
          
c         Also compute differene b/w UTC and TDT required by iers2003_etide.f
             utctdt =  taiutc(jd) + 32.184d0
                    
c         Compute IERS2003 Earth tide by calling IERS2003_etide.f
             call iers2003_etide( evecf,date(1),date(2),date(3),fhr,
     .                         sunf,moonf,dstotalf,utctdt)

          elseif (etidemod(1:6).eq.'IERS10') then

c         Also compute differene b/w UTC and TAI required by iers2010_etide.f
             utctai =  taiutc(jd)

c.        Compute IERS2010 Earth Tide by calling IERS2010_etide.f
             call iers2010_etide( evecf,date(1),date(2),date(3),fhr,
     .                         sunf,moonf,dstotalf,utctai)
        
          endif
          
c         OK, rotate results back to inertial frame
          idir = 1
          call rotsnp( idir,jd,t,tdtgpst,isptide
     .               , frame,precmod,iut1,ipole,inut
     .               , xrot,xdtrot,sidtm1,xpole,ypole ) 
          call matmpy (xrot,dstotalf,dstotal,3,3,1)

c         OK, convert results and evecf back to kilometers 
          do i = 1,3
            dstotal(i) = dstotal(i)/1000.0d0
          enddo

          
c---------Debug for IERS2003 / IERS2010 solid-Earth tide
          if( debug ) then
             write(iscrn,'(2a)') '**Debug for tides and loading:'
     .             ,' Coords and totals in km; components in mm'
             write(iscrn,'(a,i8,f12.3)') '  JD, sod: ',jd,t
             write(iscrn,'(a,3f16.8)') 'Inertial stn vector: '
     .             ,(evec(i),i=1,3)
             write(iscrn,'(a,3f16.8)') '  in spherical units:'  
     .             ,evecsph(1)/convd,evecsph(2)/convd, evecsph(3)   
             write(iscrn,'(a,a)') 'Solid Earth Tide Model :     '
     .                   ,etidemod
             write(iscrn,'(a,2x,3f8.1)') 'Inertial solid-E dXYZ : '
     .             ,(dstotal(i)*1.d6,i=1,3)    
     
c            rotate to Earth-fixed
             idir = -1
             print *,'ETIDE rotations: jd t tdtg0st ',jd,t,tdtgpst
             print *,'                 isptide sidtm1 ',isptide,sidtm1  
             print *,'           frame,precmod,xrot ',frame,precmod,xrot
             call rotsnp( idir,jd,t,tdtgpst,isptide
     .                  , frame,precmod,iut1,ipole,inut
     .                  , xrot,xdtrot,sidtm1,xpole,ypole ) 
             call matmpy (xrot,dstotal,etidefixed,3,3,1)
             write(iscrn,'(a,3x,3f8.1)') 
     .                   'E-fixed solid-E dXYZ : '
     .                   ,(etidefixed(i)*1.d6,i=1,3)
     
c            convert earth fixed x,y,z corrections to NEU components
c            kf routine rotate_geod expects units of meters: convert from km
             do i = 1,3
                etidefixed(i) = etidefixed(i)*1000.0d0
             enddo
             call rotate_geod (etidefixed, etideneu,'XYZ','NEU',evecf
     .                        ,geocn,rotmat)
             print *,'evecf geocn rotmat ',evecf,geocn,rotmat
             do i = 1,3
                etidefixed(i) = etidefixed(i)/1000.d0 
                etideneu(i) = etideneu(i)/1000.d0
             enddo
             write(iscrn,'(a,12x,3f8.1)') 'Solid-E dNEU: '
     .           ,(etideneu(i)*1.d6,i=1,3) 
          endif
c---------End debug for IERS2003 solid-Earth tide

        endif
c       End IERS2003 code 

       endif
c      End of solid-Earth tide code



c     Compute the pole tide  

      if (kbit(ietide,3)) then
           
c       compute the pole tide as per IERS Tech. Note 3 or 21
c       computation is in millimeters,dptotal returned in mm's
* MOD TAH 200219: Generate name to be passed.  Work backwards through models                             
C       lmptid = .false.
C       if( kbit(ietide,5) ) lmptid = .true.
        if( kbit(ietide,7) ) then  ! IERS20
            ptd_model = 'IERS20' 
        elseif( kbit(ietide,5) ) then  ! IERS10
            ptd_model = 'IERS10' 
        elseif( kbit(ietide,3) ) then  ! IERS10
            ptd_model = 'ZERO'  
        else 
            ptd_model = 'NONE'
        endif
        call ptide(ipole,ptd_model,jd,t,geoc,dptotal)
        do j=1,3
          dptotal(j)=dptotal(j)/1000000.0d0
        enddo       
 
c ------Debug for the pole tide component
        if (debug) then
          write(iscrn,'(a,2x,3f8.1)') 'E-fixed pole tide dXYZ :'
     .             ,(dptotal(i)*1.d6,i=1,3)
c         convert to NEU    
c         kf routine rotate_geod expects units of meters: convert from km
          do i = 1,3
            dptotal(i) = dptotal(i)*1000.0d0
          enddo
          call rotate_geod (dptotal, ptideneu,'XYZ','NEU',evecf
     .                     ,geocn,rotmat)
          do i = 1,3
            dptotal(i) = dptotal(i)/1000.0d0
            ptideneu(i) = ptideneu(i)/1000.d0
          enddo
          write(iscrn,'(a,10x,3f8.1)') 'Pole tide dNEU :'
     .               , (ptideneu(i)*1.d6,i=1,3)
        endif     
c-------End debug    

      endif


c       Compute the ocean-loading component 

c         Tidal amplitudes and phases for 11 frequencies read from table in SETUP
c         Fomulation follows IERS Tech. Note 21 (July 1996)
               
      if (kbit(ietide,4)) then
                  
c       Time for ocean-tide arguments should be UT1, but use UTC instead
c       (max error 1 sec = 10 ppm)
        utc = t - taiutc(jd) + 19.d0   
                                         
c       if # tidal components = 11, assume OSO 
c        if( otide_source.eq.'OSO' ) then
        if( notl.eq.11 ) then 

c         Scherneck routine expects true Julian date (= PEP_JD -0.5)
          xjd = jd - 0.5d0
          call ocearg ( utc,xjd,oangle ) 
                             
c       if# tidal components = 54, assume NAO
c        elseif (otide_source.eq.'NAO' ) then
        elseif (notl.eq.54 ) then

c         Matsumoto routine expects Modified Julian date 
c             (= JD - 2400000 = PEP_JD - 2400001 )  
          xmjd = jd - 2400001  + utc/86400.d0
          call ocearg2 ( xmjd, oangle )  
        
        endif
              
c        print *,'ETIDE xjd utc xmjd notl',xjd,utc,xmjd,notl 
c        print *,'oangle ',oangle
        
c       compute NEU from all components 
        do i=1,3
          otideneu(i) = 0.0d0
        enddo       
c       the standard table gives up, west, south, so reverse order and signs to get N, E, U 
c       units of table are meters and degrees, convert to km and rad for use, mm for debug 
        do i=1,notl         
c       north (1) is negative of south (3)
          otideneu(1) = otideneu(1) - otides(i,3) *
     .        ( dcos( oangle(i) - otides(i,6)*pi/180.d0 ) )
c       east (2) is negative of west (2) 
          otideneu(2) = otideneu(2) - otides(i,2) * 
     .       ( dcos( oangle(i) - otides(i,5)*pi/180.d0 ) )
c       up (3) is positive of up (1)
          otideneu(3) = otideneu(3) + otides(i,1) * 
     .       ( dcos( oangle(i) - otides(i,4)*pi/180.d0 ) ) 
        enddo

c       kf routine rotate_geod expects units of meters: convert from km
        call rotate_geod (otideneu, dototal,'NEU','XYZ',evecf
     .                   ,geocn,rotmat) 
        do i = 1,3
          dototal(i) = dototal(i)/1000.0d0
        enddo    
c       correct from CE to CM frame if necessary
        if( otidemod(8:8).eq.'E') then 
          call otlcmc( jd,t,otidemod,notl,2,docmc )
          if( debug ) then  
            write(iscrn,'(a,9x,3f8.1)') 'Otide CE dXYZ :  '
     .          , (dototal(i)*1.d6,i=1,3)
            write(iscrn,'(a,9x,3f8.1)') 'Ocean tide CMC:  '
     .          , (docmc(i)*1.d3,i=1,3)
          endif      
c         this corrects from CE to CM, but the sign that works (+)  seems
c         to be the opposite of what Scherneck says on his web page. 
c         Gerd Gendt (GFZ) confirms this.
          do i=1,3
            dototal(i) = dototal(i) + docmc(i)/1.d3
          enddo
        endif


           
c----------Debug for ocean-tide component
        if (debug) then
          write(iscrn,'(a,3f8.1)') 'E-fixed ocean tide dXYZ : '
     .             ,(dototal(i)*1.d6,i=1,3)
c         convert to NEU  
c         kf routine rotate_geod expects units of meters: convert from km
          do i = 1,3
            dototal(i) = dototal(i)*1000.0d0
          enddo   
          call rotate_geod (dototal, otideneu,'XYZ','NEU',evecf
     .                     ,geocn,rotmat)
          do i = 1,3
            dototal(i) = dototal(i)/1000.0d0
            otideneu(i) = otideneu(i)/1000.d0
          enddo
          write(iscrn,'(a,9x,3f8.1)') 'Ocean tide dNEU: '
     .             , (otideneu(i)*1.d6,i=1,3) 
        endif     
c---------End debug

      endif 

c ------------- ATM TIDES -------------------
c
c PT050628: bit 6 (not bit 5) is used for atm tides. Bit 5 in globk refers to the pole tide
      if( kbit(ietide,6) .and. latl ) then

c  PT050509: the information passed through to here concerning the atm tidal displacements
c            is the cosine and sine amplitudes of the deformation at time 1200 UT on
c            day 001 for the S1 and S2 tides. These can be converted into deformation
c            at any time according to:
c
c   displacement = C*cos ( freq*t + phas ) + S*cos ( freq*t + phas )
c
c   where  C = cosine amplitude
c          S = sine amplitude
c          freq = frequency (in radians/sec) of the tide
c                 freq S1 = 7.27220521664304e-05 rad/sec
c                 freq S2 = 0.000145444104332861 rad/sec
c          t = time in seconds
c
        freq(1) = 7.27220521664304d-05
        freq(2) = 0.000145444104332861d0

c  compute time argument
         dtatm = t - taiutc(jd) + 19.d0

c  now compute and sum the NEU tidal loadings for S1 and S2 at this epoch.
c  The values are in mm at this stage, so we divide here by 1000 to get metres
        do i=1,3
          atmtidneu(i) = 0.0d0
        enddo

        do i = 1,3
          do j=1,2         
            atmtidneu(i) =   atmtidneu(i)
     .                   + atides(j,i*2-2+1)*dcos(freq(j)*dtatm)
     .                   + atides(j,i*2-2+2)*dsin(freq(j)*dtatm)
          enddo
        enddo                        
        if( debug ) then
         print *,'S1 atides  ',(atides(1,i),i=1,6)   
         print *,'S1 freq t phi ',freq(1),dtatm,freq(1)*dtatm 
         print *,'S2 atides ',(atides(2,i),i=1,6)    
         print *,'S2 freq t phi ',freq(2),dtatm,freq(2)*dtatm 
         write(iscrn,'(a,8x,3f8.3)') 'Atmos tide dNEU : '
     .             , (atmtidneu(i),i=1,3)   
        endif
c       kf routine rotate_geod expects units of meters for site coordinate
c       convert from km
        do i = 1,3
          datotal(i) = 0.0d0
        enddo
        call rotate_geod (atmtidneu, datotal,'NEU','XYZ',evecf
     .                   ,geocn,rotmat)
        do i = 1,3
          datotal(i) = datotal(i)/1.d6
        enddo
      endif

          
c        Sum up the tidal terms
       
c     sum the E-fixed corrections, if any
      if (kbit(ietide,2).or.kbit(ietide,3).or.kbit(ietide,4)) then
        do i=1,3
          dttotal(i)=dftotal(i)+dptotal(i)+dototal(i)+datotal(i)
        enddo  
c       rotate to inertial
        idir = 1
        call rotsnp( idir,jd,t,tdtgpst,isptide
     .             , frame,precmod,iut1,ipole,inut
     .             , xrot,xdtrot,sidtm1,xpole,ypole ) 
        call matmpy (xrot,dttotal,dttotali,3,3,1)   
      endif
c     add these to the inertial (freq-indep) component
      do i=1,3
         tidevec(i)=dttotali(i)+dstotal(i)  
      enddo


c------Debug for totals

      if( debug ) then
       write(iscrn,'(a,3f11.7)') 'Sum of E-fixed terms dXYZ :  '
     .      ,(dttotal(i),i=1,3)
       write(iscrn,'(a,3f11.7)') 'E-fixed rotated to inertial: '
     .      ,(dttotali(i),i=1,3)
       write(iscrn,'(a,14x,3f11.7)')'Total inertial:'
     .      ,(tidevec(i),i=1,3)
       do i=1,3
         tidevecf(i) = etidefixed(i) + dttotal(i) 
       enddo                             
       write(iscrn,'(a,15x,3f11.7)')'Total E-fixed:'
     .      ,(dttotal(i),i=1,3)
       do i=1,3
         totalneu(i) = etideneu(i)+ftideneu(i)+ptideneu(i)+otideneu(i)
     .                +atmtidneu(i)
       enddo  
       write(iscrn,'(a,18x,3f11.7)') 'Total dNEU:',(totalneu(i),i=1,3)
       
c      if debug flag is on, create a plot file of terms in NEU (mm)
c      write only 1 set per epoch
c      CW 050118 added extra output for plot file
       if ((dnint(t)-told).gt.0.1d0.or.(dnint(t)-told).lt.-0.1d0) then

c       Write out station co-ords
cd        write(79,'(a,a,3f16.8)') 'Stn co-ords in spherical units: '
cd     .       ,sitecd,evecfsph(1)/convd,evecfsph(2)/convd, evecfsph(3)
     
cd       Write site name, tide model and time variables
cd        write(79,'(a,a4,i8,f20.8,i8,i8,i8,f16.8,f16.4,f16.4,1x,
cd     .  a8,1x,15(f6.1,1x))') 
cd     .  'Site, JD, xJD, Yr, Month, Day, FHR, utctdt, t, etide_model: '
cd     .  ,sitecd,jd,xjd,date(1),date(2),date(3),FHR,utctdt,t,etidemod
cd     .  ,(1.d6*etideneu(i),i=1,3),(1.d6*ftideneu(i),i=1,3)
cd     .  ,(1.d6*ptideneu(i),i=1,3),(1.d6*otideneu(i),i=1,3)
cd     .  ,(1.d6*totalneu(i),i=1,3)
         told = dnint(t)

       endif
      endif
c-------end debug

                                       
c--------------Non-tidal atmospheric loading ------------------------

c  convert NEU displacements (in mm) into earth-fixed XYZ displacements (in km)
      do i = 1,3
        atmneu(i) = dble(atmlod(i))
      enddo   
      call rotate_geod (atmneu, atmtotal,'NEU','XYZ',evecf
     .                   ,geocn,rotmat)
      do i = 1,3
        atmtotal(i) = atmtotal(i)/1.d6
      enddo
           
c       rotate to inertial
      idir = 1
      call rotsnp( idir,jd,t,tdtgpst,isptide
     .           , frame,precmod,iut1,ipole,inut
     .           , xrot,xdtrot,sidtm1,xpole,ypole ) 
      call matmpy (xrot,atmtotal,atmvec,3,3,1)   

      if(debug)then      
        write(iscrn,'(a,3f8.1)') 'Non-tidal atm load  dNEU :'
     .             , (atmneu(i),i=1,3)   
        write(iscrn,'(a,5x,3f11.7)')'Atm load E-fixed dXYZ : '
     .              ,(atmtotal(i),i=1,3)
        write(iscrn,'(a,4x,3f11.7)')'Atm load inertial dXYZ : '
     .              ,(atmvec(i),i=1,3)
      endif
c        stop ' stopped at end of etide'

      return
      end


