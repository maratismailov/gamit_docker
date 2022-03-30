      Subroutine kouba_gps( ttag,antbody,aybias,yrate
     .                    , xsv,vsvc,santxyz,betadg,svbcos
     .                    , week,sow,isat,iprn,isvn
     .                    , yaw_angle,ievent,ybias                        
c      this for debug
     .                    , xsun,iepoch )

c     GAMIT version of the GPS portion of Jan Kouba's eclips_May2017 routine.  
c     The original routine has in its input arguments a flag indicating the 
c     direction of the processing (forward/backward) [ idir ], a flag for 
c     eclipsing/noon-turning at this epoch [ ieclips ], the eclipse or noon-turn
c     counter for the SV [ neclips ], the start/stop times of the events
c     [ eclstm(neclips), ecletm(neclips), and the beta-angle limits for an 
c     eclipse or noon-term to be occuring. However, all of these variables 
c     are initialized, computed, or not used by the routine, so are omitted 
c     from the calling arguments for this version.  The required inputs are 
c     then just the time-tag of the observations [ ttag ], the index in the
c     y-file array [ isat ], the SV body-type [ antbody], the yaw-bias type 
c     [ aybias ] and maximum rate [yrate] as read from svnav.dat, the SV 
c     position and velocity [ xsv, vsvc ], and angles between the  Sun and SV 
c     orbit plane [ beta ] and SV vector [ svbcos ].  The return arguments  
c     are the yaw angle, an event flag [ ievent ].  Also passed for reporting 
c     purposes only are the GPS time [ week, sow ] and the PRN and SVN numbers.
c     Finally we have reformattted the code for easier reading and removed 
c     the body-axis rotations, which are done in model/svbody_coords.f. 

c     R. King January 2014 / May 2017
              
c     Routine eclips_dec2013 copyright NRCAN Geodetic Survey Division 2011,2013   
c     Contact kouba@rogers.com 
                       
c     Update History
c       Kouba's eclips.f, eclips_dec2013, eclips_jan2014, and eclips_feb2017
C       Aug. 23, 2011: Change maximum yaw rates for BLK II and BLK IIA to averages 
c                      of JPL reprocessing 1996-2008.
c       Sep 26, 2011: 1. Corrected bug causing Block IIF shadow-crossing to yaw
c                        at 0.06 deg/s even for debeadg > 8 deg
c                      2. Corected/improved IIA recovery logic
c        Dec 18, 2013  1. Corrected IIF night turns (USAF Doc.)
C                      2. Small neg beta IIF and small pos IIA noon turn
c                         (wrong) directions  for |beta| < 0.9deg
c                      3. PRN/SVN 23 IIA YBIAS= -0.5 deg up to Feb 2004 
c        Jan 24, 2014: Noon recovery corrected if IIA/IIF's have ybias=0.
c                      (not applicable currently, possible for future IIF?)
c                      and small IIA/IIF beta if statement simplied.        
c        Jan 16, 2015: the IIA body-x orientation (+x=>Sun) and logic is implemented 
c                      also for IIR/IIRM eclipsing (no SANTXYZ reversal )  
c        Mar 9, 2016-Jan 20, 2017: Add Beidou and Galileo (not in this routine). 
c        Feb 27, 2017: 1. Added code to detect a possible beta sign change during
c                         noon or shadow turns. 
c                      2. Set the IIF yas bias to -0.7 deg as recommended by Kuang
c                         et al. GPS Solut., doi:10.1007,s10291-016-0562-9 (2016).
c                      3. Correct IIA sign of yaw during noon tunrs when beta 
c                         between 0 and +0.5 (due to the yaw bias of 0.5 deg).

c         GAMIT kouba_yaw:
c          Jan 24, 2014: Initial routine checked against tables from Kouba for
c                        IIA PRN06/SV36 eclipse of 26 January 2013, IIR PRN12/SV58
c                        noon turn and night turn of 17 Jun 2011, and IIF PRN25/SV62
c                        noon turn and night turn of 17 Jun 2011. 
c         GAMIT kouba_gps: Removed all but GPS code and updated with additions noted
c                       above for changes between eclips_jan2014 and eclips_feb2017.
c                       Changed argument passed from yawtab from the yaw bias value 
c                       (ybias) to the sense of the yaw (aybias), now setting the
c                       values inside this routine. 
                            

c     Input:    
c       Required for computations:                     
c       ttag       Time argument, seconds from t/y-file start 
c       antbody    Body type (works for all GNSS)
c       aybias     svnav.dat flag for sense of the bias--values set in this routine   
c       yrate      Nominal/maximum yaw rate (deg/s) for the SV (from svnav.dat)          
c       xsv(3)     Earth-fixed SV position and velocity              
c       vsvc(3)    SV inertial velocity vector 
c       santxyz(3) Body X-axis unit vector (computed by gamit/lib/attit_yaw.f)  
c                 (Warning: the IIA body-x orientation reversed for the IIR)   
c       betadg     Angle between Sun and SV plane (deg)
c       svbcos     Cosine of angle between the SV vector and the Sun vector
c       Used only for reporting purposes
c       week, sow  GPS week, seconds-of-week
c       isat,iprn,isvn  SV array index, SV PRN and SVN 

c     Output:  
c       Required for y-file                                         
c       yaw_angle  Yaw angle, nominal or event (deg)
c       ievent  Saved values (in common/yaw/ in calling routines) of flags indicating 
c                the current status of SV events. Used internally to report start/end of 
c                events and passed out to be written on the y-file for use by MODEL.
c                 0 = no event, nominal yaw 
c                 1 = night turn or eclipse 
c                -1 = recovery phase of eclipse
c                 2 = noon turn              
c                 3 = orbit-normal mode (not used for GPS)
c       Used only for reporting in yawtab.out.DDD 
c       ybias     Yaw bias (deg) 
         
c Internal
c       ybias    IIA yaw bias=0.5 deg since Nov 1995 except PRN/SVN23?
c       eclstm  Start time for the current eciipses or noon turn (units as ttag) 
c       ecletm  End times for eciipses and noon turns for the SV (units as ttag_ 
c       pi         3.1415926...
C       anoon    SV beta angle limit (deg) for a noon turn maneuver
c       anight   SV beta angle limit (deg) for a night shadow crossing
c       svbcos   Cosine of angle between the sun and SV vectors (i.e. 
c                 the conventional 'beta', not the input beta      
c                (Warning: the IIA body-x orientation reversed for the IIR)  
c       yangle   Nominal yaw angle (deg)
c       phi      Actual yaw angle if during an event (deg)   
c       svbcos_last(maxsat),betaini(maxsat) : saved variables for each SV 
c         to precent discontinuities in events and yawing direction

c Kouba remarks:
c   SVBCOS, the COS of the angle between the sv radius vector and the sun 
c   radius vector, (the dot product of the above respective unit vectors), 
c   SVBCOS is used to test against CNOON=COS(ANOON), if sv is entering 
c   a noon maneuverer (i.e., SVBCOS > CNOON). The ANOON limit 
c   (typically < 5.7 deg), is determined within the subroutine and depends 
c   on the  sv yaw rate (YRATE)and GNSS (GPS or GLONASS). 
c
c   ANIGHT , the shadow limit is the "input " in the subroutine 
c   call statement, it can be hard coded as it is constant for a GNSS type,
c   e.g., 
c   180.D0+-13.25D0 for GPS (IPRN.GE. 32) and 180.D0+-14.20D0 for 
c   GLONASS (IPRN.GT.32), resp.). CNIGHT=COS(ANIGHT) is used for testing.
c   When SVBCOS < CNIGHT (CNIGHT is negative and close to -1), 
c   the SV enters, or it is in the shadow. 
c
c Notes: 1) Kouba uses GPS second-of-week for his time argument, but without 
c           some accommodation, this can produce a discontinuity. We keep the
c           GPS week, seconds-of-week for easy comparison with his test cases
c           but use the continuous orbital argument (satarg) as the primary 
c           argument for the yaw model    
c        2) Mapping of Kouba's 'iblk' variable into GAMIT antenna body type:
c             1  BLOCK I              4 BLOCK IIR 
c             2  BLOCK II             5 BLOCK IIR-M
c             3  BLOCK IIA            6 BLOCK IIF
c        3 Yaw biases 
c              Blocks I, II, and IIR = 0 
c              Block IIA     = +0.5
c              SV 23 (Block IIA) = -0.5 (decommisioned 2004/2/23)
c              Block IIF = -0.5 for small beta angles only
c              Block IIIA : use IIF values for now

      implicit none             

      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      
      integer*4 ievent,isat,isvn,iprn,ieclips,week,yr,doy,hr,min,i,j
      
      real*8 ttag,xsv(3),vsvc(3),yangle,santxyz(3),betadg,svbcos
     .     , eclstm,ecletm,anoon,anight,cnoon,cnight,dir,dttag
     .     , murate,det,yrate,phi,santx,santy,v(3),r(3)
     .     , yawend,twohr,halfhr,pi,dtr
     .     , yaw_angle,sec,utcoff,ybias,sow
     .     , svbcos_start(maxsat),svbcos_check,betaini(maxsat) 

                                                    
      character*1 aybias
      character*20 antbody
      character*256 message

      logical noon,night,blockII,blockIIR,blockIIF,blockIIIA
     .       ,nominal,bias_warning(maxsat)
                                                
      data bias_warning/maxsat*.false./ 
 
c** rwk temporary for debug:
      integer*4 iepoch,iprndb/2/
      logical debug /.false./,mu_code/.false./
      real*8 temp1, temp2, temp3,temp4 
      real*8 E,mu,yangle1,yangle2,yangle3,Sx,Sy,Sz,rad,cosbetamu 
      real*8 xsun(3),sinmu,cosmu,dot,sunmag,svmag,crossmu(3)
      real*8 sinalpha,cosalpha,svvmag,egmag,eamag,ermag,enmag 
      real*8 er(3),ev(3),en(3),eg(3),ea(3),es(3)

c These not passed, so save for checking
      data betaini/maxsat*0.d0/ 
      save svbcos_start,betaini

**DEBUG 
cd     if( iepoch.ge.737.and.iepoch.lt.761.and.iprn.eq.30 ) then
cd      if( iepoch.eq.1 ) then
cd      if( iprn.eq.12 .or .iprn.eq.25 )  then
cd      if( iprn.eq.30 ) then 
cd      if( iprn.eq.32  ) then  
* MOD TAH 200211: Commented out debug being turned on
*       if( iprn.eq.iprndb.and.iepoch.lt.2 ) then 
*          debug = .true.
*       else
*           debug = .false.
*       endif          


c Set the block flags - simplifies slightly the logical tests
c vs testing svantbody - 
* MOD TAH 201111: Fixed block II name, added block IIIA
      blockII  = .false.
      blockIIR = .false.
      blockIIF = .false.
      blockIIIA = .false.
      if( antbody(1:9).eq.'BLOCK I  '.or.antbody(1:9).eq.'BLOCK II '.or.
     .    antbody(1:9).eq.'BLOCK IIA' ) then
        blockII = .true.
      elseif( antbody(1:9).eq.'BLOCK IIR' ) then
        blockIIR = .true.
      elseif( antbody(1:9).eq.'BLOCK IIF' ) then
        blockIIF = .true.
      elseif( antbody(1:10).eq.'BLOCK IIIA' ) then
        blockIIIA = .true.
      else     
        if(debug) print *,'KOUBA_GPS  iepoch isat antbody '
     .        ,iepoch,isat,antbody
        call report_stat('FATAL','YAWTAB','orbits/kouba_gps'
     .               ,' ','SV body type not recognized',0)
      endif
                      
c Set the yaw bias from the svnav.dat code

      if( aybias.eq.'U') then
          ybias = 0.d0
      elseif( aybias.eq.'P') then
         ybias = 0.5d0
      elseif( aybias.eq.'N' ) then
         ybias = -0.5d0               
c        increase this from the nominal -0.5 for IIF to avoid anomalous noon turns 
***   rwk 190109: For now use Block IIF logic for Block IIIA
c         if( blockIIF ) ybias = -0.7d0
c         Kuang et al., GPS Solut.,doi:10.1007/s10291-016-0562-0, 2016
c         recommend increasing the bias for IIFs to avoid an erroneous 
c         turn for small beta
          if( blockIIF .or. blockIIIA ) ybias = -0.7d0
      elseif( aybias.eq.'Y' ) then
         ybias = dsign(ybias,betadg)
      elseif( aybias.eq.'A' ) then
         ybias = -dsign(ybias,betadg)
      else
        write(message,'(a,a1,a,i3)') 'Unknown ybias type '
     .            ,aybias,' for PRN ',iprn
        call report_stat('FATAL','YAWTAB','orbits/kouba_gps',' '
     .                       ,message,0)
      endif             
c     Kuang et al., GPS Solut.,doi:10.1007/s10291-016-0562-0, 2016
c     recommend increasing the bias for IIFs to avoid an erroneous turn for small beta
*** rwk 190109: For now use Block IIF logic for Block IIIA
* MOD TAH 201111: Removed lines below because above sets the y-boas angle.
*     Fixed indent was well
*     if( antbody(1:9).eq.'BLOCK IIF' ) ybias = -0.7d0 
*     if( antbody(1:10).eq.'BLOCK IIIA' ) ybias = -0.7d0

c     warn if not consistent with Kouba defaults
* MOD TAH 201111: Fixed indenting so code flow is more obvious.  Changed
*     tests to be consistent with logicals. (NOTE: Block I satellites
*     are grouped with block II in logical.
      if( blockII ) then
          if( ybias.ne.0.5d0.and. .not.bias_warning(isat)) then
              write(message,'(a,i2)') 
     .          'ybias not +0.5 for Block I,II, or IIA PRN ',iprn
              call report_stat('WARNING','YAWTAB','orbits/kouba_gps',' '
     .                         ,message,0) 
c             set this to avoid getting warnings at every epoch
              bias_warning(isat) = .true.
          endif 
      endif
      if( blockIIR ) then
        if( ybias.ne.0d0. and. .not.bias_warning(isat) ) then
          write(message,'(a,i2)') 
     .      'ybias not zero for Block IIR PRN ',iprn
          call report_stat('WARNING','YAWTAB','orbits/kouba_gps',' '
     .                     ,message,0)                     
          bias_warning(isat) = .true.
        endif 
      endif
      if( blockIIF .or. blockIIIA ) then    
*** rwk 190109: Temporarily use IIF logic for IIIA 
C    .    antbody(1:10).eq.'BLOCK IIIA' ) then
          if( ybias.ne.-0.7d0 .and.
     .         .not.bias_warning(isat) ) then
             write(message,'(a,i2)') 
     .         'ybias not -0.5 for Block IIF PRN ',iprn
             call report_stat('WARNING','YAWTAB','orbits/kouba_gps',' '
     .                     ,message,0) 
             bias_warning(isat) = .true.
          endif 
      endif
* MOD TAH 201111: Finish fixing indent.
       
c Set some constants
      twohr = 7200.D0
      halfhr= 1800.D0 
      pi=4.d0*datan(1.d0) 
      dtr=pi/180.D0    

c Initialize the eclipse information
      noon = .false.
      night = .false.
c     initialize start/stop time be way in the future to avoid tripping 
c     the eclipse check before an eclipse or noon turn is detected   
      eclstm = 999999.d0
      nominal = .false. 
c     assume no event from the start      
      nominal = .true.                         

C Compute the noon beta angle limit (beta zero) for a noon turn from the 
c maximum yaw rate for the SV and the actual orbit rate  (~0.00836 deg/s).
       murate= dsqrt((vsvc(1)**2+vsvc(2)**2+vsvc(3)**2)/
     & (xsv(1)**2+xsv(2)**2+xsv(3)**2))/dtr
      anoon=datan(murate/yrate)/dtr
c       rwk:  This routine computes critical beta angle (beta0) for the noon 
c       and midnight turns from the orbit rate and yaw-rate, as shown above, 
c       but does not explictly compute the beta angle for night shadowing.
c       The paper assumes that this is 13.25 deg, which I'll use here for
c       now, to be modified possibly later from the GAMIT shadow routine.
        anight = 180.d0 + 13.25d0 
      cnoon=dcos(anoon*dtr)
      cnight=dcos(anight*dtr)
c     IIRs use nominal yaw during shadowing 
c     (yaw=0) when beta < beta0 (same as noon for IIR)
      if( blockIIR ) then
         cnight = dcos((anoon+180.d0)*dtr)
         if( debug ) print *,'Blk IIR update cnight ',cnight
      endif
                            
c Set the night/noon event flags from the beta angle and the beta-angle limits
c     use the beta angle (cos) from the beginning of the event (ievent and
c     svbcos_start are saved variables)
      svbcos_check = svbcos
      if( ievent.ne.0 ) then 
         svbcos_check = svbcos_start(isat)
       else
         svbcos_check = svbcos                              
         svbcos_start(isat) = svbcos
c        this added from eclips_Feb2017 to guard against an artificial 
c        sign change for small beta
         if( dabs(betadg).le..07d0.and.betaini(isat).eq.0.d0 ) 
     .        betaini(isat)=betadg
      endif 
      if( svbcos_check .lt. cnight ) then
        night = .true.                  
        if( ievent.eq.0 ) svbcos_start(isat) = svbcos
      endif
      if( svbcos_check .gt. cnoon ) then
        noon = .true.   
        if( ievent.eq.0 ) svbcos_start(isat) = svbcos
      endif                       
  
      if( debug ) then        
        print *,' '     
        call debugtime(week,sow)   
        print *,' iepoch sow iprn antbody yrate betadg murate'
     .           ,iepoch,sow,iprn,antbody,yrate,betadg,murate
        print *,'  xsv vsvc santxyz ',xsv,vsvc,santxyz   
       if(debug) print *,
     .     'ievent svbcos svbcos_start svbcos_check night noon '
     .     ,ievent,svbcos,svbcos_start(isat),svbcos_check
     ,     ,night,noon
        print *,'anoon cnoon anight cnight yrate '
     .       ,anoon,cnoon,anight,cnight,yrate
      endif 
        
c Compute the nominal yaw angle from the orbit geometry
      if( debug ) print *,'santxyz vsvc ',santxyz,vsvc 
      yangle= dacos( (santxyz(1)*vsvc(1)+santxyz(2)*vsvc(2)+
     .               santxyz(3)*vsvc(3) ) / 
     .                 sqrt(vsvc(1)**2+vsvc(2)**2+vsvc(3)**2) )/dtr  
c     For all SVs, with IGS standard, yaw angle has opposite sign from beta 
cc** new code, no IIR sign reversal
      if( betadg.gt.0.d0 ) yangle = -yangle        
cc** end new code
c** old code with IIR sign reversal (keep for testing purposes)
cc      if( betadg.lt.0.d0 .and. blockIIR ) yangle = -yangle
cc      if( betadg.gt.0.d0 .and. .not.blockIIR ) yangle = -yangle
c** end old code 
      if( debug ) print *,' svbcos betadg yangle '
     .                   , svbcos,betadg,yangle  
c RWK 170520; This code to compute and use mu here to validate its use for kouba_galileo
c              Invoke it again only if needed

      if( mu_code.and.iprn.eq.iprndb ) then 
cd      if( iepoch.gt.1 ) stop 1 
cd      if( iprn.eq.6 ) then           
        print *,'Epoch Kouba beta yangle ',iepoch,betadg,yangle
        E = dacos(svbcos)/dtr                   
cd      cosbetamu = dcos(E*dtr)*dcos(betadg*dtr)
        mu = dacos(dcos(E*dtr)/dcos(betadg*dtr))/dtr      
        Sx = -dsin(mu*dtr)*dcos(betadg*dtr)
        Sy = -dsin(betadg*dtr)
        Sz = dcos(mu*dtr)*dcos(betadg*dtr)
        rad = dsqrt(1.d0- Sz*Sz)          
        yangle = datan2(Sy/rad,Sx/rad)/dtr 
        if( betadg.gt.0.d0 ) yangle = -yangle  
        print *,'E mu svbcos yangle ',svbcos,E,mu,yangle 
        call mu_angle( xsv,vsvc,xsun,mu,debug,iprndb )
c         do i=1,3
c           er(i) = xsv(i)/svmag 
c           ev(i) = vsvc(i)/svvmag
c           es(i) = xsun(i)/sunmag 
c         enddo       
c         ermag = dsqrt(dot(er,er))
c         print *,'er norm ',er
c         call cross(er,ev,en)
c         print *,'er ',er
c         print *,'ev ',ev 
c         print *,'en ',en
c         enmag = dsqrt(dot(en,en))  
c         do i=1,3
c           en(i) = en(i)/enmag
c         enddo 
c         print *,'en norm ',en 
c            --> orbit normal 
c         call cross(en,es,eg)                                 
c         print *,'es ',es
c         print *,'eg ',eg  
c         egmag = dsqrt(dot(eg,eg))
c         do i=1,3
c           eg(i) = eg(i)/egmag
c         enddo   
c         print *,'egmag norm ',eg
c            --> eg points to zero - 90 deg in the orbit plane
c         call cross(eg,er,ea)  
c         eamag = dsqrt(dot(ea,ea))
c         print *,'ea eamag ',ea,eamag
c            -->  angle between eg and the SV 
c         print *,'mag eg er ea ',egmag,ermag,eamag
c         sinalpha = eamag/(egmag*ermag)
c         cosalpha = dot(eg,er)/(egmag*ermag)
c         mu = datan2(sinalpha,cosalpha)/dtr - 90.d0  
        Sx = dsin(mu*dtr)*dcos(betadg*dtr)
        Sy = -dsin(betadg*dtr)
        Sz = dcos(mu*dtr)*dcos(betadg*dtr)
        rad = dsqrt(1.d0- Sz*Sz)          
        yangle = datan2(Sy/rad,Sx/rad)/dtr 
c        no sign reversal for the vector/mu method
c        if( betadg.gt.0.d0 ) yangle = -yangle
        print *,'N Sx Sy Sz rad  ',Sx,Sy,Sz,rad 
        print *,'N sinal cosal mu yangle ',sinalpha,cosalpha,mu,yangle
cd      stop           
      endif 
     
c If in a night or noon maneuver, get the start/stop times and the actual yaw

      if( night.or.noon ) then

c       indicate a potential event, reset later if non-nominal indicated
        nominal = .false.
     
c       get the orbital angle for the start (Kouba 2009 Eqs 8-9)
        det = dsqrt( (180.d0 - dacos(svbcos)/dtr)**2 - betadg**2 )
        if(debug) print *,'first det ',det 
c       initialize the yaw angle at 90 degrees (rwk: not sure why)
        phi = pi/2
C       check if already after a midnight or noon
        if( night ) then 
          if( dabs(yangle).lt.90.d0 ) det = -det
          if( det.ne.0.d0) 
     .       phi = datan2( -dtan(betadg*dtr),dsin(-det*dtr) )/dtr
           if( debug ) 
     .         print *,' NIGHT det phi ',det,phi  
        endif
      
        if( noon ) then 
          det = sqrt((dacos(svbcos)*180.d0/pi)**2-betadg**2)      
          if( dabs(yangle).gt.90.d0) det =-det
          if( det.ne.0.d0 ) 
     .        phi = datan2(dtan(betadg*dtr),-dsin(pi-det*det))/dtr
             if( debug ) print *,' NOON det phi ',det,phi 
        endif
                              
        eclstm = ttag + det/murate   
        yawend = datan(murate/yrate)/dtr  
        if(debug) print *,' ttag sow det murate eclstm yawend '
     .                    , ttag,sow,det,murate,eclstm,yawend

c       For IIR night turn only makes sense when beta < beta-zero
c       For IIA it gets here only when noon is true and that happens only 
c         when beta < beta-zero (rwk: I don't see that in the logic)
        if( (blockIIR.or.noon) .and. dabs(betadg).lt.yawend) then 
c         IIA/IIR/IIF noon or IIR midnight turn
          eclstm= eclstm -
     .       dabs(betadg)*sqrt(anoon/dabs(betadg)-1.d0)/murate
          ecletm = eclstm +
     .       2.d0*dabs(betadg)*sqrt(anoon/dabs(betadg)-1.d0)/murate
          if(debug) print *,'NOON OR NIGHT  eclstm ecletm'
     .                      , eclstm,ecletm
        endif   
c       II/IIA/IIFshadow start/end times (IIR keep nominal yaw during shadow)
        if( .not.blockIIR .and.night ) then 
          if( debug ) print *,'eclstm(mid) ',eclstm
          eclstm= eclstm -
     .           dsqrt((anight-180.d0)**2-betadg**2)/murate
          ecletm= eclstm +
     .           2.d0*dsqrt((anight-180d0)**2-betadg**2)/murate
          if( debug ) then
            temp1= dsqrt((anight-180.d0)**2-betadg*2)/murate
            temp2=2.d0*dsqrt((anight-180d0)**2-betadg**2)/murate 
            print *,'temp1 temp2 ',temp1,temp2
          endif
          if( debug) print *,'II SHADOW  eclstm ecletm ',eclstm,ecletm
        endif
        if( debug ) 
     .     print *,'aft setting times flag ttag sow eclstm ecletm '
     .              ,ttag,sow,eclstm,ecletm

c Times now set, see if within the event times
        if( ttag.le.eclstm .or. ttag.ge.(ecletm+halfhr) ) then
c         if not within the event, set nominal=true and exit
          nominal = .true.
          go to 1 
        else
c         otherwise calculate the modeled yaw
          if(debug) print *,'COMPUTING modeled yaw'
c         for small beta, use initial value throughout the eclipse 
c         to unsure no beta sign-change
          if(dabs(betadg).le..07d0.and.betaini(isat).ne.0.d0) 
     .        betadg = betaini(isat)
c         S/C position and velocity unit vectors
          do j=1,3
            v(j)=vsvc(j)/dsqrt(vsvc(1)**2+vsvc(2)**2+vsvc(3)**2)
            r(j)=xsv(j)/dsqrt(xsv(1)**2+xsv(2)**3+xsv(3)**2)
          enddo
c         orbit angle at start of event
          det = murate*(ecletm-eclstm)/2.d0 
c         IIF shadow yaw rate  (rwk: ugly to set this w/o a block check)
          yawend = ( datan2(-dtan(betadg*dtr),dsin( det*dtr)) - 
     .             datan2(-dtan(betadg*dtr),dsin(-det*dtr)) )/dtr/
     .                (ecletm-eclstm)
          if( debug ) print *,'IIF shadow yaw rate yawend ',yawend

          if( svbcos.lt.0.d0 ) then
c           SHADOW CROSSING OR NIGHT TURN    
            nominal = .false.
            if( debug ) print *,'NIGHT antbody ttag sow ecletm ievent '
     .                         ,antbody,ttag,sow,ecletm,ievent
            if( blockII.or.blockIIF.or.
*** This temporary until we know how IIIA behaves 
     .          blockIIIA ) then
              if( ttag.le.ecletm ) then
                if( ievent.eq.0 ) then
                  ievent = 1 
                  call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
                  write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .             antbody,' PRN',iprn,' eclipsing    ',doy,hr,min
     .              ,'  beta=',betadg
                  call report_stat('STATUS','YAWTAB','orbits/kouba_gps'
     .                   ,' ',message,0)   
                endif
                if( blockII ) then    
c                 IIA night turn   (rwk: why not a simply a shadow?) 
                  phi=datan2(-dtan(betadg*dtr),dsin(-det*dtr))/dtr
     .                  +sign(yrate,0.5d0)*(ttag-eclstm)
                  if( debug ) print *,'IIA night turn phi ',phi    
                endif
*** Temporarily use IIF logic for IIIA
c                if( blockIIF ) then  
                if ( blockIIF.or.blockIIIA ) then
c                 IIF night turn 
                  phi=datan2(-dtan(betadg*dtr),dsin(-det*dtr))/dtr
c                   correct IIF night crossing using computed yaw 
c                   rate (yawend) according to the USAF IIF document
     .                  + yawend*(ttag-eclstm)
                  if(debug) then
                    temp1 = dtan(betadg*dtr)
                    temp2 = dsin(-det*dtr)
                    temp3 = datan2(-temp1,temp2)/dtr
                    temp4 = yawend*(ttag-eclstm)
                    print *,'temp1,temp2,temp3,temp4 '
     .                    ,temp1,temp2,temp3,temp4
                  endif
                  if(debug) print *
     .           ,'IIF night det betadg ttag sow eclstm yawend phi'
     .                      ,det,betadg,ttag,sow,eclstm,yawend,phi   
                endif
              else
c               exiting from shadow (warning: not recommened to use for IIA)
                if( ievent.ne.-1 ) then
                  ievent = -1 
                  call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
                  write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .              antbody,' PRN',iprn,' post eclipse ',doy,hr,min
     .             ,'  beta=',betadg
                  call report_stat('STATUS','YAWTAB','orbits/kouba_gps'
     .                   ,' ',message,0)   
                endif
                if( blockII ) then
                   phi=datan2(-tan(betadg*dtr),dsin(-det*dtr))/dtr
     .                 + sign(yrate,ybias)*(ecletm-eclstm)
                   if(debug) print *,'IIA shadow exit phi ',phi 
                endif    
*** rwk 190109: Temporarily use IIF logic for IIIA
*               if( blockIIF ) then               
                if( blockIIF.or.blockIIIA ) then
c                 IIF has immediate recover, resume nominal  
                  nominal = .true.
                  go to 1          
                endif
c               compute the actual yaw difference at shadow exit
                yawend = yangle -phi 
                yawend = dmod(yawend,360.d0)
                if(dabs(yawend).gt.180.d0) 
     .             yawend = yawend - 360.d0*yawend/dabs(yawend)
                phi = phi + dsign(yrate,yawend)*(ttag-ecletm)
c               see if the actual yaw has reached the nominal yaw
                santx = yangle - phi
                santx = dmod(santx,360.d0)  
                if( debug ) print *,'yaw diff at shadow exit '
     .                ,'yangle phi santx ',yangle,phi,santx
                if(dabs(santx).gt.180d0) santx = 
     .                santx -360d0*santx/dabs(santx)
                if(debug) print *,'nominal yaw check phi santx yawend '
     .                 , phi,santx,yawend
                if( dabs(santx).gt.dabs(yawend) .or. 
     .              (yawend.ne.0.d0.and.santx/yawend.lt.0.d0) ) then
c                 the nominal yaw angle is reached, exit
                  nominal = .true.
                  go to 1 
                else
c                 still yawing, set between +- 180 deg
                  phi = dmod(phi,360.d0)
                  if(dabs(phi).gt.180.d0) phi= phi-360.d0*phi/dabs(phi)
                  if(debug) print *,'still yawing, phi +-180 ',phi
                endif 
c             endif on t < t end             
              endif
c           endif on block II or IIF 
            endif
            if( blockIIR ) then
c             IIR shadow (midnight turn) crossing
              phi= datan2(dtan(betadg*dtr),-dsin(-det*dtr))/dtr
     .             + dsign(yrate,betadg)*(ttag-eclstm)
              if(debug) print *,'IIR midnight turn phi ',phi  
              if( ievent.eq.0 ) then
                ievent = 1 
                call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
                write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .             antbody,' PRN',iprn,' night turn   ',doy,hr,min
     .            ,'  beta=',betadg
                call report_stat('STATUS','YAWTAB','orbits/kouba_gps'
     .                   ,' ',message,0)   
              endif  
              if( phi/yangle.ge.1.d0.or.phi/yangle.lt.0.d0 ) then
                nominal = .true.  
                if(debug) print *,'phi/yangle >1 or <0 use nominal yaw'
                goto 1
              endif
            endif     
c
          else
c           NOON TURN 
            nominal = .false.
            if( ievent.eq.0 ) then
              ievent = 2 
              call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
              write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .            antbody,' PRN',iprn,' noon turn    ',doy,hr,min
     .            ,'  beta=',betadg
              call report_stat('STATUS','YAWTAB','orbits/kouba_gps'
     .                   ,' ',message,0)   
              if(debug) print *,'start ievent ',ievent
            endif   
c           phi for standard noon turn for II/IIA or IIF; overwritten below 
c           if a II/IIA/IIF small-angle problem or exit is detected 
            phi= datan2(-dtan(betadg*dtr),dsin(pi-det*dtr))/dtr
     .                 -dsign(yrate,betadg)*(ttag-eclstm)
            if(debug) then 
              print *,'standard noon turn phi ',phi          
              print *,'betadg ybias ',betadg,ybias 
              print *,'betadg*dsign(1,ybias)' ,betadg*dsign(1.d0,ybias)
              print *,'betadg*ybias ',betadg*ybias 
            endif                
** rwk 190109: Temporarily use IIF logic for IIIA 
**          if( (blockII.or.blockIIF .and.       
            if( (blockII.or.blockIIF.or.blockIIIA).and.
     .         (betadg*dsign(1.d0,ybias)).le.0.9d0 .and.
     .         (betadg*ybias).gt.0.d0 ) then  
c             small negative IIF beta or small posit. IIA noon turn problem
c             use ybias instead of beta 
c             Kouba comment: in theory the above limit of 0.9 deg should abs(ybias)
              phi= datan2(-dtan(betadg*dtr),dsin(pi-det*dtr))/dtr
     .                + dsign(yrate,ybias)*(ttag-eclstm)
              if(debug)print *,'small beta betadg phi',betadg,phi
            endif
            if( blockIIR ) then 
c             IIR uses the opposite sign
              phi= datan2(dtan(betadg*dtr),-dsin(pi-det*dtr))/dtr
     .            - dsign(yrate,betadg)*(ttag-eclstm)
              if(debug) print *,'IIR noon turn yangle phi '
     .                           ,yangle,phi
c             see if a IIR exit is detected 
              if(yangle/phi.ge.1.d0 .or. phi/yangle.lt.0.d0) then
                 nominal = .true.
                 go to 1
              endif
*** rwk 190109: Temporarily use IIF logic for IIIA 
**           elseif ( blockII .or. blockIIR. or. blockIIF) then
            elseif (blockII.or.blockIIR.or.blockIIF.or.blockIIIA) then

c               see if a IIA, IIR, or IIF exit is detected
c               (the IIR check appears redundant)   
              if( debug ) then
                print *,'check small-beta noon exit nominal ',nominal
                temp1 = betadg*dsign(1.d0,ybias)
                temp2 = betadg*ybias 
                temp3 = (phi-dsign(1.d0,ybias)*360.d0)/yangle
                temp4 = phi/yangle
                print *, 'betadg ybias phi yangle '
     .             ,betadg,ybias,phi,yangle
                 print *,'temp_1-4 ',temp1,temp2,temp3,temp4 
              endif
c** rwk : The following consolidated test seems correct but produces a false result
c**       I suspect a compiler problem, so break it up
c             if( ( betadg*dsign(1.d0,ybias).le.0.9d0 .and.       
c    .              betadg*ybias.gt.0.d0 .and.
c    .              ((phi-dsign(1.d0,ybias)*360.d0)/yangle).le.1.d0.or.
c    .              ((phi-dsign(1.d0,ybias)*360.d0)/yangle).lt.0.d0 ) 
c    .                   .or.
c    .            ( (betadg*dsign(1.d0,ybias).gt.0.9d0 .or.
c    .                 betadg*ybias.le.0.d0) .and.
c    .            (((phi/yangle).ge.1.d0.or.(phi/yangle).lt.0.d0)) ) )
c**  .                   then
                
              if(  betadg*dsign(1.d0,ybias).le.0.9d0 .and.       
     .             betadg*ybias.gt.0.d0 .and.
     .             (((phi-dsign(1.d0,ybias)*360.d0)/yangle).le.1.d0.or.
     .              ((phi-dsign(1.d0,ybias)*360.d0)/yangle).lt.0.d0 ) )
     .              then
                 if( debug ) print *,'1st test T goto 1 '
                 nominal = .true.
                 goto 1  
              endif
              if(  (betadg*dsign(1.d0,ybias).gt.0.9d0 .or.
     .                 betadg*ybias.le.0.d0) .and.
     .            ((phi/yangle).ge.1.d0.or.(phi/yangle).lt.0.d0) )
     .                   then
                if( debug ) print *,'2nd test T goto 1 '
                nominal = .true.                                     
                go to 1  
              endif     
c           endif on block type for noon-turn event
            endif               
c         endif on night or noon for event 
          endif
c       endif on whether within event times
        endif
c     endif on event indicated from beta angle (night/noon)
      endif 

c     all jumps from the logic come to here
1     continue 
      if( nominal ) then 
        yaw_angle = yangle
        if( ievent.ne.0 ) then
          if(debug) print *,'end ievent ',ievent
          call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
          if( ievent.eq.1 ) then  
            write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .         antbody,' PRN',iprn,' night exit   ',doy,hr,min
     .           ,'  beta=',betadg                               
          elseif( ievent.eq.-1 ) then
            call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
            write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .         antbody,' PRN',iprn,' post exit    ',doy,hr,min
     .          ,'  beta=',betadg
          elseif( ievent.eq.2 ) then 
             call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
             write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .          antbody,' PRN',iprn,' noon exit    ',doy,hr,min
     .              ,'  beta=',betadg 
          endif  
          call report_stat('STATUS','YAWTAB','orbits/kouba_gps',' '
     .                         ,message,0)   
          ievent = 0                            
          if( debug ) print *,'reset ievent ',ievent
        endif
      else
        yaw_angle = phi 
      endif
      if( debug ) print *,'prn antbody betadg nominal xhat'
     .                    ,iprn,antbody,betadg,nominal,santxyz
      if(debug) print *,'yangle phi yaw_angle ',yangle,phi,yaw_angle
cd      if( debug.and.iepoch.gt.315 ) stop
      return
      end

c********

      Subroutine debugtime( week,sow )
      implicit none  
      integer*4 yr,doy,hr,min,week,itflag
      real*8 t,sow,sec,utcoff
      itflag = 4
      call timcon(itflag,week,sow,yr,doy,hr,min,sec,utcoff)
cd      print *,'TIME sow doy hr min ',sow,doy,hr,min
      return      
      end

      







