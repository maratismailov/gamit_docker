      Subroutine kouba_yaw( week,sow,ttag,isat,iprn,isvn,antbody,yrate
     .                    , ybias,xsv,santxyz,vsvc,beta,svbcos
     .                    , yaw_angle,ievent,svbcos_start
c      this for debug
     .                    , iepoch )

c     GAMIT version of Jan Kouba's ECLIPS_JAN2014 routine.  The original
c     routine has in its input arguments a flag indicating the direction
c     of the processing (forward/backward) [ idir ], a flag or eclipsing/
c     noon-turning at this epoch [ ieclips ], the eclipse or noon-turn)
c     counter for the SV [ neclips ], the start/stop times of the events
c     [ eclstm(neclips), ecletm(neclips), and the beta-angle limits
c     for an eclipse or noon-term to be occuring. However, all of these
c     variables appear to be initialized, computed, or not used by the
c     routine, so are omitted in this version.  The required inputs are 
c     then just the time-tag of the observations [ ttag ], the SV identifiers 
c     [ svn or prn and block ], the SV position and velocity [ xsv ], and 
c     the beta angle.   The Kouba routine has the yaw rates in a data statement, 
c     but we'll input them as read from svnav.dat.  We allow recomputing of
c     the event time-tags at each call, so there is no need to index them 
c     with the neclips variable.  Finally reformattted the code for easier
c     reading, made more explicit the block selection, and removed the body-
c     axis rotations, which are done in 'model' for the GAMIT implementation.  

c     R. King January 2014
              
c     Routine eclips_dec2013 copyright NRCAN Geodetic Survey Division 2011,2013   
c     Contact kouba@geod.nrcan.gc.ca 
                       
c       Update History: 
c         Kouba's eclips.f, eclips_dec2013, and eclips_jan2014:
C          Aug. 23, 2011: Change maximum yaw rates for BLK II and BLK IIA to averages 
c                         of JPL reprocessing 1996-2008.
C          Sep 26, 2011: 1. Corrected bug causing Block IIF shadow-crossing to yaw
c                           at 0.06 deg/s even for debeadg > 8 deg
C                        2. Corected/improved IIA recovery logic
C          Dec 18, 2013: 1. Corrected IIF night turns (USAF Doc.)
C                        2. Small neg beta IIF and small pos IIA noon turn
C                           (wrong) directions  for |beta| < 0.9deg
C                        3. PRN/SVN 23 IIA YBIAS= -0.5 deg up to Feb 2004 
c          Jan 24, 2014: Noon recovery corrected if IIA/IIF's have ybias=0.
c                        (not applicable currently, possible for future IIF?)
c                        and small IIA/IIF beta if statement simplied.
c          Jan 16, 2015: Remove the sign reversal for GPS IIR/IIR-M SVs, now
c                        IGS standard and the same as all other GNSS.
c          Jan 10, 2017: Beidou eclipsing GEO (orbit normal) yaw; IGEO & MEDO 
c                        orbit normal for |beta| < 2 (corrected below to 4 deg )
c          Jan 10, 2017: Galileo eclipsing (setting betay=beta0=0) disables 
c                        eclipsing code for Galileo and Beidou   
c          Feb 27, 2017: For GPS, Galileo, a possible beta sign change during noon
c                        or shadow turns is fixed by setting betaini if beta < 0.07 deg.
c                        

c         GAMIT kouba_yaw:
c          Jan 24, 2014: Initial routine checked against tables from Kouba for
c                        IIA PRN06/SV36 eclipse of 26 January 2013, IIR PRN12/SV58
c                        noon turn and night turn of 17 Jun 2011, and IIF PRN25/SV62
c                        noon turn and night turn of 17 Jun 2011. 
c     Input:                         
c       week, sow  GPS week, seconds-of-week for reporting purposes
c       ttag       Argument for the yaw models (see Note 1 below)
c       isat       index in SV array, used only for to status printouts
c       iprn       PRN 
c       isvn       SV number (currently not used may be needed for certain circumstnaces
c       antbody    Body type (works for all GNSS)
c       yrate      Nominal/maximum yaw rate (deg/s) for the SV (from svnav.dat)             
c       ybias      Bias in degrees of the yaw 
c       beta       SV beta angle (between SV NORMAL and Sun) (input in radians)
c       svbcos     cosine of angle between the SV vector and the Sun vector
c       xsv(6)     Earth-fixed SV position and velocity              
c       vsvc       SV inertial velocity vector 
c       santxyz    Body X-axis unit vector (computed by gamit/lib/attit_yaw.f)  
c                 (Warning: the IIA body-x orientation reversed for the IIR)


c     Output:
c       yaw_angle  Yaw angle, nominal or event (deg)
c       ievent  Saved values (in common/yaw/ in calling routines) of flags indicating 
c                the current status of SV events. Used internally to report start/end of 
c                events and passed out to be written on the y-file for use by MODEL.
c                 0 = no event, nominal yaw 
c                 1 = night turn or eclipse for GPS/Glonass or beta < 4 deg for Beidou
c                -1 = recovery phase of eclipse
c                 2 = noon turn 
         
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
c        2) Kouba makes no distinction between Block IIR-A and ``Block IIR-B, so 
c           his Block IIF is iblk=6 whereas the GAMIT Block IIF is iblk=7 

      implicit none             

      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      
      integer*4 ievent,isat,isvn,iprn,ieclips
     .        , week,yr,doy,hr,min,i,j
      
      real*8 sow,ttag,xsv(6),beta,yangle,santxyz(3),ybias,vsvc(3)
     .     , svbcos,eclstm,ecletm,anoon,anight
     .     , cnoon,cnight,dir,dttag,murate,det,yrate,betadg,phi
     .     , santx,santy,v(3),r(3),yawend,twohr,halfhr,pi,dtr
     .     , yaw_angle,sec,utcoff 
     .     , svbcos_start,svbcos_check    

                                  
      character*20 antbody
      character*256 message

      logical noon,night,blockII,blockIIR,blockIIF,glonass,beidou
     .       ,nominal
                 
                     
c Set the yaw biases: Blocks I, II, and IIR = 0 
c                     Block IIA     = +0.5
c                     SV 23 (Block IIA) = -0.5 (decommisioned 2004/2/23)
c                     BLock IIF = -0.5 for small beta angles only

c** rwk temporary for debug:
      integer*4 iepoch
      logical debug /.false./
      real*8 temp1, temp2, temp3, temp4
              

c Set the block flags - simplifies slightly the logical tests
c vs testing svantbody - 
      block II = .false.
      blockIIR = .false.
      blockIIF = .false.
      glonass = .false.
      beidou = .false. 
      if( antbody(1:9).eq.'BLOCK I  '.or.antbody(1:9).eq.'BLOCK II '.or.
     .    antbody(1:9).eq.'BLOCK IIA' ) then
        blockII = .true.
      elseif( antbody(1:9).eq.'BLOCK IIR' ) then
        blockIIR = .true.
      elseif( antbody(1:9).eq.'BLOCK IIF' ) then
        blockIIF = .true.
      elseif( antbody(1:9).eq.'GLONASS-M' ) then
        glonass = .true.
      elseif( antbody(1:6).eq.'BEIDOU' ) then 
        beidou = .true.
      else     
        if(debug) print *,'KOUBA_YAW  iepoch isat antbody '
     .        ,iepoch,isat,antbody
        call report_stat('FATAL','YAWTAB','orbits/kouba_yaw'
     .               ,' ','SV body type not recognized',0)
      endif
c     Kouba makes no distinction between Block IIR-A and ``Block IIR-B, so 
c     his Block IIF is iblk=6 whereas the GAMIT Block IIF is iblk=7 
                      
**DEBUG 
cd     if( iepoch.ge.737.and.iepoch.lt.761.and.iprn.eq.30 ) then
cd      if( iepoch.eq.1 ) then
cd      if( iprn.eq.12 .or .iprn.eq.25 )  then
cd      if( iprn.eq.30 ) then 
cd      if( iprn.eq.32  ) then
cd          debug = .true.
cd        else
cd          debug = .false.
cd       endif    


c Set the yaw bias for each system --- this now done in yawtab
c      ybias = 0.d0                                  
c      if( blockII ) ybias = 0.5d0 
c      if( isvn.eq.23 .or. blockIIF ) ybias = -0.5d0
       
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

c Redefine beta: Kouba has the input as the angle of the sun wrt to the orbit 
c normal, but all the tests require the angle of the sun wrt to the SV radius
      betadg = beta/dtr - 90.d0 
 
C Compute the noon beta angle limit (beta zero) for a noon turn from the 
c maximum yaw rate for the SV and the actual orbit rate  (~0.00836 deg/s
c for GPS, ~ 0.00888 deg/s for Glonass)  
       murate= dsqrt((vsvc(1)**2+vsvc(2)**2+vsvc(3)**2)/
     & (xsv(1)**2+xsv(2)**2+xsv(3)**2))/dtr
      anoon=datan(murate/yrate)/dtr
c       rwk:  This routine computes critical beta angle (beta0) for the noon 
c       and midnight turns from the orbit rate and yaw-rate, as shown above, 
c       but does not explictly compute the beta angle for night shadowing.
c       The paper assumes that this is 13.25 deg, which I'll use here for
c       now, to be modified possibly later from the GAMIT shadow routine.
        anight = 180.d0 + 13.25d0 
      if( glonass ) anight = 180.d0 + 14.20d0
      if( beidou ) anight = 180.d0 + 4.d0
      cnoon=dcos(anoon*dtr)
      cnight=dcos(anight*dtr)
c     IIRs use nominal yaw during shadowing and Beidou swtiches to orbit-normal
c     (yaw=0) when beta < beta0 (same as noon for IIR and Beidou)
      if( blockIIR.or.beidou ) then
         cnight = dcos((anoon+180.d0)*dtr)
         if( debug ) print *,'Blk IIR or Beidou update cnight ',cnight
      endif
                            
c Set the Glonass noon-turn yaw limit according to Dilssner 2010
      if( glonass ) then 
        yawend = 75.d0
c       interate the computation based on the orbital geometry
        do j=1,3
          yawend=dabs(datan2(-dtan(betadg*dtr),dsin(pi-
     .      dtr*murate*yawend/yrate))/dtr -      
     .      datan2(-dtan(betadg*dtr),dsin(pi+
     .      dtr*murate*yawend/yrate))/dtr)/2.D0 
        enddo 
c       update the limit using the new yawend        
        anoon = murate*yawend/yrate
        cnoon = dcos(anoon*dtr)
      endif

c Set the night/noon event flags from the beta angle and the beta-angle limits
c     use the beta angle (cos) from the beginning of the event (ievent and
c     svbcos_start are saved variables)
      svbcos_check = svbcos
      if( ievent.ne.0 ) then 
         svbcos_check = svbcos_start
       else
         svbcos_check = svbcos                              
         svbcos_start = svbcos
      endif 
      if( svbcos_check .lt. cnight ) then
        night = .true.                  
        if( ievent.eq.0 ) svbcos_start = svbcos
      endif
      if( svbcos_check .gt. cnoon ) then
        noon = .true.   
        if( ievent.eq.0 ) svbcos_start = svbcos
      endif                       
  
      if( debug ) then        
        print *,' '     
        call debugtime(week,sow)   
        print *,' iepoch sow iprn antbody yrate beta betadg murate'
     .           ,iepoch,sow,iprn,antbody,yrate,beta,betadg,murate
        print *,'  xsv vsvc santxyz ',xsv,vsvc,santxyz   
       if(debug) print *,
     .     'ievent svbcos svbcos_start svbcos_check night noon '
     .     ,ievent,svbcos,svbcos_start,svbcos_check
     ,     ,night,noon
        print *,'anoon cnoon anight cnight yrate '
     .       ,anoon,cnoon,anight,cnight,yrate
      endif 
        
c Compute the nominal yaw angle from the orbit geometry
      yangle= dacos( (santxyz(1)*vsvc(1)+santxyz(2)*vsvc(2)+
     .               santxyz(3)*vsvc(3) ) / 
     .                 sqrt(vsvc(1)**2+vsvc(2)**2+vsvc(3)**2) )/dtr  
c     For IIR, same sign as beta; for II/IIA, Glonass, and Beidou, opposite sign                  
      if( betadg.lt.0.d0 .and. blockIIR ) yangle = -yangle
      if( betadg.gt.0.d0 .and. .not.blockIIR ) yangle = -yangle
      if( debug ) print *,'betar svbcos betadg yangle '
     .                   , beta,svbcos,betadg,yangle
cd      if( iepoch.gt.1 ) stop 1 
     
c If in a night or noon maneuver, get the start/stop times 
c and the actual yaw

      if( night.or.noon ) then

c       indicate a potential event, reset later if non-nominal indicated
        nominal = .false.
     
c       Beidou switches to orbit nominal when beta < beta0, so set
c       yaw = 0 and exit
        if( beidou ) then
          phi = 0.d0
          if( debug ) print *,'Beidou beta=',beta/dtr,', set phi=0'
          goto 1
        endif  
 
c       get the orbital angle for the start (Kouba 2009 Eqs 8-9)
        det = dsqrt( (180.d0 - dacos(svbcos)/dtr)**2 - betadg**2 )
        if(debug) print *,'first det ',det 
c       initialize the yaw angle at 90 degrees (rwk: not sure why)
        phi = pi/2
C       check if already after a midnight or noon
        if( night ) then 
          if( blockIIR ) then 
            if( dabs(yangle).gt.90d0 ) det = -det
            if( det.ne.0.d0 ) 
     .        phi = datan2( dtan(betadg*dtr),-dsin(-det*dtr) )/dtr
            if( debug ) print *,'BLK IIR NIGHT det  phi ',det,phi 
          else 
c           II, IIA, IIF, Glonass
            if( dabs(yangle).lt.90.d0 ) det = -det
            if( det.ne.0.d0) 
     .         phi = datan2( -dtan(betadg*dtr),dsin(-det*dtr) )/dtr
            if( debug ) 
     .         print *,'II,IIA, IIF & Glonass NIGHT det phi ',det,phi  
          endif
        endif
      
        if( noon ) then 
          det = sqrt((dacos(svbcos)*180.d0/pi)**2-betadg**2)
          if( blockIIR ) then
            if( dabs(yangle).lt.90.d0 ) det = -det
            if( det.ne.0.d0 ) 
     .        phi = datan2(dtan(betadg*dtr),-dsin(pi-det*det))/dtr
            if( debug ) print *,'BLKIIR NOON det phi ',det,phi 
          elseif ( blockII .or. blockIIF .or. glonass  ) then
            if( dabs(yangle).gt.90.d0) det =-det
            if( det.ne.0d0 ) 
     .        phi = datan2(-dtan(betadg*dtr),dsin(pi-det*dtr))/dtr    
             if( debug ) print *,'not BLKIIR NOON det phi ',det,phi 

          endif
        endif
                              
        eclstm = ttag + det/murate   
        yawend = datan(murate/yrate)/dtr  
        if(debug) print *,' ttag sow det murate eclstm yawend '
     .                    , ttag,sow,det,murate,eclstm,yawend

c       For IIR/Glonass night turn only makes sense when beta < beta-zero
c       For IIA it gets here only when noon is true and that happens only 
c         when beta < beta-zero (rwk: I don't see that in the logic)
        if( (blockIIR.or.noon) .and. dabs(betadg).lt.yawend) then 
           if( glonass ) then 
             eclstm= eclstm- anoon/murate
             ecletm= eclstm + 2.d0*anoon/murate
           else
c            IIA/IIR/IIF noon or IIR midnight turn
             eclstm= eclstm -
     .          dabs(betadg)*sqrt(anoon/dabs(betadg)-1.d0)/murate
             ecletm = eclstm +
     .          2.d0*dabs(betadg)*sqrt(anoon/dabs(betadg)-1.d0)/murate
             if(debug) print *,'NOON OR NIGHT  eclstm ecletm'
     .                      , eclstm,ecletm
            endif
          endif   
c         II/IIA/IIF/Glonass shadow start/end times (IIR/Beioud keep nominal yaw during shadow)
          if( .not.blockIIR .and.night ) then 
            if( debug ) print *,'eclstm(mid) ',eclstm
            eclstm= eclstm -
     .             dsqrt((anight-180.d0)**2-betadg**2)/murate
            ecletm= eclstm +
     .             2.d0*dsqrt((anight-180d0)**2-betadg**2)/murate
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
            if( blockII.or.blockIIF ) then
              if( ttag.le.ecletm ) then
                if( ievent.eq.0 ) then
                  ievent = 1 
                  call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
                  write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .             antbody,' PRN',iprn,' eclipsing    ',doy,hr,min
     .              ,'  beta=',betadg
                  call report_stat('STATUS','YAWTAB','orbits/kouba_yaw'
     .                   ,' ',message,0)   
                endif
                if( blockII ) then    
c                 IIA night turn   (rwk: why not a simply a shadow?) 
                  phi=datan2(-dtan(betadg*dtr),dsin(-det*dtr))/dtr
     .                  +sign(yrate,0.5d0)*(ttag-eclstm)
                  if( debug ) print *,'IIA night turn phi ',phi    
                endif
                if( blockIIF ) then
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
                  call report_stat('STATUS','YAWTAB','orbits/kouba_yaw'
     .                   ,' ',message,0)   
                endif
                if( blockII ) then
                   phi=datan2(-tan(betadg*dtr),dsin(-det*dtr))/dtr
     .                 + sign(yrate,ybias)*(ecletm-eclstm)
                   if(debug) print *,'IIA shadow exit phi ',phi 
                endif
                if( blockIIF ) then 
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
            if( glonass ) then    
c             if yaw completed, set to nominal and exit
              if( ttag.gt.ecletm ) then
                nominal = .true.
                go to 1
              else
                yawend = yrate
                phi= datan2(-dtan(betadg*dtr),dsin(-det*dtr))/dtr
     .                    + sign(yawend,betadg)*(ttag-eclstm)            
c               yaw angle at Glonass shadow exit
                yawend = datan2(-dtan(betadg*dtr),dsin(det*dtr))/dtr
                if(yawend/phi.gt.1.d0 .or. phi/yawend.lt.0.d0) 
     .             phi = yawend
              endif 
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
                call report_stat('STATUS','YAWTAB','orbits/kouba_yaw'
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
              call report_stat('STATUS','YAWTAB','orbits/kouba_yaw'
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
            if( (blockII.or.blockIIF) .and.
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
            elseif( glonass ) then
c               see if a GLONASS or Beidou  exit is detected
                if( ttag.gt.ecletm) then
                   nominal = .true.
                   goto 1
                endif
            elseif ( blockII .or. blockIIR. or. blockIIF) then
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
          call report_stat('STATUS','YAWTAB','orbits/kouba_yaw',' '
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
      if( debug.and.iepoch.gt.315 ) stop 
      return
      end

      







