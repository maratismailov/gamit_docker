      Subroutine kouba_glonass( ttag,antbody,yrate
     .                    , xsv,vsvc,santxyz,betadg,svbcos  
     .                    , week,sow,isat,iprn,isvn
     .                    , yaw_angle,ievent,ybias
c      this for debug
     .                    , iepoch )

c     GAMIT version of the GPS portion of Jan Kouba's elips_feb2017 routine.  
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
c     Contact kouba@geod.nrcan.gc.ca 
                       
c       Update History: see kouba_gps.f
                

c     Input:    
c       Required for computations:        
c       ttag       Time argument, seconds from t/y-file start 
c       antbody    Body type (works for all GNSS)
c       yrate      Nominal/maximum yaw rate (deg/s) for the SV (from svnav.dat)          
c       xsv(3)     Earth-fixed SV position and velocity              
c       vsvc(3)    SV inertial velocity vector 
c       santxyz(3) Body X-axis unit vector (computed by gamit/lib/attit_yaw.f)  
c                 (Warning: the IIA body-x orientation reversed for the IIR)   
c       betadg     Angle between Sun and SV plane (deg)
c       svbcos     Cosine of angle between the SV vector and the Sun vector
c       Used only for reporting purposes
c       week, sow  GPS week, seconds-of-week
c       isat,iprn,isvn  SV array index, PRN and SVN 

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
c                 3 = orbit-normal mode (not used for Glonass)
c       Used only for reporting in yawtab.out.DDD 
c       ybias     Yaw bias (deg)

         
c Internal
c       ybias   y-bias value in degrees
c       eclstm  Start time for the current eciipses or noon turn (units as ttag) 
c       ecletm  End times for eciipses and noon turns for the SV (units as ttag_ 
c       pi         3.1415926...
C       anoon    SV beta angle limit (deg) for a noon turn maneuver
c       anight   SV beta angle limit (deg) for a night shadow crossing
c       yangle   Nominal yaw angle (deg)
c       phi      Actual yaw angle if during an event (deg)
c       svbcos_start(maxsat),betaini(maxsat) : saved variables for each SV 
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
c   180.D0+-14.20D0 for GLONASS (IPRN.GT.32). CNIGHT=COS(ANIGHT) is used for 
c   testing. c   When SVBCOS < CNIGHT (CNIGHT is negative and close to -1), 
c   the SV enters, or it is in the shadow. 
c
c Notes: Kouba uses GPS second-of-week for his time argument, but without 
c        some accommodation, this can produce a discontinuity. We keep the
c        GPS week, seconds-of-week for easy comparison with his test cases
c        but use the continuous orbital argument (satarg) as the primary 
c        argument for the yaw model

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
     .     , svbcos_start(maxsat),svbcos_check,betaini(maxsat)

                                  
      character*20 antbody
      character*256 message

      logical noon,night,nominal
                 
c** rwk temporary for debug:
      integer*4 iepoch
      logical debug /.false./
      real*8 temp1, temp2, temp3, temp4
         
c These not passed so save for checking
      data betaini/maxsat*0.d0/ 
      save svbcos_start,betaini
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
      ybias = 0.d0 

c Redefine beta: Kouba has the input as the angle of the sun wrt to the orbit 
c normal, but all the tests require the angle of the sun wrt to the SV radius
      betadg = beta/dtr - 90.d0 
 
C Compute the noon beta angle limit (beta zero) for a noon turn from the 
c maximum yaw rate for the SV and the actual orbit rate ~ 0.00888 deg/s.
       murate= dsqrt((vsvc(1)**2+vsvc(2)**2+vsvc(3)**2)/
     & (xsv(1)**2+xsv(2)**2+xsv(3)**2))/dtr
      anoon=datan(murate/yrate)/dtr
c       rwk:  This routine computes critical beta angle (beta0) for the noon 
c       and midnight turns from the orbit rate and yaw-rate, as shown above, 
c       but does not explictly compute the beta angle for night shadowing.
c       The paper assumes that this is 13.25 deg, which I'll use here for
c       now, to be modified possibly later from the GAMIT shadow routine.
      anight = 180.d0 + 14.20d0
      cnoon=dcos(anoon*dtr)
      cnight=dcos(anight*dtr)

c Set the noon-turn yaw limit according to Dilssner 2010
      if( dabs(betadg).lt.anoon ) yawend = 75.d0
c     iterate the computation based on the orbital geometry
      do j=1,3
        yawend=dabs(datan2(-dtan(betadg*dtr),dsin(pi-
     .      dtr*murate*yawend/yrate))/dtr -      
     .      datan2(-dtan(betadg*dtr),dsin(pi+
     .      dtr*murate*yawend/yrate))/dtr)/2.D0 
      enddo 
c     update the limit using the new yawend        
      anoon = murate*yawend/yrate
      cnoon = dcos(anoon*dtr)

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
     .         betaini(isat)=betadg
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
        print *,' iepoch sow iprn antbody yrate beta betadg murate'
     .           ,iepoch,sow,iprn,antbody,yrate,beta,betadg,murate
        print *,'  xsv vsvc santxyz ',xsv,vsvc,santxyz   
       if(debug) print *,
     .     'ievent svbcos svbcos_start svbcos_check night noon '
     .     ,ievent,svbcos,svbcos_start(isat),svbcos_check
     ,     ,night,noon
        print *,'anoon cnoon anight cnight yrate '
     .       ,anoon,cnoon,anight,cnight,yrate
      endif 
        
c Compute the nominal yaw angle from the orbit geometry
      yangle= dacos( (santxyz(1)*vsvc(1)+santxyz(2)*vsvc(2)+
     .               santxyz(3)*vsvc(3) ) / 
     .                 sqrt(vsvc(1)**2+vsvc(2)**2+vsvc(3)**2) )/dtr  
c     the yaw angle has the opposite sign as beta
      if( betadg.gt.0.d0 ) yangle = -yangle
      if( debug ) print *,'betar svbcos betadg yangle '
     .                   , beta,svbcos,betadg,yangle
cd      if( iepoch.gt.1 ) stop 1 
     
c If in a night or noon maneuver, get the start/stop times 
c and the actual yaw

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
     .       print *,'II,IIA, IIF & Glonass NIGHT det phi ',det,phi  
        endif
      
        if( noon ) then 
          det = sqrt((dacos(svbcos)*180.d0/pi)**2-betadg**2)
          if( dabs(yangle).gt.90.d0) det =-det
          if( det.ne.0d0 ) 
     .      phi = datan2(-dtan(betadg*dtr),dsin(pi-det*dtr))/dtr
        endif
                              
        eclstm = ttag + det/murate   
        yawend = datan(murate/yrate)/dtr  
        if(debug) print *,' ttag sow det murate eclstm yawend '
     .                    , ttag,sow,murate,eclstm,yawend

c       For Glonass night turn only makes sense when beta < beta-zero
        if( (noon) .and. dabs(betadg).lt.yawend) then 
           eclstm= eclstm- anoon/murate
           ecletm= eclstm + 2.d0*anoon/murate
        endif   
c       shadow start/end times 
          if( night ) then 
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
          ievent = 0 
          go to 1 
        else
c         otherwise calculate the modeled yaw
          if(debug) print *,'COMPUTING modeled yaw'
c         for small beta, use initial value throughout the eclipse 
c         to unsure no beta sign-change
          if(dabs(betadg).le..07d0.and.betaini(isat).ne.0.d0) 
     .       betadg = betaini(isat)
c         S/C position and velocity unit vectors
          do j=1,3
            v(j)=vsvc(j)/dsqrt(vsvc(1)**2+vsvc(2)**2+vsvc(3)**2)
            r(j)=xsv(j)/dsqrt(xsv(1)**2+xsv(2)**3+xsv(3)**2)
          enddo
c         orbit angle at start of event
          det = murate*(ecletm-eclstm)/2.d0 

          if( svbcos.lt.0.d0 ) then
c           SHADOW CROSSING OR NIGHT TURN    
            nominal = .false.
            if( debug ) print *,'NIGHT antbody ttag sow ecletm ievent '
     .                         ,antbody,ttag,sow,ecletm,ievent
c           if yaw completed, set to nominal and exit
            if( ttag.gt.ecletm ) then
              nominal = .true.
              go to 1
            else  
              yawend = yrate
              phi= datan2(-dtan(betadg*dtr),dsin(-det*dtr))/dtr
     .                 + sign(yawend,betadg)*(ttag-eclstm)            
c             yaw angle at Glonass shadow exit
              yawend = datan2(-dtan(betadg*dtr),dsin(det*dtr))/dtr
              if(yawend/phi.gt.1.d0 .or. phi/yawend.lt.0.d0) 
     .            phi = yawend
              if( ievent.eq.0 ) then  
                 ievent = 1
                 call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
                 write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .              antbody,' PRN',iprn,' eclipsing    ',doy,hr,min
     .              ,'  beta=',betadg
                 call report_stat('STATUS','YAWTAB','orbits/kouba_gps'
     .                   ,' ',message,0)   
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
c** TEMP: does this apply for Glonass?
c           if a II/IIA/IIF small-angle problem or exit is detected 
            phi= datan2(-dtan(betadg*dtr),dsin(pi-det*dtr))/dtr
     .                 -dsign(yrate,betadg)*(ttag-eclstm)
            if(debug) then 
              print *,'standard noon turn phi ',phi          
              print *,'betadg ybias ',betadg,ybias 
              print *,'betadg*dsign(1,ybias)' ,betadg*dsign(1.d0,ybias)
              print *,'betadg*ybias ',betadg*ybias 
            endif
c           see if a GLONASS or Beidou  exit is detected
            if( ttag.gt.ecletm) then
              nominal = .true.
              goto 1
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
      return
      end

      







