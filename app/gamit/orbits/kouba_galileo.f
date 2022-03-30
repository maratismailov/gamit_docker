      Subroutine kouba_galileo( isat,yrate,betadg,svbcos,xsv,vsvc,xsun
     .                        , ttag,week,sow,iprn,isvn,antbody
     .                        , yaw_angle,ievent
c                             this for debug
     .                        , iepoch )
                           

c     Galileo yaw based on Jan Kouba's 'Notes on Feb 2017 version of the eclips
c     subroutine'. The naming and structure in this version follows the memo
c     rather than the eclips routine. 
c
c     R. King January 2014 / May 2017
              
c     Routine eclips_dec2013 copyright NRCAN Geodetic Survey Division 2011,2013   
c     Contact kouba@rogers.com 
c    
c     Input:    
c       Required for computations:                     
c       isat       Index in SV array for ievent and saved angles 
c       yrate      Nominal/maximum yaw rate (deg/s) for the SV (from svnav.dat)     
c       betadg     Angle between Sun and SV plane (deg)
c       svbcos     Cosine of angle (mu) between the SV vector and the Sun vector
c       xsv(3)     Radius vector of the satellite
c       vsvc(3)    Velocity vector of the satellite
c       xsun(3)    Vector of the Sun wrt the Earth 
c       Used only for reporting purposes
c       ttag       Time argument, seconds from t/y-file start
c       week, sow  GPS week, seconds-of-week
c       isat,iprn,isvn  SV array index, SV PRN and SVN 
c       antbody    Body type (works for all GNSS)

c     Output:  
c       Required for y-file                                         
c       yaw_angle  Yaw angle, nominal or event (deg)
c       ievent  Saved values (in common/yaw/ in calling routines) of flags indicating 
c                the current status of SV events. Used internally to report start/end of 
c                events and passed out to be written on the y-file for use by MODEL.
c                 0 = no event, nominal yaw 
c                 1 = night turn or eclipse 
c                -1 = recovery phase of eclipse (not used for Galileo)
c                 2 = noon turn              
c                 3 = orbit-normal mode (not used for Galileo)
c       Used only for reporting in yawtab.out.DDD 
c       ybias     Yaw bias (deg)

c Internal
c       cnoon    cos(15 deg; if svbcos > cnoon, SV begins a noon turn
c       cnight   cos(195 deg); if svbcos > cnight SV begins a night turn  
c       beta     beta angle in radians, input or retained turning turn 
c       betaini(maxsat)  saved values of beta at turn beginning to avoid sign change 
c                        for very small beta ( < 0.07 deg)
c       mu       Angle of the SV in the orbital plane, measured from midnight
c       Sx, Sy, Sz  Sun angle x, y, and z components
c       Shy      Sun angle y-component when within a turn
c       Shz      Sun angle z-component when within a turn 
                           

c All angles internally in radians except input betadg and output yaw_angle 

      implicit none             

      include '../includes/dimpar.h'
      
      integer*4 ievent,isat,isvn,iprn,week,yr,doy,hr,min
      
      real*8 sow,ttag,yaw_angle,ybias,vsvc(3),svbcos,cnoon,cnight
     .     , yrate,betadg,beta,betax,betay,mu,Sx,Sy,Sz,Shy,Shz
     .     , gamma(maxsat),rad,betaini(maxsat),pi,dtr,sec,utcoff
     .     , xsv(3),xsun(3) 

      character*20 antbody
      character*256 message

c** rwk temporary for debug:
      integer*4 iepoch,iprndb/11/
      logical debug/.false./ 

c Save the initial beta angle for abs(beta) < 0.07 deg
c     (ievent also saved via its being passed to yawtab) 
      save betaini

c Conversion constants, dtr = deg-to-rad
      data  pi/3.14159265359d0/,dtr/1.745329252d-2/,betaini/maxsat*0.d0/
     .   , gamma/maxsat*1.d0/

c**DEBUG 
c      if( iepoch.ge.2153.and.iepoch.lt.2286
c     .    .and.iprn.eq.iprndb ) then
c         debug = .true.
c       else
c         debug = .false.
c       endif                                                     
c                                     
c Set the constants and convert betadg 

       betax = 15.d0*dtr
       betay = 2.d0*dtr    
       cnoon = dcos(betax)
       cnight = dcos(betax+pi)        
       beta = betadg*dtr

c Get the yr/doy/hr/min/sec for reporting

      call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
 
c Get the orbital angle from midnight ( mu ) 
                                      
      call mu_angle( xsv,vsvc,xsun,mu,debug,iprndb ) 
      mu = mu*dtr              
                   
c Nominal yaw for beta > 2 deg, otherwise use a modified Sy to restrict yaw rate
                      
      if( debug ) 
     .   print *,'iepoch betay betadg mu,svbcos cnight cnoon '
     .          , iepoch,betay,betadg,mu,svbcos,cnight,cnoon 
     .       
      if(  dabs(beta).le.betay .and.
     .   (svbcos.gt.cnoon.or.svbcos.lt.cnight ) ) then  
c       get the model for turns - same for night and noon 
        if( debug ) 
     .     print *,'turn old ievent  Sx Sy Sz ',ievent,Sx,Sy,Sz
c       see if need to save values for the start of a turn
        if( ievent.eq.0 ) then
           gamma(isat) = dsign(1.d0,Sy)
           betaini(isat) = beta
           if( debug ) print *,'start save gamma beta '
     .                          ,gamma(isat),betaini(isat)
        endif
c       for very small beta used saved value to avoid a beta sign change
        if( dabs(betadg).lt.0.07d0 ) then
          beta = betaini(isat)
        endif
        Sx = dsin(mu)*dcos(beta)
        Sy = -dsin(beta)
        Sz = dcos(mu)*dcos(beta)
        Shy = 0.5d0 * ( dsin(betay)*gamma(isat) + Sy) 
     .        + 0.5d0 * ( dsin(betay)*gamma(isat) - Sy) 
     .              *  dcos(pi*dabs(Sx)/dsin(betax))   
        Shz = dsqrt( 1.d0 - Sx**2 -Shy**2*dsign(1.d0,Sz) )
        rad = dsqrt(1.d0-Shz**2)          
        yaw_angle = datan2(Shy/rad,Sx/rad)/dtr  
        if( debug ) then                                       
          print *,'betadg beta(dg) gamma ',betadg,beta/dtr,gamma(isat)
          print *,'Sx  Sy Sz  rad ',Sx,Sy,Sz,rad 
          print *,'Sx Shy Shz psi ',Sx,Shy,Shz,yaw_angle
        endif
c       report if entering the turn (may be at first epoch)
        if(ievent.eq.0) then 
          if( svbcos.lt.cnight ) then 
            write(message,'(a10,a,i3,a,3i3,a,f6.2)')
     .        antbody(1:10),' PRN',iprn,' night turn ',doy,hr,min
     .         ,'  beta=',betadg          
            ievent = 1
            if(debug) print *,'reset ievent ',ievent
          elseif( svbcos.gt.cnoon ) then 
            write(message,'(a10,a,i3,a,3i3,a,f6.2)')
     .         antbody(1:10),' PRN',iprn,' night turn ',doy,hr,min
     .         ,'  beta=',betadg                                              
            ievent = 2                              
            if(debug) print *,'reset ievent ',ievent
          endif                            
          call report_stat('STATUS','YAWTAB','orbits/kouba_galileo'
     .                     ,' ',message,0) 
c       endif on need to report start of turn  
        endif 
  
      else
c       nominal yaw 
c       Sun unit vector wrt the SV 
        Sx = dsin(mu)*dcos(beta)
        Sy = -dsin(beta)
        Sz = dcos(mu)*dcos(beta)
        rad = dsqrt(1.d0- Sz*Sz)          
        yaw_angle = datan2(Sy/rad,Sx/rad)/dtr 
ccc       no sign reversal for the vector/mu method
ccc       if( betadg.gt.0.d0 ) yangle = -yangle
        if( debug ) then 
          print *,'ievent Sx Sy Sz rad yangle '
     .      ,ievent,Sx,Sy,Sz,rad,yaw_angle
        endif                              
c       report if first epoch after exit
        if( ievent.ne.0 ) then
          if( ievent.eq.1 ) then 
            write(message,'(a10,a,i3,a,3i3,a,f6.2)')
     .          antbody(1:10),' PRN',iprn,' night exit ',doy,hr,min
     .           ,'  beta=',betadg                      
          elseif( ievent.eq.2 ) then
            write(message,'(a10,a,i3,a,3i3,a,f6.2)')
     .          antbody(1:10),' PRN',iprn,' noon  exit ',doy,hr,min
     .           ,'  beta=',betadg              
          endif 
          call report_stat('STATUS','YAWTAB','orbits/kouba_galileo',' '
     .                     ,message,0)                      
          if(debug) print *,'reset ievent ',ievent
          ievent = 0
c       endif on need to report end of turn 
        endif  
        
c     endif on whether within a turn event
      endif 

      if( debug ) print *
     .  ,'leaving KOUBA_GALILEO iepoch isat ievent yaw_angle '
     .                   , iepoch,isat,ievent,yaw_angle 
      return
      end


