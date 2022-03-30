      Subroutine kouba_beidou( antbody,vsvc,santxyz,betadg,svbcos
     .                       , week,sow,isat,iprn,isvn
     .                       , yaw_angle,ievent,ybias 
c                             this for debug
     .                       , iepoch )

c     GAMIT version of the GPS portion of Jan Kouba's elips_feb2017 routine. 
c     See kouba_gps for an explanation of Kouba's original routine.  Beidou
c     employs nominal yaw in the ideal yaw-steering frame except for converting 
c     to orbit-normal, used for MEO and IGSO for small beta-angle, and always
c     for GEO.  Eclipsing events are detected here and reported in GAMIT.status,
c     but not used in the SV attitude computation.  The yaw-steering or orbit-
c     normal frame is implemented in model/svbody_coords.f based on the 'ievent'
c     flag computed here and written on the y-file. 

c     R. King January 2014 / May 2017
              
c     Based on Jan Kouba's routine eclips_May2017    
c     Contact kouba@rogers.com
c     Copyright Geodetic Survey Division, 2011. All rights reserved. 
c     All terms and conditions apply as detailed in "Terms and conditions for software:
                       
c     Update History: 
c       Kouba's eclips.f, eclips_dec2013, and eclips_jan2014, eclips_Feb2017, eclips_May2017
c       See ftp/kouba/eclips_May2017.f 

c Input:    
c   Required for computing the nominal yaw:           
c   vsvc(3)    SV inertial velocity vector    
c   santxyz(3) Body X-axis unit vector (computed by gamit/lib/attit_yaw.f)  
c   betadg     Angle between the Sun vector and the SV orbital frame (deg)
c   svbcos     Cosine of angle between the SV vector and the Sun vector  
c Passed for reporting purposes only:
c   week, sow  GPS week, seconds-of-week
c   isat       Index of SV for event tracking and reporting
c   iprn,isvn,antbody: PRN, SVN, and antenna-body-type

c Ouput:    
c   Required for y-file
c   yaw_ang    Yaw angle (deg) to be written on the y-file
c   ievent     Flag to indicate the attitude status       
c                 0 = nominal yaw 
c                 3 = orbit normal
c    Used only for reporting in yawtab.out.DDD 
c    ybias     Yaw bias (deg)

        
c Internal  
c     beta0        Beta angle below which orbit normal mode is invoked (4 deg)
c     pi           3.1415926...

      implicit none             

      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      
      integer*4 ievent,isat
     .        , iprn,isvn,week,yr,doy,hr,min
      
      real*8 beta,santxyz(3),vsvc(3),svbcos,anoon,anight
     .     , cnoon,cnight,betadg,pi,dtr,yaw_angle,beta0,ybias 
     .     , sow,sec,utcoff
                          
      character*20 antbody
      character*256 message

c** rwk temporary for debug:
      integer*4 iepoch
      logical debug /.false./

c Set some constants
      pi=4.d0*datan(1.d0) 
      dtr=pi/180.D0    
c     beta angle for transition to orbit-normal mode
      beta0 =4.d0
        
c  GEO always orbit-normal (yaw=0)
c  IGSO and MEO orbit-normal if beta < 4 deg; otherwise nominal

      if( antbody(9:9).eq.'G' ) then
c        GEO always orbit-normal, no need to report
         yaw_angle = 0.d0
         ievent = 3                
      else
c        IGSO, MEO orbit normal only when beta < 4 deg 
        if( dabs(betadg).le.beta0 ) then
          ievent = 3
          yaw_angle = 0.d0
          if( ievent.eq.0 ) then 
            call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
            write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .        antbody,' PRN',iprn,' start orbit normal ',doy,hr,min
     .               ,'  beta=',betadg     
              call report_stat('STATUS','YAWTAB','orbits/kouba_beidou'
     .               ,' ',message,0)             
          endif 
        else
c         yaw-steering 
          if( ievent.eq.3 ) then            
          call timcon(4,week,sow,yr,doy,hr,min,sec,utcoff)
            write(message,'(a,a,i3,a,3i3,a,f6.2)')
     .        antbody,' PRN',iprn,' exit orbit normal ',doy,hr,min
     .               ,'  beta=',betadg     
              call report_stat('STATUS','YAWTAB','orbits/kouba_beidou'
     .               ,' ',message,0)             
          endif
          ievent = 0 
          yaw_angle = dacos( (santxyz(1)*vsvc(1)+santxyz(2)*vsvc(2)+
     .                       santxyz(3)*vsvc(3) ) / 
     .                  sqrt(vsvc(1)**2+vsvc(2)**2+vsvc(3)**2) )/dtr
c         IGS standard is yaw angle opposite sign from beta
          if( betadg.gt.0.d0 ) yaw_angle = -yaw_angle
        endif
        if( debug ) print *,'betadg beta0 yaw_aangle '
     .                     , betadg,beta0,yaw_angle
      endif
                                                                                 
      return
      end

      







