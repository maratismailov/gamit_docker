      Subroutine irnss_yaw( vsvc,santxyz,betadg
     .                    , yaw_angle )
             
c     Yaw computation for IRNSS.  The SV axis orientation is the yaw-steering
c     frame with no accommodation for eclipses except that all three axes
c     are offset to direct the antenna to a ground point that will maximize
c     signal strength for India. This routine simply computes the nominal
c     yaw; the offsets are applied in model/svbody_coords.f.  However, shadowing
c     and "noon" events are detected here and reported as information in 
c     GAMIT.status.
                                                                      
c     R. King January 2014 / May 2017
              

c Input:              
c   Required for computing the nominal yaw:           
c   vsvc(3)    SV inertial velocity vector    
c   santxyz    Body X-axis unit vector (computed by gamit/lib/attit_yaw.f)  
c   betadg     Angle between the Sun vector and the SV orbital frame (deg)

c Ouput: 
c   yaw_ang    Yaw angle (deg) to be written on the y-file


c Internal                       

      implicit none             

      include '../includes/dimpar.h'
      include '../includes/orbits.h'
      
      real*8 betadg,santxyz(3),vsvc(3),svbcos,yaw_angle,pi,dtr

      character*20 antbody
      character*256 message

      logical debug/.false./
                                                              
      pi=4.d0*datan(1.d0) 
      dtr=pi/180.D0    
     
c  IRNSS always in yaw-steering mode with no shadowing events        

      yaw_angle = dacos( (santxyz(1)*vsvc(1)+santxyz(2)*vsvc(2)+
     .               santxyz(3)*vsvc(3) ) / 
     .             sqrt(vsvc(1)**2+vsvc(2)**2+vsvc(3)**2) )/dtr
c     IGS standard is yaw angle opposite sign from beta
        if( betadg.gt.0.d0 ) yaw_angle = -yaw_angle
      if( debug ) print *,'betadg yaw_aangle '
     .                   , betadg,yaw_angle
                                                                                       
      return
      end

      







