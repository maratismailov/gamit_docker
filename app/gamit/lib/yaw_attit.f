       subroutine yaw_attit( svec,sun,y_att,beta,sc_x_t,sc_y_t
     .                     , sc_z,sc_W,antbody)

c
c  Purpose: To compute the nominal attitude of a satellite at any given time. This 
c           routine is used to compute the attitude and beta angles of the satellites 
c           when entering eclipse, and also used each epoch to compute the unit vectors
c           of each satellite (and beta angle) to be used in gps_coords.f and
c           subsequently passed back to model.f
c
c PARAMETERS:
C         IN:
c             svec         : satellite coordinates wrt earth                      R*8(6)
c             sun          : sun coordinates wrt earth                            R*8(6)  
c             antbody      : body type (all GNSS)                                 C*20

c        OUT:
c              y_att       :  current yaw angle of satellite                      R*8
c              beta        :  beta angle at time satellite starts yawing          R*8
c              sc_x_t      :  unit vector of satellite x axis                     R*8(3)
c              sc_y_t      :  unit vector of satellite y axis                     R*8(3)
c              sc_z        :  unit vector of satellite z axis                     R*8(3)
c              sc_W        :  unit vector of satellite z axis                     R*8(3)

c SUBROUTINES CALLED:
c
c CREATED: 5th February 1996                 LAST MODIFIED:  21 May, 1998
c
c AUTHOR: P. Tregoning
c
c COPYRIGHT: DEPARTMENT OF EARTH AND PLANETRY SCIENCES
c            M.I.T. 1995
c 
c PT970821: model block IIR satellites with a left handed coord system on sv (ie the
c           x-axis points away from the sun)    
c         : routine also called by get_dusk to determine the beta angle at dusk times

c SCM980521: model block IIR satellites with a RH coord system x-axis now points to SUN.

c RWK140125: restore block IIR X-axis reversal to be consistent with the yaw model and
c            y-file values); code model/gps_coords.f to undo the reversal for consistency
c            with the other satellites
c            
       implicit none
c
       real*8 svec(6),sun(6)
     .       ,sc_pos(3),sc_vel(3),sc_sun_dir(3)
     .       ,earth_sc_dist,sc_Z(3),amag3,eps_angle,dot,rad2deg
     .       ,sc_x_t(3),sc_y_t(3),sc_vel_dir(3)
     .       ,sc_w(3),orbit_normal(3),beta,earth_sun_dir(3)
     .       ,a,y_att
      integer*4 i
      character*20 antbody       
   
       data rad2deg/57.29577951308232d0/
              


c PT960205: The following code has been moved from gps_coords.f into this routine
              
c RWK980113: Shift satellite coordinates into 3-vectors to use old code 
       do i=1,3
         sc_pos(i) = svec(i)
         sc_vel(i) = svec(i+3)
       enddo 

c** Create the s/c-Sun vector and normalize it.
       sc_sun_dir(1)=sun(1)-sc_pos(1)
       sc_sun_dir(2)=sun(2)-sc_pos(2)
       sc_sun_dir(3)=sun(3)-sc_pos(3)
       call normalise(sc_sun_dir,sc_sun_dir)

c** Normalize the Sun vector.
       call normalise(sun,earth_sun_dir)

c Normalize the velocity vector.
       call normalise(sc_vel,sc_vel_dir)

c** Define the s/c Z coordinate (pointing toward geocenter).
       earth_sc_dist=amag3(sc_pos)
       sc_Z(1)=-sc_pos(1)/earth_sc_dist
       sc_Z(2)=-sc_pos(2)/earth_sc_dist
       sc_Z(3)=-sc_pos(3)/earth_sc_dist

c
c** Compute the Earth-Probe-Sun (EPS angle).
       EPS_angle=acos(dot(sc_Z,sc_sun_dir))*rad2deg
c
c** The computation of the X and Y axis may be singular. When it is avoid
c   the computation and take saved values (from last step) instead.
       if (EPS_angle .ge. 1.0d-6 .or.
     &     EPS_angle .le. 180.0d0-1.0d-6) then   ! not a singularity
c
c Compute the s/c Y coordinate. It is the normal vector to the
c Earth-Probe-Sun plane.
          call cross(sc_Z,sc_sun_dir,sc_Y_t)  ! get an orthogonal vector

c  Normalize the Y vector
          call normalise(sc_Y_t,sc_Y_t)


c Compute the X coordinate by completing the system. RHS for blocks 1-3, LHS for block 4 
c The x-axis lies in the Earth_sv_sun plane. 
c BLOCK 4 satellites are now computed in RHS for consistancy with model. McClusky 980518 
c ** RWK 140108: Reverse this change to be consistent with the Kouba yaw model (do the
c                same thing in MODEL and ARC, i.e. different rotation for different blocks)
c ** RWK 170523: Remove this with the new yawtab---IIRs use IGS convention
c*            if( block.ge.4.and.block.le.6 ) then
cc            if( antbody(1:9).eq.'BLOCK IIR' ) then
cc              call cross(sc_Z,sc_Y_t,sc_X_t)    
cc            else
              call cross(sc_Y_t,sc_Z,sc_X_t) 
cc            endif  
c*** Need to code this for GLONASS, Beidou, etc.              

c      end test on singularity
       endif          
       
           

c** Define the nominal yaw angle. This is the angle between the s/c X axis
c   and a vector normal to sc_Z, on the orbit plane.
c
c
c Compute the unit vector sc_W that spans the orbit plane together with sc_Z
c (the orbit plan is defined by the s/c velocity and position vectors) and
c has the same general direction as the s/c velocity.
c       cosa=dot(sc_Z,sc_vel_dir)
c
c       sina=dsqrt(1.0d0-cosa*cosa)
c       sc_W(1)=(sc_vel_dir(1)-cosa*sc_Z(1))/sina
c       sc_W(2)=(sc_vel_dir(2)-cosa*sc_Z(2))/sina
c       sc_W(3)=(sc_vel_dir(3)-cosa*sc_Z(3))/sina  
c
c Compute the normal to the orbit plane (=pos x vel)
c       call cross(sc_W,sc_Z,orbit_normal)
c
c PT971211: this can be simplified significantly. The orbit normal is found by
c           the cross product of the sc_vel_dir and sc_Z (since these two define
c           the orbital plane. The vector in the orbital plane which is orthogonal
c           to sc_Z is then the cross product of sc_Z and orbit_normal
c           (I checked this against the numbers from the above code before changing it)

       call cross(sc_vel_dir,sc_Z,orbit_normal)
       call cross(sc_Z,orbit_normal,sc_W)

c
c** Compute the beta angle.
       beta=90.0d0-dacos(dot(orbit_normal,earth_sun_dir))*rad2deg
c
c Now compute the (ideal) yaw angle. (couter-clockwise rotation of W around Z)
       a=dot(sc_X_t,sc_W) 
       if (a .gt. 1.d0) a=sign(1.0d0,a)
       y_att = dacos(a)*rad2deg
       y_att = -sign(y_att,beta)

     
       return
       end
