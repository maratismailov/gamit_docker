CTITLE EARTH_TO_INERT

      subroutine earth_to_inert( epoch, in_xyz, out_xyz, 
     .                                  in_sys, out_sys )

      implicit none 
      
*     Routine to compute Earth coordinates to inertial coordinates
*     and visa-versa depending on in and out sys.
*     This first version is just a rough approximation taking into
*     account nominal Earth Spin only 

      include '../includes/const_param.h'
      include '../includes/cfile_def.h'
      
* PASSED VARIABLES

* epoch  - Time either MJD or JD
* in_xyz(3) - Input Earth fixed coordinates
* out_xyz(3) - Output inertial coordinates

      real*8 epoch, in_xyz(3), out_xyz(3)

* in_sys and out_sys -  Input and output systems

      character*4 in_sys, out_sys 
      
* LOCAL VARIABLES

* JD  - Full Julian date for computing angles
* Spin(3,3)  - Spin matrix from Earth fixed to Inertial (-GAST)
* gast       - Greenwich apparent sidereal time
* xrot(3,3)  - Xpolar motion rotation
* yrot(3,3)  - Ypolar motion rotation
* etoi(3,3)  - Earth to Iertial rotation matrix
* xyrot(3,3) - Product of X and Y polar motion 
* fract      - Fraction of day

      real*8 MJD, Spin(3,3), gast, xrot(3,3), yrot(3,3), etoi(3,3),
     .       xyrot(3,3), fract
      
* i,j        - Loop counters

      integer*4 i,j 

***** Check if in and out are the same
      if( in_sys(1:1).eq.out_sys(1:1) ) then
          do i = 1,3
             out_xyz(i) = in_xyz(i)
          end do
          RETURN
       end if      
      
****  check epoch to see if JD
      fract = epoch - int(epoch) + cf_ut1(1)/86400.d0 
      mjd   = int(epoch) + (-94554.d0 + 142350.d0) 
      
***** Get GST
      call gst_mjd( mjd ,  fract, gast )
C     write(*,100) mjd, fract, gast
C100  format('MJD/GAST: ',f16.5,F10.7, f10.6)

****  Construst X and Y wobble matrices
      call make_rot( 2, cf_xp(1)/3600.d0*pi/180.d0, xrot)
      call make_rot( 1, cf_yp(1)/3600.d0*pi/180.d0, yrot)
      call mult_rot( xrot, yrot, xyrot)
            
****  Construct the Spin matrix
      call make_rot( 3, -gast, spin)
      call mult_rot( spin, xyrot, etoi )
      
      
****  Now do the transformation
      do i = 1, 3
         out_xyz(i) = 0.d0
         if( out_sys(1:1).eq.'I' ) then
            do j = 1,3
               out_xyz(i) = out_xyz(i) + etoi(i,j)*in_xyz(j)
            end do
         else if ( out_sys(1:1).eq.'E' ) then

*           Use transpose of matrix etoi
            do j = 1,3
               out_xyz(i) = out_xyz(i) + etoi(j,i)*in_xyz(j)
            end do
         else
            write(*,200) in_sys(1:1), out_sys(1:1)
 200        format('EARTH_TO_INERT: Uknown output system ',2a2)
            stop 'EARTH_TO_INERT'
         end if
      end do
      
****  Thats all
      return
      end
           
CTITLE make_rot

      subroutine make_rot( type, angle, rot )

      implicit none 
      
*     Routine to make primary rotation matrices


* type  - Type of rotation  1 = x-axis, 2 = y-axis 3 = z-axis

      integer*4 type

* angle - Angle for rotation (radians)
* rot(3,3) - Rotation matrix
      
      real*8 angle, rot(3,3)
      
****  Test each of the basic types
      if( type.eq.1 ) then
          rot(1,1) = 1
          rot(1,2) = 0
          rot(1,3) = 0
      
          rot(2,1) = 0
          rot(2,2) =  cos(angle)
          rot(2,3) =  sin(angle)
     
          rot(3,1) = 0
          rot(3,2) = -sin(angle)
          rot(3,3) =  cos(angle)
           
      else if( type.eq.2 ) then
          rot(1,1) =  cos(angle)
          rot(1,2) = 0
          rot(1,3) = -sin(angle)
      
          rot(2,1) = 0
          rot(2,2) = 1
          rot(2,3) = 0
     
          rot(3,1) =  sin(angle)
          rot(3,2) =  0
          rot(3,3) =  cos(angle) 
           
      else if( type.eq.3 ) then
          rot(1,1) =  cos(angle)
          rot(1,2) =  sin(angle)
          rot(1,3) = 0
      
          rot(2,1) = -sin(angle)
          rot(2,2) =  cos(angle)
          rot(2,3) = 0
     
          rot(3,1) = 0
          rot(3,2) = 0
          rot(3,3) = 1
      else
          write(*,*) 'Unknown rotation type in make_rot',
     .               ' Type is ',type
      end if
      
***** Thats all
      return
      end
      
CTITLE MULT_ROT

      subroutine mult_rot( A, B, C )

      implicit none 
      
*     Routine to multiple C = A*B

      real*8 a(3,3), b(3,3), c(3,3)
      
      integer*4 i,j,k
      
      do i = 1, 3
         do j = 1,3
             c(i,j) = 0.d0
             do k = 1,3
                 c(i,j) = c(i,j) + a(i,k)*b(k,j)
             end do
         end do
      end do
      
*     Thats all
      return
      end
             
CTITLE GST_MJD
 
      SUBROUTINE GST_mJD(mjd, fract, GST)

      implicit none 
C
C     T.HERRING                                 4 MARCH 1981
C
C     SUBROUTINE TO COMPUTE GST (rads) GIVEN Full julian date
C
C# LAST COMPC'ED  810311:14:23                  #
C
 
      include '../includes/const_param.h'
C
 
      real*8 fract, gst
 
*       diurnv       - Diurnal sidereal velocity
*       cent         - Number of centuries since j2000.
*       mJd_0hr       - Julian date a zero hours UT
 
      real*8 mjd, t, t_0hr, gstd, diurnv, cent, mjd_0hr
 
*     Remove the fractional part of the julian date
*     Get mjd at 0:00 UT
      mjd_0hr = mjd
*                         ! Days since J2000.0
      t_0hr = mjd_0hr - 51545.d0 
*                         ! 0:00 hrs at start of day
      cent = t_0hr / 36525.d0
 
      diurnv = ( 1.002737909350795d0 + 5.9006d-11*cent
     .                               - 5.9d-15*cent**2 )
C
C**** COMPUTE GST
      gstd = ( 24110.54841d0  + 8640184.812866d0*cent
     .                        + 0.093104d0*cent**2
*                                                            ! Cycles
     .                        - 6.2d-6*cent**3 ) /86400.d0
 
      gstd = mod(gstd,1.d0)
 
*                                             ! Rads
      gst = (gstd + diurnv*fract) * 2*pi
 
***** Thats all
      return
      end
 
