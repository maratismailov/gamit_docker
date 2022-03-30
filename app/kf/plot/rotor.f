CTITLE ROTOR
 
      subroutine rotor
 
 
*     This routine will "rotate" the orientation, line feed and
*     justification for current values of sign_x and sign_y.
*     The interpretaion of sign_x and sign_y are:
*     sign_x  sign_y  Rotation
*         1     1      0.0  deg
*        -1     1    -90.0  deg
*         1    -1    +90.0  deg
*        -1    -1    180.0  deg
*
*     The justification is switched depending on the values of sign_x
*     sign_y
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   i           - Loop counter
 
      integer*4 i
 
*   ct, st      - Cos and sin of theta set to nearest integer values
*               - (used to avoid rounding error)
*   rot(2,2)    - The rotation matrix to be used
*   temc(2)     - Temporary vector used during rotation of orientation
*   teml(2)     - Temporary vector used during rotation of line feed
*   theta       - The rotation angle in rads.
 
 
      real*4 ct, st, rot(2,2), temc(2), teml(2), theta
 
***** Find out rotation angle
      theta = 0.0
*                                                             ! -pi/2
      if( sign_x.eq.-1 .and. sign_y.eq. 1 ) theta = -1.570796
*                                                             !  pi/2
      if( sign_x.eq. 1 .and. sign_y.eq.-1 ) theta =  1.570796
*                                                             !  pi
      if( sign_x.eq.-1 .and. sign_y.eq.-1 ) theta =  3.141593
 
*     Construct rotation matrix
      ct = anint( cos(theta) )
      st = anint( sin(theta) )
 
      rot(1,1) =  ct
      rot(1,2) =  st
      rot(2,1) = -st
      rot(2,2) =  ct
 
*     Now rotate orientation
 
      do i = 1,2
          temc(i) = rot(i,1)*cori(1) + rot(i,2)*cori(2)
          teml(i) = rot(i,1)*lfc(1)  + rot(i,2)*lfc(2)
      end do
 
*     Save values
      do i = 1,2
          cori(i) = temc(i)
          lfc(i)   = teml(i)
      end do
 
*     Convert character linefeed into world coordinates
 
      call comlf( lfc, lf)
 
****  Now "rotate" justification
 
      if( sign_x*sign_y.lt.0 ) then
          temc(1) = just(2)
          just(2) = just(1)
          just(1) = temc(1)
      end if
 
***** Thats all
      return
      end
 
