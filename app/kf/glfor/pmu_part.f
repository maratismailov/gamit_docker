CTITLE PMU_PART
 
      subroutine pmu_part( site, deriv, pos, comp, ut1_val )

      implicit none  
 
*     Routine to compute the pmu partials with respect to a site
*     component
 
      include '../includes/const_param.h'
 
*   site    - Global site number
*   comp    - Component of position (1=x, 2=y, 3=z)
 
      integer*4 site, comp
 
*   deriv(3)    - Derivative with respect to X,Y wobble and UT1
*   pos(3)      - Position of the site
*   ut1_val     - Value of UT1-AT (needed to alling the coordinate
*               - systems approximately (mas))
 
 
      real*8 deriv(3), pos(3), ut1_val

* dut        - Angle change in radians due to UT1 (rads)
* xy_rot(2)  - X and Y coordinates are rotation by UT1-AT

      real*8 dut, xy_rot(2)

****  First get the change in XY coordinates for the UT1-AT

      dut = ut1_val / rad_to_mas
      xy_rot(1) =  pos(1)*cos(dut) + pos(2)*sin(dut)
      xy_rot(2) = -pos(1)*sin(dut) + pos(2)*cos(dut)
 
****  Use a GOTO for the component
 
      goto ( 100, 200, 300 ) comp
 
*     X-site coordinate
  100 continue
          deriv(1) = -pos(3)/rad_to_mas
          deriv(2) = 0.d0
*         deriv(3) = (-pos(2) + ut1_val/rad_to_mas*pos(1))/rad_to_mas
          deriv(3) = -xy_rot(2)/rad_to_mas
          return
 
 
*     Y_site coordinate
  200 continue
          deriv(1) = 0.d0
          deriv(2) = pos(3)/rad_to_mas
*         deriv(3) = (pos(1) + ut1_val/rad_to_mas*pos(2))/rad_to_mas
          deriv(3) = xy_rot(1)/rad_to_mas
          return
 
*     Z_site coordinate
  300 continue
*         deriv(1) = ( pos(1) - ut1_val/rad_to_mas*pos(2))/rad_to_mas
*         deriv(2) = (-pos(2) - ut1_val/rad_to_mas*pos(1))/rad_to_mas
          deriv(1) =  xy_rot(1)/rad_to_mas
          deriv(2) = -xy_rot(2)/rad_to_mas
          deriv(3) =  0.d0
          return
 
****  Thats all
      end
 
