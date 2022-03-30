CTITLE PMU_MAIN_PART
 
      subroutine pmu_main_part( site, deriv, pos, comp, ut1_val )
 
 
*     Routine to compute the pmu partials with respect to a site
*     component.  For this routine pos is in main memory.
 
      include '../includes/const_param.h'
 
*   site    - Global site number
*   comp    - Component of position (1=x, 2=y, 3=z)
 
      integer*4 site, comp
 
*   deriv(3)    - Derivative with respect to X,Y wobble and UT1
*   pos(3)      - Position of the site
*   ut1_val     - Value of UT1-AT (needed to alling the coordinate
*               - systems approximately (mas))
 
 
      real*8 deriv(3), pos(3), ut1_val
 
****  Use a GOTO for the component
 
      goto ( 100, 200, 300 ) comp
 
*     X-site coordinate
  100 continue
          deriv(1) = -pos(3)/rad_to_mas
          deriv(2) = 0.d0
          deriv(3) = (-pos(2) + ut1_val/rad_to_mas*pos(1))/rad_to_mas
          return
 
 
*     Y_site coordinate
  200 continue
          deriv(1) = 0.d0
          deriv(2) = pos(3)/rad_to_mas
          deriv(3) = (pos(1) + ut1_val/rad_to_mas*pos(2))/rad_to_mas
          return
 
*     Z_site coordinate
  300 continue
          deriv(1) = ( pos(1) - ut1_val/rad_to_mas*pos(2))/rad_to_mas
          deriv(2) = (-pos(2) - ut1_val/rad_to_mas*pos(1))/rad_to_mas
          deriv(3) =  0.d0
          return
 
****  Thats all
      end
 
