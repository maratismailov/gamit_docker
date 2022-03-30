CTITLE SET_OCE_AVAIL
 
      subroutine set_oce_avail( avail, rad, horz )

      implicit none 
 
*     This routine will check the magnitude ofthe radial and
*     horizontal ocean loading contributions and see if
*     they are non-zero.  If they are it will set the avail bit
*     for the ocean loading.
 
*   avail   - Availability flag.  Bit 5 is radial, Bit 6 is
*           -                    horizontal.
 
      integer*4 avail
 
*   rad      - Radial contribution (delay only)
*   horz    - Horizontal contribution (delay only)
 
 
      real*4 rad, horz
 
****  Check the magnitude and set bit if non-zero
      if( abs(rad) .gt.0.d0 ) call sbit( avail, 5, 1)
      if( abs(horz).gt.0.d0 ) call sbit( avail, 6, 1)
 
****  Thats all
      return
      end
 
