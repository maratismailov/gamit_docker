CTITLE 'MAS_TO_DMS'
 
      subroutine mas_to_dms( value, deg,min,sec)

      implicit none 
 
 
c
c     routine to convert value in milliarcseconds to deg, min,sec
c
c Variables
c ---------
c value -- angle to be converted from mas to deg,min,sec
c deg, min, sec -- the converted values
c
      real*8 value, deg, min, sec
 
c
c Local variables
c ---------------
c abs_value -- absolute value of value
c sgn -- sign of value
c ideg -- integer deg
c imin -- integer min
c
      real*8 abs_value
 
c
      integer*4 sgn, ideg, imin
 
c
c.... Save the sign of value
      if( value.lt.0 ) then
         sgn = -1
      else
         sgn = 1
      end if
c
      abs_value = abs(value)
c
c.... Divide value by 3600.d3 to get deg
      ideg = abs_value / 3600.d3
      imin = (abs_value - ideg*3600.d3)/60.d3
      sec  = (abs_value - ideg*3600.d3 - imin*60.d3)/1.d3
c
c.... Save in return variables
      deg = sgn*(ideg + 1.d-6)
      min = imin
****  make sure that deg is not zero (and therefore will not carry
*     the sign.  (Because of this min should be writen with I3 
*     format).
      if( ideg.eq.0 ) then
          min = sgn*imin
      end if
c
      return
      end
 
 
