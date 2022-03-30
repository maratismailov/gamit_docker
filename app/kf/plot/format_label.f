CTITLE    ...................................................................
 
      subroutine format_label(vmin, vmax, ref_val, labform)
c
c     Routine to evaluate and appripriate format format for the labels
c
c Variables
c ---------
c vmin, vmax -- the mimumin and maximum values of the axis
c ref_val -- reference value removed from data to avoid rounding error
c labform -- the  format most appropriate
c
      real*4 vmin, vmax
 
c
      real*8 ref_val
 
c
      character*(*) labform
 
c
c Local variables
c ---------------
c ndigitb -- number of digits before decimal point
c ndigita -- number of digits after deciminal point
c edigit  -- number of deciminal places needed for E format
c value -- maximum value  of label
c range -- the range of the scales
c
      integer*4 ndigitb, ndigita, edigit
 
c
      real*8 value, range
 
c
c.... Get maximum label value
      value = max(abs(vmin+ref_val), abs(vmax+ref_val))
*                                   ! get range for spacing between tics
      range = abs(vmax-vmin)/10.0
c
c.... See if we have value numbers
*                                   ! stops log error
      if( value.eq.0 ) value = 0.1
*                                   ! stops log error
      if( range.eq.0 ) range = 0.1
c
      if( value.gt.1 ) then
         ndigitb = log10(value)+1.
         if( range.gt.1 ) then
            ndigita = 0
         else
            ndigita = -(log10(range)-1.)
         end if
      else
         ndigitb = -(log10(value)-2.)
         if( range.gt.1 ) then
            ndigita = 0
         else
            ndigita = -(log10(range)-1.)
         end if
      end if
*
*     Do a final check on format
      if( ndigita.lt.0 ) ndigita = 5
      if( ndigitb.lt.0 ) ndigitb = 5
c
*                                                ! use F format
      if( ndigitb+ndigita.le.10 ) then
         write(labform,'( "(F",i2,".",i2,")" )')
     .                ndigitb+ndigita+2,ndigita
*                                                ! use E format
      else
         ndigita = log10(range)
         ndigitb = log10(value)
         if( ndigita.gt.20 ) ndigita = 20
         if( ndigitb.gt.20 ) ndigitb = 20
 
         edigit = abs(ndigitb-ndigita) + 1
         write(labform,'( "(E",i2,".",i2,")" )')
     .                edigit+7, edigit
      end if
c
      return
      end
 
