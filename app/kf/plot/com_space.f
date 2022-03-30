CTITLE    ..................................................................
 
      subroutine com_space(vmin, vmax, tic_space)
c
c     routine to compute the 'optimum' spacing bewteen the tic marks
c     This translates as between 6 and 10 tic marks per axis.
c
c Variables
c ---------
c vmin, vmax -- the minimum and maximum values on the axis
c tic_space  -- the computed tic spacing in world coordinates
c
      real*4 vmin, vmax, tic_space
 
c
c Local variables
c ---------------
c span   -- the interval spanned by the data
c power  -- the power of the span
c iexp   -- power or power - 1 if power is less than zero
c ival   -- the number of ticmarks
c
      real*4 span, power
 
c
      integer*4 iexp, ival

*     Check size of vmin and vmax
*                                                     ! No scale are set, set unit values
      if( abs(vmin).gt.1.d20 .or. abs(vmax).gt.1.d20 ) then
          vmin = -1.0
          vmax = 1.0
      end if
 
c
c.... Compute span and number of decades spanned
      span = abs(vmax-vmin)
c
c.... See if span is zero
      if( span.eq.0 ) then
*                           ! convenient number to stop log error
          span = 1
      end if
c
      power = log10(span)
      if( power.lt.0 ) then
         iexp = power - 1
      else
         iexp = power
      end if
c
c.... Compute number of tics
      ival = span/(10.**iexp) + 0.5
      if( ival.eq.3 .or. ival.eq.4 ) ival = 5
      if( ival.gt.5 ) ival = 10
c
      tic_space = ival*10.**(iexp-1)
c
      return
      end
 
