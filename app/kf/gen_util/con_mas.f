ctitle
 
      subroutine con_mas( val, pos_mas)

      implicit none 
 
 
c
c     routine to convert deg, min, sec into milli-arc-sec
c
c Variables
c ---------
c val  -- real*8 three element array containing deg, min, sec. The
c     sign is assumed to be on the first non-zero element
c pos_mas -- val converted to milli-arc-seconds
c signum -- the sign of the position
c
      real*8 val(3), pos_mas
 
c
      real*8 signum
 
c
c.... get the sign of the position
*                               ! get sign from first value
      if( val(1).ne.0 ) then
         signum = sign(1.d0, val(1))
c
*                               ! try next entry
      else
*                               ! get sign from second value
         if( val(2).ne.0 ) then
            signum = sign(1.d0, val(2))
*                               ! get sign from last value
         else
            signum = sign(1.d0, val(3))
         end if
      end if
c
c.... convert position to mas
      pos_mas = signum* ( abs(val(1))*3600.d0 +
     .                    abs(val(2))*  60.d0 +
     .                    abs(val(3))          ) * 1.d3
c
      return
      end
 
c.........................................................................
