CTITLE SWITCH_EMA8
 
      subroutine switch_ema8(a,b,temp)
 
      implicit none
 
c
c     Switches ema real*8 values
c     Only A and B are real*8
 
 
      real*8 a, b, temp
 
c
      temp = a
      a    = b
      b    = temp
c
      end
 
