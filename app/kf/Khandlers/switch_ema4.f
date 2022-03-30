CTITLE SWITCH_EMA4
 
      subroutine switch_ema4(a,b,temp)

      implicit none
 
 
c
c     Switches real*4 values
c     Only A and B are EMA
 
      real*4 a, b, temp
 
c
      temp = a
      a    = b
      b    = temp
c
      end
 
