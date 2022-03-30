      real*8 function frac(x)

c     return the fractional part of a number

      real*8 x

      frac = x - dint(x)

      return
      end
