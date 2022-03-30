      real*8 function rnhalf(x)

c     round to nearest half integer
c     e.g. rnhalf (1.67) = 1.5

      real*8 x

      rnhalf = dnint(2.0d0 * x)/2.0d0

      return
      end

