      real*8 function rnqrtr(x)

c     round to nearest quarter integer
c     e.g. rnqrtr (1.71) = 1.75

      real*8 x

      rnqrtr = dnint(4.0d0 * x)/4.0d0

      return
      end

