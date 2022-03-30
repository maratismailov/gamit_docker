c
      subroutine avgvar(ndata,data,amean,var)
c
c     Sequential calculate the mean and variance of a time series
c
      implicit none

      integer*4 ndata

      real*8 data,amean,var,dev,dsp,varr

c     input:
c     dsp --- increment of arc span(length)
c     data --- one data point of the arc
c
c     output: (also input)
c     amean --- updated mean value for this time series
c     var --- updated variance for this time series
c
c

      if (ndata.eq.1) then
         amean = data
         var = 0.0d0
         goto 30
      endif
      dev = data - amean
      dsp = 1.0d0/dble(ndata)
      amean = amean + dev*dsp
      varr = var**2*(1.0d0 - dsp)
      var = dsqrt(varr + (dev*(1.0d0 + dsp))**2*dsp)

 30   continue

      return
      end

