      subroutine allan (f0,xx,tt,aa,nn,nper)

c     given a series of phase measurements, calculate
c     the Allan Standard Deviation


c     input:
c        nominal frequency in Hz
         real*8 f0
c        phase (ignored if less than -1.d99)
         real*8 xx(*)
c        time in seconds
c        assumed correct, even if xx is not good
         real*8 tt(*)
c        number of points
         integer*4 nn
c     output
c        allan standard deviation (dimensionless)
         real*8 aa(*)
c        number of different periods
c        return zero in case of error
         integer nper
c     local
         real*8 dv1,dv2
         integer ngood,i,j,k
         real*8 allansd
         real*8 period
         real*8 sum,dt1,dt2

      nper = 0

c     write (6,*) 'ALLAN: nn ',nn

c     deduce sampling period
      period = tt(2) - tt(1)
      if (period .lt. 0.0d0) then
         print *,'ALLAN: bad period ',period
         return
      endif

c     OK, we want to make measurements out to a sampling period
c     of span/2, or nn/2 * dt.  k is sampling period
      do k = 1,nn/2
         ngood = 0
         sum = 0.d0
c        Take k trips through the data
c        We should be able to interleave k times
         do j = 1,k
c           if we have 3 good points, then take two derivatives
            do i=j,nn-k-k,k
               if (xx(i)     .gt. -1.d99 .and.
     .             xx(i+k)   .gt. -1.d99 .and.
     .             xx(i+k+k) .gt. -1.d99) then
                  dt1 = tt(i+k)   - tt(i)
                  dt2 = tt(i+k+k) - tt(i+k)
                  dv1 = (xx(i+k)  - xx(i)  - f0*dt1)/(f0*dt1)
                  dv2 = (xx(i+k+k)- xx(i+k)- f0*dt2)/(f0*dt2)
                  sum = sum + (dv2 - dv1)**2
                  ngood = ngood + 1
                endif
            enddo
         enddo

         if (ngood .gt. 2) then
            nper = nper + 1
            allansd = dsqrt(sum/(2.d0*ngood))
c           write (6,'(1x,i4,1x,2(1pe16.4))')
c    .      ngood,nper*period,allansd
            aa(nper) = allansd
         else
c           ran out of points, stop here
            return
         endif
      enddo

      return

      end

