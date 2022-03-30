      program inter

c  test program to check that my linear interpolation works!
c  P. Tregoning
c  8 January 2004

      implicit none

      integer*4 a,b,ncol,nrow,i,j
      parameter (a=3,b=10)
 
      real*4 times(b),values(a,b)
      real*8 obstime,interpolated(a)

      ncol = 2
      nrow = 5

c  put some values in the array to be interpolated
      do i=1,nrow
        times(i) = i*1.0
        do j=1,3
          values(j,i) = i*2.0 + float(i**2) + (j-1)*10.0
        enddo
      enddo

c  print out the values matrix
      print*,'input values'
      do i=1,nrow
        print*,(values(j,i),j=1,ncol)
      enddo      

c  now call lininterp and see what comes back
      obstime = 4.4
      print*,'calling lininterp for time',obstime
      call lininterp(a,b,times,values,obstime,ncol,nrow,interpolated)
      print*,'obstime and interp',obstime,(interpolated(i),i=1,ncol)

      end

