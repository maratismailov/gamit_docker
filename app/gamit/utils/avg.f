      program avg

c     average a column of numbers

      IMPLICIT REAL*8(A-H,O-Z)
      parameter (maxbsl = 5000)

c     read a 132 char line
      character*80   file1,file2

      dimension value(maxbsl)


      in = 0
      file1 = 'standard input'
      do i = 1,maxbsl

         read (unit = 5,
     .      fmt     = *,
     .      iostat  = ios,
     .      err     = 1010,
     .      end     = 1020) value(i)

         in = in + 1
      enddo


 1010 continue
      call ferror (ios,6)

 1020 continue

      sum = 0.d0
      do 200 i = 1,in
         sum = sum + value(i)
 200  continue

      if (in .gt. 0) then
         ravg = sum/in
      endif

c     compute sample variance
      devsum = 0.d0
      do 210 i = 1,in
         devsum = devsum + (value(i) - ravg)**2
 210  continue

c     make standard deviation
      if (in .gt. 1) then
         sigma = dsqrt(devsum/(in-1))
      endif

      write (6,700) 'N      ',in
 700  format (1x,a20,1x,i20)
      write (6,800) 'Sum    ',sum
      write (6,800) 'Average',ravg
      write (6,800) 'S.D.   ',sigma
 800  format (1x,a20,1x,1pe20.13)

      end




