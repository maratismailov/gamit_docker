      program allan

c     Take the 2-point Allan Variance
c     of two columns of numbers

      IMPLICIT REAL*8(A-H,O-Z)
      parameter (maxbsl = 5000)

c     read a 132 char line
      character*80   file1,file2

      dimension time(maxbsl),value(maxbsl),deriv(maxbsl)


c     read in time and phase
      in = 0
      file1 = 'standard input'
      do i = 1,maxbsl

         read (unit = 5,
     .      fmt     = *,
     .      iostat  = ios,
     .      err     = 1010,
     .      end     = 1020) time(i),value(i)

c           convert time from minutes to seconds
            time(i) = time(i) * 60.0d0

         in = in + 1
      enddo


 1010 continue
      call ferror (ios,6)

 1020 continue

c     form derivative, assuming that this is L1
      freql1 = 1.57542D9
      do 200 i = 1,in - 1
         dt = time(i+1)-time(i)
c         deriv(i) = ((value(i+1) - value(i)) - freql1*dt)/(freql1*dt)

         deriv(i) = (value(i+1) - value(i))/(freql1*dt)
         write (6,*) time(i),value(i),deriv(i)
 200  continue

c     sum of squares
      sum = 0.d0
      do 250 i = 1,in - 2
         sum = sum + (deriv(i+1)-deriv(i))**2
250   continue

      if (in .gt. 2) then
         alsd = dsqrt(sum/(2*(in -1)))
      endif

      write (6,700) 'N        ',in
 700  format (1x,a20,1x,i12)
      write (6,800) 'Allan SD ',alsd
 800  format (1x,a20,1x,1pe12.2)

      end




