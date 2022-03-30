      program histgm

c     produce a histogram


      IMPLICIT REAL*8(A-H,O-Z)
      parameter (maxbsl = 5000)

c     read an 80 char line
      character*80   line

      dimension value(maxbsl)
      dimension vhi(maxbsl/10),vlo(maxbsl/10),nn(maxbsl/10+1)



c     first command line argument is number of bins
c     if zero, assume 11 bins
      call rcpar (1,line)
      if (nblen(line) .gt. 0) then
         read (line(1:nblen(line)),*) nbins
         if (nbins .eq. 0) nbins = 11
      else
         nbins = 11
      endif

c     second command line argument is size of bins
c     if zero, assume 1/3 sigma
      call rcpar (2,line)
      if (nblen(line) .gt. 0) then
         read (line(1:nblen(line)),*) wide
      else
         wide = 0.0d0
      endif



c     third command line argument is option
c     1 for regular histogram
c     2 for plottable output
      call rcpar (3,line)
      if (nblen(line) .gt. 0) then
         read (line(1:nblen(line)),*) iopt
         if (iopt .eq. 0) iopt = 1
      else
         iopt = 1
      endif


c     initialize
      do i = 1,nbins
         nn(i) = 0
      enddo

      in = 0
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

c     sum and extrema
      sum = 0.d0
      vmax = -1.0d100
      vmin =  1.0d100
      do 200 i = 1,in
         sum = sum + value(i)
         vmax = max(value(i),vmax)
         vmin = min(value(i),vmin)
 200  continue

c     mean value
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

      if (wide .eq. 0.0d0) wide = sigma/3.d0

c     round the value to the nearest multiple
c     of the bin size
      ravg = dnint(ravg/wide)*wide

c     create the lower and upper bounds of the bins
      do i = 1,nbins
         vlo(i) = ravg + (i - 1 - nbins/2.) * wide
         vhi(i) = vlo(i) + wide
      enddo
c     first and last bins go to extreme values
      vlo(1)     = min(vmin,vlo(1))
      vhi(nbins) = max(vmax,vhi(nbins))

c     make stuff fall in the bins
      do i = 1,in
         do j = 1,nbins
            if (value(i).lt.vhi(j).and.value(i).ge.vlo(j)) then
               nn(j) = nn(j) + 1
            endif
         enddo
      enddo

c     print it out
      do i = 1,nbins
         if (iopt .eq. 2) then
            write (6,800) vlo(i),0
            write (6,800) vlo(i),nn(i)
            write (6,800) vhi(i),nn(i)
            write (6,800) vhi(i),0
         else
            write (6,850) vlo(i),vhi(i),nn(i)
         endif
      enddo

 800  format (1x,1pg10.2,1x,i4)
 850  format (2(1x,1pg12.4),1x,i4)

      end




