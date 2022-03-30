      program extrap

c     extrapolate X,Y,Z coordinates to a new epoch.

      IMPLICIT REAL*8(A-H,O-Z)

      character*80   file1
      character*8    code
      character*240  line
      integer        ios,in
      integer        rcpar
      external       rcpar

      in = 0
      file1 = 'standard input'
      ios = 0

c     first command line arg is epoch in decimal years.
      ios = rcpar(1,line)
      if (ios .eq. 0) then
         read (line,*) t0
      else
         print *, 'Need reference epoch as first argument...'
         stop
      endif

      print *,'XYZ propagated to epoch ',t0
      print *,'(a8,1x,f20.4,1x,f20.4,1x,f20.4,1x)'

      do 1010 while (ios .eq. 0)
         read (unit = 5,
     .      fmt     = '(a)',
     .      iostat  = ios,
     .      err     = 1010,
     .      end     = 1020) line


         if (line(1:1) .eq. ' ') then
            read (line(2:9),*,err=1010,iostat=ios) code
            read (line(10:nblen(line)),*,iostat=ios,err=1010)
     .            x,y,z,xd,yd,zd,t1
            x1 = x + xd * (t1 - t0)
            y1 = y + yd * (t1 - t0)
            z1 = z + zd * (t1 - t0)
            write (6,'(a8,1x,f20.4,1x,f20.4,1x,f20.4,1x)',
     .           iostat=ios,err=1010)
     .           code,x1,y1,z1
         endif

         in = in + 1

 1010 continue
      if (ios .ne. -1) then
        call ferror (ios,6)
      endif

 1020 continue

      end




