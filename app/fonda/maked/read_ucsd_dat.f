      subroutine read_ucsd_dat(ifil1,ifil2)
c
c     read UCSD data file and create FONDA format data file
c     error model:  sigma**2 = a**2 + (b*L)**2
c     a = 3mm;   b = 0.2ppm
c     (Savage & Prescott 1973, JGR, 78, pp.6001-6008)
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*8 name1,name2,sitnam(120)
      character*30 note
      integer iyear,iday,month,id1,id2
      integer*4 julday
      integer ifil1,ifil2,i,j,jobs,itp,igp
      logical old1,old2
c
      small = 1.0d-5
c
c     default format
      if (infmt(1:1).eq.'*')
     . infmt(1:30) = '(6x,a8,1x,a8,1x,3i2,12x,f11.4)'
c     
c     First, sort out all working sites and experiment number.
      nsit = 0
      iexp = 0
      do 20 i = 1,10000
         read (ifil1,fmt=infmt,end=50) name1,name2,iyear,month,iday,obs
c        skip unrealistic data
         if (obs.lt.small) goto 20
         if (nsit.eq.0) then
            nsit = nsit+1
            sitnam(nsit) = name1
            nsit = nsit+1
            sitnam(nsit) = name2
            iexp = iexp+1
            goto 20
         endif 
c        check if it is a new site
         old1 = .false.
         old2 = .false.
         do 40 j = 1,nsit
            if (name1.eq.sitnam(j)) old1 = .true.
            if (name2.eq.sitnam(j)) old2 = .true.
            if (old1.and.old2) goto 60
 40      continue
         if (.not.old1) then
            nsit = nsit+1
            sitnam(nsit) = name1
         endif
         if (.not.old2) then
            nsit = nsit+1
            sitnam(nsit) = name2
         endif
c        check observation time
 60      iexp = iexp+1
 20   continue
 50   print*,'nsit,iexp=',nsit,iexp
      if (nsit.lt.1) goto 120
      jobs = iexp
c
c     site id
      note = '{network site number          '
      write (ifil2,30) nsit,note
      write (ifil2,'(a8)') (sitnam(i),i=1,nsit)
c
c     experiment number (always 1)
      i = 1
      itp = 4
      note = '{experiment number            '
      write (ifil2,30) i,note
      note = '{exp. index, obs. type, number'
      write (ifil2,140) i,itp,jobs,note

c     second, shift the format one by one
      rewind (ifil1)
      igp = 1
      aterm = 3.0d0
      bterm = 0.2d-3
      do 90 i = 1,iexp
 85      read (ifil1,fmt=infmt,end=120) 
     .      name1,name2,iyear,month,iday,obs
c        skip unrealistic data
         if (obs.lt.small) goto 85
         err = dsqrt(aterm**2+(bterm*obs)**2)
         time1 = 1900.0d0+julday(month,iday,iyear,1)/365.2422d0
c        identify two sites
         old1 = .false.
         old2 = .false.
         do 70 j = 1,nsit
            if (name1.eq.sitnam(j)) then
               old1 = .true.
               id1 = j
            endif
            if (name2.eq.sitnam(j)) then
               old2 = .true.
               id2 = j
            endif
            if (old1.and.old2) goto 100
 70      continue

 100     write (ifil2,150) time1,obs,err,name1,name2

 90   continue

c
c10   format (6x,a8,1x,a8,1x,3i2,12x,f11.4)
 30   format (i5,25x,a30)
 140  format (3i5,15x,a30)
 150  format (1x,f10.4,f14.5,f14.8,2(2x,a8))
c
 120  continue
      return
      end
