      subroutine read_usgs_new_dat(ifil1,ifil2,nnet,mode,infmt)
c
c     read USGS new data file and create FONDA format data file
c     error model:  sigma**2 = a**2 + (b*L)**2
c     a = 3mm;   b = 0.2ppm
c     (Savage & Prescott 1973, JGR, 78, pp.6001-6008)
c     mode = 1: distinguesh network
c     mode = 2: merge same site name from different network
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*4 netid(120),netido,netid1
      character*8 name1,name2,sitnam(120)
      character*30 note
      character*50 infmt
      character*256 string
      integer iyear,iday,month,id1,id2
      integer*4 julday,id
      logical old1,old2
c
      small = 1.0d-5
c
c     default format
      if (infmt(1:1).eq.'*')
     . infmt(1:33) = '(a4,2x,a8,1x,a8,1x,3i2,12x,f11.4)'
c     
c     First, sort out all working sites and experiment number.
      nsit = 0
      iexp = 0
      netido = 'xxxx'
      do i = 1,nnet+2
         read (ifil1,'(2x)') 
      enddo
      do 20 i = 1,10000
         read (ifil1,fmt=infmt,err=20,end=50) 
     .      netid1,name1,name2,year,obs
c        skip unrealistic data
         if (obs.lt.small) goto 20
         if (nsit.eq.0.or.netid1.ne.netido) then
            nsit = nsit+1
            sitnam(nsit) = name1
            netid(nsit) = netid1
            nsit = nsit+1
            sitnam(nsit) = name2
            netid(nsit) = netid1
            netido = netid1
            i1 = iyear
            i2 = month
            i3 = iday
            iexp = iexp+1
            goto 20
         endif 
c        check if it is a new site
         old1 = .false.
         old2 = .false.
         do 40 j = 1,nsit
            if (name1.eq.sitnam(j).and.netid1.eq.netid(j))
     .         old1 = .true.
            if (name2.eq.sitnam(j).and.netid1.eq.netid(j))
     .         old2 = .true.
            if (old1.and.old2) goto 60
 40      continue
         if (.not.old1) then
            nsit = nsit+1
            sitnam(nsit) = name1
            netid(nsit) = netid1
         endif
         if (.not.old2) then
            nsit = nsit+1
            sitnam(nsit) = name2
            netid(nsit) = netid1
         endif
c        
 60      iexp = iexp+1
 20   continue
 50   print*,'nsit,iexp=',nsit,iexp
      if (nsit.lt.1) goto 120
      jobs = iexp
c
c     site id
      note = '{network site number          '
      write (ifil2,30) nsit,note
      write (ifil2,'(a4,2x,a8)') (netid(i),sitnam(i),i=1,nsit)
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
     .      netid1,name1,name2,iyear,month,iday,obs
c        skip unrealistic data
         if (obs.lt.small) goto 85
         err = dsqrt(aterm**2+(bterm*obs)**2)
         iyear = iyear+1900
         id = julday(month,iday,iyear,3)
c        identify two sites
         old1 = .false.
         old2 = .false.
         do 70 j = 1,nsit
            if (mode.le.1) then
               if (name1.eq.sitnam(j).and.netid1.eq.netid(j)) then
                  old1 = .true.
                  id1 = j
               endif
               if (name2.eq.sitnam(j).and.netid1.eq.netid(j)) then
                  old2 = .true.
                  id2 = j
               endif
            endif
            if (mode.eq.2) then
               if (name1.eq.sitnam(j)) then
                  old1 = .true.
                  id1 = j
               endif
               if (name2.eq.sitnam(j)) then
                  old2 = .true.
                  id2 = j
               endif
            endif
            if (old1.and.old2) goto 100
 70      continue

 100     write (ifil2,150) iyear,id,obs*1.0d-3,err,id1,id2

 90   continue

c
c10   format (a4,2x,a8,1x,a8,1x,3i2,12x,f11.4)
 30   format (i5,25x,a30)
 140  format (3i5,15x,a30)
 150  format (1x,2i5,2f14.8,2i5)
c
 120  continue
      return
      end
