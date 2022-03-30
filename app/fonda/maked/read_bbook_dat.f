      subroutine read_bbook_dat(ifil1,ifil2,mode,type)
c
c     read BLUE BOOK data file and create FONDA format data file
c     mode = 1: distinguesh network
c     mode = 2: merge same site name from different network
c     type = 1: mark to mark distance
c     type = 2: Clarke-arc distance
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*8 name1,name2
      character*30 note
      character*120 line
      character*1  idir1,idir2
      integer iyear,iday,month,mode,type
      integer*4 julday
      integer i,j,k,it,ioer,ioei
      integer ifil1,ifil2,ii,lcode,isit1,isit2,iod,iom,ios
      integer nblen
      logical find1,find2
c
      small = 1.0d-5
c     
c     First, sort out all working sites and experiment number.
      iexp = 0
      do i = 1,31
         igroup(i) = 0
      enddo
      rewind (ifil1)
      do 20 i = 1,30000
         read (ifil1,'(a120)', end=50, err=20)
     .     line
         read (line,'(7x,i2)',err=20) it
c        pick up observations of directions
         if (it.eq.20.or.it.eq.22) igroup(3) = igroup(3)+1
c        pick up observations of azimuth
         if (it.eq.60) igroup(1) = igroup(1)+1
c        pick up observations of baseline length
         if (it.eq.52) igroup(4) = igroup(4)+1
c        pick up observations of deflection
         if (it.eq.85) igroup(31) = igroup(31)+1
c        pick up observations of geoid height (not used for the time being)
c         if (it.eq.84) igroup(19) = igroup(19)+1
 20   continue

c     get total experiment number
 50   do i = 1,31
         if (igroup(i).gt.0) iexp = iexp+1
         if (igroup(i).gt.0) itype(iexp) = i
      enddo
c
c     network site information
      note = '{network site number          '
      write (ifil2,30) nnet,note
      write (ifil2,'(a8)') (sname(i),i=1,nnet)
c
c     write down experiment number 
      note = '{experiment number            '
      write (ifil2,30) iexp,note

      note = '{exp. index, obs. type, number'
c     second, shift the format one by one
      do 90 i = 1,iexp
         rewind (ifil1)
         it = itype(i)
         write (ifil2,40) i,it,igroup(it),note
c
c        various formats depend on observation types
         call blank(infmt)
         if (it.eq.1)
     .   infmt(1:51) = 
     .     '(7x,i2,1x,i3,24x,i4,i2,i2,5x,i3,10x,i3,i2,i3,6x,i3)'
         if (it.eq.3)
     .   infmt(1:51) = 
     .     '(7x,i2,1x,i3,24x,i4,i2,i2,5x,i3,10x,i3,i2,i4,i4,i4)'
         if (it.eq.4)
     .   infmt(1:51) = 
     .     '(7x,i2,1x,i3,21x,i2,i2,i2,5x,i3,14x,i6,i4,1x,i3,i4)'
         if (it.eq.31)
     .   infmt(1:32) = 
     .     '(7x,i2,1x,i3,50x,i4,a1,4x,i4,a1)'

         do 60 j = 1,30000
            read (ifil1,'(a120)', end=90, err=60) line
            read (line,'(7x,i2)',err=60) ii
            if (ii.ne.60.and.ii.ne.20.and.ii.ne.22.and.ii.ne.52
     .         .and.ii.ne.85) goto 60
            ioei = 0
            ioer = 0
            lcode = nblen(line)
            if (it.eq.1) read (line,fmt=infmt,end=90,err=60) 
     .         ii,isit1,iyear,month,iday,isit2,iod,iom,ios,ioer
            if (ii.eq.20) read (line,fmt=infmt,end=90,err=60) 
     .         ii,isit1,iyear,month,iday,isit2,iod,iom,ios,ioei,ioer
            if (ii.eq.22) read (line,'(7x,i2,1x,i3,37x,i3,10x,
     .         i3,i2,i4,i4,i4)',end=90,err=60) 
     .         ii,isit1,isit2,iod,iom,ios,ioei,ioer
            if (it.eq.4) read (line,fmt=infmt,end=90,err=60) 
     .         ii,isit1,iyear,month,iday,isit2,iod,ios,ioei,ioer
            if (it.eq.31) read (line,fmt=infmt,end=90,err=60) 
     .         ii,isit1,iod,idir1,ios,idir2
c           skip type-mismatch data
            if (it.eq.1.and.ii.ne.60) goto 60 
            if (it.eq.3.and.(ii.ne.20.and.ii.ne.22)) goto 60 
            if (it.eq.4.and.ii.ne.52) goto 60 
            if (it.eq.31.and.ii.ne.85) goto 60 
c
            if (ioer.le.0) ioer = ioei
c           avoid JULDAY stop
            if (month.le.0) then
               print*,' Month is missed for record ',j
               month = 6
            endif   
            if (iday.le.0) then
               print*,' Day is missed for record ',j
               iday = 1 
            endif
c           if the site name mismatch, reject the data
               find1 = .false.
               find2 = .false.
               do k = 1,nnet
                  if (isit1.eq.iesit(k)) name1 = sname(k)
                  if (isit1.eq.iesit(k)) find1 = .true.
                  if (it.ge.2.and.it.le.10) then
                     if (isit2.eq.iesit(k)) name2 = sname(k)
                     if (isit2.eq.iesit(k)) find2 = .true.
                  endif
               enddo
               if (.not.find1) goto 60
               if (it.ge.2.and.it.le.10.and..not.find2) goto 60

c           Compute Errors and Output data

c           Azimuth
            if (it.eq.1) then
               sec = dble(ios)*1.0d-2
               er  = dble(ioer)*1.0d-2
               time1 = 1900.0d0+julday(month,iday,iyear,1)/365.2422d0
               write (ifil2,110) 
     *            time1,iod,iom,sec,er,name1
            endif

c           Direction
            if (it.eq.3) then
               sec = dble(ios)*1.0d-2
               er  = dsqrt(dble(ioer)**2+dble(ioei)**2)*1.0d-2
               if (ii.eq.20) then
                  time1 = 1900.0d0+
     .                    julday(month,iday,iyear,1)/365.2422d0
                  time0 = time1
               else
                  time1 = time0
               endif
               write (ifil2,110) 
     *            time1,iod,iom,sec,er,name1,name2
            endif

c           Reduced Distance
            if (it.eq.4) then
               blen = dble(iod)+dble(ios)*1.0d-4
c              units need to be metres
c              error model:  sigma**2 = a**2 + (b*L)**2
               er  = dsqrt((dble(ioer)*(1.0d-7)*blen)**2+
     .         (dble(ioei)*1.0d-4)**2)*1.0d+3
               if (iyear.lt.1900) iyear = iyear+1900
               time1 = 1900.0d0+julday(month,iday,iyear,1)/365.2422d0
               write (ifil2,150) 
     *            time1,blen,er,name1,name2
            endif

c           Deflection
            if (it.eq.31) then
               defe = dble(ios)*1.0d-2
               defn = dble(iod)*1.0d-2
               if (idir1.eq.'S'.or.idir1.eq.'s') defn = -defn
               if (idir2.eq.'W'.or.idir2.eq.'w') defe = -defe
               time1 = rtime
               write (ifil2,140) 
     *            time1,name1,defe,defn
            endif
 60      continue

 90   continue

c
 30   format (i5,25x,a30)
 40   format (1x,3i5,14x,a30)
 110  format (1x,f10.4,2i5,f14.8,f9.4,2(2x,a8))
 140  format (1x,f10.4,3x,a8,2(2x,f8.3))
 150  format (1x,f10.4,f14.5,f14.5,2(2x,a8))
c
 120  continue
      return
      end
