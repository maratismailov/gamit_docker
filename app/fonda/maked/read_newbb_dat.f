      subroutine read_newbb_dat(ifil1,ifil2,mode,type)
c
c     read BLUE BOOK data file and create FONDA format data file
c     mode = 1: distinguesh network
c     mode = 2: merge same site name from different network
c     type = 1: mark to mark distance (currently hardwired in 
c               call to this subroutine!)
c     type = 2: Clarke-arc distance (currently never called)
c               see implementation below.
c
      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*8 name1,name2
      character*30 note
      character*120 line
      character*1  idir1,idir2,code
      integer iyear,iday,ihh,imm,month,mode,type
      integer*4 julday
      integer i,j,it,ioer,ioei,iod,iom,ios
      integer ifil1,ifil2,isit1,isit2,ii,lcode
      integer nblen
      real*8 dtoyr,htoyr,mtoyr
      logical find1,find2
c
      small = 1.0d-5
      
c     define date/time conversion constants
      dtoyr = 365.2422d0
      htoyr = dtoyr*24.0d0
      mtoyr = htoyr*60.0d0
c     
c     First, sort out all working sites and experiment number.
      iexp = 0
      do i = 1,31
         igroup(i) = 0
      enddo
      rewind (ifil1)
      do 20 i = 1,30000
         
c        initialise variables
         ihh = 0
         imm = 0

         read (ifil1,'(a120)', end=50, err=20)
     .     line
         read (line,'(7x,i2,1x,i4)',err=20) it,isit1
c        pick up observations of directions
         if (it.eq.20.or.it.eq.22) then
            read (line,'(50x,i4)',err=20) isit2
            call chk_sname(3,isit1,isit2,find1,find2,name1,name2)
            if (.not.find1.or..not.find2) goto 20
            igroup(3) = igroup(3)+1
         endif
c        pick up observations of azimuth
         if (it.eq.60) then
            read (line,'(50x,i4)',err=20) isit2
            call chk_sname(1,isit1,isit2,find1,find2,name1,name2)
            if (.not.find1.or..not.find2) goto 20
            igroup(1) = igroup(1)+1
         endif
c        pick up observations of baseline length
         if (it.eq.52) then
            read (line,'(45x,i4)',err=20) isit2
            call chk_sname(4,isit1,isit2,find1,find2,name1,name2)
            if (.not.find1.or..not.find2) goto 20
            igroup(4) = igroup(4)+1
         endif
c        pick up observations of deflection
         if (it.eq.85) then
            call chk_sname(31,isit1,isit2,find1,find2,name1,name2)
            if (.not.find1) goto 20
            igroup(31) = igroup(31)+1
         endif
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
     .   infmt(1:56) = 
     .     '(7x,i2,1x,i4,25x,i2,i2,i2,i2,i2,1x,i4,9x,i3,i2,i3,6x,i3)'
         if (it.eq.3)
     .   infmt(1:56) = 
c    .     '(7x,i2,1x,i4,25x,i2,i2,i2,5x,i4,9x,i3,i2,i4,i4,i4)'
     .     '(7x,i2,1x,i4,25x,i2,i2,i2,i2,i2,1x,i4,9x,i3,i2,i4,i4,i4)'
         if (it.eq.4)
     .   infmt(1:57) = 
     .     '(7x,i2,1x,i4,20x,i2,i2,i2,i2,i2,1x,i4,14x,i5,i4,a1,i3,i4)'
         if (it.eq.31)
     .   infmt(1:32) = 
     .     '(7x,i2,1x,i4,49x,i4,a1,4x,i4,a1)'

         do 60 j = 1,30000
            read (ifil1,'(a120)', end=90, err=60) line
            read (line,'(7x,i2)',err=60) ii
            if (ii.ne.60.and.ii.ne.20.and.ii.ne.22.and.ii.ne.52
     .         .and.ii.ne.85) goto 60
c           skip type-mismatch data
            if (it.eq.1.and.ii.ne.60) goto 60 
            if (it.eq.3.and.(ii.ne.20.and.ii.ne.22)) goto 60 
            if (it.eq.4.and.ii.ne.52) goto 60 
            if (it.eq.31.and.ii.ne.85) goto 60 
            ioei = 0
            ioer = 0
            lcode = nblen(line)

c           Azimuth
            if (it.eq.1) read (line,fmt=infmt,end=90,err=60) 
     .         ii,isit1,iyear,month,iday,ihh,imm,isit2,iod,iom,ios,
     .         ioer

c           Direction (1st)
            if (it.eq.3.and.ii.eq.20) 
     .         read (line,fmt=infmt,end=90,err=60) 
     .         ii,isit1,iyear,month,iday,ihh,imm,isit2,iod,iom,ios,
     .         ioei,ioer

c           Direction (subsequent)
            if (it.eq.3.and.ii.eq.22) 
     .         read (line,'(7x,i2,1x,i4,36x,i4,9x,
     .         i3,i2,i4,i4,i4)',end=90,err=60) 
     .         ii,isit1,isit2,iod,iom,ios,ioei,ioer

c           Reduced distance
            if (it.eq.4) read (line,fmt=infmt,end=90,err=60) 
     .         ii,isit1,iyear,month,iday,ihh,imm,isit2,iod,ios,
     .         code,ioei,ioer

c           Deflection
c           don't need to read hours and minutes for deflection
            if (it.eq.31) read (line,fmt=infmt,end=90,err=60) 
     .         ii,isit1,iod,idir1,ios,idir2

c           if the site name mismatch, reject the data
            call chk_sname(it,isit1,isit2,find1,find2,name1,name2)
            if (.not.find1) goto 60
            if (it.ge.1.and.it.le.10.and..not.find2) goto 60
c
            if (it.eq.31) goto 31
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

c           convert two digit years to 1900s (note hardwired here!)

            if (iyear.lt.1900) iyear = iyear+1900

c           Compute Errors and Output data
c           error unit is tenth of second -- based on Andrea's statement

c           Azimuth
            if (it.eq.1) then
c              convert to units of seconds
               sec = dble(ios)*1.0d-2
c              er  = dble(ioer)*1.0d-2
c              need to propogate error here!??
               er  = dsqrt(dble(ioer)**2+dble(ioei)**2)*1.0d-2
c              time1 = 1900.0d0+julday(month,iday,iyear,1)/dtoyr
               time1 = iyear + (julday(month,iday,iyear,3)-1)/dtoyr
c              add hours and minutes
               time1 = time1 + ihh/htoyr+imm/mtoyr
               write (ifil2,110) 
     *            time1,iod,iom,sec,er,name1
            endif

c           Direction
            if (it.eq.3) then

               sec = dble(ios)*1.0d-2
c              convert to units of 1 second and propogate errors
               er  = dsqrt(dble(ioer)**2+dble(ioei)**2)*1.0d-2
               if (ii.eq.20) then
c                 time1 = 1900.0d0+ 
c    .                    julday(month,iday,iyear,1)/dtoyr
                  time1 = iyear + (julday(month,iday,iyear,3)-1)
     .                    /dtoyr
c                 add hours and minutes
                  time1 = time1 + ihh/htoyr+imm/mtoyr
                  time0 = time1
               else
                  time1 = time0
               endif
               write (ifil2,110) 
     *            time1,iod,iom,sec,er,name1,name2
            endif

c           Reduced distance
            if (it.eq.4) then

c              make sure distances are mark-to-mark for now!
               if (code.ne.'X') then
                  print*,'Warning! Distance is not mark-to-mark'
                  goto 60
               endif

               blen = dble(iod)+dble(ios)*1.0d-4
c              units need to be metres
c              error model:  sigma**2 = a**2 + (b*L)**2
               er  = dsqrt((dble(ioer)*(1.0d-7)*blen)**2+
     .         (dble(ioei)*1.0d-4)**2)*1.0d+3
c              time1 = 1900.0d0+julday(month,iday,iyear,1)/dtoyr
               time1 = iyear + (julday(month,iday,iyear,3)-1)/dtoyr
c              add hours and minutes
               time1 = time1 + ihh/htoyr+imm/mtoyr
               write (ifil2,150) 
     *            time1,blen,er,name1,name2
            endif

c           Deflection
 31         if (it.eq.31) then
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
c---------------------------------------------------------------
      subroutine chk_sname(it,isit1,isit2,find1,find2,name1,name2)

      implicit real*8(a-h,o-z)
      include 'maked.fti'

      character*8 name1,name2
      logical find1,find2
      integer it,isit1,isit2,k
c
c     match the site index and assign names
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
c
      return
      end
