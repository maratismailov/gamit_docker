      program ficaspan


c     dump the time tags from a FICA file so we can plot them
c     Kurt Feigl 5 March 87
c
      character*80 fname1,fname2
      integer*4 ioerr1,i
C     block should start with this
      character*4 blkstart
      integer*4 iweek,j,iwkn_start
C     decimal hour since start of day
      real*8 dechr
      real*8 sow_start
      character*40 version
      integer itset


      include '../includes/makex.h'

c     passed values
C     .true. at end of file
      logical     fend,
C     .true. to print things
     &            debug,
C     .true. on error
     &            ferr,
C     .true. once we have a week num
     &            week_set

C     1 for data, 2 for ephemeris
      integer*4   iflag,
C     for broadcast ephem
     &            nprn,
C     sat id num PRN
     &            svid(maxchn),
C     tracking mode
     &            tmode(maxchn),
C     gps week number for epoch
     &            igpswk(maxchn),
C     quality vector L1,L2
     &            iqvec(maxchn,2),
     .            iwkn_file

C     gps sec of week for epoch
      real*8      gpssec(maxchn),
C     L1 doppler phase (DATA!)
     &            dofl1(maxchn),
C     L2 doppler phase (DATA!)
     &            dofl2(maxchn),
C     L1 pseudorange
     &            prgl1(maxchn),
     &            prgl2(maxchn),
C     signal to noise ratio L1,L2
     &            denrat(maxchn,2),
C     broadcast ephemeris params
     &            bc_ephem(16),
C     broadcast clock params
     &            bc_clock(6),
C     additional SV quantites from sub-frame 1
     .            subfr1(7)



      itset = 0
      iweek = 0
      version = 'FICASPAN'
      print *,version

C     block should start with this
      blkstart = 'BLK '

      print *, 'FICA file? '
      read '(a80)', fname1

      uficaf = 10
      open(unit = uficaf,
     +   file   = fname1,
     +   status = 'old',
     +   form   = 'formatted',
     +   iostat = ioerr1)
      if (ioerr1.ne.0) then
         print *, 'Error opening file: ', fname1
         call ferror(ioerr1,6)
         stop
      else
         print *, 'Opened file: ', fname1
      endif

      i = index(fname1,' ') -1
      fname2 = fname1(1:i)//'.span'
      open(unit = 20,
     +   file   = fname2,
     +   status = 'unknown',
     +   form   = 'formatted',
     +   iostat = ioerr1)
      if (ioerr1.ne.0) then
         print *, 'Error opening file: ', fname2
         call ferror(ioerr1,6)
         stop
      else
         print *, 'Opened file: ', fname2
      endif

      fend = .false.
      ferr = .false.
      debug = .false.
      uscren = 6
      qinfor = .false.
      week_set = .false.

c     temporary kluge!

 15   if (.not. fend .and. .not. ferr) then
         call dofica
     &      (debug,iflag,fend,ferr,nprn,svid,tmode,igpswk,
     &      gpssec,dofl1,dofl2,prgl1,prgl2,denrat,iqvec,
     &      bc_ephem,bc_clock,subfr1,iwkn_file,week_set,
     .      iwkn_start,sow_start)

C        here we have PHASE data
         if (iflag.eq.1 .and. week_set) then
C          loop over channels
           do j=1,maxchn
              if (svid(j).ne.0) then
c                hours since the beginning of the day
                 dechr = dmod(gpssec(j),86400.0d0)/3600.
                 write(20,205) dechr,svid(j)
 205             format (f20.1,5x,i3)
               endif
            enddo
         endif

         goto 15
      endif
c

      if (fend) then
         print *, ' '
         print *, 'Reached end of file: ',fname1
      endif

      close (20)
      close (10)


      stop
      END



