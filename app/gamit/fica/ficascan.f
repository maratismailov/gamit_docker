      program ficascan


c     scan a FICA file and count blocks
c     Kurt Feigl 5 March 87
c
      include '../includes/makex.h'
      real*8 f1f(maxfic)
      integer*4 f1i(maxfic)
      character*8 f1c(maxfic)
      character*80 fname1,buff80
C     block number, seen, read, written
      integer  tallyb(4,100)
      integer*4 ioerr1,i
C     number of elements in array
      integer n1f,n1i,n1c
C     block should start with this
      character*4 blkstart
C     number of the block
      integer iblkid
      integer*4 iweek,iwk1,j
      real*8 sow1
c     collection interval
      real*8 asampl
      character*40 version
      integer itset,nblock,iflag

c     values from receiver operator
      character*16 auser
C     user-input site code
     .,            asite
C     reciever serial number
     .,            arcvr
C     operator input antenna ht.
     .,            antht
C     atmospheric pressure
      real*8       apress
C     atmospheric temperature
     .,            atemp
C     atmospheric pressure
     .,            ahumid
C     software version
     .,            swveru

      itset = 0
      iweek = 0
      uscren = 6
      version = 'FICASCAN  v. 1.4 of 89/08/11 14:28:02'
      write (uscren,*) version

C     block should start with this
      blkstart = 'BLK '
      call uptaly(tallyb, 0, 0, 100, nblock)

      write (uscren,*) 'FICA file? '
      read '(a80)', fname1

      open(unit = 10,
     +   file   = fname1,
     +   status = 'old',
     +   form   = 'formatted',
     +   iostat = ioerr1)
      if (ioerr1.ne.0) then
         write (uscren,*) 'Error opening file: ', fname1
         call ferror(ioerr1,6)
         stop
      else
         write (uscren,*)  'Opened file: ', fname1
      endif

c     read until end of file
      iblkid = -1
 10   if (ioerr1.eq.0) then
         call rfica (10,6,iblkid,f1f,f1i,f1c,n1f,n1i,n1c,ioerr1)
         if (iblkid .eq. 0) then
            write (uscren,210) (f1c(i),i=1,8)
 210        format (1x,'FICA header information:',/,1x,8a8)
         endif

         call uptaly(tallyb, iblkid, 1, 100, nblock)
c        guess what time it is

         if (iblkid .eq. 6 .or. iblkid .eq. 55) then
            sow1 = f1f(3)
            if (iweek .gt. 0) then
               iwk1 = iweek
               itset = itset+1
            endif
         else if (iblkid .eq. 9) then
            iweek = int(f1f(6))
         else if (iblkid .eq. 101) then
            call blk101
     .      ( .false.,iflag,apress,atemp,ahumid,asampl,swveru
     .      , auser,asite,arcvr,antht
     .      , f1f,f1i,f1c,n1f,n1i,n1c )
         else if (iblkid .eq. 70) then
            iwk1 = f1i(2)
            sow1   = f1f(1)
            itset = itset+1
         else if (iblkid .eq. 80 .or.
     .            iblkid .eq. 1180) then
            iwk1 = f1i(2)
            sow1   = f1f(1)
            itset = itset+1
         else if (iblkid .eq. 401) then
            itset = itset+1
            iwk1 = int(f1f(4))
c           this should end in .08 and be the time at the end of
c           the phase measurement
            sow1 = f1f(3) + (f1f(5) - f1f(2)) + f1f(8)/2.d0
         endif

c        tell us the time of the first few blocks
         if (itset.gt.0 .and. itset .lt. 5) then
            write (buff80,310) itset
 310        format (1x,'Good epoch # ',i2)
            call wtime (6,iwk1,sow1,buff80)
         endif

         goto 10
      endif
c

      if (ioerr1.eq.-1) then
         write (uscren,*)  ' '
         write (uscren,*) 'Reached end of file: ',fname1
         call wtime (uscren,iwk1,sow1,'Ending epoch   ')
      endif

      write (uscren,*)  '        '
      write (uscren,*)  ' ** BLOCK COUNTS ** '
      write (uscren,*)  'block  saw'
      do i = 1,nblock
         if (tallyb(2,i).ne.0) then
               write (uscren, '(1x,i4,1(2x,i4))') (tallyb(j,i),j=1,2)
         endif
      enddo

      stop
      END

      subroutine uptaly(tallyb, iblock, itype, maxtal, nblock)
c     update the read, write tally
C     -1 to return the index of iblock (overwrites itype)
      integer itype
C      0 to initialize the counts

C      1 to increment seen count
c
C      2 to increment read count
c
C      3 to increment written count
c
      integer iblock,i,nblock,j,maxtal
      integer tallyb(4,maxtal)

      nblock = 53

      if (itype.eq.-1) then
         do i = 1,nblock
            if (tallyb(1,i) .eq. iblock) then
               itype = i
               return
            endif
         enddo
C     initialize everything
      else if (itype.eq.0) then
         tallyb(1,1) = 1
         tallyb(1,2) = 2
         tallyb(1,3) = 3
         tallyb(1,4) = 6
         tallyb(1,5) = 7
         tallyb(1,6) = 8
         tallyb(1,7) = 9
         tallyb(1,8) = 10
         tallyb(1,9) = 11
         tallyb(1,10) = 12
         tallyb(1,11) = 13
         tallyb(1,12) = 50
         tallyb(1,13) = 51
         tallyb(1,14) = 52
         tallyb(1,15) = 53
         tallyb(1,16) = 54
         tallyb(1,17) = 55
         tallyb(1,18) = 56
         tallyb(1,19) = 57
         tallyb(1,20) = 58
         tallyb(1,21) = 59
         tallyb(1,22) = 62
         tallyb(1,23) = 63
         tallyb(1,24) = 1001
         tallyb(1,25) = 1002
         tallyb(1,26) = 1003
         tallyb(1,27) = 1004
         tallyb(1,28) = 1005
         tallyb(1,29) = 1006
         tallyb(1,30) = 1007
         tallyb(1,31) = 1008
         tallyb(1,32) = 1009
         tallyb(1,33) = 101
         tallyb(1,34) = 102
         tallyb(1,35) = 109
         tallyb(1,36) = 162
         tallyb(1,37) = 400
         tallyb(1,38) = 401
         tallyb(1,39) = 402
         tallyb(1,40) = 403
         tallyb(1,41) = 404
         tallyb(1,42) = 405
         tallyb(1,43) = 410
         tallyb(1,44) = 411
         tallyb(1,45) = 420
         tallyb(1,46) = 421
         tallyb(1,47) = 423
         tallyb(1,48) = 424
         tallyb(1,49) = 425
         tallyb(1,50) = 426
         tallyb(1,51) = 70
         tallyb(1,52) = 80
         tallyb(1,53) = 1180
         do i = 2,4
             do j = 1,nblock
                tallyb(i,j) = 0
             enddo
         enddo
      else
         do i = 1,nblock
            if (tallyb(1,i) .eq. iblock) then
               tallyb(itype+1,i) = tallyb(itype+1,i) + 1
               return
            endif
         enddo
      endif
      return
      end


