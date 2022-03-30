      program ficacomp


c     compare 2 FICA files
c     Kurt Feigl 5 March 87
c
      real*8 f1f(500)
      real*8 f2f(500)
      integer*4 f1i(500)
      integer*4 f2i(500)
      character*8 f1c(500)
      character*8 f2c(500)
      character*80 fname1,fname2
C     block number, seen, read, written
      integer*2 tallyb(4,36)
      integer*4 ioerr1,ioerr2,i
C     number of elements in array
      integer n1f,n1i,n1c
C     number of elements in array
      integer n2f,n2i,n2c
C     block should start with this
      character*4 blkstart
C     number of the block
      integer blkid1
C     number of the block
      integer blkid2
C     only deal with this one
      integer blknew
C     pointer to tally for blknew
      integer talptr
C     fractional difference
      real*8 diff
C     allowable difference 1 part in thresh
      real*8 thresh
      character*40 version
      integer*4 j

      common /tally/ tallyb

      version = 'FICACOMP v. 1.2 of 88/08/30 11:29:15'
      print *,version

C     block should start with this
      blkstart = 'BLK '
      call uptaly(0,0)

      print *, 'First file? '
      read '(a80)', fname1
      print *, 'second file? '
      read '(a80)', fname2

      print *, 'Which block? '
      read *, blknew
      talptr = -1
      call uptaly(blknew,talptr)

c      print *, 'Minimum item number? '
c      read *, itemmin
c      print *, 'Maximum item number? '
c      read *, itemmax

c      print *, 'Threshold? (1 part in t) Enter t. '
c      read *, thresh

       thresh = 1e14

      open(unit = 10,
     +   file   = fname1,
     +   status = 'old',
     +   form   = 'formatted',
     +   iostat = ioerr1)
      if (ioerr1.ne.0) then
         print *, 'Error opening file: ', fname1
         goto 999
      else
         print *, 'Opened file: ', fname1
      endif

      open(unit = 11,
     +   file   = fname2,
     +   status = 'old',
     +   form   = 'formatted',
     +   iostat = ioerr2)
      if (ioerr2.ne.0) then
          print *, 'Error opening file: ', fname2
          goto 999
      else
         print *, 'Opened file: ', fname2
      endif


      print *, ' '
      print *, 'Found the following differences:'

  5   if (ioerr1 .eq. 0 .or. ioerr2 .eq. 0) then
c
c        read until you get the block you want
         blkid1 = -1
 10      if (ioerr1.eq.0 .and. blkid1 .ne. blknew) then
            call rfica (10,6,blkid1,f1f,f1i,
     +                 f1c,n1f,n1i,n1c,ioerr1)
            call uptaly(blkid1,1)
            goto 10
         endif
c
c        read until you get what you want
         blkid2 = -1
 15      if (ioerr2.eq.0 .and. blkid2 .ne. blknew) then
            call rfica (11,6,blkid2,f2f,f2i,
     +                  f2c,n2f,n2i,n2c,ioerr2)
            call uptaly(blkid2,2)
            goto 15
         endif

c

         if (ioerr1.eq.-1) then
            print *, ' '
            print *, 'Reached end of file: ',fname1
         endif
         if (ioerr2.eq.-1) then
            print *, 'Reached end of file: ',fname2
         endif


         if (ioerr1.ne.0 .and. ioerr1.ne.-1) then
            CALL ferror(ioerr1,6)
            ioerr1 = 0
         else if (ioerr2.ne.0 .and. ioerr2.ne.-1) then
            CALL ferror(ioerr2,6)
            ioerr2 = 0
         else if (blkid1.eq.blknew.and.blkid2.eq.blknew) then
            if (max(n1f,n2f) .gt. 0) then
               do i = 1,max(n1f,n2f)
                  diff = f1f(i) - f2f(i)
                  if (dabs(diff) .gt.
     +                dabs(max(f1f(i), f2f(i))) / thresh) then
                      print 110, blkid1,tallyb(2,talptr),
     +                          tallyb(3,talptr),i,
     +                          f1f(i),f2f(i),diff
  110                format (' Block: ',3(i3,1x),' F item: ',
     +                        i3,3(1x,1PD20.13))
                  endif
               enddo
            endif
            if (max(n1i,n2i) .gt. 0) then
               do i = 1,max(n1i,n2i)
                  if (f1i(i) .ne. f2i(i))
     +               print 120,blkid1,tallyb(2,talptr),
     +                         tallyb(3,talptr),
     +                         i,f1i(i),f2i(i),
     +                         f1i(i)-f2i(i)
  120              format (' Block: ',3(i3,1x),' I item: ',i3,3(1x,i12))
               enddo
            endif
            if (max(n1c,n2c) .gt. 0) then
               do i = 1,max(n1c,n2c)
                  if (f1c(i) .ne. f2c(i))
     +               print 130, blkid1,tallyb(2,talptr),
     +                         tallyb(3,talptr),
     +                         i,f1c(i),f2c(i)
  130               format (' Block: ',3(i3,1x),' C item: ',i3,2(1x,a8))
               enddo
            endif
         endif
         goto 5
      endif

      write (6,*) '        '
      write (6,*) ' ** BLOCK COUNTS ** '
      write (6,*) 'block  f1   f2 '
      do i = 1,36
         if (tallyb(2,i).ne.0 .or. tallyb(3,i).ne.0) then
               write (6, '(1x,i4,2(2x,i3))') (tallyb(j,i),j=1,3)
         endif
      enddo

  999 stop

      END


      subroutine uptaly(blocknum, itype)
c     update the read, write tally
C     -1 to return the index of blocknum (overwrites itype)
      integer itype
C      0 to initialize the counts

C      1 to increment seen count
c
C      2 to increment read count
c
C      3 to increment written count
c
      integer blocknum,i,j
      integer*2 tallyb(4,36)
      common /tally/ tallyb


      if (itype.eq.-1) then
         do i = 1,36
            if (tallyb(1,i) .eq. blocknum) then
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
         do i = 2,4
             do j = 1,36
                tallyb(i,j) = 0
             enddo
         enddo
      else
         do i = 1,36
            if (tallyb(1,i) .eq. blocknum) then
               tallyb(itype+1,i) = tallyb(itype+1,i) + 1
               return
            endif
         enddo
      endif
      return
      end





