      program ficalabl


c     Read a fica file and translate it into another fica file
c
c     Convert block 70s into block 80s per Clynch's specs

c     Kurt Feigl February, 1989

c
C     file names
      character*80 fname1,fname2

C     FICA float array
      real*8 f1f(500), f2f(500)
C     FICA integer array
      integer*4 f1i(500),f2i(500)
C     FICA character array
      character*8 f1c(500),f2c(500)
C     number of elements in array
      integer n1f,n1i,n1c
C     number of elements in array
      integer n2f,n2i,n2c
C     number of the block
      integer iblkid1
C     number of the block
      integer iblkid2
C     block number, seen, read, written
      integer*2 tallyb(3,60)
C     I/O error codes
      integer*4 ioerr1,ioerr2
      integer i,j
C     logical unit numbers
      integer lficin,lficout,lscren
      character*40 version
C     username
      character*16 uname
      character*3 buff3
C     run time yr,mo,dy,hr,mn,sc
      integer*4 irtime(6)

      integer*4 idumb

      common /tally/ tallyb

      lscren = 6
      call uptaly(0,0)

      version = 'FICALABL  v. 1.1 of 89/02/18 '
      print *,version

      write (6,'(a,$)') 'Input FICA file name? '
      read '(a80)', fname1

      write (6,'(a,$)') 'Output FICA file name? '
      read '(a80)', fname2

      lficin = 10
      open(unit = lficin,
     +   file   = fname1,
     +   status = 'old',
     +   form   = 'formatted',
     +   iostat = ioerr1)
      if (ioerr1.ne.0) then
         print *, 'Error opening file: ', fname1
         call ferror(ioerr1,lscren)
         stop
      else
         print *, 'Opened file: ', fname1
      endif

      lficout = 12
      open(unit = lficout,
     +   file   = fname2,
     +   status = 'old',
     +   form   = 'formatted',
     +   iostat = ioerr2)
      if (ioerr2.ne.0) then
         print *, 'Error opening file: ', fname2
         call ferror(ioerr2,lscren)
         stop
      else
         print *, 'Opened file: ', fname2
      endif

C     WRITE A HEADER IN BLOCK 0

c     get the run time and version number
      call getdat(irtime(1),irtime(2),irtime(3))
      call gettim(irtime(4),irtime(5),irtime(6),idumb)
      write (lscren, 2000) version,(irtime(i),i=1,6)
 2000 format (/,1x,a32,1x,i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,
     .   ':',i2.2,//)

c     get the current time
      call getdat(irtime(1),irtime(2),irtime(3))
      call gettim(irtime(4),irtime(5),irtime(6),idumb)

c     program name
      write(f2c(1),1)
   1  format('FICALABL')

c     program version
      f2c(2) = '1.1'

c     institution
      write(f2c(3),2)
   2  format('MIT GPS')

c     date of conversion
      write(f2c(4),3) irtime(1)-1900,irtime(2),irtime(3)
   3  format (i2.2,'/',i2.2,'/',i2.2)

c     time of conversion
      write(f2c(5),13) irtime(4),irtime(5)
  13  format (1x,i2.2,':',i2.2,1x)

c     name of this file
c     assume GAMIT file naming convention, e.g. BLHL6.350
      write(f2c(6),'(a5,a3)') fname1(1:5),fname1(7:9)

c     name of the user
      call getusr(uname)
      write(f2c(7),'(1x,a7)') uname

c     source format
      write (f2c(8),'(a8)') ' FICA   '

c     write these goodies to the FICA header in block 0 and screen
      n2f = 0
      n2i = 0
      n2c = 8
      iblkid2 = 0
      call wfica(lficout,iblkid2,f2f,f2i,f2c,n2f,n2i,n2c)
      write (lscren,210) (f2c(i),i=1,8)
 210  format (10a8)

c     if this is a file with a naked header, then
c     read first line and copy it into a block 0
      read(lficin,210) (f2c(i),i=1,10)
      write (buff3,'(a3)') f2c(1)
      if (buff3 .ne. 'BLK') then
         n2f = 0
         n2i = 0
         n2c = 10
         iblkid2 = 0
         call wfica(lficout,iblkid2,f2f,f2i,f2c,n2f,n2i,n2c)
         write (lscren,210) (f2c(i),i=1,10)
      endif

      rewind (lficin)

c     loop until end of file or too many errors
 10   call rfica (lficin,lscren,iblkid1,f1f,f1i,f1c,n1f,n1i,n1c,ioerr1)

      if (ioerr1 .eq. 0) then
         call uptaly(iblkid1,1)
         if (iblkid1 .eq. 70) then
            iblkid1 = 80
         endif
         call wfica (lficout,iblkid1,f1f,f1i,f1c,n1f,n1i,n1c)
         call uptaly(iblkid1,2)

c        loop for another read
         goto 10
      else
         if (ioerr1 .gt. 0) then
            call ferror (ioerr1,lscren)
         else
            write (lscren,220) fname1
 220        format (1x,'Reached end of file:',/,a)
         endif
      endif
c


      write (6,*) '        '
      write (6,*) ' ** BLOCK COUNTS ** '
      write (6,*) 'block  read  wrote'
      do i = 1,60
         if (tallyb(2,i).ne.0 .or. tallyb(3,i).ne.0) then
               write (6, '(1x,i4,2(2x,i4))') (tallyb(j,i),j=1,3)
         endif
      enddo

      close (10)
      close (12)


      stop
      END


      subroutine uptaly (blocknum, itype)
c     update the read, write tally
      integer itype
C      0 to initialize the counts

C      1 to increment read count
c
C      2 to increment written count
c
      integer blocknum,i,j,nmax
      integer*2 tallyb(3,60)
      common /tally/ tallyb

      nmax = 52

      if (itype.eq.-1) then
         do i = 1,nmax
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
         do i = 2,3
             do j = 1,nmax
                tallyb(i,j) = 0
             enddo
         enddo
      else
         do i = 1,nmax
            if (tallyb(1,i) .eq. blocknum) then
               tallyb(itype+1,i) = tallyb(itype+1,i) + 1
               return
            endif
         enddo
      endif
      return
      end





