      program ficamerg


c     merge blocks from 2 different 2 FICA files
c     Kurt Feigl July 88
c
      real*8 ff(500)
      real*8 f1f(500)
      real*8 f2f(500)
      integer*4 fi(500)
      integer*4 f1i(500)
      integer*4 f2i(500)
      character*8 fc(500)
      character*8 f1c(500)
      character*8 f2c(500)
      character*80 filename1,filename2,filenameo,uname
      integer*4 ios,ios1,ios2,i,j
C     number of elements in array
      integer nf,ni,nc
C     number of elements in array
      integer n1f,n1i,n1c
C     number of elements in array
      integer n2f,n2i,n2c
C     number of the block
      integer iblkid1,iblkid2,iblkid
C     GPS week
      integer*4 iwk1,iwk2,iweek
C     only deal with this one
      integer iblkf1,iblkf2
C     pointer to tally for iblkf1
      integer tallyptr
C     yr,mo,day,hr,min,sec
      integer*4 irtime(6)
C     return time diference in seconds
      real*8 secdif
C     second of GPS week
      real*8 sow1,sow2
      logical get1,get2

      integer maxtal,idumb
      integer*2 tallyb(4,50)
      common /tally/ tallyb

      maxtal = 50

c     flags to read another block from file 1, file 2
      get1 = .true.
      get2 = .true.

      call uptaly(0,0)
      tallyptr = -1
      call uptaly(iblkf1,tallyptr)

      print *, 'First input file? '
      read '(a80)', filename1
      print *, 'Second input file? '
      read '(a80)', filename2

      print *, 'Output file? '
      read '(a80)', filenameo


      print *, 'Which block from first file? '
      read *, iblkf1
      print *, 'Which block from second file? '
      read *, iblkf2

      print *, 'Week number? '
      read *, iweek


      open(unit = 10,
     +   file   = filename1,
     +   status = 'old',
     +   form   = 'formatted',
     +   iostat = ios1)
      if (ios1.ne.0) then
         print *, 'Error opening file: ', filename1
         call ferror(ios1,6)
         stop
      else
         print *, 'Opened file: ', filename1
      endif

      open(unit = 11,
     +   file   = filename2,
     +   status = 'old',
     +   form   = 'formatted',
     +   iostat = ios2)
      if (ios2.ne.0) then
          print *, 'Error opening file: ', filename2
          call ferror(ios2,6)
          stop
      else
         print *, 'Opened file: ', filename2
      endif

      open(unit = 12,
     +   file   = filenameo,
     +   status = 'new',
     +   form   = 'formatted',
     +   iostat = ios)
      if (ios.ne.0) then
          print *, 'Error opening file: ', filenameo
          call ferror (ios,6)
          stop
      else
         print *, 'Opened file: ', filenameo
      endif

c     prepare to write a header
c     get the current time
      call getdat(irtime(1),irtime(2),irtime(3))
      call gettim(irtime(4),irtime(5),irtime(6),idumb)

c     program name
      write(fc(1),1)
   1  format('FICAMERG')

c     program version
      fc(2) = ' 1.7'

c     institution
      write(fc(3),2)
   2  format('MIT GPS')

c     date of conversion
      write(fc(4),3) irtime(1)-1900,irtime(2),irtime(3)
   3  format (i2.2,'/',i2.2,'/',i2.2)

c     time of conversion
      write(fc(5),13) irtime(4),irtime(5)
  13  format (1x,i2.2,':',i2.2,1x)


c     name of this file
c     assume GAMIT file naming convention, e.g. BLHL6.350
      write(fc(6),'(a8)') filenameo(1:5)//filenameo(7:9)

c     name of the user
      call getusr(uname)
      write(fc(7),'(1x,a7)') uname

c     source format
      write (fc(8),7)
   7  format ('FICA',4x)

c     write these goodies to the FICA header in block 0 and screen
      nf = 0
      ni = 0
      nc = 8
      iblkid = 0
      call wfica(12,iblkid,ff,fi,fc,nf,ni,nc)
      call wfica(6, iblkid,ff,fi,fc,nf,ni,nc)


c     loop until the end of both files
 5    continue
      if (ios1+ios2 .gt. -2 .and. ios1+ios2 .le. 0) then
c        keep on reading
c        read file 1 until you get the block you want
         if (get1) then
            iblkid1 = -1
 10         if (ios1.eq.0 .and. iblkid1 .ne. iblkf1) then
               call rfica (10,6,iblkid1,f1f,f1i,f1c,n1f,n1i,n1c,ios1)
               call uptaly(iblkid1,1)
c              pass through block zero and block 101
               if (iblkid1 .eq. 0 .or. iblkid1 .eq. 101) then
c                 write the block from file 1
                  call wfica (12,iblkid1,f1f,f1i,f1c,n1f,n1i,n1c)
                  call uptaly(iblkid1,3)
               endif
               goto 10
            endif
            get1 = .false.
            if (ios1.eq.-1) then
               print *, 'Reached end of file: ',filename1
            endif
         endif
c
c        read file 2 until you get what you want
         if (get2) then
            iblkid2 = -1
 15         if (ios2.eq.0 .and. iblkid2 .ne. iblkf2) then
               call rfica (11,6,iblkid2,f2f,f2i,f2c,n2f,n2i,n2c,ios2)
               call uptaly(iblkid2,2)
c              pass through block zero and block 101
               if (iblkid2 .eq. 0 .or. iblkid2 .eq. 101) then
c                 write the block from file 2
                  call wfica (12,iblkid2,f2f,f2i,f2c,n2f,n2i,n2c)
                  call uptaly(iblkid2,3)
               endif
               goto 15
            endif
            get2 = .false.
            if (ios2.eq.-1) then
               print *, 'Reached end of file: ',filename2
            endif
         endif


         if (ios1.gt.0) then
            call ferror(ios1,6)
         endif
         if (ios2.gt.0) then
            call ferror(ios2,6)
         endif

         if (ios1 .eq. 0) then
c           try to figure out what time it is
            if (iblkf1 .eq. 6 .or. iblkf1 .eq. 55) then
               iwk1 = iweek
               sow1 = f1f(3)
            else if (iblkf1 .eq. 9) then
               iwk1 = int(f1f(6))
               sow1 = f1f(13)
            else if (iblkf1 .eq. 401) then
               iwk1 = int(f1f(4))
               sow1 = f1f(3) + (f1f(5) - f1f(6))
            endif
         endif
         if (ios2 .eq. 0) then
c           try to figure out what time it is
            if (iblkf2 .eq. 6 .or. iblkf2 .eq. 55) then
               iwk2 = iweek
               sow2 = f2f(3)
            else if (iblkf2 .eq. 9) then
               iwk2 = int(f2f(6))
               sow2 = f2f(13)
            else if (iblkf2 .eq. 401) then
               iwk2 = int(f2f(4))
               sow2 = f2f(3) + (f2f(5) - f2f(6))
            endif
         endif

CD        print *,'time tags ',iwk1,iwk2,sow1,sow2

         if (ios1 .eq. 0 .and. ios2 .eq. 0) then
            if (secdif(iwk2,sow2,iwk1,sow1) .ge. 0.) then
c              write the block from file 1
               call wfica (12,iblkid1,f1f,f1i,f1c,n1f,n1i,n1c)
               call uptaly(iblkid1,3)
c              set up to read another block from file 1
               if (ios1.ne.-1) get1 = .true.
            else
c              write the block from file 2
               call wfica (12,iblkid2,f2f,f2i,f2c,n2f,n2i,n2c)
               call uptaly(iblkid2,3)
c              set up to read another block from file 2
               if (ios2.ne.-1) get2 = .true.
            endif
         else
            if (ios1 .eq. -1) then
c              write the block from file 2
               call wfica (12,iblkid2,f2f,f2i,f2c,n2f,n2i,n2c)
               call uptaly(iblkid2,3)
c              set up to read another block from file 2
               get1 = .false.
               if (ios2.ne.-1) get2 = .true.
            else if (ios2 .eq. -1) then
c              write the block from file 1
               call wfica (12,iblkid1,f1f,f1i,f1c,n1f,n1i,n1c)
               call uptaly(iblkid1,3)
c              set up to read another block from file 1
               get2 = .false.
               if (ios1.ne.-1) get1 = .true.
            endif
         endif
         goto 5
      endif


      write (6,*) '        '
      write (6,*) ' ** BLOCK COUNTS ** '
      write (6,*) 'block file1 file2 wrote'
      do i = 1,maxtal
         if (tallyb(2,i).ne.0 .or. tallyb(3,i).ne.0) then
               write (6, '(2x,i4,3(2x,i4))') (tallyb(j,i),j=1,4)
         endif
      enddo


      stop
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
      integer blocknum,i,nmax,j
      integer*2 tallyb(4,50)
      common /tally/ tallyb

      nmax = 50

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
         tallyb(1,nmax) = 426
         do i = 2,4
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


