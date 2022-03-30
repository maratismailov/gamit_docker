      program ficachop

c     read a FICA file and cut out various blocks

c     Kurt Feigl June 88



      integer*4        ibkids(30),nbkid,i,j,
     .                 uscren,ukeybd,uficin,uficot

      character*80     version,fficin,fficot
C     username
      character*16     uname
C     number of the block
      integer*4 iblkid
      integer*4 ios
C     run time yr,mo,dy,hr,mn,sc
      integer*4 irtime(6)
c     the 3 FICA arrays
      real*8 ff(500)
      integer*4 fi(500)
      character*8 fc(500)
c     the number of elements in them
      integer*4 nf,ni,nc
      integer*4 idumb

      logical fend,ferr
      fend = .false.
      ferr = .false.

      data version/'FICACHOP v 1.6 88/07/20 11:24:44'/

C     unit numbers
      ukeybd = 5
      uscren = 6
      uficin = 10
      uficot = 11

c     get the run time and version number
      call getdat(irtime(1),irtime(2),irtime(3))
      call gettim(irtime(4),irtime(5),irtime(6),idumb)
      write (6, 2000) version,(irtime(i),i=1,6)
 2000 format (/,1x,a32,1x,i4,'-',i2.2,'-',i2.2,1x,i2.2,':',i2.2,
     .   ':',i2.2,//)

      write (6, *)' FICACHOP is a simple utility to select given blocks'
      write (6, *)' from an input FICA file and write them to a FICA'
      write (6, *)' file.  In addition to passing the selected blocks,'
      write (6, *)' FICACHOP also passes through the header blocks '
      write (6, *)' (id = 0), and writes its own header block.'
      write (6, *)' '


      write (uscren,5)
  5   format (1x,'Enter input FICA file name')
      read  (ukeybd,'(a)') fficin

      write (uscren,7)
  7   format (1x,'Enter output FICA file name')
      read  (ukeybd,'(a)') fficot

      write (uscren,10)
  10  format (1x,'Enter FICA block id numbers to extract.',/,
     .        1x,'Please enter one id per line, and finish with -1')

c     loop until you get a -1
      i = 0
  12  continue
      read  (ukeybd,*) j
      if (j .ge. 0) then
         i = i+1
         ibkids(i) = j
         goto 12
      endif
      nbkid = i

CD     write (uscren,30) (ibkids(i),i=1,nbkid)
CD 30   format (1x,'Extracting blocks: ',20(1x,i3))

      open (unit = uficin,
     .      file = fficin,
     .      status = 'old',
     .      iostat = ios)

      if (ios.ne.0) then
         call ferror(ios,6)
         stop
      endif

      open (unit = uficot,
     .      file = fficot,
     .      status = 'new',
     .      iostat = ios)

      if (ios.ne.0) then
         call ferror(ios,6)
         stop
      else
c        prepare to write a header
c        get the current time

c        program name
         write(fc(1),1)
   1     format('FICACHOP')

c        program version
         fc(2) = ' 1.6'

c        institution
         write(fc(3),2)
   2     format('MIT GPS')

c        date of conversion
         write(fc(4),3) irtime(1)-1900,irtime(2),irtime(3)
   3     format (i2.2,'/',i2.2,'/',i2.2)

c        time of conversion
         write(fc(5),13) irtime(4),irtime(5)
  13     format (1x,i2.2,':',i2.2,1x)

c        name of this file
c        assume GAMIT file naming convention, e.g. BLHL6.350
         write(fc(6),'(a5,a3)') fficin(1:5),fficin(7:9)

c        name of the user
         call getusr(uname)
         write(fc(7),'(1x,a7)') uname

c        source format
         write (fc(8),'(a8)') ' FICA   '

c        write these goodies to the FICA header in block 0 and screen
         nf = 0
         ni = 0
         nc = 8
         iblkid = 0
         call wfica(uficot,iblkid,ff,fi,fc,nf,ni,nc)
         write (uscren,210) (fc(i),i=1,8)

c        begin loop to reading FICA blocks
 500     continue
         call rfica (uficin,6,iblkid,ff,fi,fc,nf,ni,nc,ios)
CD        write(uscren,510) iblkid
CD 510     format(1x,'Read block ',i4)

C        end of file
         if (ios .eq. -1) then
            fend = .true.
            stop
         else if (ios .ne. 0 ) then
            write (uscren,110)
 110        format (1x,'REPEATED ERROR IN FICA FILE: ')
            call ferror(ios,6)
            ferr = .true.
            stop
         endif

c        if the block is in the list, pass it through
c        pass all block 0s through anyway
         if (.not. ferr) then
            if (iblkid .eq. 0) then
               write (uscren,210) (fc(i),i=1,8)
 210           format (1x,'FICA header information',/,8(1x,a8,/))
            endif
            do i = 1,nbkid
               if (iblkid .eq. ibkids(i) .or. iblkid .eq. 0) then
                  call wfica(uficot,iblkid,ff,fi,fc,nf,ni,nc)
               endif
            enddo
            goto 500
         endif
      endif

      close (uficot)
      close (uficin)

      stop
      end


