      subroutine rfica (luin,luout,iblkid,ff,fi,fc,nf,ni,nc,ioerror)
c     read a FICA block

      implicit none

C     logical unit of input FICA file (previously opened)
      integer luin
C     logical unit of error output file (previously opened)
      integer luout
C     blocks should start with this
      character*4 blkstr
C     they do start with this
      character*4 blkmrk
C     number of elements in array
      integer nf,ni,nc
C     number of the block
      integer iblkid
      integer ioerror,i
c     the FICA arrays
      integer maxfic
      parameter (maxfic = 500)
      real*8 ff(maxfic)
      integer*4 fi(maxfic)
      character*8 fc(maxfic)
c     variables for status reporting
      integer*4 len,rcpar
      character*16 fname
      character*80 buf80,prog_name
      character*256 message

      integer*4 ierct
      save ierct

      data ierct/0/
      data blkstr/'BLK '/

c     get calling program name and file name for report_stat
      len = rcpar(0,prog_name)
      inquire ( unit=luin, name=fname, iostat=ioerror )
      if( ioerror.ne.0 ) goto 1000
               
c     initialize the arrays to avoid leaving non-ascii characters
c     that will cause read errors later
      do i=1,maxfic
        ff(i) = 0.
        fi(i) = 0
        fc(i) = ' '
      enddo

c     This character read necessary mainly to protect against files
c     with a header rather than a Block 0

   10 read(unit=luin,fmt='(a80)',err=1000,end=500,iostat=ioerror) buf80
      if( buf80(1:3) .ne. 'BLK' ) goto 10

      read(unit   = buf80,
     +   fmt    = '(A4,1x,4I5)',
     +   err    = 1000,
     +   end    = 500,
     +   iostat = ioerror)
     +   blkmrk, iblkid, nf, ni, nc

CD     print *,'Block: ',iblkid,nf,ni,nc
      if (blkmrk .eq. blkstr) then
         if (nf .gt. 0) then
            read(unit   = luin,
     +         fmt    = '(4(1PD20.13))',
     +         err    = 1000,
     +         end    = 500,
     +         iostat = ioerror) (ff(i),i=1,nf)
         endif
         if (ni .gt. 0) then
            read(unit   = luin,
     +      fmt    = '(6I12)',
     +      err    = 1000,
     +      end    = 500,
     +      iostat = ioerror) (fi(i),i=1,ni)
         endif
         if (nc .gt. 0) then
            read(unit   = luin,
     +         fmt    = '(10a8)',
     +         err    = 1000,
     +         end    = 500,
     +         iostat = ioerror) (fc(i),i=1,nc)
         endif
      else
         write(message,'(a,z8)') blkmrk
         call report_stat('WARNING',prog_name,'lib/rfica',' '
     .          , 'Unrecognized block delimiter (hex)',0)
      endif

 500  continue
 1000 continue

c     Handle errors other than EOF here.  Only pass nonzero ioerror back
c     if EOF or if error count is too high.
      if (ioerror .gt. 0) then
         ierct = ierct+1
         write(message,'(a,i5)')
     .      'Read error in FICA file. Error count is ',ierct
         write(luout,'(a)') message
c**         call report_stat('WARNING',prog_name,'lib/rfica',' ',message,0)
         if (ierct .gt. 100) then
            write(message,'(a)') 'Over 100 read errors in FICA file: '
            write (luout,'(a)') message
            call report_stat('WARNING',prog_name,'lib/rfica',fname
     .                      ,message,0)
         else
            ioerror = 0
         endif
      endif

      return
      end

