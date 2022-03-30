c  Subroutines to interface the fortran cview program with the
c  X11 Xlib wrappers which are written in C.
c  Written by Andrea Donnellan 7/24/1990
c  E-mail andrea@seismo.gps.caltech.edu

      subroutine  gtsize(msg,size,igerr)
c     return x and y extents of text
      character*(*) msg
      integer*2 size(2)
      integer*4 igerr

      character*320 string

      string = msg
      call gtsizec(string,size,len(msg))

      return
      end

      subroutine gmsg(iwin,msg)
      character*(*) msg
      integer*2 iwin(4)
      integer*2 size(2)

      character*320 string

      string = msg
      call gtsizec(string,size,len(msg))
      call gmsgc(iwin,string,size)

      return
      end


      subroutine gdump (ibmx0,ibmy0,nx,ny,lcolor,idump)
c     dump the screen to the printer
c     origin
      integer*2 ibmx0,ibmy0
c     height and width
      integer nx,ny
c     status code = 0 if OK
      integer idump
c     .true. if color
      logical lcolor

      integer nblen
      character*80 command

      integer istat,istat2,istat3,system,unlink,ftell
      integer ioff5,ioff6

      EXTERNAL lunit,unlink,fcheck

      istat = 0
      istat2 = 0
      istat3 = 0

      if (lcolor) then
         write (command,*)
     .   'screendump | rasfilter8to1 | psraster -i | lpr &'
      else
         write (command,*)
     .   'screendump | psraster -i | lpr &'
      endif

c     remember where we are in unit 5
      ioff5 = ftell(5)
      if (ioff5.lt.0) then
         write (6,*) 'GDUMP: error on ftell. IOFF5 = ',ioff5
         call ferror(-ioff5,6)
      endif

c     remember where we are in unit 6
      ioff6 = ftell(6)
      if (ioff6.lt.0) then
         write (6,*) 'GDUMP: error on ftell. IOFF6 = ',ioff6
         call ferror(-ioff6,6)
      endif

c     execute the shell command
      istat = system (command(1:nblen(command)))
      if (istat .ne. 0) then
         write (6,*) 'GDUMP: istat = ',istat
         call ferror(istat,6)
      endif

c     restore file pointer on unit 5, because system screws it up.
c** function      istat2 = fseek(5,ioff5,0)
c** intrinsic     call fseek(5,ioff5,0)
c      call fseekg(5, ioff5, 0, istat2)
      if (istat2.ne.0) then
         write (6,*) 'GDUMP: error on seek on unit 5.'
         call ferror(istat2,6)
      endif

c     restore file pointer on unit 6, because system screws it up.
c** function      istat3 = fseek(6,ioff6,0)
c** intrinsic     call fseek(6,ioff6,0)
c      call fseekg(6, ioff6, 0, istat3)
      if (istat3.ne.0) then
         write (6,*) 'GDUMP: error on seek on unit 6.'
         call ferror(istat2,6)
      endif

      if (istat .eq. 0 .and. istat2 .eq. 0 .and. istat3 .eq. 0) then
         idump = 0
      else
         idump = -1
      endif

      return
      end


