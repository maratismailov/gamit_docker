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


      subroutine gdump (ibmx0,ibmy0,nx,ny,idump)
c     dump the screen to the printer
c     origin
      integer*2 ibmx0,ibmy0
c     height and width
      integer*2 bmwdth,bmhght
c     status code = 0 if OK
      integer idump
      integer nx,ny
      integer system
c      external system

      integer nblen
      character*80 command

      command = '~/gu/com/sd 1'
      idump = system (command(1:nblen(command)))

      return
      end


