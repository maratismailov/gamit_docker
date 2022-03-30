c     GRAPHICS SUBROUTINES
c     SUN version
c     written under CGI

      subroutine gsetup (ibmx0,ibmy0,bmwdth,bmhght,lcolor,csize,ipos)

c     set up the graphics package
c     if already intialized, don't bother
c     obviously a very system dependent routine
c
c     returns:
c        bmwdth bit map width
c        bmhght bit map height
c        lcolor .true. if this is a color node
c        ipos   x,y of cursor when we started this
c        csize  width of a character in screen units
c       *name   name of this thing, so that it may be closed

      implicit real (a-h,o-z)
      implicit integer (i-n)

c     passed values
      integer*2 csize(2),ipos(2)
      logical lcolor
      integer*2 bmwdth,bmhght,ibmx0,ibmy0

      INCLUDE "/usr/include/f77/cgidefs77.h"

      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

      parameter(ibig=256)

c     local values
      integer*4 ix0,iy0,sx,sy

      integer*4 ifont,index,icfont,prec
      integer*4 icolor,hgt,ipath,lignh,lignv
      real*4 efac,space,bx,by,ux,uy,hfac,cfac

      character screenname*(ibig)
      integer screenlen
      character windowname*(ibig)
      integer windowlen,widowfd,retained,dd,cmapsize
      character cmapname*(ibig)
      integer cmaplen
      integer cfopenvws,cftextalign
      integer flags
      character ptr*(ibig),devid*(ibig)
      integer noargs

      integer halign,valign
      real hcalind,valignd

      call cfopencgi()

c     this works, but not under suntools
c     dd = 2
c     This seems to work under suntools
      dd = 9
      igerr = cfopenvws(name,screenname,windowname,
     .     windowfd,retained,dd,cmapsize,cmapname,flags,
     .     ptr, noargs)

      if (igerr .eq. 4) then
         write (*,*) 'CVIEW does not run under SUNTOOLS'
         stop
      endif

      call cfqphyscsys(name,ix0,iy0,nx,ny,sx,sy)
      ibmx0 = ix0
      ibmy0 = iy0
      bmwdth = nx
      bmhght = ny

c     character height
      hgt = ny/100
      call cfcharheight(hgt)

c     character alignment
c     normal
      halign = 3
c     top
      valign = 1
      igerr = cftextalign(halign,valign,hcalignd,vcalignd)

c     interior fill style is SOLIDI, with perimeter visibility OFF
      call cfintstyle(1,0)

c     fill color is black
      call cfflcolor(1)

c     size is as big as the screen
      call cfvdcext(ix0,iy0,ix0+nx,iy0+ny)

c     clipping is ON
      call cfclipind(1)

c     text attributes
      call cfqtextatts(ifont,index,icfont,prec,efac,space,icolor,
     .     hgt,bx,by,ux,uy,ipath,lignh,lignv,hfac,cfac)

c     character size
      csize(1) = hgt/efac
      csize(2) = hgt

c     color
      lcolor = .false.

c     original cursor position
      ipos(1) = ix0
      ipos(2) = iy0

      return

      end

      subroutine gmsg (iwin,msg)
c
c     prints a message in centered in the designated window
c
      implicit real (a-h,o-z)
      implicit integer (i-n)
      integer*2 iwin(4)
      integer*2 size(2)
      integer*4 igerr,cfrectangle
      integer x,y
      integer conx,cony,llpx,llpy,ulpx,ulpy,urpy,urpx
      character*(*) msg
      character*1 nextch

c     keep in mind that in CGI, Y increases upward

      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

      llpx = iwin(1)
      llpy = ny - (iwin(2) + iwin(4))
      urpx = iwin(1) + iwin(3)
      urpy = ny - iwin(2)

c     make it white
      call cfflcolor(0)
      igerr = cfrectangle(llpx,llpy,urpx,urpy)

c     fill color is black
      call cfflcolor(1)

c     size of text
      call gtsize(msg,size,igerr)

      x = iwin(1) + iwin(3)/2 - size(1)/2
      y = ny - (iwin(2) + iwin(4)/2 - size(2)/2)
      call cftext(x,y,msg)

      return
      end


      subroutine  gtsize (msg,size,igerr)
c     return x and y extents of text
      character*(*) msg
      integer*2 size(2)
      integer*4 igerr,i,len,nch,cfqtextext
      integer conx,cony,llpx,llpy,ulpx,ulpy,urpy,urpx
      character*1 nextch

      do i=len(msg),1,-1
         if (msg(i:i) .ne. ' ') goto 10
      enddo

 10   continue

c     protect against emptiness
      if (i .eq. 1) i = len(msg)

      nextch = ' '
      igerr = cfqtextext
     .(msg(1:i),nextch,conx,cony,llpx,llpy,ulpx,ulpy,urpx,urpy)

      size(1) = urpx - ulpx
      size(2) = ulpy - llpy

      return
      end



      subroutine gmvcur (ipos,igerr)
c     move the cursor to position x=ipos(1),y=ipos(2)
c     I don't think we need this one on the sun

      integer*2 ipos(2)
      integer*4 igerr

      return
      end


      subroutine gclip (iwin,lon)

c     turn graphics clipping on or off
c     lon is .true. to turn it on
c     the iwin is defined as x,y,dx,dy

      logical lon
      integer*2 iwin(4)
      integer*4 igerr,ixmin,ixmax,iymin,iymax,iflag,cfcliprect,cfclipind

c     keep in mind that in CGI, Y increases upward

      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

      ixmin = iwin(1)
      ixmax = iwin(1) + iwin(3)
      iymin = ny - (iwin(2) + iwin(4))
      iymax = ny - iwin(2)

      if (lon) then
         iflag = 2
      else
         iflag = 0
      endif

      igerr = cfcliprect(ixmin,ixmax,iymin,iymax)
      igerr = cfclipind(iflag)

      return
      end

      subroutine gmouse (key,pos,wait)
c     Get input from mouse
c     input:
c         wait .true. to go into a wait state
c     output:
c         pos  x,y of cursor
c         key  mouse button

c     passed values
      integer*2 pos(2)
      character key*1
      logical wait
      integer igerr

      include '../includes/macdep.h'

c     initialize input device
      integer devclass,devnum,x,y,xlist(2),ylist(2),n,choice
      real val
      character*256 string
      integer segid
      integer pickid
      integer trigger,trigger2
      integer valid
      integer timeout
      integer echotype
      integer exlow,eylow,exup,eyup
      integer cfinitlid,cftrackon,cfassoc,cfassoc,cfassoc,cfreqinp

c     keep in mind that in CGI, Y increases upward


      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

      logical done
      save done
      data done /.false./

      if (.not. done) then
c        initialize mouse
c        IC_STROKE
c        devclass = 1
c        IC_LOCATOR
         devclass = 0
c        why not?
         devnum = 1
         igerr = cfinitlid(devclass,devnum,x,y,xlist,ylist,n,val,
     .           choice,string,segid,pickid)

c        turn on tracking

c        default cursor
         echotype = 0
         igerr = cftrackon(devclass,devnum,echotype,
     .           exlow,eylow,exup,eyup,x,y,xlist,ylist,n,val,
     .           choice, string, segid, pickid)

c        associate trigger
c        left mouse button
         igerr = cfassoc(2,devclass,devnum)
c        middle mouse button
         igerr = cfassoc(3,devclass,devnum)
c        right mouse button
         igerr = cfassoc(4,devclass,devnum)
         done = .true.
      endif

c     wait forever or 1 microsecond
      if (wait) then
         timeout = -1
      else
         timeout = 1
      endif

      igerr = cfreqinp(devclass,devnum,timeout,valid,trigger,x,y,
     .        xlist,ylist,n,val,choice,string,segid,trigger2,pickid)

      pos(1) = x
      pos(2) = ny - y

      if (trigger .eq. 2) then
         key = MOUSE1
      else if (trigger .eq. 3) then
         key = MOUSE2
      else if (trigger .eq. 4) then
         key = MOUSE3
      else
         key = '.'
      endif

c        dissasociate
c         igerr = cfdissoc(trigger,devclass,devnum)

c        release input device
c         igerr = cfrelidev(devclass,devnum)


      return
      end



      subroutine gline (ix1,iy1,ix2,iy2,igerr)
c     draw a line from (ix1,iy1) to (ix2,iy2)
      integer*2 ix1,iy1,ix2,iy2
      integer*4 igerr
      integer*4 ix4(2),iy4(2),cfpolyline

      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

c     keep in mind that in CGI, Y increases upward

      integer*4 igerr,n4
      parameter (n4 = 2)

      ix4(1) = ix1
      iy4(1) = ny - iy1
      ix4(2) = ix2
      iy4(2) = ny - iy2
      igerr = cfpolyline(ix4,iy4,n4)

      return
      end


      subroutine gbox (iwin,igerr)
c     fill in window iwin
      integer*2 iwin(4)
      integer*4 ixmin,ixmax,iymin,iymax,cfrectangle
      integer*4 igerr


c     keep in mind that in CGI, Y increases upward

      integer*4 nx,ny,name
      common /gcom/ nx,ny,name


      ixmin = iwin(1)
      ixmax = iwin(1) + iwin(3)
      iymin = ny - (iwin(2) + iwin(4))
      iymax = ny - iwin(2)

c     interior fill style is SOLIDI, with perimeter visibility ON
      call cfintstyle(1,1)

      if (ixmin .lt. ixmax .and. iymin .lt. iymax) then
         igerr =  cfrectangle(ixmin,iymin,ixmax,iymax)
      else
         igerr = 0
      endif

c     interior fill style is SOLIDI, with perimeter visibility OFF
      call cfintstyle(1,0)

      return
      end


      subroutine gcirc (ix,iy, iradi, igerr)
c     draw an open circle of radius iradi at ix,iy
      integer*2 ix,iy,iradi
      integer*4 igerr,cfcircle,ix4,iy4,irad4

c     keep in mind that in CGI, Y increases upward

      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

c     interior fill style is HOLLOW, with perimeter visibility ON
      call cfintstyle(0,1)

      ix4 = ix
      iy4 = ny - iy
      irad4 = iradi

      igerr = cfcircle(ix4, iy4, irad4)

c     interior fill style is SOLIDI, with perimeter visibility OFF
      call cfintstyle(1,0)

      return
      end

      subroutine gdot (ix,iy, iradi, igerr)
c     draw an filled circle of radius iradi at ix,iy
      integer*2 ix,iy,iradi
      integer*4 igerr,cfcircle,ix4,iy4,irad4

c     keep in mind that in CGI, Y increases upward

      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

      ix4 = ix
      iy4 = ny - iy
      irad4 = iradi

      igerr = cfcircle(ix4, iy4, irad4)

      return
      end

      subroutine gplus (ix,iy, iradi, igerr)
c     draw a plus sign of radius iradi at ptcent(1),ptcent(2)
      integer*2 ix,iy,iradi
      integer*4 igerr

      call gline (ix-iradi,iy,ix+iradi,iy,igerr)
      call gline (ix,iy-iradi,ix,iy+iradi,igerr)

      return
      end

      subroutine gend (igerr)

c     terminate graphics package
      integer*4 igerr
      integer defflag,index

      integer*4 nx,ny,name
      common /gcom/ nx,ny,name

      defflag = 0
      index = 0
      call cfclrvws(name,defflag,index)
      call cfclosevws(name)
      call cfclosecgi()

      return
      end



      subroutine gdump (ibmx0,ibmy0,nx,ny,idump)
c     dump the screen to the printer
c     origin
      integer*2 ibmx0,ibmy0
c     height and width
      integer*2 bmwdth,bmhght
c     status code = 0 if OK
      integer idump, system
      integer nx,ny

      integer nblen
      character*80 command

      command = '~/gu/com/sd 1'
      idump = system (command(1:nblen(command)))

      return
      end


