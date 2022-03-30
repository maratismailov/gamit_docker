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
c
      INCLUDE '/sys/ins/base.ins.ftn'
      INCLUDE '/sys/ins/gpr.ins.ftn'
      INCLUDE '/sys/ins/time.ins.ftn'
      INCLUDE '/sys/ins/gmf.ins.ftn'
      INCLUDE '/sys/ins/kbd.ins.ftn'

c     origin
      integer*2 ibmx0,ibmy0
c     height and width
      integer*2 bmwdth,bmhght
c     color?
      logical color
c     character size
      integer*2 csize(2)
c     original cursor position
      integer*2 ipos(2)

      integer*2 ibmsiz(2), font_id, nch, size(2)
     .,         disp_bitmap_size(2), xtext, ytext,xhld
     .,         xl0(2), yl0(2), key_set(16), event_type,pos(2)
     .,         ytext1, xtext1, config, disp(31), disp_len
     .,         disp_len_ret, spacew

      integer*4 curbmdesc,attr,curatdesc
      integer*2 rastop(8),opos(2),siz(2),hpid,idir


      integer*4 init_bitmap, status
      integer*4 value(0:15)

      integer cplanes
      integer nblen                                { function returns # in string

      logical lcolor,active
      character*80 fname



      data value / 0, 16777215, 65280, 57600, 51200, 44800, 38400,
     &             32000, 16777000, 14803200, 13158400, 11513600,
     &             9868800, 8224000, 9868950, 16711700/
c

c     get screen size
      disp_len=31
      call gpr_$inq_disp_characteristics
     .   (gpr_$borrow,int2(1),disp_len,disp,disp_len_ret,status)
      cplanes = disp(15)
      bmwdth = disp(5)
      bmhght = disp(6)
      if (cplanes .gt. 1) then
         lcolor = .true.
      else
         lcolor = .false.
      endif

c     origin
      ibmx0 = disp(3)
      ibmy0 = disp(4)

      ibmsiz(1) = bmwdth
      ibmsiz(2) = bmhght
CD     print *,'Width, height',bmwdth,bmhght

      call gpr_$init
     .   (gpr_$borrow,int2(1),ibmsiz,int2(0),init_bitmap,status)

c     how big is a space in default font?
      call gpr_$inq_text(font_id,idir,status)
      call gpr_$inq_space_size(font_id,spacew,status)
CD     print *,'GSETUP: space width',spacew

c
c     load color table if color node
      if (lcolor) then
         call gpr_$set_color_map (int4(0), int2(16), value, status)
         fname = '/sys/dm/fonts/f7x13'
      else
         fname = '/sys/dm/fonts/legend'
      end if
      call gpr_$load_font_file
     .    (fname,int2(nblen(fname)),font_id,status)

c     make the cursor pattern a cross by drawing one in a little bit map
      siz(1)=15
      siz(2)=15
      hpid = int2(0)
      call gpr_$allocate_attribute_block(curatdesc,status)
      if (status.ne.0) call error_$print(status)
      call gpr_$allocate_bitmap(siz,hpid,curatdesc,curbmdesc,status)
      if (status.ne.0) call error_$print(status)
      call gpr_$set_bitmap(curbmdesc,status)
      if (status.ne.0) call error_$print(status)
      call gpr_$move(int2(0 ),int2( 7),status)
      if (status.ne.0) call error_$print(status)
      call gpr_$line(int2(14),int2( 7),status)
      if (status.ne.0) call error_$print(status)
      call gpr_$move(int2( 7),int2( 1),status)
      if (status.ne.0) call error_$print(status)
      call gpr_$line(int2( 7),int2(14),status)

c     switch back to screen bit map
      call gpr_$set_bitmap(init_bitmap,status)
      if (status.ne.0) call error_$print(status)

c     make the cursor pattern what we just drew
      call gpr_$set_cursor_pattern(curbmdesc,status)
      if (status.ne.0) call error_$print(status)

c     set the cursor origin
      opos(1) = int2(7)
      opos(2) = int2(7)
      call gpr_$set_cursor_origin(opos,status)
      if (status.ne.0) call error_$print(status)

c     draw in white
      call gpr_$set_draw_value(int4(1),status)
      call gpr_$set_fill_value(int4(1),status)

c     white letters on black background
      if (lcolor) then
         call gpr_$set_text_value(int4(1),status)
         call gpr_$set_text_background_value(int4(0),status)
        if (status.ne.0) call error_$print(status)
      else
         call gpr_$set_text_value(int4(0),status)
         call gpr_$set_text_background_value(int4(1),status)
         if (status.ne.0) call error_$print(status)
      endif

c     find out how big characters are
      call gpr_$set_text_font (font_id,status)
      if (status.ne.0) call error_$print(status)
      call gpr_$inq_text_extent ('M', int2(1), csize, status)
      if (status.ne.0) call error_$print(status)


c     where is the cursor?
      call gpr_$inq_cursor(curbmdesc,rastop,active,ipos,opos, status)
      if (status.ne.0) call error_$print(status)
CD     print *,'GSETUP: cursor origin',opos(1),opos(2)
CD     print *,'GSETUP: cursor position',ipos(1),ipos(2)

CD     call gpr_$inq_bitmap_dimensions(curbmdesc,siz,hpid,status)
      if (status.ne.0) call error_$print(status)
CD     print *,'GSETUP: cursor pattern bitmap size',siz(1),siz(2)
      if (status.ne.0) call error_$print(status)

      call gpr_$set_cursor_active (.false., status)

      return
      end

      subroutine gmsg (window,msg)
c
c     prints a message in centered in the designated window
c
      integer*2 x,y,window(4), nch, csize(2),start(2),xend
      integer nblen
      integer*4 status
      character*(*) msg
c
      nch = int2(nblen(msg))

c     size of text
      call gpr_$inq_text_extent (msg, nch, csize, status)

c     offset
      call gpr_$inq_text_offset (msg, nch, start, xend, status)

c     make it black
      call gpr_$set_fill_value (int4(0),status)
      call gpr_$rectangle(window,status)
      call gpr_$set_fill_value (int4(1),status)

      x = window(1) + int2((window(3)-csize(1))/2) + start(1)
      y = window(2) + int2((window(4)-csize(2))/2) + start(2)
      call gpr_$move (x, y, status)
      call gpr_$text (msg,nch,status)

      return
      end

      subroutine gmouse (key,pos,wait)
c     Get input from mouse
c     input:
c         wait .true. to go into a wait state
c     output:
c         pos  x,y of cursor
c         key  mouse button

      INCLUDE '/sys/ins/gpr.ins.ftn'
      INCLUDE '/sys/ins/kbd.ins.ftn'

      character key*1
      logical wait
      integer*2 key_set(16), event_type,pos(2)
      integer*4 status
      logical unobs

      data key_set / 16 * 16#ffff/

      call gpr_$set_cursor_active (.true., status)
      call gpr_$enable_input (gpr_$buttons, key_set, status)
      call gpr_$enable_input (gpr_$locator, key_set, status)

      event_type = gpr_$no_event
      key = '.'

      if (wait) then
         do while (event_type .ne. gpr_$buttons)
            unobs = gpr_$event_wait (event_type, key, pos, status)
            call gpr_$set_cursor_position(pos,status)
         enddo
      else
 10      continue
         key = '.'
         event_type = gpr_$no_event
         unobs = gpr_$cond_event_wait (event_type, key, pos, status)
         call gpr_$set_cursor_position(pos,status)
c        Ignore upstrokes and updates
         if (key .eq. KBD_$M1U) goto 10
         if (key .eq. KBD_$M2U) goto 10
         if (key .eq. KBD_$M3U) goto 10
         if (event_type .eq. gpr_$locator)        goto 10
         if (event_type .eq. gpr_$locator_update) goto 10
         if (event_type .eq. gpr_$locator_stop)   goto 10
      endif

      call gpr_$set_cursor_active (.false., status)

      return
      end


      subroutine gclip (window,lon)

c     turn graphics clipping on or off
c     lon is .true. to turn it on
c     the window is defined as x,y,dx,dy

      logical lon
      integer*2 window(4)
      integer*4 status


      call gpr_$set_clip_window(window,status)
      call gpr_$set_clipping_active(lon,status)

      return
      end


      subroutine gbox (window,igerr)
c     fill in window
      integer*2 window(4)
      integer*4 igerr

      call gpr_$rectangle(window,igerr)

      return
      end

      subroutine  gtsize (string,size,igerr)
c     return x and y extents of text
      character*(*) string
      integer*2 size(2),nch
      integer*4 igerr
      integer nblen

      nch = int2(nblen(string))
      call gpr_$inq_text_extent (string, nch, size, igerr)

      return
      end

      subroutine gline (ix1,iy1,ix2,iy2,igerr)
c     draw a line from (ix1,iy1) to (ix2,iy2)
      integer*2 ix1,iy1,ix2,iy2
      integer*4 igerr

      call gpr_$move (ix1,iy1,igerr)
      call gpr_$line (ix2,iy2,igerr)

      return
      end

      subroutine gdot (ix,iy, iradi, igerr)
c     draw an filled circle of radius iradi at ix,iy
      integer*2 ix,iy,iradi
      integer*4 igerr
      integer*2 ptcent(2)

      ptcent(1) = ix
      ptcent(2) = iy

      call gpr_$circle_filled(ptcent, iradi, igerr)

      return
      end

      subroutine gcirc (ix,iy, iradi, igerr)
c     draw an open circle of radius iradi at ix,iy
      integer*2 ix,iy,iradi
      integer*4 igerr
      integer*2 ptcent(2)

      ptcent(1) = ix
      ptcent(2) = iy

      call gpr_$circle(ptcent, iradi, igerr)

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

      subroutine gmvcur (ipos,igerr)
c     move the cursor to position x=ipos(1),y=ipos(2)

      integer*2 ipos(2)
      integer*4 igerr

      call gpr_$set_cursor_position (ipos, igerr)

      return
      end


      subroutine gend (igerr)
c     terminate graphics package
      integer*4 igerr

      call gpr_$set_cursor_active (.true., igerr)
      call gpr_$terminate (.true.,igerr)

      return
      end


      subroutine gdump (ibmx0,ibmy0,nx,ny,lcolor,idump)
c     dump the screen to the printer
c     origin
      integer*2 ibmx0,ibmy0
c     height and width
      integer*2 bmwdth,bmhght
c     status code = 0 if OK
      integer idump
c     .true. if color
      logical lcolor

      print *,'GDUMP: PRINT not enabled on Apollo! Use FILE.'

      idump = 1

      return
      end







