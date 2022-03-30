CTITLE    ..................................................................
 
      subroutine put_label

      implicit none
c
c     Routine to get the label information from the buffer
c     and put the label on the screen.  The user may give the
c     position and orientation of the label or select these
c     using the graphics device.
c
c Include files
c -------------
*                          ! the parameter file
      include 'plot_param.h'
c
*                          ! the common block
      include 'plot_com.h'
c
c
c Local variables
c ---------------
c st_x , st_y -- the start position of the label in world coordinates
c st_z        -- z-coordinate of label (dummy)
c or_x , or_y -- orientation of the label
c pl_x , pl_y -- the vector rotated 90 deg from the orientation
c
      real*4 st_x, st_y, st_z, or_x, or_y, pl_x, pl_y
      real*4 pc_x, pc_y   ! Percent values for label location
 
c
c label -- the label to be written
c ilabel -- integer equivalent to label used by Graphics 1000/II
c xv1 , yv1 -- virtual coordinates of a piont
c xv2 , yv2 -- virtual coordinates of a second point
c inqm -- index for finding question mark
c inqu -- index for finding double quote
c
      character*80 label
 
c
      integer*4 ilabel(40)
 
c
      equivalence (label,ilabel)
c
      real*4 xv1, yv1, xv2, yv2
 
c
 
      integer*4 inqm, inqu
 
c Functions
c ---------
c
      integer*4 i, ierr, ichar, iperc

      logical percent_pos  ! Set true if position is % value
 
c
c Scratch common
c --------------
c new_buffer -- a copy of the buffer for concatination
c
      character*80 new_buffer
 
c
      common new_buffer
 
c
c.... Before reading the command buffer, set 2 occurrances of '?'
c     with -999.0 0.0
      do i = 1,2
c
c....    See where " is so that we dont change a ? inside a quote
         inqu = index(buffer,'"')
*                               ! no " so set to large value
         if( inqu.eq.0 ) then
            inqu = len(buffer)
         end if
c
         inqm = index(buffer,'?')
*                                                 ! yes there is a ? in string
         if( inqm.gt.0 .and. inqm.lt.inqu ) then
            new_buffer = buffer(:inqm-1) // ' -999.0 0.0 ' //
     .         buffer(inqm+1:)
            buffer = new_buffer
         end if
      end do

*     See if there if the first value is % number
      iperc = index(buffer,'%')
* MOD TAH 210629: Added test that % index > 0
      if( iperc.gt.0 .and. iperc.lt.inqu ) then
          buffer(iperc:iperc) = ' '
          read(buffer(9:),*,iostat=ierr,err=1000, end=1000)
     .         pc_x, pc_y, or_x, or_y
          st_x = pc_x/100.0*(scale_size(2)-scale_size(1))
          st_y = pc_y/100.0*(scale_size(4)-scale_size(3))

      else   ! Regular read
          
c....     Get the absolute position from the buffer
          read(buffer(9:),*,iostat=ierr,err=1000, end=1000)
     .        st_x, st_y, or_x, or_y
      endif
c
c.... Now get the label
      call read_label(buffer,label)
c
c.... See if we need to determine the position of the label
*                                                    ! get position from
      if( st_x.eq.-999. .or. st_y.eq. -999. ) then
c                                                      graphics device
*                                          ! get virtual cooridinates
         call jwloc(1,1,ichar, xv1, yv1 )
c                                            point selected
c
*                                            ! convert to world coords.
         call jvtow(xv1,yv1, st_x,st_y,st_z)
c
c....    Referr the start to the lower left hand connor of plot
         st_x = st_x - scale_size(1)
         st_y = st_y - scale_size(3)
c
         write(termlu,'(" Start coordinate of label ",2f10.3)')
     .      st_x, st_y
c
      end if

*     See if label is to be output in map mode
      if( map_mode ) then
          call maptrn( st_y, st_x, xv1, yv1 )
          st_x = xv1 - scale_size(1)
          st_y = yv1 - scale_size(3)
      end if
c
c.... Move the cursor
C     call j2mov(st_x+scale_size(1), st_y+scale_size(3))
      call s2mov(st_x+scale_size(1), st_y+scale_size(3))
c
c.... See if we need to get orientation
      if( or_x.eq.-999. .or. or_y.eq.-999. ) then
c
         call jwloc(1,1,ichar, xv2, yv2 )
c
         or_x = xv2 - xv1
         or_y = yv2 - yv1
c
         write(termlu,'(" Orientation ",2f10.3)') or_x, or_y
c
      end if
c
c.... Determine the vector a 90 deg to or_x,or_y
      pl_x = -or_y
      pl_y =  or_x
c
c.... Set the orientation
C     call jcori(or_x,or_y, 0.0,  pl_x,pl_y, 0.0 )
c
c.... Set the justifivation
C     call jjust(0.0, 0.0)
      call scj(or_x,or_y, 0.0, 0.0)
c
c.... Now write out label
C     call write_label(label)
      call save_label(label)
c
      return
c
c.... Error return from read
 1000 continue
      call report_error('IOSTAT',ierr,'decod',buffer,0,'READ_LABEL')
c
      end
 
