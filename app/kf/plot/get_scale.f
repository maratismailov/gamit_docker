CTITLE    ..................................................................
 
      subroutine get_scale
c
c     Routine to get the scale information
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                         ! the common block
      include 'plot_com.h'
c
c Local variables
c ---------------
c values -- the absolute values of the scales
c
      real*8 values(4)
 
      integer*4 i, ierr, ichar
*   xv1, xv2    - Virtual x coordinates of window
*   yv1, yv2    - Virtual y coordinates of window
*   st_x, st_y, st_z    - World coordinates of lower left cornor of
*               - window
*   end_x, end_y, end_z - World coordinates of upper right cornor of
*               - window
 
      real*4 xv1, xv2, yv1, yv2, st_x, st_y, st_z, end_x, end_y, end_z

 
c
c.... Get the scale values from the buffer

*     See if we will sense the scales
      if( index(buffer,'?').gt.0 ) then
          write(termlu,100)
 100      format(' Position cursor to lower left hand cornor'
     .          ,' of window')
          call jwloc( 1,1, ichar, xv1, yv1)
          call jvtow(xv1,yv1,  st_x, st_y, st_z)
          scale_size(1) = st_x
          scale_size(3) = st_y
 
          write(termlu, 110)
  110     format(' Position cursor to upper right hand cornor',
     .           ' of window')
          call jwloc( 1,1, ichar, xv2, yv2)
          call jvtow(xv2,yv2, end_x, end_y, end_z)
          scale_size(2) = end_x
          scale_size(4) = end_y
      else
          read(buffer(9:),*,iostat=ierr, err=1000,end=1000) values
c
c....     Now remove the reference values to match the plot
          scale_size(1) = values(1) - ref_valx
          scale_size(2) = values(2) - ref_valx
          scale_size(3) = values(3) - ref_valy
          scale_size(4) = values(4) - ref_valy
      end if
c
c.... Set the use_def_scale false
      use_def_scale = .false.
c
c.... Now set the scales. Force scale setting
      scale_set = .false.
c
      call set_scale
c
*                                 ! Update window as well
      if( reset_scales ) then
          do i = 1,4
              poly_window(i) = scale_size(i)
          end do
      end if
 
      return
c
c.... Error label
 1000 continue
      call report_error('IOSTAT',ierr,'decod',buffer,0,'GET_SCALE')
      return
      end
 
