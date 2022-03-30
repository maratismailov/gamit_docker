CTITLE    .................................................................
 
      subroutine get_view(gbuffer)
c
c     Routine to get the view information.  The view size we use
c     deterimines what amounts of the total viewing area will be
c     used for the plot i.e., a view size of 0,1,0,1 will totally
c     fill the screen. The view size is implemented by scaling the
c     aspect ratio of the device being used.
c
c Include files
c -------------
*                          ! the parameter file
      include 'plot_param.h'
c
*                          ! the control file
      include 'plot_com.h'
c
c Variables
c ---------
c gbuffer -- user command buffer
c ierr    -- error flag
c
      character*(*) gbuffer
 
c
      integer*4 ierr
 
c
c Functions
c ---------
c trimlen -- HP utility
c
c
c.... Read the view information from the user buffer
      read(gbuffer(9:),*, iostat=ierr, err=1000 ) view_size
c
c.... Now set the view, set the view_set false to ensure that the
c     view will be set
      view_set = .false.
c
      call set_view
c
c.... Check for error
 1000 continue
      if( ierr.ne.0 ) then
          call report_error('IOSTAT',ierr,'decod',gbuffer,0,
     .        'GET_VIEW')
      end if
c
      return
      end
 
