CTITLE    ...................................................................
 
      subroutine set_view
c
c     Routine the set the view area for the GRAPHICS 1000/II
c     package.
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
c view_arg -- the arguments for the view size
c
      real*4 view_arg(4)
 
c
c scratch common
c --------------
c
      common view_arg
 
c
c.... See if the call to jview has already been made
      if( view_set ) return
c
c.... Compute the arguments for the call to view
      view_arg(1) = 1.0*view_size(1)
      view_arg(2) = 1.0*view_size(2)
      view_arg(3) = aspect_ratio*view_size(3)
      view_arg(4) = aspect_ratio*view_size(4)
c
c.... Set the view size
      call jview(view_arg(1), view_arg(2) ,view_arg(3), view_arg(4))
c
c.... Indicate that view has been set
      view_set = .true.
c
      return
      end
 
