CTITLE pltsp
 
      subroutine pltsp
c
c     This Segment will fit polynomials to the data and
c     compute rms scatters.  The results are output to users
c     terminal and stored in the poly_labels so that they can
c     be output to the plot itself.
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                          ! the common block
      include 'plot_com.h'
c
*                          ! the ema declaration for pltsl
      include 'plot_ema.h'
c
 
*     Get the polynomial window
      if( pel.eq.34 ) then
          call get_poly_window
      end if
 
*     Only compute polynomial if we have the correct command number
 
*                             ! Fit polynomial
      if( pel.eq.35 ) then
          call poly_fit( ema_data( ix_array), ema_data( iy_array),
     .                   ema_data(ipt_array) )
      end if
 
*                             ! Draw the polynomial
      if( pel.eq.37 ) then
          call pdraw
      end if
 
*                             ! Mark window for polynomial
      if( pel.eq.38 ) then
          call mark_window
      end if
 
*                             ! Identify point
      if( pel.eq.39 ) then
          call identify( ema_data( ix_array), ema_data( iy_array),
     .                   ema_data(ipt_array) )
      end if
 
c.... Thats all for this segment
*                   ! return to main segment
      return
c
      end
 
 
 
