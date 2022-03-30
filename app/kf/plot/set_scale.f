CTITLE    ................................................................
 
      subroutine set_scale

      implicit none 
c
c     Routine to associate user coordinates with cornor of view area
c
c Include files
c -------------
*                          ! the parameter file
      include 'plot_param.h'
c
*                          ! the common block
      include 'plot_com.h'
c
c Local variables
c ---------------
c xmin, xmax -- min and max values of scales in x direction (may have
c     sign reverse of sign_x is -1.
c ymin, ymax -- min and max values of scales in y direction (see above)
c
      real*4 xmin, xmax, ymin, ymax
      integer*4 i

      real*4 prev_scale(4,2)   ! Previous value of scale which we can
              ! recover (used after fit to overlaid data).
      real*8 prev_ref_valx(2), prev_ref_valy(2)
      save prev_scale, prev_ref_valx, prev_ref_valy
 
c
c.... See if scale already set
      if( scale_set ) return
c
c.... Set the scale; see if we should use default values
      if( use_def_scale ) then
         do i = 1,4
            scale_size(i) = default_scale(i)
         end do
      end if

*     See if scale POP is to be done
      if( pop_scale ) then
*         Adjust for change in reference values
          scale_size(1) = prev_scale(1,2)-(ref_valx-prev_ref_valx(2))
          scale_size(2) = prev_scale(2,2)-(ref_valx-prev_ref_valx(2))
          scale_size(3) = prev_scale(3,2)-(ref_valy-prev_ref_valy(2))
          scale_size(4) = prev_scale(4,2)-(ref_valy-prev_ref_valy(2))
         pop_scale = .false.
      end if

*     Make sure scales make sense
      if( scale_size(2).eq.scale_size(1) ) scale_size(2) = 
     .                              0.01*scale_size(1)+1.d-6
      if( scale_size(4).eq.scale_size(3) ) scale_size(4) = 
     .                              0.01*scale_size(3)+1.d-6

c.... Now set the scale.  If we are plotting with axes reversed
c     change the sign of scales
      if ( sign_x.eq. -1 ) then
         xmin = -scale_size(2)
         xmax = -scale_size(1)
      else
         xmin = scale_size(1)
         xmax = scale_size(2)
      end if
c
      if( sign_y.eq. -1) then
         ymin = -scale_size(4)
         ymax = -scale_size(3)
      else
         ymin = scale_size(3)
         ymax = scale_size(4)
      end if
c
      call jwind(xmin, xmax, ymin, ymax)
c
c.... Indicate scales are set
      scale_set = .true.

*     Save the scale we have just set and move the previous
*     one up the list
      prev_scale(:,2) = prev_scale(:,1)
      prev_ref_valx(2) = prev_ref_valx(1)
      prev_ref_valy(2) = prev_ref_valy(1)
      prev_scale(:,1) = scale_size(:)
      prev_ref_valx(1) = ref_valx
      prev_ref_valy(1) = ref_valy
c
      return
      end
 
