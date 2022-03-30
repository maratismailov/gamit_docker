CTITLE GET_POLY_WINDOW
 
      subroutine get_poly_window
 
      implicit none

*     This routine will get the window to be used for fitting the
*     poynomial to the data.  There are several ways the data may
*     be entered.
*     ? WINDOW <?>    ! Interactive select of window (Lower Left hand
*                     ! and then upper right hand cornors)
*     ? WINDOW <lower left> <upper right>
*     where lower left is a pair of values for the lower left hand cornor
*          relative the the lower left hand cornor of the plot
*     and   upper right is a pair of values for the upper right hand cornor
*          relative to the lower left hand cornor. (For FIELD type 0 ie.
*          calender date fields the values are entered in days.)
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   ichar       - Character used to end JWLOC (not used)
*   ierr        - IOSTAT error reading data
 
      integer*4 ichar, ierr
 
*   xv1, xv2    - Virtual x coordinates of window
*   yv1, yv2    - Virtual y coordinates of window
*   st_x, st_y, st_z    - World coordinates of lower left cornor of
*               - window
*   end_x, end_y, end_z - World coordinates of upper right cornor of
*               - window
 
      real*4 xv1, xv2, yv1, yv2, st_x, st_y, st_z, end_x, end_y, end_z
 
*   values(4)   - total values for the window positions
 
      real*8 values(4)
 
*     See if the ? option has been given, if so that sense the window
 
*                                         ! Sense the position
      if( index(buffer,'?').gt.0 ) then
 
          write(termlu, 100)
  100     format(' Position cursor to lower left hand cornor'
     .           ' of window')
          call jwloc( 1,1, ichar, xv1, yv1)
          call jvtow(xv1,yv1,  st_x, st_y, st_z)
 
          write(termlu, 110)
  110     format(' Position cursor to upper right hand cornor',
     .           ' of window')
          call jwloc( 1,1, ichar, xv2, yv2)
          call jvtow(xv2,yv2, end_x, end_y, end_z)
 
*                                         ! Read from buffer
      ELSE
 
          read( buffer(9:), * ,iostat=ierr) st_x, st_y, end_x, end_y
          call report_error('IOSTAT',ierr,'decod', buffer, 0,
     .                      'GET_POLY_WINDOW')
          st_x = st_x   + scale_size(1)
          st_y = st_y   + scale_size(3)
          end_x = end_x + scale_size(1)
          end_y = end_y + scale_size(3)
 
      END IF
 
***** Now compute the window
 
      poly_window(1) = st_x
      poly_window(2) = end_x
      poly_window(3) = st_y
      poly_window(4) = end_y
 
***** Tell user
 
      write( termlu,200 ) poly_window(1) - scale_size(1),
     .                    poly_window(3) - scale_size(3),
     .                    poly_window(2) - scale_size(1),
     .                    poly_window(4) - scale_size(3)
 
  200 format(' Polynomial window LF cornor ',2(f12.4,1x),/,
     .       '                   UR cornor ',2(f12.4,1x)   )
 
***** Thats all
      return
      end
   
 
