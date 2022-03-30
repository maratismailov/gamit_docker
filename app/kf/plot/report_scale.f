CTITLE    ..................................................................
 
      subroutine report_scales(scale)

      implicit none
c
c     Routine to report the scale of the plots to the user.  The
c     reporting will convert julian dates back to year, month, day
c     hour  minute second if need be.  Either the default scales
c     or the set scales can be reported by this routine
c
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
      include 'plot_com.h'
c
c
c Variables
c ---------
c scale -- the scales (relative to the reference values) for the plot
c     (Order: min x, max x, min y, max y)
c true_scale -- the absolute scale of the plot.
c dates      -- the dates corresponding to bounds
c i,j        -- loop counters
c output_line -- buffer used to form the output line
c lab_format  -- the format to be used with real*8 output
c
 
      real*4 scale(4)
 
      real*8 true_scale(4)
  
      integer*4 dates(5) , i, j
      integer*4 absf3  ! Absolute value of field(3).
      real*8 sectag
 
 
      character*80 output_line
 
 
      character*10 lab_format
 
c
c.... See if we have read any data (and hence have valid scales)
*                                  ! no scales yet
      if( num_data.eq.0 ) return
 
c.... START, compute absoute scale
      true_scale(1) = ref_valx + scale(1)
      true_scale(2) = ref_valx + scale(2)
      true_scale(3) = ref_valy + scale(3)
      true_scale(4) = ref_valy + scale(4)
 
c
c.... Do the x_axis scales
      output_line = ' X axis scale '
 
*                                   ! time axis
* MOD TAH 200331: Added field 3 date type
      if( x_field(1).eq.0 .or. x_field(1).eq.3 ) then
* MOD TAH 140338: Allow input of YY DoY but will be reported as Y M D
          if( x_field(3).lt.0 ) then
              absf3 = abs(x_field(3))+1
          else
              absf3 = x_field(3)
          end if

          do i = 1,2
              call mjd_to_ymdhms(true_scale(i),dates, sectag)

              if( absf3.lt.6 ) then
                  write(output_line(16+(i-1)*20:),100)
     .                (dates(j),j=1, 5)
              else
                  write(output_line(16+(i-1)*27:),100)
     .                (dates(j),j=1, 5), sectag
 100              format(i4,"/",i2,"/",i2,1x,i2,":",i2,1x,F6.2)
              endif

          end do
 
*                                   ! Normal data
      Else
 
          call format_label(scale(1), scale(2), ref_valx, lab_format)
          do i = 1,2
              write(output_line(20*i:),lab_format) true_scale(i)
          end do
c
      end if
c
c.... Write line
      call report(output_line)
c
c.... Do the y_axis scales
      output_line = ' Y axis scale '
 
*                                   ! time axis
      if( y_field(1).eq.0 .or. y_field(1).eq.3 ) then
* MOD TAH 140338: Allow input of YY DoY but will be reported as Y M D
          if( y_field(3).lt.0 ) then
              absf3 = abs(y_field(3))+1
          else
              absf3 = y_field(3)
          end if
          do i = 1,2
              call mjd_to_ymdhms(true_scale(i+2),dates, sectag)

              if( absf3.lt.6 ) then
                  write(output_line(16+(i-1)*20:),100)
     .                (dates(j),j=1, 5)
              else
                  write(output_line(16+(i-1)*27:),100)
     .                (dates(j),j=1, 5), sectag
              endif
          end do
 
*                                   ! Normal data
      Else
 
          call format_label(scale(3), scale(4), ref_valy, lab_format)
          do i = 1,2
              write(output_line(20*i:),lab_format) true_scale(i+2)
          end do
c
      end if
c
c.... Write line
      call report(output_line)
c
c.... Thats all
      return
      end
 
