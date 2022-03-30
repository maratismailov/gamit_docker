ctitle
 
      subroutine report_break(idev)

      implicit none 
 
c
c     routine to report the occurrence of clock breaks to the
c     logical unit idev.
c
*                                  ! the SOLVK parameter file
      include '../includes/kalman_param.h'
c
*                                  ! VALUES block of data file
      include '../includes/obs_values.h'
*                                  ! NAMES block of data file
      include '../includes/obs_names.h'
c
c Variables
c ---------
c idev  -- the output device to which clock breaks are reported.
c date(5) -- Date of the clock break
c sectag   -- Seconds tag of break
 
      integer*4 idev, date(5)
      real*8 sectag 
 
*   i,j     - Loop counters
*   num_out - number of clock break values to output.  This
*           - will either the number found or the maximum
*           - number allowed.  At most, max_clk_brk values
*           - are saved.
 
      integer*4 i, num_out
 
c.... Output the break information, if there are any breaks
      if( num_clk_brk.eq.0 ) return
c
      write(idev,100) num_clk_brk, data_base, version
  100 format(/' There are ',i4,' clock breaks in experiment ',a,
     .   ' Vers. ',i3)
c
c.... Check to see if maximum number of breaks exceeded
      if( num_clk_brk.gt. max_clk_brk ) then
         write(idev,150) num_clk_brk, max_clk_brk
  150    format(' *** WARNING *** ',i2,' clock breaks exceedes maximum',
     .      ' number allowed',/,18x,
     .      ' Only the largest ',i2,' breaks will be processed')
         num_out = max_clk_brk
*                 ! set number of breaks to be output to number observed
      else
         num_out = num_clk_brk
      end if
c
c.... Now write the values which were allowed
      do i = 1, num_out
         call jd_to_ymdhms(clk_brk_epoch(i), date, sectag )
         write(idev,200) site_names(clk_brk_site(i)),  date, sectag,
     .      clk_brk_mag(i)/1000.d0
  200    format(' Site ',a8,' break at ',i4,'/',i2,'/',i2,1x,i2,':',
     .       i2,'.',f10.2, '  Magnitude ', f12.1,' ns')
      end do
c
      return
      end
c .....................................................................
