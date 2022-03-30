CTITLE    ...............................................................
 
      subroutine report_window(window)
c
c     Routine to report the window for polynomial fitting relative to
c     the lower left hand cornor of the plot.
c
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
 
*                         ! the control common block
      include 'plot_com.h'
c
c
*   ierr        - IOSTAT error
*   outlen      - Length of output line
*   trimlen     - HP function for length of string
 
      integer*4 ierr, outlen, trimlen
 
*   window(4)    - Values of the window to be presented.
*   outwind(4)   - Values to be output
 
      real*4 window(4), outwind(4)
 
*   xlab_format - Format to be used for x labels
*   ylab_format - Format to be used for y labels
 
      character*10 xlab_format, ylab_format
 
*   outline      - Line to be output
 
      character*80 outline
 
      common outline
 
c
c.... See if we have read any data (and hence have valid window)
*                                  ! no scales yet
      if( num_data.eq.0 ) return
 
*     Get the window relative to left hand cornor
      outwind(1) = window(1) - scale_size(1)
      outwind(2) = window(2) - scale_size(1)
      outwind(3) = window(3) - scale_size(3)
      outwind(4) = window(4) - scale_size(3)
 
*     Now get the format needed
      call format_label(outwind(1),outwind(2), 0.d0, xlab_format)
      call format_label(outwind(3),outwind(4), 0.d0, ylab_format)
 
      write(termlu,100)
  100 format(' Polynomial Window (WRT Lower Left Cornor [LLC]'
     .      ,' of plot)')
 
      outline = ' LLC '
      write(outline(6:),xlab_format, iostat=ierr) outwind(1)
      outlen = trimlen(outline)
      write(outline(outlen+2:),ylab_format, iostat=ierr) outwind(3)
      outlen = trimlen(outline)
 
*     See if line will be too long
*                               ! Write this line and set up for next
      if( outlen.gt.40 ) then
          write(termlu,'(a)', iostat=ierr) outline(1:outlen)
          outlen = 0
      end if
 
*     Now add URC values
      outline(outlen+2:) = ' URC '
      write(outline(outlen+6:),xlab_format, iostat=ierr) outwind(2)
      outlen = trimlen(outline)
      write(outline(outlen+2:),ylab_format, iostat=ierr) outwind(4)
 
      outlen = trimlen(outline)
      write( termlu,'(a)', iostat=ierr) outline(1:outlen)
 
***** Thats all
      return
      end
