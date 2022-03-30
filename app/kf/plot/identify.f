CTITLE IDENTIFY
 
      subroutine identify( x_data, y_data, pt_data )
 
 
*     Routine to identify a data point or points.  The output can optional
*     be sent to a LU or file. (The results will be appended to the end of
*     of the file.  All points with the current character size are listed
*     to the output device. (NOTE: By making the character size the size of
*     plot in mm then all points can be listed)
*
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*   i,j,k   - Loop counters
*   ierr    - IOSTAT error for file work
*   indx    - Pointer for READ_LINE
*   lenline - Length of output line
*   out_unit    - Unit number for output. If .ne.200 then LU
*   nout        - Number of values found
*   point(2)    - Point type and edit flag for current data point
*   pt_data(2,1)    - Point type and edit flag of data
*   trimlen     - HP function for length of string
 
      integer*4 i, ierr, indx, lenline, out_unit, nout,
     .    point(2), pt_data(2,1), trimlen

*   rtchr   - Charcater from sense charcater

      character*4 rtchr
 
*   x_data(2,1) - X data with sigma
*   y_data(2,1) - Y data with sigma
*   xv,yv       - Virtual coordinates of point
*   xw,yw,zw    - World coordinates of point relative to reference
*               - values
*   xtolv,ytolv - Virtual tolerance of match of point
*   xtolw,ytolw - World tolereance of match to point
 
 
      real*4 x_data(2,1), y_data(2,1), xv,yv, xw,yw,zw, xtolv,ytolv,
     .    xtolw,ytolw
 
*   xvalue(2)   - X value and sigma
*   yvalue(2)   - Y value and sigma
 
      real*8 xvalue(2), yvalue(2)
 
*   out_file    - Name of output file (defaults to '1' if not given)
 
      character*64 out_file
 
*   outline    - Line to be output to the file
 
      character*80 outline 
 
***** First try to get output name
      out_unit = 200
      indx = 9
      if( trimlen(buffer(9:)).gt.0 ) then
          call read_line(buffer,indx, 'CH', ierr, i, out_file)
      else
          ierr = -1
      end if
      if( ierr.ne.0 .or. trimlen(out_file).eq.0 ) then
          out_file = '6'
      end if

      call open_lu(out_unit, out_file, ierr, 'append')
      call report_error('IOSTAT',ierr,'append', out_file, 0,
     .                  'IDENTIFY/PLOT')
      if( ierr.ne.0 ) out_unit = 6
 
*     Get the tolerance for match in x and y directions using
*     current character size
      xtolv = charsz_x/x_size_mm
      ytolv = charsz_y/(y_size_mm*aspect_ratio)
 
*     Now convert to world coordinates
      xtolw = xtolv*(scale_size(2)-scale_size(1))/
     .              (view_size(2) -view_size(1))
      ytolw = ytolv*(scale_size(4)-scale_size(3))/
     .              (view_size(4) -view_size(3))
 
*     Now cursor point and keep looping while rtchr is 'R'
 
*                     !  'R'
      rtchr = 'R'
*                                                  ! 'R' or 'r'
      do while ( rtchr(1:1).eq.'R' .or.  rtchr(1:1).eq.'r' )
 
*         Get cursor position (in virtual coordinates)
          call jwloc( 1, 1, rtchr, xv, yv)
*         Convert to world coordinates
          call jvtow( xv,yv, xw,yw,zw )
 
*****     Now loop over data and see who is close
*                         ! Number of values output
          nout = 0
          do i = 1, num_data
              if( abs(x_data(1,i)-xw).le.xtolw .and.
*                                                      ! We found a point
     .            abs(y_data(1,i)-yw).le.ytolw ) then
 
                  nout = nout + 1
 
*                 Get the values
 
                  xvalue(1) = ref_valx + x_data(1,i)
                  xvalue(2) = x_data(2,i)
                  yvalue(1) = ref_valy + y_data(1,i)
                  yvalue(2) = y_data(2,i)
 
                  point(1) = pt_data(1,i)
                  point(2) = pt_data(2,i)
 
******            Now output
                  outline = ' '
                  call make_line( xvalue, x_field, outline)
                  call make_line( yvalue, y_field, outline)
 
*                 Now add point data
                  lenline = trimlen(outline)
                  write(outline(lenline+2:),200, iostat=ierr)
     .                    point, point(2)
  200             format(i3,2x,i6,' (',o6,'b)')
 
*                 Now write the line
                  write(out_unit,'(A)',iostat=ierr)
     .                outline(1:trimlen(outline))
                  call report_error('IOSTAT',ierr,'writ',outline,0,
     .                              'IDENTIFY')
*                         ! Point with tolerance of data
              end if
*                         ! Looping over data
          end do
 
          write(*,220) nout
  220     format(1x,i4,' points listed')
 
*                         ! Looping until 'R' not hit
      end do
 
***** Close the output
      if( out_unit.eq.200 ) call fmpclose(out_unit)
      return
      end
 
