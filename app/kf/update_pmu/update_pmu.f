 
CTITLE UPDATE_PMU
 
      program update_pmu
 
*     J.L. Davis                   2:20 PM  WED., 27  MAY , 1987
*
C                                                           in
*     Program to update earth-orientation values (PM & UT1)
*     a KalObs file.  Delay and delay rate theoreticals as well
*     as the table values are changed.  Note that only four points
*     are used for the new tabular values.
 
* MOD TAH 880304  Added code to remove a linear trend (weighted by
*     the time difference from mid_epoch) from the file values of
*     polar motion and UT1.  This is particularly important for
*     UT1 which has a high slope. Changes made in GET_NEW_TABLES.
 
* MOD TAH 880317 Introduced new UT1 type.  Bit 2 of type is set if
*     series is not regularized. ie. if short period UT1 terms have
*     not been removed. (GSFC and NGS are usually like this)
*
* MOD TAH 881228 Added output of corrections to a print device so that
*     we can keep a record of the changes.  Also fixed leap second problem
*     when we are crossing a boundary and using UT1-UTC tables.

* MOD TAH 910115: Added check on date at which update_pmu was run so that
*     runs which did not save the pmu_calc values could be detected.  These
*     runs all occurred before 90/11/25

      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
 
*       num_header          - Number of header lines in EOR file
*                           -   (Default = 10)
*   ,   ut1_def             - Definition of values in UT1 column.
*                           -   1=UT1-UTC (default);2=UT1-TAI
*   ,   printer             - LU for print device
 
      integer*4 num_header, ut1_def, printer
 
*       new_ut1_inf(4)      - New UT1 info array
*   ,   new_ut1_pts(max_pts)      - New UT1 tabular values
*   ,   new_wob_inf(3)      - New X-Y PM info array
*   ,   new_wob_pts(2,max_pts)    - New X-Y PM tabular values (mas)
 
      real*8 new_ut1_inf(4), new_ut1_pts(max_pts), new_wob_inf(3),
     .    new_wob_pts(2,max_pts)

*  date_bug - Date of not saving pmu_calc bug fixed.

      integer*4 date_bug(5)

*  jd_bug  - Julian date at which the bug of not saving pmu_calc was
*            fixed
*  sectag  - Seconds tag

      real*8 jd_bug, sectag

 
*       pmu_file            - PM-UT1 file name
*   ,   KalObs_file         - KalObs file name
 
      character*256 pmu_file, KalObs_file

      data date_bug / 1990, 11, 25, 0, 0 /
 
***** Get the runstring parameters
      call get_runstring(pmu_file,KalObs_file,num_header,ut1_def,
     .                   printer)

*     Set the bug fix jd

      sectag = 0.d0
      call ymdhms_to_jd( date_bug, sectag, jd_bug )
 
***** Open the KalObs file, get the header info, and read the "a priori"
*     block
      call get_initial_info(KalObs_file)

* MOD TAH 900920: Check to see the number of data points in the tables.
*     If it is greater than 5 then reset the value.
      if( fut1_inf(3).gt.max_pts .or. fwob_inf(3).gt.max_pts ) then
          write(printer,100) fut1_inf, fwob_inf
 100      format('*** WARNING *** UT1/Wobble tables are too large.',
     .           ' Current INF array ares:',/,
     .           ' UT1:    ',f12.2,1x,f5.2,1x,f2.0,1x,f5.0,
     .           ' Wobble: ',f12.2,1x,f5.2,1x,f2.0,/,
     .           ' Resetting values to aviod problems.')      
          fut1_inf(3) = max_pts
          fwob_inf(3) = max_pts
      end if
 
***** Open the new PM-UT1 file and get the new tabular values
      call get_new_tables(pmu_file,num_header,new_ut1_inf,new_ut1_pts,
     .    new_wob_inf,new_wob_pts,ut1_def)
 
***** Update theoretical values
      call pmu_update_theo(new_ut1_inf,new_ut1_pts,
     .    new_wob_inf,new_wob_pts, jd_bug)
 
***** Update values for PM and UT1 at mid_epoch
      call update_mid_pmu(new_ut1_inf,new_ut1_pts,
     .    new_wob_inf,new_wob_pts,printer, KalObs_file)
 
***** Update apriori info in KalObs file
      call update_apriori_info(pmu_file, ut1_def)
 
***** Close all files
      call finish
 
      end
 
CTITLE GET_RUNSTRING
 
      subroutine get_runstring(pmu_file,KalObs_file,num_header,ut1_def,
     .                         printer)
 
*     J.L. Davis                   4:42 PM  WED., 27  MAY , 1987
*
*     Get the runstring parameters and call proper_runstring if
*     necessary
 
 
*       DecimalToInt        - Character to integer conversion
*   ,   error1              - Error flag
*   ,   error2              - Error flag
*   ,   error3              - Error flag
*   ,   len_KalObs_file     - Length of string descriptor
*   ,   len_pmu_file        -   "    "     "        "
*   ,   len_num_header      - Length of string with integer value
*   ,   len_ut1_def         - Length of string with integer value
*   ,   len_printer         - Length of printer string
*   ,   num_header          - Number of header records
*                           - (default=10)
*   ,   printer             - Print device
*   ,   rcpar               - HP runstring utility
*   ,   ut1_def             - Flag for definition of UT1
 
      integer*4 DecimalToInt, error1, error2, error3, len_KalObs_file,
     .    len_pmu_file, len_num_header, len_ut1_def, len_printer,
     .    num_header, printer, rcpar, ut1_def
 
*       pmu_file            - File with PM-UT table
*   ,   KalObs_file         - KalObs file name
 
      character*(*) pmu_file, KalObs_file
 
*       num_header_str*5    - Input string with num_header
*   ,   help_file*16        - Help file descriptor
 
      character num_header_str*5, help_file*20, ut1_def_str*5,
     .    printer_str*128
 
      data help_file / 'update_pmu.hlp' /
 
***** Initialize input strings
      pmu_file       = ' '
      KalObs_file    = ' '
      num_header_str = ' '
      ut1_def_str    = ' '
      printer_str    = ' '
 
***** Read in the runstrings
      len_pmu_file    = rcpar(1,pmu_file      )
      len_KalObs_file = rcpar(2,KalObs_file   )
      len_num_header  = rcpar(3,num_header_str)
      len_ut1_def     = rcpar(4,ut1_def_str   )
      len_printer     = rcpar(5,printer_str   )
 
***** Initialize the numeric quantities with default values
      num_header =  3
      ut1_def    =  1
      printer    =  6
 
***** Convert character "num_header" and ut1_def to integer
      if (len_num_header .ne. 0)
     .    num_header = DecimalToInt(num_header_str,error1)
      if (len_ut1_def    .ne. 0)
     .    ut1_def    = DecimalToInt(ut1_def   _str,error2)
      if (len_printer    .ne. 0) then
          printer = 201
          call open_lu(printer, printer_str, error3, 'append')
          if ( error3.ne.0 ) printer = 6
      end if
 
***** Do we need to call proper_runstring?
      if (error1 .ne. 0)          call proper_runstring(help_file,
     .                                 'UPDATE_PMU',1)
      if (error2 .ne. 0)          call proper_runstring(help_file,
     .                                 'UPDATE_PMU',1)
      if (error3 .ne. 0)          call proper_runstring(help_file,
     .                                 'UPDATE_PMU',1)
      if (len_pmu_file .eq. 0)    call proper_runstring(help_file,
     .                                 'UPDATE_PMU',1)
      if (len_KalObs_file .eq. 0) call proper_runstring(help_file,
     .                                 'UPDATE_PMU',1)
 
      end
 
CTITLE GET_INITIAL_INFO
 
      subroutine get_initial_info(KalObs_file)
 
*     J.L. Davis                   4:58 PM  WED., 27  MAY , 1987
*
*     Do some initial opening and reading of the KalObs file
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       error                   - Error flag
 
      integer*4 error
 
*       KalObs_file             - KalObs file descriptor
 
      character*(*) KalObs_file
 
 
      character*16 prog
 
      data prog / 'GET_INITIAL_INFO' /
 
***** Open KalObs file
      call open_KalObs(KalObs_file,error)
 
***** Error?
      call report_error('FMGR',error,'open',KalObs_file,1,prog)
 
***** Read the KalObs header
      call rw_KalObs_header('R',error)
 
***** Error?
      call report_error('FMGR',error,'read','KalObs header',1,prog)
 
***** Output the header
      call out_header(6,KalObs_file)
 
      end
 
CTITLE GET_NEW_TABLES
 
      subroutine get_new_tables(pmu_file,num_header,new_ut1_inf,
     .    new_ut1_pts,new_wob_inf,new_wob_pts,ut1_def)
 
*     J.L. Davis                   5:15 PM  WED., 27  MAY , 1987
*
*
*     Routine to get new table values from the PM-UT1 file.  The
*     following procedure is used:
*
*     1.  The PM-UT1 file containing the new values is opened (see
*         subroutine SCAN_PMU_FILE for format of this file)
*
*     2.  For each of PMX, PMY, and UT1-TAI, the following pro-
*         cedure is performed:
*
*     2.a The PM-UT1 file is scanned using subroutine SCAN_PMU_FILE.
*         All values of the requested component of PM-UT1 are returned
*         that have epoch-tags within NUM_SPAN*FWOB_INF(2) days of
*         MID_EPOCH (NUM_SPAN*FUT1_INF(2) for UT1), up to a maximum
*         number MAX_RETURN_VALUES.  For example, if NUM_SPAN=4, and
*         FWOB_INF(2)=5 (indicating a spacing in the KalObs table
*         of 5 days), then all values from the PM-UT1 file with
*         epoch tags within 20 days of MID_EPOCH will be returned.
*
* MOD TAH 880304 2.a.1 Remove linear trend first before smoothing
*         with Gaussian filter.
*
*     2.b The returned values, which may be spaced unevenly, are
*         now interpolated onto an evenly spaced table using
*         weighted Gaussian interpolation.  These new tables differ
*         from the original tables in the data base and in the
*         KalObs file in that there will only be four points
*         in the new table.  An "INF" array is created for these points.
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       max_return_val              - See above
*   ,   num_span                    - See above
 
      integer*4 max_return_val, num_span
 
      parameter (max_return_val = 55)
      parameter (num_span       =  max_pts)
 
*       error                       - Error flag
*   ,   i, j                        - Loop counter
*   ,   num_header                  - Number of header records
*                                   -   at top of PM-UT1 file
*   ,   num_return                  - Number of values returned
*                                   -   by file scan
*   ,   ut1_def                     - UT1 definition;  see help file
*       num_grid                    - Number of points in the grid
*                                     returned.  Should be fut1_inf(3)
 
      integer*4 error, i, j, num_header, num_return, ut1_def, num_grid
 
*       day_limit                   - Range (days) for table return
*   ,   new_ut1_inf(4)              - New info array for UT1
*   ,   new_ut1_pts(max_pts)              - New array of points for UT1
*   ,   new_wob_inf(3)              - New info array for X-Y PM
*   ,   new_wob_pts(2,max_pts)            - New array of points for X-Y PM
*   ,   return_dates(max_return_val)    - Array of Julian dates
*   ,   return_values(max_return_val)   - Table of X, Y, or UT1
*                                       -   values
*   ,   sigmas(max_return_val)      - Sigmas corresponding to
*                                   -   RETURN_VALUES
*   ,   temp_pts(max_pts)                 - Temporary storage of points
*       trend(2)                    - Offset and rate removed
*                                   - before smoothing and then
*                                   - added back into new tabular
*                                   - values
*       tut1_inf(4), twob_inf(3)  - Temporay values of fut1_inf and fwob_inf 
*                      in case we need to change them,
 
      real*8 day_limit, new_ut1_inf(4), new_ut1_pts(max_pts), 
     .    new_wob_inf(3), new_wob_pts(2,max_pts),
     .    return_dates(max_return_val),
     .    return_values(max_return_val), sigmas(max_return_val),
     .    temp_pts(max_pts), tut1_inf(4), twob_inf(3)
 
*       pmu_file                    - File with PM-UT1 values
 
      character*(*) pmu_file
 
*       prog*14                     - Routine name
 
      character prog*14
 
      data prog / 'GET_NEW_TABLES' /
 
***** Open PM-UT1 file
      open (105, file = pmu_file, iostat = error, status = 'old')
 
***** Error?
      call report_error('IOSTAT',error,'open',pmu_file,1,prog)

*     Copy the fwob_inf and fut1_inf arrays
      do i = 1,3
         twob_inf(i) = fwob_inf(i)
         tut1_inf(i) = fut1_inf(i)
      end do
      tut1_inf(4) = fut1_inf(4)
      num_grid = nint(fut1_inf(3))
 
***** How many days should be included for X-Y PM?
      day_limit = twob_inf(2) * num_span
 
***** Loop over X and Y wobble
      do i = 1, 2
 
*****     Get nearest points for X/Y wobble.  Units are mas.
          call scan_pmu_file(105,day_limit,max_return_val,i,num_header,
     .        return_values,return_dates,sigmas,num_return)
 
*****     Interpolate onto the new table. Units are mas.
          call grid(return_values,return_dates,sigmas,num_return,
     .        twob_inf,num_grid,temp_pts)
 
*****     Copy temporary array into new wobble array, changing units from
*         arcseconds to milliarcseconds
          do j = 1, num_grid
              new_wob_pts(i,j) = 1.0D+03 *  temp_pts(j)
          end do
 
      end do
 
***** Make new "information" array for wobble
      call get_new_wob_inf(twob_inf,new_wob_inf)
 
***** How many days should be included for UT1?
      day_limit = tut1_inf(2) * num_span
 
***** Get nearest points for UT1.  Units are seconds.
      call scan_pmu_file(105,day_limit,max_return_val,3,num_header,
     .    return_values,return_dates,sigmas,num_return)
 
***** Change UT1 values from whatever they were to UT1-TAI
      call change_ut1(return_values,return_dates, num_return,ut1_def)
 
***** Interpolate onto grid.  Units are seconds.
      call grid(return_values,return_dates,sigmas,num_return,tut1_inf,
     .          num_grid, new_ut1_pts)
 
***** Make new "information array for UT1
      call get_new_ut1_inf(tut1_inf,new_ut1_inf)
 
***** Close PM-UT1 file
      close (105)
 
      end
 
CTITLE SCAN_PMU_FILE
 
      subroutine scan_pmu_file(lu,day_limit,max_return_val,icol,
     .    num_header,return_values,return_dates,sigmas,num_return)
 
*     J.L. Davis                   3:24 PM  THU., 28  MAY , 1987
*
*     Routine to read through the PM-UT1 file and return all
*     values within DAY_LIMIT days of mid_epoch
*
*     Format of file is:
*
*     HEADER  |
*       .     |
*       .     } NUM_HEADER header lines (default=10)
*       .     |
*     HEADER  |
*     YR MO DY HR MN   X   SX    Y    SY    UT1-UTC SUT1
*                          .
*                          .
*                          .
*                        (EOF)
*
*     where YR is the year (e.g. 1984), MO is the month (1-12),
*     DY is the day of month, HR is the hour of day (0-23), MN is
*     the minute of the hour, X is the X-pole (arcsec), SX is the
*     X-pole uncertainty (arcsec), Y is the Y-pole (arcsec), SY is
*     the Y-pole uncertainty (arcsec), UT1 is in seconds, and
*     SUT1 is the uncertainty of UT1 (seconds).  Note that the definition
*     of the UT1 column is variable (e.g., this may be a column of
*     UT1-UTC values, or UT1-TAI values), but the units are still
*     seconds.  Any line with nonnumeric values in the values slots
*     are ignored.  Any line with nonnumeric values in the sigmas
*     slot will have an arbitrary sigma value of 2 mas or .13 msec
*     assigned.  For example, suppose we have the entry
*
*     1981 10 01 00 00 -.0471 .0026 .2007 ----- .22315 .00019
*
*     The X-pole and UT1 are read normally.  However, the Y-pole
*     will have a sigma of 2 mas.
 
*     MOD JLD 880713 Checked that sigma read is greater than zero.
*                    If not, change to 2 mas.
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       use_point       - .true. if all READ_LINE calls go OK
 
      logical use_point
 
*       date(5)         - YMDHM
*   ,   dum             - Dummy for READ_LINE
*   ,   i               - Loop counter
*   ,   icol            - 1 for X, 2 for Y, 3 for UT1
*   ,   indx            - Index of line for READ_LINE
*   ,   line_err        - IOSTAT returned from READ_LINE
*   ,   lu              - Lu of PM-UT1 file
*   ,   max_return_val  - Maximum number of points to return
*   ,   num_header      - Number of header records to be skipped
*   ,   num_return      - Number of values in returned arrays
*   ,   rerr            - Error flag for file read
 
      integer*4 date(5), dum, i, icol, indx, line_err, lu,
     .    max_return_val, num_header, num_return, rerr
 
*       day_limit       - Range for accepting point (days)
*   ,   jul_date        - Julian date
*   ,   return_dates(1) - Array of JDS's
*   ,   return_values(1)  - Returned array of values
*   ,   sig             - Sigma returned from READ_LINE
*   ,   sigmas(1)       - Sigmas from file
*   ,   val             - Value returned from READ_LINE
 
      real*8 day_limit, jul_date, return_dates(1), return_values(1),
     .    sig, sigmas(1), val
 
*       cdum            - Dummy for READ_LINE
 
      character*2 cdum
 
*       line            - Line read from file
 
      character*80 line
 
***** Format for read of line
  100 format(A)
 
***** Rewind PM-UT1 file
      rewind (lu)
 
***** Skip over headers
      do i = 1, num_header
          read (lu,100) line
      end do
 
***** Initialize values for read through file
      num_return = 0
      rerr       = 0
 
***** Loop until a read error occurs
      do while (rerr .eq. 0)
 
*****     Read a line
          read (lu, 100, iostat = rerr) line
 
*****     Continue if no read error
          if (rerr .eq. 0) then
 
*****         Initially, we will use this point
              use_point = .true.
 
*****         Look for the date at the start of the line
              indx = 1
 
*****         Loop over the first five values in the line:  the date
              do i = 1, 5
 
*****             Get the next value
                  call read_line(line,indx,'I4',line_err,date(i),cdum)
 
*****             If there was an error, do not use this point
                  if (line_err .ne. 0) use_point = .false.
 
              end do
 
*****         Skip over columns we will not use
              do i = 1, icol - 1
 
*****             Read a value
                  call read_line(line,indx,'CH',line_err,dum,cdum)
 
*****             Read a sigma
                  call read_line(line,indx,'CH',line_err,dum,cdum)
 
              end do
 
*****         Now read the value
              call read_line(line,indx,'R8',line_err,val,cdum)
 
*****         Error?
              if (line_err .ne. 0) use_point = .false.
 
*****         Now read the sigma
              call read_line(line,indx,'R8',line_err,sig,cdum)
 
*****         Error on sigma--still use point, but set sig to 2 mas
*             NOTE:  Since USE_POINT is not set, point is still not used.
*                    JLD 880713
              if (line_err .ne. 0) then
 
*****             Set to 2 mas
                  sig = 0.002
 
*****             Convert to msec if UT1
                  if (icol .eq. 3) sig = sig / 15.0
 
              end if
 
*****         MOD JLD 880713:  If the sigma is zero, change to 2 mas
              if (sig .le. 1.0d-08) then
 
*****             Set to 2 mas
                  sig = 0.002
 
*****             Convert to msec if UT1
                  if (icol .eq. 3) sig = sig / 15.0
 
              end if
 
*****         Has the reading of this point gone OK?
              if (use_point) then
 
*****             Determine the JD of this point
                  call ymdhms_to_jd(date,0.0D0,jul_date)
 
*****             Is this point within DAY_LIMIT of MID_EPOCH?
                  if (abs(jul_date-mid_epoch) .le. day_limit) then
 
*****                 Add another point to the array
                      num_return = num_return + 1
 
*****                 But check that we have room
                      if (num_return .gt. max_return_val) then
 
*****                     If we do not have room, print a message
                          write(*,150)
  150                     format(' *** Max number of return values ',
     .                           ' exceeded.')
 
*****                     Do not add to table
                          num_return = num_return - 1
 
                      else
 
*****                     We have room
                          return_values(num_return) = val
                          return_dates(num_return)  = jul_date
                          sigmas(num_return)        = sig
 
*                                         ! Room enough in table
                      end if
 
*                                         ! Near enough to MID_EPOCH
                  end if
 
*                                         ! OK readin line
              end if
 
*                                         ! OK readin file
          end if
 
*                                         ! Loop over file
      end do
 
      end
 
CTITLE GRID
 
      subroutine grid(values,dates,sigmas,num_values,grid_info,num_grid,
     .                grid_array)
 
*     J.L. Davis                   2:55 PM  FRI., 29  MAY , 1987
*
*     Routine to take the (presumed) unevenly spaced values read
*     from the PM-UT1 file and grid them onto and evenly spaced
*     array.  The beginning epoch for the array is stored in
*     GRID_INFO(1) and the separation of points in days is stored
*     in GRID_INFO(2).  The number of points to be gridded is
*     NUM_GRID.
 
 
*       i,j                     - Loop counter
*   ,   num_grid                - Number of points to be gridded
*   ,   num_values              - Number of values read from PM-UT1
*                               -   file
      integer*4 i,j, num_grid, num_values
 
*       dates(1)                - Array of JD's
*   ,   FWHM                    - FWHM used to smooth data
*   ,   grid_array(1)           - The output array of gridded points
*   ,   grid_info(2)            - The "info" array from the database
*   ,   sigmas(1)               - The sigmas read from PM-UT1 file
*   ,   sig_dum                 - Dummy sigma for WT_GAUSS_FILTER
*   ,   time                    - epoch of smoothed value
*       trend(2)                - Trend removed locally about
*                               - grid point
*   ,   values(1)               - The points read from PM-UT1 file
 
      real*8 dates(1), FWHM, grid_array(1), grid_info(2), sigmas(1),
     .    sig_dum, time, trend(2), values(1)
 
***** Sort the points into ascending time order
      call order_array(values,dates,sigmas,num_values)
 
***** Determine the appropriate FWHM (days)
      call wt_filter_fwhm(dates,sigmas,num_values,FWHM)

* MOD TAH 901128:  shortned the full width at half max value

      fwhm = fwhm / 2.d0
 
***** Report
      write(*,100) FWHM
  100 format(' Using a FWHM of ',F4.1,' days.')
 
***** Loop over points in table
      do i = 1, num_grid
 
*****     Find the time for this point
          time = grid_info(1) + dble(i-1) * grid_info(2)
 
*         Remove the local trend of the data
          call remove_trend( values, dates, sigmas, num_values, time,
     .                       time, trend )
 
*****     Find the smooth value at this time
          call wt_gauss_filter(dates,values,sigmas,
     .        num_values,time,grid_array(i),sig_dum,FWHM)
 
*         Add back in the trend which was removed (only offset needed
*         since we are at the epoch of the polynomial)
          grid_array(i) = grid_array(i) + trend(1)
 
*         Reform the values array
          do j = 1, num_values
              values(j) = values(j) + trend(1) +
     .                    trend(2)*(dates(j)-time)
          end do
      end do
 
      end
 
CTITLE ORDER_ARRAY
 
      subroutine order_array(values,dates,sigmas,num_values)
 
*     J.L. Davis                   3:30 PM  FRI., 29  MAY , 1987
*
*     Routine to time sort the arrays read from the PM-UT1 file
 
 
*       i, j            - Loop index
*   ,   num_values      - Number of values read
 
      integer*4 i, j, num_values
 
*       dates(1)        - Dates corresponding to the values
*   ,   sigmas(1)       - Sigmas corresponding to the values
*   ,   temp            - Temporary storage for exchange
*   ,   values(1)     - Times and values from PM-UT1 file
 
      real*8 dates(1), sigmas(1), temp, values(1)
 
***** Bubble sort
      do i = 1, num_values - 1
 
          do j = 1, num_values - i
 
              if (dates(j) .gt. dates(j+1)) then
 
*****             Exchange times, values, and sigmas
                  temp       = dates(j)
                  dates(j)   = dates(j+1)
                  dates(j+1) = temp
 
                  temp        = values(j)
                  values(j)   = values(j+1)
                  values(j+1) = temp
 
                  temp        = sigmas(j)
                  sigmas(j)   = sigmas(j+1)
                  sigmas(j+1) = temp
 
              end if
 
          end do
 
      end do
 
      end
 
CTITLE GET_NEW_WOB_INF
 
      subroutine get_new_wob_inf(fwob_inf,new_wob_inf)
 
*     J.L. Davis                   3:43 PM  FRI., 29  MAY , 1987
*
*     Routine to create new wobble information array
*
 
*       fwob_inf(3)             - Wobble info array from KalObs
*   ,   new_wob_inf(3)          - The new wobble info array
 
      real*8 fwob_inf(3), new_wob_inf(3)
 
***** Assign values
*                                     ! Start JD
      new_wob_inf(1) = fwob_inf(1)
*                                     ! Table spacing (days)
      new_wob_inf(2) = fwob_inf(2)
*                                     ! Number of table values
      new_wob_inf(3) = fwob_inf(3)
 
      end
 
CTITLE GET_NEW_UT1_INF
 
      subroutine get_new_ut1_inf(fut1_inf,new_ut1_inf)
 
*     J.L. Davis                   3:43 PM  FRI., 29  MAY , 1987
*
*     Routine to create new UT1-TAI information array
*
 
*       fut1_inf(4)             - UT1 info array from KalObs
*   ,   new_ut1_inf(4)          - The new UT1 info array
 
      real*8 fut1_inf(4), new_ut1_inf(4)
 
***** Assign values
*                                     ! Start JD
      new_ut1_inf(1) = fut1_inf(1)
*                                     ! Table spacing (days)
      new_ut1_inf(2) = fut1_inf(2)
*                                     ! Number of table values
      new_ut1_inf(3) = fut1_inf(3)
*                                     ! Units of table (sec)
      new_ut1_inf(4) = 1
 
      end
 
CTITLE FINISH
 
      subroutine finish
 
*     J.L. Davis                   4:06 PM  MON.,  1  JUNE, 1987
*
*     Routine to close files
*
 
 
      integer*4 error
 
      call close_KalObs(error)
 
      end
 
CTITLE CHANGE_UT1
 
      subroutine change_ut1(UT1_values,return_dates, num_values,ut1_def)
 
*     J.L. Davis                   4:25 PM  MON.,  1  JUNE, 1987
*
*     Routine to change values read in UT1 column of PM-UT1 file
*     into UT1-TAI.  The flag UT1_DEF tells the program the
*     definition of the values of this column.  The allowed
*     values are:
*
*               Bit
*     UT1_DEF = Bit 1 off Values are UT1-UTC
*     UT1_DEF = Bit 1 on  Values are UT1-TAI
*               Bit 2 off Values are regularized
*               Bit 2 on  Short period terms should be removed
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_apr.h'
 
*       i                       - Loop counter
*   ,   num_values              - Number of values in UT1 array
*   ,   num_center              - Value at center of epoch
*   ,   ut1_def                 - Definition of values; see above
 
      integer*4 i, num_values, num_center, ut1_def
 
*       dut1                    - Change in UT1 values for tides (tsec)
*       ut1_values(1)           - Values of UT1
*       return_dates(1)         - JD's of values
*       ut1_center              - Value of ut1 at center.  Needed
*                               - to correct for leap seconds
      real*8 dut1, ut1_values(1), return_dates(1), ut1_center
 
*       kbit                    - Bit checking function
 
      logical kbit
 
***** UT1_def = 1:  Values are UT1-UTC
*     Check size of ut1_value.  This will allow us to see if it is
*     UT1-UTC or UT1-TAI
      if( abs(ut1_values(1)).lt.2.d0 ) then
          call sbit( ut1_def, 1, 0 )
          write(*,100)
 100      format(' UT1 definition for UT1-UTC being used')
      else
          call sbit( ut1_def, 1, 1 )
          write(*,110)
 110      format(' UT1 definition for UT1-TAI being used')
      end if
 
*                                         ! Convert to UT1-AT
      if (.not. kbit(ut1_def,1) ) then
 
*****     Loop over values read from PM-UT1 file
          num_center = num_values/2 + 1
          ut1_center = ut1_values(num_center) - tai_utc
          do i = 1, num_values
 
*****         [UT1-TAI] = [UT1-UTC] - [TAI-UTC]
              ut1_values(i) = ut1_values(i) - tai_utc
              if( (ut1_values(i)-ut1_center).lt.-0.5 ) then
                  ut1_values(i) = ut1_values(i) + 1.0d0
              end if
              if( (ut1_values(i)-ut1_center).gt. 0.5 ) then
                  ut1_values(i) = ut1_values(i) - 1.d0
              end if
          end do
 
      end if
 
***** UT1_def = 2;  Values are UT1-TAI
*                                     ! Remove short period terms
      if (kbit(Ut1_def,2) ) then
 
*****     Compute short period ut1 corrections
          do i = 1, num_values
              call short_period_ut1( return_dates(i), dut1)
              ut1_values(i) = ut1_values(i) - dut1
          end do
 
      end if
 
      end
 
CTITLE PMU_UPDATE_THEO
 
      subroutine pmu_update_theo(new_ut1_inf,new_ut1_pts,
     .    new_wob_inf,new_wob_pts, jd_bug)
 
*     J.L. Davis                   2:10 PM  TUE.,  2  JUNE, 1987
*
*     Routine to loop through KalObs file and update the theoretical
*     values for the delays using the new PM-UT1 tables and the
*     partials.
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
 
*       error                   - Error flag
*   ,   i                       - Loop counter
 
      integer*4 error, i, j

*   kbit    - Logical function to check bits

      logical kbit
 
*       delta_pmx               - Change to X pole
*   ,   delta_pmy               - Change to Y pole
*   ,   delta_ut1               - Change to UT1-TAI
*   ,   delta_tau               - Change to group delay
*   ,   delta_tau_dot           - Change to phase-delay rate
*   ,   new_ut1_inf(1)          - New UT1 info array
*   ,   new_ut1_pts(1)          - New UT1-TAI array (seconds)
*   ,   new_wob(2,2)            - New wobble interpolated values
*   ,   new_wob_inf(1)          - New wobble info array
*   ,   new_wob_pts(1)          - New wobble array (mas)
*   ,   new_ut1(2)              - New UT1 interpolated values
*   ,   old_wob(2,2)            - "Old" wobble values determined
*                               -   from interpolation using old
*                               -   table.
*   ,   old_ut1(2)              - "Old" UT1 value
*       dut1                    - Tidal contribution to UT1 (tsec) (needed if
*                                 ut1_tide is 0.  We will set value
*                                 to -1 when we are finished)
*       calc_adj                - Tidal contribution to UT1 needed when comparing
*                                 the full UT1 in PMU_CALC with the interpolated
*                                 values from the tables (mas) 
*       jd_bug                  - Date at which bug in not saving pmu fixed
 
      real*8 delta_pmx, delta_pmy, delta_ut1, delta_tau, delta_tau_dot,
     .    new_ut1_inf(1), new_ut1_pts(1), new_wob(2,2), new_wob_inf(1),
     .    new_wob_pts(1), new_ut1(2), old_wob(2,2), old_ut1(2), dut1,
     .    calc_adj, jd_bug, tol

*       jd_last                 - Last julian date at which we ran update_pmu

      real*8 jd_last

*       use_pmu_calc            - Indicates that we should use pmu_calc values

      logical use_pmu_calc

***** Get the last time we ran update_pmu.  Jd_last is returned -1 if
*     update_pmu has not been run before.

      call get_last_pmu( user_pmu_dsc, jd_last )

*     Set the logicals for use_pmu_calc.  We should use for all calc_versions
*     greater than 7.05 unless pmu_update has been run before the bug in 
*     saving pmu_calc was fixed.
      use_pmu_calc = .false.
      if( Calc_ver.gt.7.05 ) use_pmu_calc = .true.
      if( jd_last.gt.0.d0 .and. jd_last.lt.jd_bug ) then
          use_pmu_calc = .false.
      end if

      if( use_pmu_calc ) then
          write(*,'('' Using CALC supplied PMU values'')')
      else
          write(*,'('' NOT USING CALC SUPPLIED PMU VALUES'')')
          write(*,'('' Bug JD, Last JD are '',2f12.2)') jd_bug,
     .                 jd_last
      end if

****  Set tolerance on output
      tol = 0.01d0
 
***** Loop over all points in KalObs
      do i = 1, num_obs
 
*****     Read in the record for this observation
          call rw_KalObs_block('R','data',site,error,i)
 
*****     Check for error
          call corrupt_KalObs('IOSTAT',error,i,'read',1)
 
*****     Get the "old" polar motion and UT1 values for this epoch by
*         interpolating using the "old" tables.  Table values are in
*         milliarcseconds: so are interpolated values
          call Lagrange_intp(fwob_inf,fwob_pts,epoch,old_wob,2)
          call Lagrange_intp(fut1_inf,fut1_pts,epoch,old_ut1,1)
 
*****     Get the "new" values corresponding to those above
          call Lagrange_intp(new_wob_inf,new_wob_pts,epoch,new_wob,2)
          call Lagrange_intp(new_ut1_inf,new_ut1_pts,epoch,new_ut1,1)
 
*****     Calculate change in PM and UT1, all in mas
*         Get the short period UT1 correction first.
          call short_period_ut1( epoch, dut1)
          if( .not.use_pmu_calc) then
              delta_pmx = new_wob(1,1) - old_wob(1,1)
              delta_pmy = new_wob(2,1) - old_wob(2,1)
              delta_ut1 = (new_ut1(1) - old_ut1(1)) * 15.0D+03
          else
              delta_pmx = new_wob(1,1) - pmu_calc(1)
              delta_pmy = new_wob(2,1) - pmu_calc(2)
*             Remove the short period UT1 from the CALC total version
              delta_ut1 = new_ut1(1)*15.d3 - 
     .                                  (pmu_calc(3)-dut1*15.d3) 
          end if

*****     Calculate change in group delay due to these changes.  No unit
*         changes for PM are necesary since derivative is in psec/mas and
*         PM values are in mas.  For UT1, partial derivative is in
*         psec/mas and values are in seconds, so a unit change is
*         necessary
*                                                 ! X
          delta_tau = delta_pmx * pmu_part(1,1)
*                                                 ! Y
     .              + delta_pmy * pmu_part(2,1)
*                                                 ! UT1
     .              + delta_ut1 * pmu_part(3,1)
 
*****     Change theoretical value of group delay, phase delay, and
*         single-band delay.
*                                                              ! Group
          db_theoretical(1) = db_theoretical(1) + delta_tau
*                                                              ! Phase
          db_theoretical(2) = db_theoretical(2) + delta_tau
*                                                              ! SB delay
          db_theoretical(3) = db_theoretical(3) + delta_tau
 
*****     Calculate change to phase delay rate, assuming that the
*         contribution due to the PM-UT1 rate is zero
          delta_tau_dot = delta_pmx * pmu_part(1,2)
     .                  + delta_pmy * pmu_part(2,2)
     .                  + delta_ut1 * pmu_part(3,2)
 
*****     Add this contribution to the theoretical phase_delay rate
          db_theoretical(4) = db_theoretical(4) + delta_tau_dot
 
*****     See if we are developing errors relatiave to calc values
          if( Calc_ver.gt.7.05 ) then
             do j = 1,2
                 if( abs(old_wob(j,1)-pmu_calc(j)).gt.0.05 .and.
     .               use_pmu_calc  ) then
                     write(*,250) epoch,j, old_wob(j,1), pmu_calc(j)
 250                 format(' *** ERROR At epoch ',f12.4,' Wobble ',i2,
     .                      ' Differs from calc (O,C): ',2f8.3)
                     tol = 1.5d0*tol

*                    If we are observation #1 then kill now before any
*                    thing is saved.
                     if( i.eq.1 .and.
     .                   abs(old_wob(j,1)-pmu_calc(j)).gt.tol*100) then
                         stop ' Update_pmu_tab: Kill on 1st observation'
                     end if
                 end if
              end do

*             If the original tables in calc did not have short period   
*             terms, then remove from pmu_calc before comparing values.
              if( kbit(data_notes,11) .and. ut1_tide.eq.0 ) then
                  calc_adj = 0.d0
              else
                  calc_adj = dut1*15.d3
              end if 

              if( abs(old_ut1(1)*15.d3-(pmu_calc(3)-calc_adj))
     .            .gt.0.10 .and. use_pmu_calc ) then
                  write(*,260) epoch, old_ut1(1)*15000.d0, pmu_calc(3)
 260              format(' *** ERROR At epoch ',f12.4,' UT1 ',
     .                   ' Differs from calc: (O,C) ',2f12.3)
                 tol = 1.50*tol
*                If we are observation #1 then kill now before any
*                thing is saved.
                 if( i.eq.1 .and.
     .               abs(old_ut1(1)*15.d3-(pmu_calc(3)-calc_adj))
     .                  .gt.tol*100 ) then
                     stop ' Update_pmu_tab: Killed on first observation'
                 end if

              end if

****	      Save values
              pmu_calc(1) = new_wob(1,1)
              pmu_calc(2) = new_wob(2,1)
*             Add the tidal term to pmu_calc(3) so that we save the total
*             value of UT1 used.
              pmu_calc(3) = (new_ut1(1)+dut1) *15.d3
          end if

*****     Write this observation back to file

          call rw_KalObs_block('W','data',site,error,i)
 
*****     Check for error
          call corrupt_KalObs('IOSTAT',error,i,'writ',1)
 
      end do
 
      end
 
CTITLE UPDATE_MID_PMU
 
      subroutine update_mid_pmu(new_ut1_inf,new_ut1_pts,
     .    new_wob_inf,new_wob_pts, printer, KalObs_file)
 
*     J.L. Davis                   5:25 PM  TUE.,  2  JUNE, 1987
*
*     This routine does the following:
*
*     1. Updates the values of UT1_APR and WOB_APR to reflect
*        the values in the new tables.
*
*     2. Copies the new tables into the location for the old tables
*
*     3. Writes the a priori block to KalObs file
*
*     4. Closes the KalObs file
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       date(5)                 - Date opf first UT1 point
*       dev                     - Current Output device
*       error                   - Error flag
*       i,j,k                   - Loop counters
*       printer                 - Print device
*       trimlen                 - Return length of string
 
      integer*4 date(5), dev, i,j,k, printer, trimlen
 
*       new_ut1_inf(4)          - New UT1 info array
*   ,   new_ut1_pts(4)          - New UT1 array
*   ,   new_wob_inf(3)          - New wobble info array
*   ,   new_wob_pts(2,4)        - New wobble array
*       Sectag                  - Seconds tag
 
      real*8 new_ut1_inf(4), new_ut1_pts(max_pts), new_wob_inf(3),
     .    new_wob_pts(2,max_pts), Sectag
 
 
      character*14 prog
 
*       KalObs_file            - Name of the KalObs file
 
      character*(*) KalObs_file
 
      data prog / 'UPDATE_MID_PMU' /
 
***** Use Lagrangian interpolation and the new tables to get the
*     "a priori" values for PM-UT1 at the epoch PMU_EPOCH
      call Lagrange_intp(new_wob_inf,new_wob_pts,pmu_epoch,wob_apr,2)
      call Lagrange_intp(new_ut1_inf,new_ut1_pts,pmu_epoch,ut1_apr,1)
 
***** Convert the units of the apriori UT1 value from seconds to mas
      ut1_apr(1) = 15.0D+03 * ut1_apr(1)
      ut1_apr(2) = 15.0D+03 * ut1_apr(2)
 
***** Tell the user the new values being used
      call jd_to_ymdhms( fut1_inf(1), date, sectag )
 
***** Output the correction,  Do LU1 and the printer (if different)
      do k = 1,2
          dev = 0
          if( k.eq.1 ) dev = 6
          if( k.eq.2 .and. printer.ne.6 ) dev = printer
*                                ! Output the values
          if( dev.ne.0 ) then
              write(dev,50) KalObs_file(1:trimlen(KalObs_file))
  50          format(/' For KalObs File ',a)
 
              write(dev,100) (fut1_inf(i),i=1,3), date,
     .                       (new_ut1_inf(i),i=1,3)
  100         format(' Old and New UT1 information array',
     .       /' Start JD ',f12.3,' Step and # ',2(f5.2,1x),' Date ',
     .         I4,'/',i2,'/',i2,1x,i2,':',i2,
     .       /' Start JD ',f12.3,' Step and # ',2(f5.2,1x) )
 
              write(dev,120) (fut1_pts(i), new_ut1_pts(i),
     .             (fut1_pts(i)-new_ut1_pts(i))*15.d3,
     .              i=1,nint(fut1_inf(3)))
  120         format(' UT1 points [old,new (tsec), delta (mas)]',/,
     .                10(2f12.7,2x,f12.3,:/) )
 
 
              write(dev,160) (fwob_inf(i),i=1,3),
     .                       (new_wob_inf(i),i=1,3)
  160         format( ' Old and New wobble information array',
     .       /' Start JD ',f12.3,' Step and # ',2(f5.2,1x),
     .       /' Start JD ',f12.3,' Step and # ',2(f5.2,1x) )
 
              write(dev,180) ((fwob_pts(j,i), new_wob_pts(j,i),
     .              (fwob_pts(j,i)-new_wob_pts(j,i)),j=1,2),
     .               i=1,nint(fwob_inf(3)))
  180         format(' Wobble points (old,new,del) x and y (mas)',/
     .               10( 3f11.5,4x,3f11.5,/) )
*                   ! Need to output
          end if
*                   ! Looping over two devices
      end do

* MOD TAH Modded lengths of moves to match integer*4 variables.
 
***** Zero out the old PM-UT1 tables
c     call ClearBuffer(fwob_pts,20)
c     call ClearBuffer(fut1_pts,10)
 
***** Transfer the info arrays
c     call MoveWords(new_wob_inf,fwob_inf, 6)
c     call MoveWords(new_ut1_inf,fut1_inf, 8)
 
***** Transfer table values
c     call MoveWords(new_wob_pts,fwob_pts,16)
c     call MoveWords(new_ut1_pts,fut1_pts, 8)

* MOD TAH: Don't be fancy just copy arrays
      do i = 1, 3
         fwob_inf(i) = new_wob_inf(i)
      end do
      do i = 1,4
         fut1_inf(i) = new_ut1_inf(i)
      end do
      do i = 1,nint(fut1_inf(3))
         fut1_pts(i) = new_ut1_pts(i)
         do j = 1,2
            fwob_pts(j,i) = new_wob_pts(j,i)
         end do
      end do

****  Thats all
      return
      end
 
CTITLE UPDATE_APRIORI_INFO
 
      subroutine update_apriori_info(pmu_file, ut1_def)
 
*     J.L. Davis                   1:47 PM  SAT.,  6  JUNE, 1987
*
*     Routine to make the final changes to the apriori block, and
*     to write the block to the KalObs file.
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       error               - Error flag
*   ,   i                   - Loop counter
 
      integer*4 error, date(5), i, len_file, TrimLen, ut1_def,
     .          len_des, start_file, jerr

      real*8 sectag
 
 
      character*20 prog
 
*       pmu_file            - PMU file name
 
      character*(*) pmu_file
 
      data prog / 'UPDATE_APRIORI_INFO' /
 
***** Get the current date
      call systime(date,sectag)
 
***** How many characters does this file name have?
      len_file = TrimLen(pmu_file)
 
***** Create the polar motion descriptor
      write(user_pmu_dsc,100,iostat=jerr) (date(i),i=1,5), ut1_def
  100 format('PMU on ',I4,'/',I2.2,'/',I2.2,1X,I2.2,':',
     .    I2.2,' UT1_def ',i2)

      len_des = trimlen(user_pmu_dsc)
      start_file = len_des + len_file - 79
      if( start_file.le.0 ) start_file = 1
      user_pmu_dsc(len_des+2:) = pmu_file(start_file:len_file)
 
***** Write this out
      write(*,150) user_pmu_dsc(1:trimlen(user_pmu_dsc))
  150 format(A)

***** Now set the tidal terms in UT1 flag.  In update_pmu we use
*     regularized UT1 in tables, therefore set the flag to indicate
*     this
      ut1_tide = -1
*     Update data_notes to say that we have set this value
      call sbit(data_notes,11,1)
 
***** Write A priori block to KalObs
      call increment_KalVer
      call rw_kalobs_header('W', error )
c     call rw_KalObs_block('W','apr',aprioris,error,0)
 
***** Check for error
      call report_error('FMP',error,'writ','a prioris',0,prog)
 
***** Close KalObs file
      call close_KalObs(error)
 
***** Check for error
      call report_error('FMP',error,'clos','KalObs',0,prog)
 
      end
 
CTITLE REMOVE_TREND
 
      subroutine remove_trend( values, dates, sigmas, num, start, mid,
     .                         trend )
 
*     Routine to remove a linear trend from values.  The slope estimate
*     is weighted by the sigmas and by the time from mid(_epoch).
 
*   i,j         - Loop counters
*   num         - Number of values in values
 
      integer*4 i, num
 
*   a(3), b(2)  - Nomral equations and solution vector
*   dates(1)    - Dates of the values
*   det         - Determinate of normal equations
*   dt          - Time difference between start and date (days)
*   mid         - Mid epoch of VLBI data
*   sigmas(1)   - Sigmas of the values
*   start       - Epoch to which trend will be referred
*   trend(2)    - Offset and rate
*   values(1)   - Values to be fitted
*   wgh         - Weight given to each point
 
      real*8 a(3), b(2), dates(1), det, dt, mid, sigmas(1), start,
     .    trend(2), values(1), wgh
 
***** First check to see if we have enough data
 
*                             ! NOT ENOUGH
      if( num.lt.2 ) then
          stop ' UPDATE_PMU Abort: Not enough data in tables'
      end if
 
*     Clear normal equations
      do i = 1,3
          a(i) = 0.d0
      end do
      do i = 1,2
          b(i) = 0.d0
      end do
 
****  Start incrementing normal equations
 
      do i = 1, num
          dt = dates(i) - start
          wgh = 1.d0/ (sigmas(i)**2 * (abs(dates(i)-mid)+.1d0)**2 )
 
          a(1) = a(1) + wgh
          a(2) = a(2) + wgh*dt
          a(3) = a(3) + wgh*dt*dt
 
          b(1) = b(1) + wgh*values(i)
          b(2) = b(2) + wgh*dt*values(i)
      end do
 
*     Invert and get the trend
      det = a(1)*a(3) - a(2)*a(2)
 
*     Chech the det
*                                     ! Not a good solution
      if( det.lt.1.d-7*a(1) ) then
          stop ' UPDATE_PMU aborted:  Trend removal singular'
      end if
 
      trend(1) = (a(3)*b(1) - a(2)*b(2) )/det
      trend(2) = (a(1)*b(2) - a(2)*b(1) )/det
 
****  Now remove trend from data
      do i = 1,num
          values(i) = values(i) - trend(1) - trend(2)*(dates(i)-start)
      end do
 
****  Thats all
      return
      end

CTITLE get_last_pmu

      subroutine get_last_pmu( user_pmu_dsc, jd_last )

*     This routine will read the user_pmu_dsc and get the last time
*     that update pmu was run.  -1 is returned if no last run

* jd_last   - Julian date of last run

      real*8 jd_last

* user_pmu_dsc - Update pmu description line

      character*(*) user_pmu_dsc

* LOCAL

* trimlen  - Length of string
* indx     - Position in string
* date(5)  - Date run from line

      integer*4 trimlen, indx, date(5), i, ierr

* sectag   - Seconds tag
      real*8 sectag

***** First see if we have anything in string

      if( trimlen(user_pmu_dsc).eq.0 .or.
     .    ichar(user_pmu_dsc(1:1)).eq.0 ) then
          jd_last = -1.d0
          return
      end if

***** Get position of time
      indx = index( user_pmu_dsc,' on') + 6

*     Kill off the / characters
      call sub_char( user_pmu_dsc, '/',' ')

      read(user_pmu_dsc(indx:),*,iostat=ierr) (date(i),i=1,3)
      call report_error('IOSTAT',ierr,'read',user_pmu_dsc,1,
     .                  'UPDATE_PMU')

      date(4) = 0
      date(5) = 0
      sectag  = 0
      call ymdhms_to_jd( date, sectag, jd_last )

****  Thats all
      return
      end

