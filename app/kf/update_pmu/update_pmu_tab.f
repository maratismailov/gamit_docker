
CTITLE UPDATE_PMU_TAB
 
      program update_pmu_tab
 
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
*     UT1 which has a high slope. Changes made in GET_NEW_TABS.
 
* MOD TAH 880317 Introduced new UT1 type.  Bit 2 of type is set if
*     series is not regularized. ie. if short period UT1 terms have
*     not been removed. (GSFC and NGS are usually like this)
*
* MOD TAH 881228 Added output of corrections to a print device so that
*     we can keep a record of the changes.  Also fixed leap second problem
*     when we are crossing a boundary and using UT1-UTC tables.

* MOD TAH 910115: Added check on date at which update_pmu_tab was run so that
*     runs which did not save the pmu_calc values could be detected.  These
*     runs all occurred before 90/11/25

* MOD TAH 910418: This is a version of update_pmu_tab which expects tabular
*     points for its input tables.  It reads the data and then selects
*     the closest group of 5 tabular points.

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
*     If it is greater than max_pts then reset the value.
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
      call get_new_tabs(pmu_file,num_header,new_ut1_inf,new_ut1_pts,
     .    new_wob_inf,new_wob_pts,ut1_def, printer)
 
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
 
      data help_file / 'update_pmu_tab.hlp' /
 
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
     .                                 'UPDATE_PMU_tab',1)
      if (error2 .ne. 0)          call proper_runstring(help_file,
     .                                 'UPDATE_PMU_tab',1)
      if (error3 .ne. 0)          call proper_runstring(help_file,
     .                                 'UPDATE_PMU_tab',1)
      if (len_pmu_file .eq. 0)    call proper_runstring(help_file,
     .                                 'UPDATE_PMU_tab',1)
      if (len_KalObs_file .eq. 0) call proper_runstring(help_file,
     .                                 'UPDATE_PMU_tab',1)
 
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
 
CTITLE GET_NEW_TABS
 
      subroutine get_new_tabs(pmu_file,num_header,new_ut1_inf,
     .    new_ut1_pts,new_wob_inf,new_wob_pts,ut1_def, printer)
 
*	This routine will scan the polar motion/UT1 table and return the 
*	information arrays for the mid_epoch of the experiment 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
*       max_return_val              - See above
 
      integer*4 max_return_val 

      parameter (max_return_val = 55)
 
*       error                       - Error flag
*   ,   i, j                        - Loop counter
*   ,   num_header                  - Number of header records
*                                   -   at top of PM-UT1 file
*   ,   num_return                  - Number of values returned
*                                   -   by file scan
*   ,   ut1_def                     - UT1 definition;  see help file
*       printer                     - Print unit number
*       line_err                    - IOSTAT error decoding line
*       im                          - Point number a center of data
*                                     span
*       num_grid                    - Number of values to be saved
*                                     in the grid.
 
      integer*4 error, i, j, num_header, num_return, ut1_def, 
     .          printer, line_err, im, indx, date(5), num_grid

      logical OK
 
*       day_limit                   - Range (days) for table return
*   ,   new_ut1_inf(4)              - New info array for UT1
*   ,   new_ut1_pts(max_pts)              - New array of points for UT1
*   ,   new_wob_inf(3)              - New info array for X-Y PM
*   ,   new_wob_pts(2,max_pts)            - New array of points for X-Y PM
*   ,   return_dates(max_return_val)    - Array of Julian dates
*   ,   xp(max_return_val), yp(max_return_val), ut(max_return_val)   
*                                    - Table of X, Y, or UT1 values
*       jd  - Generic JD
*       sectag - Seconds tag for conversion of dates
*       spacing - Spacing between data.
 
      real*8 day_limit, new_ut1_inf(4), new_ut1_pts(max_pts), 
     .    new_wob_inf(3), new_wob_pts(2,max_pts), 
     .    return_dates(max_return_val),
     .    xp(max_return_val), yp(max_return_val), ut(max_return_val),
     .    jd, spacing, sectag, vals(10)
 
*       pmu_file                    - File with PM-UT1 values
 
      character*(*) pmu_file

*       line    - line read from input

      character*256 line
 
*       prog*14                     - Routine name
 
      character prog*14
 
      data prog / 'GET_NEW_TABS' /
 
***** Open PM-UT1 file
      open (105, file = pmu_file, iostat = error, status = 'old')
 
***** Error?
      call report_error('IOSTAT',error,'open',pmu_file,1,prog)

 
***** How many days should be included for X-Y PM/UT1:  Set to 25 days
*     before and after mid-epoch
      day_limit = 25

***** Skip over headers
      Do i = 1, num_header
          read (105,'(a)') line
      end do
 
*     Counter for number of values
      j = 0
      OK = .true.

*     Now loop finding the data we need
      do while ( error.eq.0 )
           read(105,'(a)',iostat=error) line

           if( error.eq.0 .and. line(1:1).eq.' ' ) then
               indx = 1
               call multiread(line,indx,'I4', line_err,date,prog,5)
               sectag = 0.d0
               call ymdhms_to_jd( date, sectag, jd )

*              See if within limits
               if( abs(mid_epoch-jd).lt.day_limit ) then
                   j = j + 1
                   num_return = j
                   return_dates(j) = jd

*                  Now decode the rest of the line
                   call multiread(line,indx,'R8', line_err, vals,
     .                             prog, 5)
                   xp(j) = vals(1)
                   yp(j) = vals(3)
                   ut(j) = vals(5)
                   if( line_err.ne.0 ) OK = .false.
               end if

*              See if we are past the useful point
               if( jd.gt. mid_epoch+day_limit ) then
                   error = -1
               end if
          end if
      end do

***** See if OK to continue
      if( .not.OK ) then
          write(printer,100)
 100      format(' *** ERROR *** decoding line from PMU file')
          stop ' UPDATE_PMU_TAB: Aborted due to error in PMU file'
      end if
      if( num_return.le.4 ) then
          write(printer,110)
 110      format(' *** ERROR *** Not enough data read from PMU file')
          stop ' UPDATE_PMU_TAB: Aborted, not enough data in PMU file'
      end if

***** Change UT1 values from whatever they were to UT1-TAI
      call change_ut1(ut,return_dates, num_return,ut1_def)

***** Now check out the spacing of the data
      spacing = return_dates(2) - return_dates(1)

***** See if spacing is constant
      do i = 2, num_return-1
          if( abs(spacing - (return_dates(i+1) - return_dates(i))).gt.
     .        0.0001d0 ) then
              write(printer,120) spacing, i,
     .                          (return_dates(j),j=1,num_return)
 120          format(' *** ERROR *** Nominal spacing is ',f3.1,' days.',
     .               ' For point ',i3,' spacing different',/,
     .               ' Dump of julian dates of data is',/,
     .               10(5F12.3,/) )
              stop ' UPDATE_PMU_TAB: Aborted, table not uniform'
           end if
      end do

***** Now find out which entries we actually want.  First find closest
*     data point
*     Set the number of points to use
      num_grid = 5 + 1.d0/spacing
      if( num_grid.gt.max_pts ) num_grid = max_pts
      im = int((start_epoch-return_dates(1))/spacing)

*     Check we have enough enough data
      if( return_dates(num_return).le. 
     .    return_dates(im+num_grid-1)     ) then
          write(printer,140) num_grid, spacing, im, num_return,
     .            return_dates(im)
 140      format(' *** ERROR *** Not enough returned values',/,
     .           I4,' points wanted with spacing of ',f5.2,
     .              ' days.  First point is index ',i4,' and ',
     .           i4,' points in returned array (JD ',f12.3,')')
          stop ' UPDATE_PMU_TAB: Aborted; Not enough return data'
      end if

*     Check to see that the number of points will span the data
*     range in the experiment
      if( end_epoch.gt.return_dates(im+num_grid)+spacing) then
          write(printer,160) num_grid, spacing, im, num_return,
     .            return_dates(im)
 160      format(' *** ERROR *** Not enough points in INF arrays',/,
     .           I4,' points wanted with spacing of ',f5.2,
     .              ' days.  First point is index ',i4,' and ',
     .           i4,' points in returned array (JD ',f12.3,')')
          stop ' UPDATE_PMU_TAB: Aborted; Not points in INF array'
      end if


*     Now use the two points priori and imediately after this point
      j = 0
      do i = im, im+num_grid-1
         j = j + 1

*        Convert wobble to mas.  The UT1 stays in seconds
         new_wob_pts(1,j) = xp(i)*1000.d0
         new_wob_pts(2,j) = yp(i)*1000.d0
         new_ut1_pts(j) = ut(i)
      end do

      new_wob_inf(1) = return_dates(im)
      new_wob_inf(2) = spacing
      new_wob_inf(3) = num_grid

*
      new_ut1_inf(1) = return_dates(im)
      new_ut1_inf(2) = spacing
      new_ut1_inf(3) = num_grid
      new_ut1_inf(4) = 1.d0
 
***** Close PM-UT1 file
      close (105)
 
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
*       ut1_values(num_values)           - Values of UT1
*       return_dates(num_values)         - JD's of values
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
      else
          call sbit( ut1_def, 1, 1 )
      end if
      if( .not.kbit( ut1_def,1) ) then
          write(*,100) 
 100      format(' UT1 definition set to UT1-UTC')
      else
          write(*,110)
 110      format(' UT1 definition set to UT1-TAI')
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
*       calc_adj                - Tidal contribution to UT1 needed when 
*                                 comparing the full UT1 in PMU_CALC with
*                                 the interpolated
*                                 values from the tables (mas) 
*       jd_bug                  - Date at which bug in not saving pmu fixed
*       tol                     - Tol set on out put of errors
 
      real*8 delta_pmx, delta_pmy, delta_ut1, delta_tau, delta_tau_dot,
     .    new_ut1_inf(1), new_ut1_pts(1), new_wob(2,2), new_wob_inf(1),
     .    new_wob_pts(1), new_ut1(2), old_wob(2,2), old_ut1(2), dut1,
     .    calc_adj, jd_bug, tol

*       jd_last                 - Last julian date at which we ran update_pmu_tab

      real*8 jd_last

*       use_pmu_calc            - Indicates that we should use pmu_calc values

      logical use_pmu_calc

***** Get the last time we ran update_pmu_tab.  Jd_last is returned -1 if
*     update_pmu_tab has not been run before.

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

***** Set tolerance on error output
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
                 if( abs(old_wob(j,1)-pmu_calc(j)).gt.tol .and.
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
     .            .gt.tol .and. use_pmu_calc ) then
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

****          Save values
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
*   ,   new_ut1_pts(max_pts)          - New UT1 array
*   ,   new_wob_inf(3)          - New wobble info array
*   ,   new_wob_pts(2,max_pts)        - New wobble array
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
     .               i=1,nint(new_ut1_inf(3)))

  120         format(' UT1 points [old,new (tsec), delta (mas)]',/,
     .        4(2f12.7,2x,f12.3,:/) )
 
 
              write(dev,160) (fwob_inf(i),i=1,3),
     .                       (new_wob_inf(i),i=1,3)
  160         format( ' Old and New wobble information array',
     .       /' Start JD ',f12.3,' Step and # ',2(f5.2,1x),
     .       /' Start JD ',f12.3,' Step and # ',2(f5.2,1x) )
 
              write(dev,180) ((fwob_pts(j,i), new_wob_pts(j,i),
     .              (fwob_pts(j,i)-new_wob_pts(j,i)),j=1,2),
     .               i=1,nint(new_wob_inf(3)))
  180         format(' Wobble points (old,new,del) x and y (mas)',/
     .              10( 3f11.5,4x,3f11.5,/) )
*                   ! Need to output
          end if
*                   ! Looping over two devices
      end do


* MOD TAH: Don't be fancy just copy arrays
      do i = 1, 3
         fwob_inf(i) = new_wob_inf(i)
      end do
      do i = 1,4
         fut1_inf(i) = new_ut1_inf(i)
      end do
      do i = 1,nint(new_ut1_inf(3))
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
  100 format('TAB on ',I4,'/',I2.2,'/',I2.2,1X,I2.2,':',
     .    I2.2,' UT1_def ',i2)

      len_des = trimlen(user_pmu_dsc)
      start_file = len_des + len_file - 79
      if( start_file.le.0 ) start_file = 1
      user_pmu_dsc(len_des+2:) = pmu_file(start_file:len_file)
 
***** Write this out
      write(*,150) user_pmu_dsc(1:trimlen(user_pmu_dsc))
  150 format(A)

***** Now set the tidal terms in UT1 flag.  In update_pmu_tab we use
*     regularized UT1 in tables, therefore set the flag to indicate
*     this
      ut1_tide = -1
*     Update data_notes to say that we have set this value
      call sbit(data_notes,11,1)
 
***** Write A priori block to KalObs
      call increment_KalVer
      call rw_kalobs_header('W', error )
 
***** Check for error
      call report_error('FMP',error,'writ','a prioris',0,prog)
 
***** Close KalObs file
      call close_KalObs(error)
 
***** Check for error
      call report_error('FMP',error,'clos','KalObs',0,prog)
 
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
     .                  'UPDATE_PMU_tab')

      date(4) = 0
      date(5) = 0
      sectag  = 0
      call ymdhms_to_jd( date, sectag, jd_last )

****  Thats all
      return
      end

