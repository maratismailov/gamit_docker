CTITLE
 
      subroutine plts8
c
c     This segment decodes the PLOT runstring information
c     and opens the user's control file if needed.
c
c Include files
c -------------
*                          ! the parameter file
      include 'plot_param.h'
c
*                          ! the common block
      include 'plot_com.h'
c
*   DecimalToInt     - Conversion of string to an integer value
*                    - ie. the string '13' will give integer value 13.
*   ierr             - IOSTAT error
 
      integer*4 DecimalToInt, ierr, i

      integer*4 k, len_run, rcpar 
 
c
c.... Now see how many parameters were passed in runstrings
*                                          ! open the device
      if( runstring(1)(1:1).ne.' ' ) then
c
*                           ! unit number for device
         userlu = 100
         call open_lu(userlu, runstring(1), ierr,'old')
         if( userlu.eq.100 ) then
             terminal = .false.
         else
             terminal = .true.
         end if
c
c....    Check the error on open
*                                              ! error ocurred, see if
         if( ierr.ne.0 ) then
*                                            ! Stop program
            call report_error('IOSTAT',ierr,'open',runstring(1),
     .                      1,'PLOT')
*                  ! open error
         end if
c
*                                            ! nothing in runstring, so
      else
c                                              use loglu function
         userlu = 6
         terminal = .true.
c
      end if
c
c.... Now get the plot lu
*                                           ! get from runstring
      if( runstring(2)(1:1).ne.' '  ) then
         meta_file = runstring(2)
      end if
*                                           ! use user terminal
      plotlu = 6
c
c.... See if file name passed, assign name. The file will be read later
c     when we have more information about which columns of file contain
c     data.
      num_files = 0
      if( runstring(3)(1:1).ne.' ' ) then
c                                             ! save the file name
         input_file = runstring(3)
*        Save this file as first in saved list
         num_files = 1 
         saved_files(1) = input_file
c
*              ! set names to blank
      else
c
         input_file = '  '
c
      end if
c
c.... See if the number of header records in the file was passed
*                                           ! headers passed
      if( runstring(4)(1:1).ne.' ' ) then
          input_num_headers = decimaltoint( runstring(4), ierr )
*                                 ! error decoding
          if( ierr.ne.0 ) then
 
c....         Report and kill progrm
              call report_error('decimalToInt',ierr,'decod',
     .            runstring(4),1,'PLOT')
 
          end if
c
c.....    Compute actual number of headers we can save
          actual_num_headers = min(input_num_headers,max_headers)
c
      end if

      ignore_col1 = .true.
      if( runstring(5)(1:1).ne.' ' ) then
          i = decimaltoint( runstring(5), ierr )
          if( i.gt.0 .and. ierr.eq.0 ) then
              ignore_col1 = .false.
          end if
      end if

*     Now see if more file names have been passed
      len_run = 1
      k = 5
      do while ( len_run.gt.0 .and. num_files.lt.max_files )
           k = k + 1
           len_run = rcpar(k, runstring(5) )
           if( len_run.gt.0 ) then
               num_files = num_files + 1
               saved_files(num_files) = runstring(5)
           endif
      end do

c
c.... Thats all
*                   ! return to main segment
      return
c
      end
 
 
