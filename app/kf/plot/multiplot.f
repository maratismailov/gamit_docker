      program multiplot
 
 
*     Program to read a values file, and to output multiple plots
*     num on a page a time.  The root passed in the runstring gives
*     the root part of the name of the baseline files to be created
*     and the name of the meta file to be saved by PLOT.

* MOD TAH 960812: Changed to run cplotx

*     Runstring:
*     Multiplot,<values>,<num>, <start y,m,d>, <end y,m,d>, y_grad,
*                       edit, <sigma_limit>, <root> <start>, <skip>
*
*     <root>        is the root for the baseline files and the meta
*                   file.
*     <values>      is name of values file
*     <start y,m,d> is start yr,mon,day. If not given plot sets values.
*     <end y,m,d>   is end yr,mon,day
*     <y_grad>      is graduated steps for y_axis (m)
*     <edit>        is name of file for editted data.  These are greater
*                   than 3 scaled sigmas or sigma > 1 meter
*     <sigma_limit> number of sigmas to throw out. Default 3.
*     <start>       is number of values in to start with. Default 1.
*     <skip>        is number of values to skip. Default 0
*
 
      integer*4 max_data
 
      parameter ( max_data = 5000 )
 
*   num_data        - Number of data read
*   used_data       - Actual number of data left after editing
*   i,j,k           - Loop control
*   ierr            - IOSTAT error
*   outlu           - Output lu.
*   irow, icol      - Current row and column of plot
*   start           - start number
*   skip            - number to skip each time
*   date(5)         - ymdhm read from file or elsewhere
*   trimlen         - Returns length of string
*   edit_flag(max_data)     - Edit flags for data
*   len_buffer       - Length of buffer read from file
*   indx            - Index for readline
*   ir              - Used for setting range
*   ifive(5)        - Five parameters returned from plot (not
*                   - used)
      integer*4 num_data, used_data, i,j,ierr, outlu, irow, icol,
     .    start, skip, date(5), trimlen, edit_flag(max_data),
     .    len_buffer, indx, ir, ifive(5)
 
*   offset_unit     - Offset for unit pointer (not really used)
*   num_per_page        - Number of plots per page
*   mxr, mxc            - Maximum rows and columns per page
*   page            - Current page number
*   edit_unit       - Unit number for edited data file.  (maybe
*                     6 if no file name given)
 
      integer*4 offset_unit, num_per_page, mxr, mxc, page,
     .          edit_unit
 
*   sectag          - Seconds tag
*   jd              - jd of data
*   epoch(max_data) - Epochs of data points
*   length(max_data)- Length read from values file
*   sigma(max_data) - Sigmas of data read from file
*   sigma_limit     - N sigma to discard
*   start_jd        - start jd if given in runstring
*   end_jd          - end jd if given in runstring
*   y_grad          - step to be used on plot (if given in
*                   - runstring)
*   nom_length      - Nominal length to remove from values
*   ref_length      - Reference length removed from all values
*   ref_epoch       - Reference epoch removed from all values
*   value           - value read from file
*   sig             - Sigma read from file
*   mean            - Mean length from file
*   nrms            - NRMS scatter
*   mean_epoch      - Mean epoch from file
*   slope           - Slope in mm/yr read-- converted to m/day
*   min_length      - minimum length
*   max_length      - maximim length
*   min_jd, max_jd  - MIn and max JD of values
*   range           - range needed for min and max length
*   avlength        - Median length for min and max length
*   ymin, ymax      - ymin and ymax needed for plot
 
      real*8 sectag, jd, epoch(max_data), length(max_data),
     .    sigma(max_data), sigma_limit, start_jd, end_jd, y_grad,
     .    nom_length, ref_epoch, value, sig, mean, nrms,
     .    mean_epoch, slope, min_length, max_length, min_jd, max_jd,
     .    range, avlength, ymin, ymax
 
*       eof         - Indicates end of file
      logical eof
 
*   cdum            - Dummy character
 
      character*4 cdum
 
*   cdate           ! Current time and data as character string
      character*32  cdate
 
*   values_file     - Name of values file
*   edit_file       - Name of edit file
*   plot_ctrl       - Name of plot control file
*   plot_data       - Name of plot data file
*   runname         - Name of plot when it runs
*   root            - Root of the files to be generated for the
*                   - plot data files, and the meta file.
*   meta_file       - Name of the metafile (if root is passed)
 
      character*128 values_file, edit_file, plot_ctrl, plot_data,
     .    runname, root, meta_file
 
*             title - Title of plots
*             runline   - Line for FmpRunProgram
 
      character*132 title, runline
 
*             baseline  - Name of baseline
*             buffer    - Buffer for reading values into
 
      character*80 baseline, buffer

*     plot_residual - Indicates that we should plot residuals
*     and not the total value

      logical plot_residual
 
      common / data_block / epoch, length, sigma, edit_flag
 
***** Set names
      plot_ctrl = 'multi_con.plt'
 
***** Start, decode runstring
      call get_runstring(outlu, values_file, num_per_page, start_jd,
     .                end_jd, y_grad, edit_file, sigma_limit,
     .                    root, start, skip)

*     See if we ploting residuals
      if( num_per_page.lt.0 ) then
          num_per_page = abs(num_per_page)
          plot_residual = .true.
      else
          plot_residual = .false.
      end if
 
*     Finish up calculation of rows and columns
      mxc = 2
      if( num_per_page.gt.9 ) then
          mxc = 3
      end if
      if( num_per_page.eq.1 ) then
          mxc = 1
      end if
      mxr = num_per_page/mxc
 
      irow = mxr
      icol = mxc
      page = 0
 
*     Generate current time
      call systime( date, sectag)
      write(cdate,100) date
 100  format('Run: ',i4,'/',i2,'/',i2,1x,i2,':',i2)
 
***** Open up the values file, the edit file, amd the plot files
 
      open(100, file=values_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',values_file,1,'multiplot')
 
*     open(200, file=edit_file, iostat=ierr, status='unknown')
      call open_lu(edit_unit, edit_file, ierr, 'append')
      call report_error('IOSTAT',ierr,'open',edit_file,1,'multiplot')
 
      open(201, file=plot_ctrl, iostat=ierr, status='unknown')
      call report_error('IOSTAT',ierr,'open',plot_ctrl,1,'multiplot')
 
*     Write out header to file
      write(201,200) cdate(1:trimlen(cdate)),
     .    values_file(1:trimlen(values_file))
 200  format('*',/,
     .       '* PLOT control file generated by MULTIPLOT on ',a,/,
     .       '* Values file : ',a)
 
****  Now start looping through the values file
      sectag = 0
      eof = .false.
      do while ( .not.eof )
 
*****     Start getting the data
          read(100,'(a)',iostat=ierr) title
          if( title(1:11).eq.'VALUES_FILE') then
              read(100,'(a)',iostat=ierr) title
          end if

          if( trimlen(title).eq.0 ) then
              title = '* SOLUTION OF UNKNOWN TYPE'
          end if
*                                         ! No more data
          if( ierr.ne.0 ) then
              eof = .true.
*                                         ! Still data in file
          else
              read(100,'(a)', iostat=ierr) baseline
*                                         ! Skip a line
              read(100,'(a)') buffer
 
*             open the output data file
              call open_plot_data( 202, baseline, root, plot_data )
 
*****         Now start reading values
              do j = 1, start-1
*                                         ! Skip to start place
                  read(100,'(a)') buffer
              end do
 
*****         Now start reading
              len_buffer = 1
*                                         ! Counter for num_data
              i = 0
              do while ( len_buffer.gt.0 )
                  read(100,'(a)',iostat=ierr) buffer
                  len_buffer = trimlen(buffer)
*                                                 ! Values here
                  if( len_buffer.gt.0 .and.
     .               buffer(1:1).eq.' '   ) then
                      i = i + 1
                      indx = 1
                      call multiread(buffer,indx,'I4', ierr, date,
     .                                cdum,5)
                      call ymdhms_to_mjd(date, sectag,jd)
                      call read_line(buffer, indx, 'R8', ierr,
     .                               value, cdum)
                      call read_line(buffer, indx, 'R8', ierr,
     .                               sig, cdum)

*                     See if we are ploting residuals. If we
*                     are get the next two values.
                      if( plot_residual ) then
                         call read_line(buffer, indx, 'R8', ierr,
     .                                  value, cdum)
                         call read_line(buffer, indx, 'R8', ierr,
     .                                  sig, cdum)
                      end if
 
*                     Save reference values if first point: DONT NEED
*                     TO DO YET
                      epoch(i) = jd
                      length(i) = value
                      sigma(i) = sig
 
*                     Delete any Weird dates
                      if( date(1).lt.70 ) then
                          i = i - 1
                      end if
                  end if
 
*****             Now skip line (if need)
                  j = skip
                  do while ( j.gt.0 .and. len_buffer.gt.0 )
                      j = j - 1
                      read(100,'(a)') buffer
                      len_buffer = trimlen(buffer)
                  end do
*                                 ! Reading the values from the file
              end do
 
*****         Save number of data
              num_data = i
 
*****         Now get the mean, and the NRMS
              read(100,'(a)') buffer
              indx = 6
              call read_line(buffer,indx,'R8',ierr, mean, cdum)
              indx = 68
              call read_line(buffer, indx,'R8', ierr,nrms, cdum)
 
*****         Get slope and mean epoch
              read(100,'(a)') buffer
              indx = 6
              call read_line(buffer,indx,'R8', ierr, slope, cdum)
 
*             Check size of slope (aviods problems with small amount
*             of data)
              if( abs(slope).gt. 200 ) then
                  slope = 0
              end if
 
*                                              ! Convert mm/yr to m/day
              slope = slope/365.d0/1000.d0
              indx = 71
              call read_line(buffer,indx,'R8', ierr, mean_epoch,
     .                       cdum)

*             If we plotinmg resiuduals, set mean and slope to zero
              if( plot_residual ) then
                  mean = 0.d0
                  slope = 0.0d0
              end if
 
****          Convert mean_epoch to refernce jd
              date(1) = mean_epoch
              date(2) = 1
              date(3) = (mean_epoch-date(1))*365.d0
              date(4) = 0
              date(5) = 0
              call ymdhms_to_mjd( date, sectag, ref_epoch)
 
*****         Skip a line
              read(100,'(a)', iostat=ierr) buffer
              if( ierr.ne.0 ) eof = .true.
 
*****         Now start setting up the plot and data file
              call edit_data(epoch,length,sigma, edit_flag, baseline,
     .             num_data, sigma_limit, nrms, mean, ref_epoch,
     .             slope, min_length, max_length,used_data,
     .             min_jd, max_jd, edit_unit)
 
*****         Get the nominal length
              nom_length = aint(mean/10)*10.d0
              avlength = (max_length+min_length)/2
              if( plot_residual ) then
                  nom_length = 0.d0
                  avlength   = 0.d0
              end if

              len_buffer = trimlen(baseline)
              write(baseline(len_buffer+1:), 300) avlength
  300         format(' +',f14.3,' m')
 
*****         If y_grad was given set scales
              if( y_grad.ne.0 ) then
                  range = (max_length - min_length)*1.2
                  ir = range/y_grad + 1
                  avlength = (max_length + min_length)/2 - nom_length
                  ymin = avlength - ir*y_grad/2
                  ymax = avlength + ir*y_grad/2
              else
                  range = max_length - min_length
                  ymax = max_length - nom_length + 0.1*range
                  ymin = min_length - nom_length - 0.1*range
              end if
 
*****         Now output the data file
              if( used_data.gt.1 ) then
                  call write_data(epoch, length, sigma, edit_flag,
     .                 num_data, nom_length, title, baseline)
                 close(202)
 
*****             Write out plot contol file
                  if( start_jd.ne.0 ) then
                      min_jd = start_jd
                      max_jd = end_jd
                  end if
 
*                Advance to next plot position
                 call advance_rc( 201, irow, icol, page, values_file,
     .                        cdate, mxr, mxc )
 
                  call write_plot( ymin, ymax, y_grad, start_jd, end_jd,
     .                  used_data, min_jd, max_jd, plot_data)
                 write(*,220) baseline(1:trimlen(baseline)),
     .                    irow, icol, page
 220              format(' Processed ',a,' Row ',i3,' Col ',i3,
     .                  ' Page ', i3)
 
 
*                     ! Any data to plot
              end if
*                     ! No eof yet
          end if
 
*****     Check to see if break issued
C         if( ifbrk() ) then
C             Eof = .true.
C         end if
*                     ! while eof not found
      end do
 
***** Thats all
      close(100)
      close(200)
      close(201)
 
***** Generate the meta file name and run plot
 
      if( root(1:3).ne.'mp_' ) then
          meta_file = root(1:trimlen(root)) // 'meta'
      else
          meta_file = ' '
      end if
 
      runline = 'cplotx ' // plot_ctrl(1:trimlen(plot_ctrl)) 
*      call fmprunprogram(runline, ifive, runname, 1,6,offset_unit)
      call execute(runline, ifive, 1,6,offset_unit)
 
***** That all
      end
 
 
CTITLE GET_RUNSTRING
 
      subroutine get_runstring(outlu, values_file, npp, start_jd,
     .           end_jd, y_grad, edit_file, sigma_limit, root,
     .           start, skip)
 
*     Routine to decode the runstring
 
*   i,j,k           - Loop control
*   ierr            - IOSTAT error
*   outlu           - Output lu.
*   start           - start number
*   skip            - number to skip each time
*   date(5)         - ymdhm read from file or elsewhere
*   rcpar           - Read runstring
*   len_run         - Length of runstring
*   npp             - Number of plots per page
 
      integer*4 i, ierr, outlu, start, skip, date(5), rcpar,
     .    len_run, npp, decimaltoint
 
*   sectag          - Seconds tag
*   jd              - jd of data
*   epoch(1)        - Epochs of data points
*   length(1)       - Length read from values file
*   sigma(1)        - Sigmas of data read from file
*   sigma_limit     - N sigma to discard
*   start_jd        - start jd if given in runstring
*   end_jd          - end jd if given in runstring
*   y_grad          - step to be used on plot (if given in
*                   - runstring)
*   nom_length      - Nominal length to remove from values
 
      real*8 sectag, sigma_limit, start_jd, end_jd, y_grad
 
*   values_file     - Name of values file
*   edit_file       - Name of edit file
*   root            - Root of the file names
 
      character*(*) values_file, edit_file, root
 
*   runstring       - Runstring
 
      character*128 runstring
 
***** Start decoding
      outlu = 6
 
      len_run = rcpar(1,values_file)
      if( len_run.le.0 ) then
          call proper_runstring('multiplot.hlp','multiplot',0)
          Stop ' MULTI_PLOT stop: Incorrect runstring'
      end if
 
****  Get the number of plots per page
      len_run = rcpar(2,runstring)
      if( len_run.gt.0 ) then
          npp = decimaltoint( runstring, ierr)
          if( ierr.ne.0 .or. npp.eq.0 ) len_run = 0
      end if
*                             ! If any problem use default
      if( len_run.le.0 ) then
          npp = 12
      end if
 
 
***** Get start date if it is threre
      start_jd = 0
      do i = 1, 3
          len_run = rcpar(2+i,runstring)
          if( len_run.gt.0 ) then
              read(runstring,*,iostat=ierr) date(i)
              call report_error('IOSTAT',ierr,'decod',runstring,1,
     .                     'Get_runstring')
          else
              date(i) = 0
          end if
      end do
      date(4) = 0
      date(5) = 0
      sectag = 0
      if( date(1).gt.0 ) then
          call ymdhms_to_mjd(date, sectag, start_jd)
      end if
 
***** Get end date if it is threre
      end_jd = 0
      do i = 1, 3
          len_run = rcpar(5+i,runstring)
          if( len_run.gt.0 ) then
              read(runstring,*,iostat=ierr) date(i)
              call report_error('IOSTAT',ierr,'decod',runstring,1,
     .                     'Get_runstring')
          else
              date(i) = 0
          end if
      end do
      date(4) = 0
      date(5) = 0
      sectag = 0
      if( date(1).gt.0 ) then
          call ymdhms_to_mjd(date, sectag, end_jd)
      end if
 
***** Get y_grad
      y_grad = 0
      len_run = rcpar(9,runstring)
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr) y_grad
          call report_error('IOSTAT',ierr,'decod',runstring,1,
     .                'Get_runstring')
      end if
 
***** Get edit file name
      len_run = rcpar(10,edit_file)
      if( len_run.le.0 ) edit_file = '6'
 
***** Get sigma_limit
      sigma_limit = 10
      len_run = rcpar(11,runstring)
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr) sigma_limit
          call report_error('IOSTAT',ierr,'decod',runstring,1,
     .                'Get_runstring')
      end if
 
***** Get the root to the file names
      len_run = rcpar(12, root)
      if( len_run.le.0 ) then
          root = 'mp_'
      end if
 
***** Get strat
      start = 0
      len_run = rcpar(13,runstring)
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr) start
          call report_error('IOSTAT',ierr,'decod',runstring,1,
     .                'Get_runstring')
      end if
 
***** Get skip
      skip = 0
      len_run = rcpar(14,runstring)
      if( len_run.gt.0 ) then
          read(runstring,*,iostat=ierr) skip
          call report_error('IOSTAT',ierr,'decod',runstring,1,
     .                'Get_runstring')
      end if
 
***** Thats all
      return
      end
 
CTITLE EDIT_DATA
 
      subroutine edit_data(epoch,length,sigma, edit_flag, baseline,
     .           num_data, sigma_limit, nrms, mean, ref_epoch, slope,
     .           min_length, max_length, used_data, min_jd, max_jd,
     .           edit_unit)
 
 
*     Routine to scan data and edit "bad points" using sigma_limit
*     and size of error bar
*
*   num_data        - Number of data read
*   used_data       - Actual number of data used
*   i,j,k           - Loop control
*   ierr            - IOSTAT error
*   edit_flag(1)    - Edit flags for data
*   date(5)         - Date needed for output
*   trimlen         - Legnth of string
*   edit_unit       - Unit number for edited data file.
 
      integer*4 num_data, used_data, i, edit_flag(1),
     .    date(5), trimlen, edit_unit
 
*   epoch(1)        - Epochs of data points
*   length(1)       - Length read from values file
*   sigma(1)        - Sigmas of data read from file
*   sigma_limit     - N sigma to discard
*   ref_epoch       - Reference epoch removed from all values
*   mean            - Mean length from file
*   nrms            - NRMS scatter
*   slope           - Slope in mm/yr read-- converted to m/day
*   est_length      - Estimated length
*   residual        - residual to length
*   min_length      - Minum length
*   max_length      - Max length
*   jd              - Julian date
*   sectag          - Seconds tag
*   min_jd, max_jd      - Min and max times
*   range           - Range of time for increasing sizes a little.
 
      real*8 epoch(1), length(1), sigma(1), sigma_limit, ref_epoch,
     .    mean, nrms, slope, est_length, residual, min_length,
     .    max_length, jd, sectag, min_jd, max_jd, range
 
*             baseline - Name of baseline
 
      character*(*) baseline
 
*       outline     - Indicates name has been output
 
 
      logical outline
 
***** Start looping
      outline = .false.
      min_length = 1.d20
      max_length = -1.d20
      used_data = 0
 
      max_jd = 0
      min_jd = 1.d20
 
      do i = 1, num_data
 
*         Get estimated length
          est_length = mean + slope*(epoch(i)-ref_epoch)
          residual = length(i) - est_length
 
          if( sigma(i).gt. 0.3 ) then
              edit_flag(i) = 1
          else
              if( abs(residual/sigma(i)).gt.sigma_limit ) then
                  edit_flag(i) = 2
              else
                  edit_flag(i) = 0
                  min_length = min(min_length, length(i))
                  max_length = max(max_length, length(i))
 
                  min_jd = min(min_jd, epoch(i))
                  max_jd = max(max_jd, epoch(i))
 
              end if
          end if
 
****      See if we should output to edit file
          if( edit_flag(i).ne.0 ) then
              if( .not.outline ) then
                  write(edit_unit,'(a)') baseline(1:trimlen(baseline))
                  outline = .true.
              end if
 
              jd = epoch(i)
              call mjd_to_ymdhms(jd, date, sectag)
              date(1) = date(1) - 1900
              write(edit_unit,100) date, length(i), sigma(i), residual,
     .                       residual/sigma(i), edit_flag(i)
  100         format(5i3,1x,f14.4,1x,f8.4,1x,f8.4,1x,f8.2,1x,i3)
          else
*                                          ! Increment
              used_data = used_data + 1
          end if
      end do
 
****  Set the start and stop time
      range = max_jd - min_jd
*                                          ! .5 days added in incase no
      min_jd = min_jd - 0.02*range - 0.5
*                                          ! data.
      max_jd = max_jd + 0.02*range + 0.5
 
****  Thats all
      return
      end
 
CTITLE WRITE_DATA
 
      subroutine write_data(epoch, length, sigma, edit_flag,
     .                      num_data, nom_length, title, baseline)
 
 
*     Routine to write out the data file for plot
 
*   num_data        - Number of data read
*   i,j,k           - Loop control
*   date(5)         - ymdhm read from file or elsewhere
*   trimlen         - Returns length of string
*   edit_flag(1)    - Edit flags for data
*   len_buffer      - Length of buffer read from file
 
      integer*4 num_data, i, date(5), trimlen, edit_flag(1)
 
*   sectag          - Seconds tag
*   jd              - jd of data
*   epoch(1)        - Epochs of data points
*   length(1)       - Length read from values file
*   sigma(1)        - Sigmas of data read from file
*   value           - value read from file
*   nom_length      - Nominal length (m)
 
      real*8 sectag, jd, epoch(1), length(1), sigma(1), value,
     .    nom_length
 
*             title - Title of plot
*   baseline        - Name of baseline
 
 
      character*(*) title, baseline
 
***** Loop over the data writing out non--edited values
 
      write(202,'(a)') title(1:max(1,trimlen(title)))
      write(202,'(a)') baseline(1:trimlen(baseline))
      write(202,'(a)') ' '
 
****  Now loop
      do i = 1, num_data
          jd = epoch(i)
          call mjd_to_ymdhms(jd, date, sectag)
          date(1) = date(1) - 1900
          if( edit_flag(i).eq.0 ) then
              value = length(i) - nom_length
              write(202,100) date, value, sigma(i)
  100         format(5i3,1x,f12.4,1x,f8.4)
          end if
      end do
 
***** Thats all
      return
      end
 
CTITLE WRITE_PLOT
 
      subroutine write_plot( ymin, ymax, y_grad, start_jd, end_jd,
     .                       num_data, min_jd, max_jd, plot_data)
 
 
*     Routine to write out the plot control file
 
*   datestart(5)    - Start date
*   dateend(5)      - End date
*   num_data        - Number of data points
 
      integer*4 datestart(5), dateend(5), num_data
 
*   start_jd        - start jd if given in runstring
*   end_jd          - end jd if given in runstring
*   y_grad          - step to be used on plot (if given in
*                   - runstring)
*   ymin, ymax      - ymin and ymax needed for plot
*   min_jd, max_jd  - Start and stop times (used to position
*                   - labels)
*   step            - Spacing between labels
*   beg             - Start of label positions
*   sectag          - Seconds tag
 
      real*8 start_jd, end_jd, y_grad, ymin, ymax, min_jd, max_jd,
     .    step, beg, sectag
 
*   plot_data       - Name of the plot data file for this
*                   - baseline
 
      character*(*) plot_data
 
* LOCAL VARIABLES
*   trimlen         - Length of string.
 
      integer*4 trimlen
 
****  Start outputing the control file
 
      write(201,*) '  File ', plot_data(1:trimlen(plot_data))
 
      write(201,*) '  X_field 0 1 5'
      write(201,*) '  Y_field 1 6 7 "Residual (m)"'
      write(201,*) '  Read'
      write(201,*) '  Y_Scale ', ymin, ymax
 
      call mjd_to_ymdhms(min_jd, datestart, sectag)
      datestart(1) = datestart(1) - 1900
      call mjd_to_ymdhms(max_jd, dateend, sectag)
      dateend(1) = dateend(1) - 1900
      write(201,*) '  X_scale ', datestart, dateend
 
      write(201,*) '  point 1'
      write(201,*) '  draw'
C     write(201,*) '  axes'
 
      write(201,*) '  Poly years mm 365.25 1000.'
      write(201,*) '  fit 0 '
      write(201,*) '  Pdraw '
      if( num_data.gt.2 ) then
          write(201,*) ' Fit 1'
          write(201,*) ' Pdraw'
      end if
 
      step = (ymax-ymin)/35
*                                  ! 5% in from edge
      beg = (max_jd - min_jd)/20
 
      write(201,*) ' Label ',beg, ymax-ymin-step,' 1 0 :h1'
      write(201,*) ' Label ',beg, ymax-ymin-2*step,' 1 0 :h2'
      write(201,*) ' Label ',beg, 4*step,' 1 0 :p1'
      write(201,*) ' Label ',beg, 3*step,' 1 0 :p2'
      write(201,*) ' Label ',beg, 2*step,' 1 0 :p3'
      write(201,*) ' Label ',beg, step  ,' 1 0 :p4'

*     Now rescale and write
      write(201,*) ' y_scale ', (ymin-ymax)/2, (ymax-ymin)/2
      write(201,*) ' axes '
 
 
***** Thats all
      return
      end
 
CTITLE OPEN_PLOT_DATA
 
      subroutine open_plot_data( unit, baseline, root, plot_data )
 
*     This routine will open the data file in which the baseline lengths
*     for the current baseline wil be saved.  (The name of the data file
*     will be generated from root and baseline

* unit          - Unit number for output
  
      integer*4  unit
 
*   baseline        - Second record from values file group, containd
*               - the baseline
*   root        - Root of the data file name
*   plot_data   - The concatinated name of the data file
 
      character*(*) baseline, root, plot_data
 
* LOCAL VARIABLES
 
*   indx        - Pointer in string
*   ierr        - IOSTAT error flag
*   trimlen     - Length of string
 
      integer*4 indx, ierr, trimlen
 
*   site1, site2    - Will contain the site names in the baselines
 
      character*8 site1, site2

*   soln            - Solution number read for this solution

      character*4 soln
 
****  Strip of the baseline name
 
      indx = 1
      call getword( baseline, site1, indx )
      call getword( baseline, site2, indx )
*     Now really get site two
      call getword( baseline, site2, indx )

*     Now get the solution number
      call getword( baseline, soln, indx)
      call getword( baseline, soln, indx)
 
*     Now construct file name
 
      plot_data = root(1:trimlen(root)) // site1(1:trimlen(site1)) //
     .    '_' // site2(1:trimlen(site2)) // '.dat' // 
     .           soln(1:max(1,trimlen(soln)))
 
      open(202, file=plot_data, iostat=ierr, status='unknown')
      call report_error('IOSTAT',ierr,'open',plot_data,
     .                   1,'open_plot_data')
 
****  Thats all
      return
      end
