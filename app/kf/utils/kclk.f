 
      program kclk
 
*     This program will read through a KalObs file and save the values
*     of any clock jumps which occur.  The program requires as a KalObs
*     file and a file containing the list of jumps.  The jumps are
*     in ns (converted internally to ps) and are incrementally added
*     to existing values.  If CLEAR is given as the jump value then
*     current values will be cleared.
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
 
*         max_kclk  - Maximum number of values allowed
 
      integer*4 max_kclk
 
      parameter ( max_kclk = 100 )
 
 
*   i,j,k       - Loop counters
*   ierr        - Fmp file error flag
*   len_runstring   - Length of runstring
*   rcpar       - HP runstring reading utility
*   unit_in     - Input unit number.  If nothing is passed than std
*                 input
*                 is used.
 
 
      integer*4 i,j, ierr, len_runstring, rcpar, unit_in
 
*   kclk_site(max_kclk) - Sites at which clock offsets have
*       - been entered.
*   num_kclk    - Number of clock offsets entered
*   is      - Site number in baseline to which jump referrs
*           - (Either 1 or 2. If -1 then not this baseline)
 
      integer*4 kclk_site(max_kclk), num_kclk, is
 
*   val_kclk(3,max_kclk)  - Values of the clock offsets entered
*           - rate, and acceleration (ps,ps/dy, ps/dy**2)
*   epoch_kclk(max_kclk)    - epochs to which values are applied.
*           - (all times greater than equal to this value)
*   dt      - Time since start time of the experiment (start_epoch).  
*             This should be the same time as used by PLTSL for fitting
*             polynomials.
 
      real*8 val_kclk(3,max_kclk), epoch_kclk(max_kclk), dt
 
*   clear_kclk(max_kclk)    - Logical indicating that existing
*           - values should cleared.(Allowed for each value so
*           - value so that all sites can be cleared at once.
*   interactive - True if this is an interactive session
*   write_out   - Indicates that we should write out KalObs
*           - data record.
 
 
      logical clear_kclk(max_kclk), interactive, write_out
 
*   KO_name     - Name of the Kalobs file to be updated
*   In_name     - Name of the input file.
 
 
      character*132 KO_name, in_name
 
***** Get the Kalobs file name
 
      len_runstring = rcpar(1, KO_name )
 
*                                     ! Name not given, STOP
      if( len_runstring.eq.0 ) then
          call proper_runstring('kclk.hlp', 'kclk', 1)
          stop ' kclk: Incorrect Runstring, add KalObs name'
      end if
 
*     Get the name of the input file (if none given assumed 5)
      len_runstring = rcpar(2, in_name)
      if( len_runstring.le.0 ) in_name = '5'
 
      unit_in = 100
      call open_lu(unit_in, in_name, ierr, 'old')
      call report_error('IOSTAT',ierr,'open',in_name,1,'kclk')
 
****  See if interactive session (by unit_in not being 100)
      if( unit_in.eq.100 ) then
          interactive = .false.
      else
          interactive = .true.
      end if
 
***** Now open the KalObs file
 
      call open_KalObs( KO_name, ierr )
      call report_error('FmpOpen',ierr,'open',KO_name, 1, 'kclk')
 
***** Read header
 
      call RW_KalObs_header('R', ierr)
      call report_error('Fmpread',ierr,'read',KO_name, 1, 'kclk')
 
***** Now start reading the command file
 
      call read_kclk( unit_in, max_kclk, num_kclk, epoch_kclk,
     .        kclk_site, val_kclk, clear_kclk)
 
      if( num_kclk.eq.0 ) then
          write(*,*) ' KCLK: No clock jumps entered'
          stop ' KCLK: Terminated'
      endif
 
***** Now loop over data file, updating as we go
 
      do i = 1, num_obs
 
          call RW_KalObs_block('R','DATA',site, ierr, i )
          call report_error('FmpRead',ierr,'read','data block',1,
*                                     ! Kill if error
     .                      'kclk')
 
*         See if this baseline should be updated
          write_out = .false.
          do j = 1, num_kclk
              is = -1
              if( site(1).eq.kclk_site(j) ) is = 1
              if( site(2).eq.kclk_site(j) ) is = 2
 
              if( epoch.ge. epoch_kclk(j) .and.is.gt.0 ) then
 
*                 See if we should clear first
                  if( clear_kclk(j) ) then
                      clock_jump(is,1) = 0
                      clock_jump(is,2) = 0
                  else
                      dt = epoch - start_epoch
                      clock_jump(is,1) = clock_jump(is,1) + 
     .                                   val_kclk(1,j) +
     .                                   val_kclk(2,j)*dt +
     .                                   val_kclk(3,j)*dt**2
                      clock_jump(is,2) = clock_jump(is,2) +
     .                                   val_kclk(2,j)/186.4d0 +
     .                                 2*val_kclk(3,j)*dt/186.4d0
                  end if
 
                  write_out = .true.
              end if
          end do
 
*         write out data record
          if( write_out ) then
              call RW_KalObs_block('W','DATA',site, ierr,i)
              call report_error('FmpWrite',ierr,'writ','data block',1,
*                                     ! kill of error
     .                'kclk')
          end if
 
      end do
 
***** Update data notes and write the header back out
*     Say that clock jumps have been added.
      call sbit( data_notes,10,1)
      call increment_KalVer 
      call RW_KalObs_header('W', ierr)
      call close_KalObs( ierr )
 
***** Thats all
      end
 
CTITLE READ_KCLK
 
      subroutine read_kclk( unit_in, max_kclk, num_kclk, epoch_kclk,
     .        kclk_site, val_kclk, clear_kclk)
 
*     This routine will read the input file and get the jump values from
*     this file.
 
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
 
*         max_kclk  - Maximum number of values allowed
 
      integer*4 max_kclk
 
*   ierr        - Fmp file error flag
*   unit_in     - Input unit number.  If nothing is passed than std
*                 input
*                 is used.
 
 
      integer*4 ierr, unit_in
 
*   kclk_site(max_kclk) - Sites at which clock offsets have
*       - been entered
*   indx    - Pointer in string
*   iel     - Site number pulled from string
*   date(5) - Date and hr min at which the offset should be
*           - applied
*   num_kclk    - Number of clock offsets entered
 
      integer*4 kclk_site(max_kclk), indx, iel, date(5), num_kclk,
     .          trimlen, jndx, jerr
 
*   val_kclk(3,max_kclk)  - Values of the clock offsets entered
*           - rate and acceration. (ps,ps/dy. ps/day**2)
*   value   - Value in ns from control file (ns)
*   rate    - Rate in ns/day
*   acc     - Acceration in ns/day**2
*   epoch_kclk(max_kclk)    - epochs to which values are applied.
*           - (all times greater than equal to this value)
*   tepoch  - Temporary epoch (used until we sure all is OK with
*           - line)
*   sec_tag - Seconds tag for the ymdhms to jd conversion.
 
      real*8 val_kclk(3,max_kclk), value, epoch_kclk(max_kclk), tepoch,
     .    sec_tag, rate, acc
 
*   clear_kclk(max_kclk)    - Logical indicating that existing
*           - values should cleared.(Allowed for each value so
*           - value so that all sites can be cleared at once.
*   tclear  - Temporary storage for clear status
*   interactive - True if this is an interactive session
 
 
      logical clear_kclk(max_kclk), tclear, interactive
 
*   line        - Line read from input file
*   cvalue      - Values number as a string (used to check for clear)
 
 
      character*132 line, cvalue
 
****  If we are interactive write header
 
      if( interactive ) then
          write(*,100)
 100      format(/,' KCLK: Program to add clock jump values to KalObs'
     .            ' file.',/,
     .            ' Input is <site name > yy mm dd hh mm [Value'
     .            ' (ns)/CLEAR]')
      end if
 
****  Now loop until end of file
      num_kclk = 0
      ierr = 0
      do while ( ierr.eq.0 )
 
          if( interactive ) then
              write(*,'(a,$)') '? '
          end if
 
          read(unit_in, '(a)', iostat=ierr ) line
          call casefold ( line )
          if( ierr.eq.0 .and. trimlen(line).gt.0 .and.
*                                         ! OK so continue
     .        line(1:1).eq.' ' ) then
 
*             Get the site name
              indx = 1
              call get_cmd( line, site_names, num_sites,
     .                        iel, indx )
 
*                                     ! Site found
              if( iel.gt.0 ) then
                  call multiread(line, indx, 'I4', ierr, date,
     .                        cvalue, 5)
 
                  sec_tag = 0.d0
                  call ymdhms_to_jd( date, sec_tag, tepoch)
 
*                 Now pull off the value.  Save current pointer in
*                 case we need to read value from line
                  jndx = indx
                  call getword( line, cvalue, indx)
                  call casefold( cvalue )
 
*****             See if it says clear
*                                                 ! OK close enough
                  if( cvalue(1:1).eq.'C' ) then
                      tclear = .true.
                      value = 0.d0
                  else
                      call read_line( line, jndx, 'R8', ierr, value,
     .                    cvalue)
                      tclear = .false.

*                     Now see if the rate and acceration have been given
                      call read_line( line, jndx, 'R8', jerr, rate,
     .                     cvalue)
                      if( jerr.eq.0 ) then
                           call read_line( line, jndx, 'R8', jerr,
     .                          acc, cvalue)
                           if( jerr.ne.0 ) acc = 0.d0
                      else
                           rate = 0.d0
                      end if
                  end if
 
*****             Now if all OK, then svae values
                  if( ierr.eq.0 ) then
                      num_kclk = num_kclk + 1
                      if( num_kclk.gt.max_kclk ) then
                          write(*,150) max_kclk
 150                      format(/' Too many clock jumps added.'
     .                        ' Maximum allowed is ',i4)
                          stop ' KCLK: Terminated, Too many jumps'
                      end if
 
*                     Now save values
                      epoch_kclk(num_kclk) = tepoch
                      val_kclk(1,num_kclk) = value*1000.d0
                      val_kclk(2,num_kclk) = rate *1000.d0
                      val_kclk(3,num_kclk) = acc  *1000.d0
                      clear_kclk(num_kclk) = tclear
                      kclk_site(num_kclk) = iel
                  end if
*                             ! Site not found
              else
                  write(*,200) line(1:trimlen(line))
 200              format(' Site not found in line: ',a)
              end if
 
*                             ! Read OK
          end if
*                             ! Looping over input
      end do
 
***** Thats all
      return
      end
 
 
 
