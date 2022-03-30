      program edit_val

      implicit none 
 
*     This program will edit a values file generated with local and
*     global h-files by removing the duplicate entry associated with
*     global file.  The tolerance on the day match is 1 day.
*
*     % edit_val [sigma limit] <inputs .... >
*
*     where [sigma limit] is the largest sigma that will be output
*         <inputs ... > is a list of file to be operated on.  These
*         files are first renamed (to orgininal name with .org extent)
*         and then lines with matching dates have a * placed in column
*         one.
*
*   ierr, jerr  - IOSTAT errors
*   rcpar   - Reads the runstring
*   len_run - Length of runstring
*   trimlen - Length of string
*   fmprename - Renames a file (Libhp1000 emulation)
*   i,j     - Loop counters
*   nr      - Sequential number in runstring
*   date(5) - Date of determination (assummed 0 hrs)
*   num_edits   - Number of edit times actually read from input.
*   ndel    - Mumber of deleted data.
*   indx    - Pointer in string
 
      integer*4 ierr, jerr, rcpar, len_run, trimlen, i, nr, date(5),
     .         ndel, indx, fmprename, kerr
 
*   prev_written - Logical to indicate previous line written
 
      logical prev_written 
 
*   sectag  - Seconds tag for JD computation (set to zero)
*   jd      - Julian date of current line.
*   vals(2) - Current values
*   prev_jd - Julian date on record
*   prev_vals(2) - Previous records values 
*   sigma_limit  - Largest sigma to be output
 
      real*8 sectag, jd, vals(2), prev_jd, prev_vals(2), sigma_limit
 
*   input_file  - Read from runstring
*   output_file - Read from runstring
 
      character*128 input_file, output_file
 
*   line        - Line read from input.
*   prev_line   - Previous line from file (pending output)
 
      character*256 line, prev_line
 
*   cd          - Dummy string from read_line and multiread.
 
      character*8 cd

      prev_vals(1) = 0.d0
      prev_vals(2) = 0.d0
      prev_jd      = 0.d0
 
****  Decode the runstring
      len_run = rcpar(1, line)
      if( len_run.le.0 ) then
          call proper_runstring('edit_val.hlp', 'edit_val',1)
      else
          read(line,*, iostat=ierr) sigma_limit
          write(*,100) sigma_limit
 100      format(' EDIT_VAL: Removing duplicates from values file',/,
     .           '           Sigma limit ',f8.4,' m')
          call report_error('IOSTAT',ierr,'decod',line,1,'edit_val')
      end if

*     Now start looping over the inputs.
      nr = 1
      len_run = 1
      do while ( len_run.gt.0 )
          nr = nr + 1
          len_run = rcpar(nr, input_file )
 
*         If we have something process
          ierr = -99
          if( len_run.gt.0 ) then
 
*             First rename input.  Generate new name.
              output_file = input_file
              input_file(len_run+1:) = '.org'
              ierr = fmprename(output_file, input_file, ' ')
              call report_error('IOSTAT',ierr,'renam',output_file,
     .            0, 'edit_val')
          end if
 
*         If still OK, open input and create the output
          if( ierr.eq.0 ) then
              open(101,file=input_file, iostat=ierr, status='old')
              open(200,file=output_file, iostat=ierr, status='new')
              call report_error('IOSTAT',ierr,'creat',output_file,
     .            0, 'edit_val')
          end if
 
****      Now loop over the input file
          ndel = 0
          prev_written = .true.
          do while ( ierr.eq.0 )
              read(101,'(a)', iostat=ierr) line
              jerr = -1
*                                     ! Process
              if( ierr.eq.0 ) then
                  if( trimlen(line).eq.0 .or.
*                                                     ! Just write
     .                line(1:1).ne.' '       ) then
                      if( .not.prev_written )
     .                write(200,'(a)', iostat=ierr)
     .                      prev_line(1:max(1,trimlen(prev_line)))
*                     write the current line
                      prev_written = .true.
C                     write(200,'(a)', iostat=ierr)
C    .                    line(1:max(1,trimlen(line)))
                      call report_error('IOSTAT',ierr,'writ',line,0,
     .                    'edit_val')
                      jerr = -1
 
*                                     ! See if we should edit.
                  else
 
*                     Get date from line
                      indx = 1
                      jerr = 0
                      do i = 1,5
                          call read_line( line, indx, 'I4', kerr, 
     .                                    date(i), cd)
                          if( kerr.ne.0 ) jerr = -1
                      end do
                  end if
              end if
 
*             If no error so far continue
              if( ierr.ne.0 ) jerr = ierr
              if( jerr.eq.0 ) then
                  call ymdhms_to_jd( date, sectag, jd )
                  call multiread(line, indx, 'R8',jerr,vals,cd,2)

                  if( vals(2).gt.sigma_limit ) then
                      if( .not.prev_written )
     .                write(200,'(a)', iostat=ierr)
     .                  prev_line(1:max(1,trimlen(prev_line)))
                      prev_written = .true.
                      ndel = ndel + 1
                  end if

*                 Now check to see if we are < 1 day away from 
*                 previous value
                  if( abs(jd-prev_jd).lt. 0.99 .and. 
     .                vals(2).le.sigma_limit ) then
*                     values are with in a day.  See if the estmates
*                     are the same
                      if( abs(vals(1)-prev_vals(1)).lt.0.001 ) then
*                         The entries seem to be same.  Write out the
*                         one with the smallest sigma
                          if( vals(2).lt.prev_vals(2)) then
                              if( vals(2).le.sigma_limit)
     .                        write(200,'(a)', iostat=ierr)
     .                                line(1:max(1,trimlen(line)))
                              prev_jd = jd
                              prev_vals(1) = vals(1)
                              prev_vals(2) = vals(2)
                              prev_line = line
                              prev_written = .true.
                              ndel = ndel + 1
                          else
                              if( .not.prev_written .and.
     .                             prev_vals(2).le.sigma_limit)
     .                        write(200,'(a)', iostat=ierr)
     .                          prev_line(1:max(1,trimlen(prev_line)))
                              prev_written = .true.
                          end if
                       else
*                         values differ.  Write the previous line if
*                         it have not been written
                          if( .not.prev_written )
     .                    write(200,'(a)', iostat=ierr)
     .                          prev_line(1:max(1,trimlen(prev_line)))
*                         Now save the current values
                          prev_jd = jd
                          prev_vals(1) = vals(1)
                          prev_vals(2) = vals(2)
                          prev_line = line
                          prev_written = .false.
                      end if
                  else
*                     The date has changed.  If previous value not
*                     written write now and save the values.
                      if( .not.prev_written )
     .                write(200,'(a)', iostat=ierr)
     .                      prev_line(1:max(1,trimlen(prev_line)))
*                     Now save the current values if sigma limit OK
                      if( vals(2).le.sigma_limit ) then
                          prev_jd = jd
                          prev_vals(1) = vals(1)
                          prev_vals(2) = vals(2)
                          prev_line = line
                          prev_written = .false.
                      end if
                  end if
              else
*                 Error decoding the data (probably a comment line)
                  if( .not.prev_written )
     .            write(200,'(a)', iostat=jerr)
     .                  prev_line(1:max(1,trimlen(prev_line)))
*                 write the current line
                  write(200,'(a)', iostat=jerr)
     .                  line(1:max(1,trimlen(line)))
                  prev_written = .true.
              end if
*                         ! Looping over data set.
          end do
 
*         Tell user what is going on
          if( len_run.gt.0 .and. ierr.eq.-1 ) then
              write(*,200) nr-1, output_file(1:trimlen(output_file)),
     .                    ndel
  200         format(I4,'. File: ',a,'. ',i4,' Data edited')
              close(200)
C             close(101, iostat=ierr, status='delete')
              close(101, iostat=ierr)
              call report_error('IOSTAT',ierr,'delet',input_file,
     .                          0,'edit_val')
          else
              if( len_run.gt.0 ) then
                  write(*,220) nr-1,
     .                output_file(1:trimlen(output_file))
  220             format(I4,'. File: ',a,'. PROBLEM WITH NEW FILE')
              end if
          end if
*                     ! Looping over the runstring
      end do
 
****  Thats all
      end
 
 
 
 
 
 
 
 
 
 
