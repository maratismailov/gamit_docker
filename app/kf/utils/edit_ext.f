      program edit_ext

      implicit none 
 
*     This program will edit date and quanitity files generated by
*     extract.  The runstring for the program is:
*
*     % edit_ext [edit file] <inputs .... >
*
*     where [edit_file] is a file containing a list of dates to be
*         "editted", acheived by placing * in column one.
*     and <inputs ... > is a list of file to be operated on.  These
*         files are first renamed (to orgininal name with .org extent)
*         and then lines with matching dates have a * placed in column
*         one.
*
*   max_edits   - Maximum number of edits allowed
 
      integer*4 max_edits
 
      parameter ( max_edits = 1000 )
 
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
 
      integer*4 ierr, jerr, rcpar, len_run, trimlen, i,j, nr, date(5),
     .    num_edits, ndel, indx, fmprename
 
*   bad         - Indicates if a data point is bad or not.
 
      logical bad
 
*   edit_jds(max_edits) - The list of julian dates to be
*           - to be edited.
*   sectag  - Seconds tag for JD computation (set to zero)
*   jd      - Julian date of current line.
 
 
      real*8 edit_jds(max_edits), sectag, jd
 
*   edit_file   - Name of file containing list of edit times.
*   input_file  - Read from runstring
*   output_file - Read from runstring
 
      character*128 edit_file, input_file, output_file
 
*   line        - Line read from input.
 
      character*256 line
 
*   cd          - Dummy string from read_line and multiread.
 
      character*8 cd
 
****  Decode the runstring
      len_run = rcpar(1, edit_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'edit_ext', 'edit_ext', 1)
      end if
      open(100, file=edit_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',edit_file,1,
     .                'edit_ext')
 
***** Now read the edit file
 
      sectag = 0.d0
      j = 0
 
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr ) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
*                                                 ! Get the time
     .        trimlen(line).gt.0 ) then
*                                             ! from the line
 
              indx = 1
              call multiread(line, indx, 'I4',jerr, date, cd, 5)
              if( jerr.eq.0 .and. j.lt.max_edits) then
                  j = j + 1
                  call ymdhms_to_jd( date, sectag, edit_jds(j))
              end if
          end if
      end do
 
****  Print out message.
      num_edits = j
      write(*,100) num_edits, edit_file(1:trimlen(edit_file))
  100 format(/' EDIT_EXT: There are ',i4,' edits epochs in ',a)

      if( num_edits.eq.0 ) stop ' No data to edit'
 
*     Now start looping over the inputs.
      nr = 1
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
     .            0, 'edit_ext')
          end if
 
*         If still OK, open input and create the output
          if( ierr.eq.0 ) then
              open(101,file=input_file, iostat=ierr, status='old')
              open(200,file=output_file, iostat=ierr, status='new')
              call report_error('IOSTAT',ierr,'creat',output_file,
     .            0, 'edit_ext')
          end if
 
****      Now loop over the input file
          ndel = 0
          do while ( ierr.eq.0 )
              read(101,'(a)', iostat=ierr) line
              jerr = -1
*                                     ! Process
              if( ierr.eq.0 ) then
                  if( trimlen(line).eq.0 .or.
*                                                     ! Just write
     .                line(1:1).ne.' '       ) then
                      write(200,'(a)', iostat=ierr)
     .                    line(1:max(1,trimlen(line)))
                      call report_error('IOSTAT',ierr,'writ',line,0,
     .                    'edit_ext')
 
*                                     ! See if we should edit.
                  else
 
*                     Get date from line
                      indx = 1
                      call multiread( line, indx, 'I4', jerr, date,
     .                    cd, 5)
                  end if
              end if
 
*             If no error so far continue
              if( jerr.eq.0 ) then
                  call ymdhms_to_jd( date, sectag, jd )
 
*                 Loop over edit times and see if match
                  bad = .false.
                  do i = 1, num_edits
*                                                             ! edit
                      if( abs(jd-edit_jds(i)).lt.1.3d-3 ) then
                          bad = .true.
                      end if
                  end do
 
*                 If point is bad, then add * to column one
                  if( bad ) then
                      line(1:1) = '*'
                      ndel = ndel + 1
                  end if
 
*                 Now write out the line
                  write(200,'(a)', iostat=ierr)
     .                    line(1:max(1,trimlen(line)))
                  call report_error('IOSTAT',ierr,'writ',line,0,
     .                            'edit_ext')
 
*                         ! No error on read
              end if
*                         ! Looping over data set.
          end do
 
*         Tell user what is going on
          if( len_run.gt.0 .and. ierr.eq.-1 ) then
              write(*,200) nr-1, output_file(1:trimlen(output_file)),
     .                    ndel
  200         format(I4,'. File: ',a,'. ',i4,' Data edited')
              close(200)
              close(101, iostat=ierr, status='delete')
              call report_error('IOSTAT',ierr,'delet',input_file,
     .                          0,'edit_ext')
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
 
 
 
 
 
 
 
 
 
 