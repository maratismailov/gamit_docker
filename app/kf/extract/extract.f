 
      Program extract
 
 
*     This program will read any solution (ASCII) file and extract
*     information from the solution based on the field information
*     provided by the user.  See EXTRACT::HELP for details of the
*     commands.
*
*                                 11:02 PM WED., 26 Aug., 1987
 
      include '../includes/extract_comm.h'
 
*   iel     - Command number
*   ierr    - IOSTAT error
*   indx    - Current position in command buffer
*   trimlen - HP function for length of string
 
      integer*4 iel, ierr, indx, trimlen
 
*   eof     - Indicates end of file or END command given
 
      logical eof
 
*   buffer  - Line read from command file
 
      character*120 buffer
 
***** Decode the runstring, and open the files if they are given
 
      call init_extract
 
      call decode_extract_run
 
*     Now loop over the command file given by the user, and
*     process each command
 
      eof = .false.
      do while ( .not. eof )
 
          read(unit_comm,'(a)',iostat=ierr) buffer
          indx = 1
 
*                                     ! Check for end of file
          if( ierr.eq.-1 ) then
              eof = .true.
*                                     ! See if error
          else
              call report_error('IOSTAT',ierr,'read',command_file,
*                                             ! Kill if errror
     .                          1,'EXTRACT')
 
*             If there was no error then decode and process command
              iel = 0
              if( buffer(1:1).eq.' ' ) then
 
                  call get_cmd( buffer, extract_commands,
     .                 num_extract_commands, iel, indx )
              end if
 
              if( iel.gt.0 ) then
                  call process_extract_com( buffer, indx, iel, eof)
              else
*                                          ! Report unknown command
                  if( iel.eq.-1 ) then
                      write(*,100) buffer(1:max(1,trimlen(buffer)))
  100                 format(' Unkown command: ',a)
                  end if
              end if
*                         ! No read error
          end if
*                         ! Loop over command file
      end do
 
***** Thats all
      end
 
CTITLE DECODE_EXTRACT_RUN
 
      subroutine decode_extract_run
 
 
*     Routine to decode runstring:
*     CI> EXTRACT command_file, <input file>, <output file>
*     where command_file is name of command file (and must be
*                        given (May be an LU)
*           input file is the name of the file to be read and the
*                        information extracted from
*           output file is the name of the output file
*
 
      include '../includes/extract_comm.h'
 
*   ierr        - IOSTAT error opening files
*   len_run     - Length of runstring
*   RCPAR       - HP function to read runstring
 
      integer*4 ierr, len_run, RCPAR
 
***** Get the command file name
 
      len_run = rcpar(1, command_file )
      if( len_run.gt.0 ) then
          unit_comm = 100
          call open_lu( unit_comm, command_file, ierr, 'old' )
          call report_error('IOSTAT',ierr,'open',command_file,1,
*                                                 ! Kill if error
     .                      'DECODE_EXTRACT_RUN')
      else

          call proper_runstring('extract.hlp','extract',-1) 
c         write(*,100)
c 100     format(/' EXTRACT: Incorrect runstring',/,
c    .            ' CI> EXTRACT, command_file, <Input file>,'
c    .            ' <output file>',/,
c    .            ' See EXTRACT::HELP for more information'/)
c         stop ' EXTRACT Terminated: Incomplete runstring'
      end if
 
***** Try to get the input file
 
      len_run = rcpar(2, input_file )
      if( len_run.gt.0 ) then
          unit_in = 101
          call open_lu( unit_in, input_file, ierr, 'old')
          call report_error('IOSTAT',ierr,'open',input_file,1,
*                                                 ! Kill if error
     .                      'DECODE_EXTRACT_RUN')
      else
          input_file = ' '
      end if
 
***** Try to get the output file
 
      len_run = rcpar(3, output_file )
      if( len_run.eq.0 ) then
*                             ! Use terminal
          output_file = '6'
      end if
 
*                                         ! Check for wild card
      if( input_file(1:1).ne.' ' ) then
          call wild_card( output_file, input_file )
      end if
      unit_out = 200
      call open_lu(unit_out, output_file, ierr, 'unknown') 
      call report_error('IOSTAT',ierr,'open',output_file,1,
*                                             ! Kill if error
     .                  'DECODE_EXTRACT_RUN')
 
***** Thats all
      return
      end
 
CTITLE PROCESS_EXTRACT_COM
 
      subroutine process_extract_com( buffer, indx, iel, eof)
 
 
*     This routine does all the work.  It will process the commands
*     as they are recieved.
 
      include '../includes/extract_comm.h'
 
*   iel     - Command number
*   ierr    - IOSTAT error
*   indx    - Current position in command buffer
*   num     - General number value (usually for field number)
*   snum    - Another general number used
 
      integer*4 iel, ierr, indx, num, snum
 
*   eof     - Indicates end of file or END command given
 
      logical eof
 
*   cdumm   - Dummy argument from READ_LINE
 
      character*1 cdumm
 
*   buffer  - Line read from command file
 
      character*(*) buffer
 
***** Use a computed goto to get commands
 
      GOTO (  100,  200,  300,  400,  500,  600,  700,  800,  900, 1000,
     .       1100  ) iel
 
***** END   :  Set EOF and exit
  100 continue
          eof = .true.
          return
 
***** INPUT : Name of input file
  200 continue
          call read_line( buffer, indx,'CH',ierr,iel, input_file)
 
*         Try to open
          unit_in = 101
          call open_lu(unit_in, input_file, ierr, 'old')
          call report_error('IOSTAT',ierr,'open',input_file,0,
     .                      'PROCESS_EXTRACT_COM')
          if( ierr.ne.0 ) then
              input_file = ' '
          end if
          return
 
***** OUTPUT: Name of output file
  300 continue
          call read_line( buffer, indx,'CH',ierr,iel, output_file)
 
*         Check for wild cards
*                                           ! Input name given so check
          if( input_file(1:1).ne.' ' ) then
*                                           ! for wild card
              call wild_card( output_file, input_file )
          end if
 
*         Try to open
          unit_out = 200
          call open_lu(unit_out, output_file, ierr, 'unknown')
          call report_error('IOSTAT',ierr,'open',output_file,0,
     .                      'PROCESS_EXTRACT_COM')
          if( ierr.ne.0 ) then
              output_file = ' '
          end if
 
          return
 
***** DESCRIPT: Description command for each field. Form is
*     DESCRIPT <field #> "label".  Label must in double quotes ".
*
  400 continue
 
          call read_line(buffer, indx, 'I4',ierr, num, cdumm)
          call report_error('IOSTAT',ierr,'read',buffer,0,
     .                      'PROCESS_EXTRACT_COM')
 
          if( ierr.eq.0 ) then
              call decode_quotes(buffer,indx,field_title(num))
          end if
          return
 
***** TITLE :title for top of output file. Form
*     TITLE "label" or TITLE number, where number is line number
*     from input file
  500 continue
 
*         Try to Read_line first
          call read_line(buffer,indx,'I4',ierr, title_num,cdumm)
*                                 ! Get string between quotes
          if( ierr.ne.0 ) then
              call decode_quotes(buffer,indx,title_label)
          end if
          return
 
***** FIELD command
  600 continue
          call decode_field_com( buffer,indx )
          return
 
***** BEGIN  command.  Sets the start fields which must be present before
*     the field is searched for
  700 continue
 
*         Get field number
          call read_line(buffer,indx,'I4', ierr, num, cdumm)
          call report_error('IOSTAT',ierr,'read',buffer,0,
     .                      'PROCESS_EXTRACT_COM')
 
*         Now get the number of the START ie. 1 for first, 2 for
*         second string to be found after the first has been found
          if( ierr.eq.0 ) then
              call read_line(buffer,indx,'I4',ierr, snum, cdumm)
              call report_error('IOSTAT',ierr,'read',buffer,0,
     .                          'PROCESS_EXTRACT_COM')
          end if
 
*         Now get the label to be found
          if( ierr.eq.0 ) then
              call decode_quotes(buffer,indx,start_label(snum, num))
          end if
          return
 
***** FINISH   command.  Gives the label which will cause start found
*     status to be set false.
  800 continue
 
*         Get field number
          call read_line(buffer,indx,'I4', ierr, num, cdumm)
          call report_error('IOSTAT',ierr,'read',buffer,0,
     .                      'PROCESS_EXTRACT_COM')
 
          if( ierr.eq.0 ) then
              call decode_quotes(buffer,indx,end_label(num))
          end if
          return
 
***** NORESET command: Tells which fields should not be set to not found
*     after and output has been done.
  900 continue
 
*         Loop getting all of the values
          ierr = 0
          snum = 0
          do while ( ierr.eq.0 )
              call read_line(buffer,indx,'I4',ierr,num, cdumm)
*                                         ! Field number found
              if( ierr.eq.0 ) then
                  snum = snum + 1
*                                         ! Say to not reset
                  reset_field(num) = 1
              end if
          end do
 
*         If no values given then clear all (If the user wants to he
*         say NORESET CLEAR
          if( snum.eq.0 ) then
              do num = 1, max_fields
                  reset_field(num) = 0
              end do
          end if
 
          return
 
***** OUTFORM  :  Output format
 1000 continue
 
          call read_line(buffer,indx,'I4',ierr,num,cdumm)
          call report_error('IOSTAT',ierr,'read',buffer,0,
     .                      'PROCESS_EXTRACT_COM')
 
          if( ierr.eq.0 ) then
              call decode_quotes( buffer,indx, outformat(num))
          end if
          return
 
***** RUN command: reads the input file, and outputs the output file
*     according to the field information
 1100 continue
          call extract_data
 
          return
 
***** Thats all
      end
 
CTITLE DECODE_QUOTES
 
      subroutine decode_quotes( buffer, indx, string)
 
 
*     Routine to extract string between "" and save the result
*     in string
 
*   i,j         - Loop counter
*   indx        - Start character position in buffer, returns
*               - with position after quote
*   len_buffer  - Length of the buffer
*   start       - start position in buffer
*   stop        - stop position in buffer
*   trimlen     - HP length function
 
      integer*4 i,j, indx, len_buffer, start, stop, trimlen
 
*   buffer      - Line to be decoded
*   string      - Returned string
 
      character*(*) buffer, string
 
***** Start at indx and find first "
 
      string = ' '
      len_buffer = trimlen(buffer)
      i = indx
      do while ( i.lt.len_buffer )
*                                         ! First quote found
          if( buffer(i:i).eq.'"' ) then
              start = i+1
 
*             Search for second quote
              j = start
              do while ( j.le.len_buffer )
*                                                  ! Second quote found
                  if( buffer(j:j).eq.'"' ) then
                      stop = j - 1
 
*                     Save string if start and stop span a string
                      if( stop.ge.start ) then
                          string = buffer(start:stop)
                          indx = j + 1
*                                             ! To cause exit
                          j = len_buffer+2
                          i = len_buffer+2
*                                 ! Double quotes found
                      end if
*                                 ! keep looking for second "
                  else
                      j = j + 1
                  end if
*                                 ! Searching for second quote
              end do
 
*             Check we found second quote
*                                         ! Didnot find
              if( j.eq.len_buffer+1) then
                  write(*,100) buffer(indx:max(indx,trimlen(buffer)))
  100             format(' ** WARNING ** did not find second " in ',
     .                   a)
                  i = len_buffer + 2
*                                 ! Did not find second quote
              end if
 
*                                 ! Keep searching for first "
          ELSE
              i = i + 1
          END IF
*                                 ! Seaching for first quote
      end do
 
****  Make sure we found first quote
      if( i.eq.len_buffer ) then
          write(*,150) buffer(indx:max(indx,trimlen(buffer)))
  150     format(' ** WARNING ** did not find first " in ',a)
      end if
 
****  Thats all
      return
      end
 
CTITLE DECODE_FIELD_COM
 
      subroutine decode_field_com( buffer,indx )
 
 
*     Routine to decode the field command.  The form of the command
*     is:
*     FIELD # "descriptor" #_args type {FORMAT 0/1 "(format)"
*                                      {READL  0/1 entry numbers
*     OR
*     FIELD # CLEAR
*
 
      include '../includes/extract_comm.h'
 
*   fnum        - Field number
*   i,j         - Loop counters and temporary values
*   ierr        - READ_LINE error
*   indx        - Current position in buffer
*   indx_copy   - Copy of indx
*   item        - Item numbers to be used with READ_LINE
*               - decoding of the fields
*   trimlen     - HP function for length of string
 
      integer*4 fnum, i,j, ierr, indx, indx_copy, item, trimlen
 
*   cdumm       - Dummy argument for READ_LINE
 
      character*1 cdumm
 
*   arg_type    - Type of argument (CH, I2 or R8)
 
      character*2 arg_type
 
*   arg_format  - Type of file decoding scheme (either FORMAT
*               - or READL
 
      character*6 arg_format
 
*   buffer      - Line to be decoded
 
      character*(*) buffer
 
****  First get the field number
 
      call read_line(buffer,indx,'I4',ierr,fnum,cdumm)
      call report_error('IOSTAT',ierr,'decod',buffer,0,
     .                  'DECODE_FIELD_COM')
      if( ierr.ne.0 ) RETURN
 
*                                 ! Clear field information
      decode_data(fnum) = 0
 
***** Now see if CLEAR option is given
      indx_copy = indx
      call read_line(buffer,indx,'CH',ierr,i,arg_format)
      call casefold( arg_format )
 
*     See if CLEAR command given
      if( arg_format(1:5).eq.'CLEAR') then
*                                     ! Clear field information
          decode_data(fnum) = 0
          RETURN
      end if
 
***** Get back indx and keep decoding the field
      indx = indx_copy
 
*     Get the string we must search for
      call decode_quotes( buffer,indx, field_label(fnum))
 
*     Get the number of arguments in this field
 
      call read_line(buffer,indx,'I4',ierr,items_per_field(fnum),
     .               cdumm)
 
      call report_error('IOSTAT',ierr,'decod',buffer,0,
     .                  'DECODE_FIELD_COMMAND')
 
      if( ierr.ne.0 ) RETURN
 
***** Get the type of data in the field
 
      call read_line(buffer,indx,'CH',ierr,i,arg_type)
 
***** See what type of field it is
      call casefold( arg_type )
*                                         ! I*2
      if( arg_type.eq.'I4') then
          decode_data(fnum) = 1
*                                         ! Real*8
      else if ( arg_type.eq.'R8' ) then
              decode_data(fnum) = 2
*                                         ! Character
      else if ( arg_type.eq.'CH' ) then
              decode_data(fnum) = 3
      end if
 
*     See if we found any thing
*                                         ! Did not find type
      if( decode_data(fnum).eq.0 ) then
          write(*,100) arg_type, buffer(1:trimlen(buffer))
 100      format(/' Unknown argument type ',a,' from ',a)
          RETURN
      end if
 
****  See how we should read line
 
      call read_line(buffer,indx,'CH',ierr,i,arg_format)
      call casefold( arg_format )
 
*     Now get the starting position.  0 for next character after
*     field label; 1 for from start of line containing field label
 
      call read_line( buffer,indx,'I4',ierr,j,cdumm)
      call sbit( decode_data(fnum), 5, j)
 
*                                             ! Get the format
      if( arg_format(1:6).eq.'FORMAT') then
          call decode_quotes( buffer,indx, informat(fnum) )
          call sbit( decode_data(fnum), 4,1)
*                                             ! Assume readl given
      else
 
*         Now loop getting the entry numbers
          ierr = 0
          do j = 1, items_per_field(fnum)
 
              call read_line(buffer,indx,'I4',ierr,item,cdumm)
              readl_ents(j,fnum) = item
          end do
      end if
 
****  Thats all
      return
      end
 
CTITLE EXTRACT_DATA
 
      subroutine extract_data
 
 
*     This now does all the work.  It will read through the input
*     file and extract the field data and write the results to the
*     output file.  When it is finished it will rewind the input
*     file in case we want to read it again.
 
      include '../includes/extract_comm.h'
 
*   i,j,k       - Loop counters
*   ierr,jerr   - IOSTAT errors
*   trimlen     - HP function for length of string
 
      integer*4 i,j, ierr,jerr, trimlen
 
*   EOF_in      - Indicates eof on input file
*   ifbrk       - HP function set true if BR typed
 
      logical EOF_in, ifbrk
 
*   line_in     - Line number in input file
*   line_out    - Line number in output file
 
      integer*4 line_in, line_out
 
*   inline      - Line read from input file
*   outline     - Line to be output to output file
 
 
      character*256 inline, outline
 
****  Firstly get the number of fields we need to process and the numbers
*     of start labels for each field
      num_fields = 0
 
      do i = 1, max_fields
          if( decode_data(i).ne.0 ) num_fields = i
 
*         See if we have any start labels
          do j = 1, max_starts
              if( trimlen(start_label(j,i)).gt.0 ) then
                  num_starts(i) = j
              end if
          end do
      end do
 
***** Now write out the descriptors at the top of the file
 
      ierr = 0
      if( trimlen(title_label).ne.0 ) then
          write(unit_out,'(''*'',1x,a)', iostat=jerr)
     .        title_label(1:trimlen(title_label))
              call report_error('IOSTAT',jerr,'writ',output_file,0,
     .                          'EXTRACT_DATA')
*                         ! Read input file to line number need
      else
          do i = 1, title_num
              read(unit_in,'(a)',iostat=ierr) inline
              call report_error('IOSTAT',ierr,'read',input_file,0,
     .                          'EXTRACT_DATA')
          end do
 
          write(unit_out,'(''*'',1x,a)', iostat=jerr)
     .        inline(1:max(1,trimlen(inline)))
              call report_error('IOSTAT',jerr,'writ',output_file,0,
     .                          'EXTRACT_DATA')
      end if
 
***** If any errors get out now
      if( jerr.ne.0 .or. ierr.ne.0 ) RETURN
 
***** Now write out the remaining headers
 
      do i = 1, num_fields
*                                             ! Field being searched
          if( decode_data(i).ne.0 ) then
*                                                     ! Output user
              if( trimlen(field_title(i)).gt.0 ) then
                  write(unit_out,'(a,1x,i3,1x,a)', iostat=jerr)
     .                '*',items_per_field(i),
     .                field_title(i)(1:trimlen(field_title(i)))
              else
                  write(unit_out,'(a,1x,i3,1x,a)', iostat=jerr)
     .                '*',items_per_field(i),
     .                field_label(i)(1:max(1,trimlen(field_label(i))))
              end if
              call report_error('IOSTAT',jerr,'writ',output_file,0,
     .                          'EXTRACT_DATA')
          end if
      end do
 
*                                 ! Get out we have errors
      if( jerr.ne.0 ) RETURN
 
***** Now start searching for fields
 
      call init_search_fields
 
      EOF_in = .false.
 
      line_in  = 0
      line_out = 0
 
*                                  ! Search until end of file
      do while ( .not. EOF_in )
 
*         Read next line
          read(101,'(a)',iostat=ierr) inline
          
          if( ierr.ne.0 ) then
              EOF_in = .true.
*                                 ! Start decoding
          else
 
              line_in = line_in + 1
 
*             For each field see if we can find the start records
*             if these have not been found
              call find_starts( inline )
 
*             Check for end field.  Test moved from below.
              call find_ends( inline )
 
*             Now check to see if we can find the fields we need
              if( trimlen(inline).gt.0 ) call find_fields( inline )
 
*             If we have found the fields set up the output
              if( trimlen(inline).gt.0 ) call make_output( outline )

*             If all is found then output the results
              if( all_found ) then
 
                  line_out = line_out + 1
 
                  write(unit_out,'(a)',iostat=jerr)
     .                outline(1:max(1,trimlen(outline)))
                  call report_error('IOSTAT',jerr,'writ',outline,0,
     .                              'EXTRACT_DATA')
 
*                 Now reset the found information
                  call reset_field_inf
 
*                 If error on output then stop
                  if( jerr.ne.0 ) EOF_in = .true.
 
              end if
 
*             Check to break issued.
 
c              if( ifbrk() ) then
c                  call extract_status( eof_in, line_in, line_out,
c     .                                 outline )
c              end if
 
*             See if we can find end field
C             call find_ends( inline )
 
*                         ! No error reading input
          end if
*                         ! Looping over input file
      end do
 
***** Now rewind the input file
      rewind ( 101 )
 
***** Thats all
      return
      end
 
CTITLE INIT_SEARCH_FIELDS
 
      subroutine init_search_fields
 
 
*     This routine will initialize all of the search field information
*     before we start reading the input file.
 
      include '../includes/extract_comm.h'
 
*   i,j,k       - Loop counters
 
      integer*4 i,j
 
***** Set the start labels found to false
 
      do i = 1, num_fields
          do j = 1, max_starts
              start_found(j,i) = .false.
          end do
      end do
 
*     Set field found false
      do i = 1, num_fields
          field_found(i) = .false.
      end do
 
      all_found = .false.
 
***** Thats all
      return
      end
 
CTITLE FIND_STARTS
 
      subroutine find_starts( inline )
 
 
*     Routine to check line to see if it has start labels
 
      include '../includes/extract_comm.h'
 
*   i,j     - Loop counters
*   pos     - position of start label in line
*   len_start   - Length of the start field
*   trimlen     - HP function for length of string
 
      integer*4 i,j, pos, len_start, trimlen
 
*   inline      - Line read from input file
 
      character*(*) inline
 
****  Loop over all of the fields and all of the starts in the
*     fields
 
      do i = 1, num_fields
          j = 1
          do while ( j.le.num_starts(i) )
 
*                                               ! We have not yet found
              if( .not.start_found(j,i) ) then
 
*                 Search for this start label
                  len_start = trimlen( start_label(j,i) )
                  if ( len_start.gt.0 ) then
                      pos = index ( inline,
     .                              start_label(j,i)(1:len_start) )
                  else
*                                 ! Say we found
                      pos = 1
                  end if
 
*                                         ! Found
                  if( pos.gt.0 ) then
                      start_found(j,i) = .true.
*                                         ! Do not search for the rest of the
                  else
*                                         ! starts (no until we have found
*                                         ! all earlier ones)
                      j = num_starts(i) + 1
                  end if
              end if
 
*                             ! Look for next start
              j = j + 1
          end do
      end do
 
***** Thats all
      return
      end
 
CTITLE FIND_FIELDS
 
      subroutine find_fields( inline )
 
 
*     This routine will search for the field labels, if all the
*     start labels have been found.
 
      include '../includes/extract_comm.h'
 
*   i,j,k       - Loop counters
*   pos         - position of start label in line
*   len_field   - Length of the field label
*   trimlen     - HP function for length of string
 
      integer*4 i,j, pos, len_field, trimlen
 
*   inline      - Line read from input file
 
 
      character*(*) inline
 
****  Loop over each of the fields
 
      do i = 1, num_fields
 
*         First see if we have found all of the start labels
          all_starts_found = .true.
          do j = 1, num_starts(i)
              all_starts_found = all_starts_found .and.
     .                           start_found(j,i)
          end do
 
*         See if we have found all starts
*                                         ! See if we can find field
          if ( all_starts_found ) then
 
*             Search for this field label
              len_field = trimlen( field_label(i) )
              if ( len_field.gt.0 ) then
                  pos = index ( inline,
     .                          field_label(i)(1:len_field) )
              else
*                             ! Say we found
                  pos = 1
              end if
 
*                                     ! Found
              if( pos.gt.0 ) then
                  field_found(i) = .true.
                  pos = pos + len_field
              end if
 
*             Now extract the fields information if we found field
              if( field_found(i) .and. pos.gt.0 ) then
                  call read_field( inline, pos, i)
              end if
*                                     ! all starts found
          end if
*                                     ! Looping over fields
      end do
 
***** Now see if we found all fields
      all_found = .true.
      do i = 1, num_fields
          all_found = all_found .and. field_found(i)
      end do
 
***** Thats all
      return
      end
 
CTITLE READ_FIELD
 
      subroutine read_field( inline, pos, fnum )
 
 
*     Routine to extract the information from the input line
*     given the field information given earlier
 
      include '../includes/extract_comm.h'
 
*   count       - Index for readlining entries
*   data_type   - Type of data to read (1=I*2, 2=R*8, 3=CH)
*   ent     - Position in Ifield when Rfield read with
*           - READ_LINE
*   fnum    - Field number being extracted
*   i,j     - Loop counters
*   idum    - Dummy integer value
*   ierr    - IOSTAT error
*   indx    - Indx for use in format or in read_line
*   next_readl  - Next entry to be extracted by readline
*   pos     - Current position in line after field label
 
      integer*4 count, data_type, ent, fnum, i,j, idum, ierr, indx,
     .    next_readl, pos, cand
 
*   kbit    - Bit checking function
 
      logical kbit
 
*   cdum    - Dummy for READ_LINE
 
      character*1 cdum
 
*   type(3) - Character string with READLINE type
 
      character*2 type(3)
 
*   inline  - Line read from input file
 
      character*(*) inline
 
      data type / 'I4', 'R8', 'CH' / 

      ent = 0
 
***** Get the starting position we should use (either pos, or 1)
 
*                                          ! Start at beginning of
      if( kbit(decode_data(fnum),5) ) then
*                                         ! line
          indx = 1
*                                         ! Start from current position
      else
          indx = pos
      end if
 
***** Now see if we should use readline or format
      data_type = cand(decode_data(fnum),7)
 
*                                           ! Use a format
      if( kbit(decode_data(fnum),4) ) then
 
*         Read depending of data type
*                                         ! I*2 read
          if( data_type.eq.1 ) then
              read(inline(indx:),informat(fnum), iostat=ierr)
     .            (Ifield(j,fnum),j=1,items_per_field(fnum))
              call report_error('IOSTAT',ierr,'read',inline(indx:),
     .            0, 'READ_FIELD')
          end if
 
*                                         ! R*8 read
          if( data_type.eq.2 ) then
              read(inline(indx:),informat(fnum), iostat=ierr)
     .            (Rfield(j,fnum),j=1,items_per_field(fnum))
              call report_error('IOSTAT',ierr,'read',inline(indx:),
     .            0, 'READ_FIELD')
          end if
 
*                                         ! CH read
          if( data_type.eq.3 ) then
              read(inline(indx:),informat(fnum), iostat=ierr)
*                                         ! Only one item for CH
     .            Cfield(fnum)
              call report_error('IOSTAT',ierr,'read',inline(indx:),
     .            0, 'READ_FIELD')
          end if
 
*                         ! Use readline to get inforation
      ELSE
 
          count = 1
          do i = 1, items_per_field(fnum)
 
              next_readl = readl_ents(i,fnum)
              if( data_type.eq.1 ) ent = i
*                                                    ! Allow for R*8
              if( data_type.eq.2 ) ent = 2*(i-1) + 1
 
*                                             ! Skip over non wanted fields
              do j = count, next_readl-1
                  call read_line(inline,indx,'CH',ierr,idum,cdum)
              end do
 
              j = next_readl
              call read_line(inline, indx, type(data_type), ierr,
     .             Ifield(ent,fnum),Cfield(fnum))
              call report_error('IOSTAT',ierr,'readlin',
     .             inline(indx:),0, 'READ_FIELD')
 
              count = next_readl + 1
          end do
      END IF
 
***** Thats all
      return
      end
 
CTITLE FIND_ENDS
 
      subroutine find_ends( inline )
 
 
*     Routine to search for the end labels.  If these are found
*     then the start_found and field_found values are set false.
 
      include '../includes/extract_comm.h'
 
*   i,j         - Loop counters
*   pos         - position of end label in line
*   len_end     - Length of the end field
*   trimlen     - HP function for length of string
 
      integer*4 i,j, pos, len_end, trimlen
 
*   inline      - Line read from input file
 
      character*(*) inline
 
****  Loop over all of the fields and see if we find the end labels
 
      do i = 1, num_fields
 
          len_end = trimlen( end_label(i) )
          if ( len_end.gt.0 ) then
              pos = index ( inline, end_label(i)(1:len_end) )
 
*                                     ! Found
              if( pos.gt.0 ) then
                  do j = 1, num_starts(i)
                      start_found(j,i) = .false.
                  end do
                  field_found(i) = .false.
              end if
          end if
      end do
 
***** Thats all
      return
      end
 
CTITLE MAKE_OUTPUT
 
      subroutine make_output( outline )
 
 
*     This routine actually makes up the output line.
 
      include '../includes/extract_comm.h'
 
*   data_type   - Type of data to be output
*   i,j,k   - Loop counters
*   ierr    - IOSTAT error
*   pos     - Current position in outline
*   trimlen - HP function for length of string
 
      integer*4 data_type, i,j,ierr, pos, trimlen, cand
 
*   outline     - Line to be output to the file
 
      character*(*) outline
 
***** See if we have found all the fields
 
      if( .not.all_found ) RETURN
 
*     Build up the output line
      outline = ' '
 
      do i = 1, num_fields
 
          pos = trimlen(outline) + 2
          if( pos.eq.2 ) pos = 1
 
          data_type = cand(decode_data(i),7)
 
*                                             ! There is field
          if( decode_data(i).ne.0 ) then
*                                             ! information
 
*             Write based on type
*                                             ! I*2
              if( data_type.eq.1 ) then
                  write(outline(pos:),outformat(i),iostat=ierr)
     .                (Ifield(j,i),j=1,items_per_field(i))
              end if
 
*                                             ! R*2
              if( data_type.eq.2 ) then
                  write(outline(pos:),outformat(i),iostat=ierr)
     .                (Rfield(j,i),j=1,items_per_field(i))
              end if
 
*                                             ! CH
              if( data_type.eq.3 ) then
                  write(outline(pos:),outformat(i),iostat=ierr)
     .                Cfield(i)
              end if
 
              call report_error('IOSTAT',ierr,'with writ',
     .                outformat(i),0,'MAKE_OUTPUT')
          end if
      end do
 
****  Thats all
      return
      end
 
CTITLE RESET_FIELD_INF
 
      subroutine RESET_field_inf
 
 
*     Routine to reset the field found values after the output line
*     has been written
 
      include '../includes/extract_comm.h'
 
*   i       - Loop counter
 
      integer*4 i
 
***** Loop over all of the fields and reset if we are supposed to
 
      do i = 1, num_fields
          if( reset_field(i).eq.0 ) then
              field_found(i) = .false.
          end if
      end do
 
***** Thats all
      return
      end
 
CTITLE INIT_EXTRACT
 
      subroutine init_extract
 
 
*     This routine initalzies the field information for EXTRACT
 
      include '../includes/extract_comm.h'
 
*   i,j,k   - Loop counters
 
      integer*4 i,j
 
***** Start off, set all found values to false
 
      num_fields = 0
      title_label = ' '
*                             ! Use first line as title
      title_num   =  1
 
      do i = 1, max_fields
          decode_data(i) = 0
          field_found(i) = .false.
 
          do j = 1, max_starts
              start_found(j,i) = .false.
              start_label(j,i) = ' '
          end do
 
          field_label(i) = ' '
          field_title(i) = ' '
          end_label(i)   = ' '
          informat(i)    = '(*)'
          outformat(i)   = '(10(f13.4,1x))'
 
          items_per_field(i) = 0
          reset_field(i)     = 0
          num_starts(i)      = 0
 
      end do
 
*     Now set up some defaults
      field_label(1) = 'EXPERIMENT date :'
      field_title(1) = 'Experiment date'
      informat(1)    = '(1x,i4,4(1x,i2))'
      outformat(1)   = '(1x,i4,4i3)'
      decode_data(1) = 9
      items_per_field(1) = 5
*                                 ! Will not be reset when output
      reset_field(1) = 1
 
 
****  Thats all
      return
      end
 
CTITLE EXTRACT_STATUS
 
      subroutine extract_status( eof_in, line_in, line_out, outline)
 
 
*     Routine to report the status of the current run, and give the
*     user the option of stopping run or program
 
      include '../includes/extract_comm.h'
 
*   i,j,k       - Loop counters
*   trimlen     - Length of string
 
      integer*4 i,j, trimlen
 
*   line_in     - Line number in input file
*   line_out    - Line number in output file
 
      integer*4 line_in, line_out
 
*   eof_in      - EOF indicator for input file (set true if user
*               - says to stop this command
 
      logical eof_in
 
*   ians        - Answer to question (A -- abort program, S -- stop
*               - this run, C -- continue (default))
 
      character*1 ians
 
*   outline     - Last line written to file
 
 
      character*(*) outline
 
***** Write out the status information.  Write to error unit (unit 0)
*     incase the output is being redirected.
 
      write(0,100) line_in, line_out
 100  format(' EXTRACT STATUS: Line ',i5,' in input, Line ',i5,
     .       ' in output')
 
      do i = 1, num_fields
          if( decode_data(i).ne.0 ) then
              write(0,150) i, field_label(i), field_found(i),
     .                     (start_found(j,i),j=1,num_starts(i))
  150         format(' Field ',i2,' Label ',a20,' Field found ',L5,
     .               ' Starts found ',10L1)
          end if
      end do
 
      write(0,200) outline(1:max(1,trimlen(outline)))
  200 format(' Last line written to output ',/,a)
 
      write(0,250)
  250 format(' Option: A-Abort EXTRACT, S-stop this search,',
     .       ' C-continue ? ',$)
      read(*,*) ians
 
      call casefold( ians )
 
      if( ians.eq.'A' ) STOP ' EXTRACT Terminated : At break'
*                                         ! Will stop this search
      if( ians.eq.'S' ) eof_in = .true.
 
***** Thats all
      return
      end
 
CTITLE EXTRACT_COMBD
 
      block data extract_combd
 
 
*     Sets up the commands for EXTRACT
 
      include '../includes/extract_comm.h'
 
*                                         ! 1--End command
      data extract_commands / 'END     '
*                                         ! 2--Input file name
     .,                       'INPUT   '
*                                         ! 3--Output file name
     .,                       'OUTPUT  '
*                                         ! 4--Description for each field
     .,                       'DESCRIPT'
*                                         ! 5--Title for output
     .,                       'TITLE   '
*                                         ! 6--Field information
     .,                       'FIELD   '
*                                         ! 7--Start label information
     .,                       'BEGIN   '
*                                         ! 8--End label information
     .,                       'FINISH  '
*                                         ! 9--No reset after write
     .,                       'NORESET '
*                                         !10--Output format
     .,                       'OUTFORM '
*                                         !11--Go get the data
     .,                       'RUN     '/
 
      end
 
 
