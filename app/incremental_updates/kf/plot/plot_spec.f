      subroutine decode_data(field, bak_array, values, point)
C     subroutine decode_data(field, inbuffer , values, point)
c
c     Routine to decode the data string read from the data
c     file.  To do this the data record is transferred from ema
c     to main memory where it is decoded using the field information
c
c
c Include files
c -------------
*                           ! the parameter file
      include 'plot_param.h'
c
*                           ! the common block
      include 'plot_com.h'
c
c Variables
c ---------
c field  -- the information on the fields of data
c     field(1) -- data type: 0 = time, 1 = normal
c     field(2) -- start column
c     field(3) -- For time fields = number of columns in field.
c        These are interpretted as year,month,day, hour, min, and
c        seconds.  For normal fields this value is the column for
c        the standard deviation.
c bak_array -- the ema file of data
c inbuffer  -- the character buffer containing the data record
c
c values   -- the returns values from the strings.  This will be
c     julian date for time fields, for normal fields these will be
c     the values.
c point -- the type of point to be plotted and its edit flag
 
      integer*4 field(3), point(2)
 
c
      real*8 values(2)
 
c
      integer*4 bak_array(1)
 
c
C     character*(*) inbuffer
c
c
c
c Local variables
c ---------------
c unw_flag  -- the unweight flag from the bak_array
c
c idummy   -- an array of integer values used to get to the
c     the correct column number in the file
c ivalues  -- integer values to be read for data arguments
c
c local_data -- the copy of the ema file record
c char_data  -- local data in charter form (equivalence to local_data)
c fjldy_8    -- real*8 function to return julian data from month, day,
c     year.
c ierr       -- error flag
c i          -- loop counter
c indx       -- Pointer in line for READLINE
 
      integer*4 unw_flag,  ivalues(6), local_data(int_recl),
     .    ierr, i, indx

* MOD TAh 140328: Added negative field 3 to allow yr doy hr min sec format
      integer*4  absf3   ! Abs value of field three to allow yr doy date


      real*8 decyrs  ! Used for JD to DecYrs
 
c
      character*(char_recl) char_data
 
c
*   icdum(30)       - Dummy space for characters
 
      character*1 icdum(30)
 
      equivalence (local_data,char_data)
c
      real*8 sectag 
 
c
c.... Intialize the values
      values(1) = 0.d0
      values(2) = 0.d0
*                       ! initialize to 'A'
      point(1)     = 1
      point(2)     = 0
c
c.... copy ema data to a local string
      do i = 1, int_recl
          local_data(i) = bak_array(i)
      end do
c
c.... Start decoding string
*                                   ! read time arguments
      if( field(1).eq.0 ) then
 
          indx = 1
          call pos_indx( char_data, field(2), indx )
          if( field(3).gt.0 ) then 
* MOD TAH 140328: Date field Y M D H M s (original form) 
             if( field(3).lt.6 ) then 
                 call multiread( char_data, indx, 'I4', ierr, ivalues,
     .                   icdum, field(3) )
                 sectag = 0.d0
             else
                 call multiread( char_data, indx, 'I4', ierr, ivalues,
     .                   icdum, 5 )
                 call read_line( char_data, indx, 'R8', ierr, sectag,
     .                   icdum )
             end if
          else
* MOD TAH 140328: New date form at Y DoY H M S
             absf3 = abs(field(3))
             if( absf3.lt.5 ) then 
                 call multiread( char_data, indx, 'I4', ierr, ivalues,
     .                   icdum, absf3 )
                 sectag = 0.d0
             else
                 call multiread( char_data, indx, 'I4', ierr, ivalues,
     .                   icdum, 4 )
                 call read_line( char_data, indx, 'R8', ierr, sectag,
     .                   icdum )
             end if
*            Now insert Jan into month
             ivalues(3:5) = ivalues(2:4)
             ivalues(2) = 1    ! Jan + DOY will give correct MJD 
          endif
 
          call ymdhms_to_mjd(ivalues,sectag,values(1))
          values(2) = 0.d0
 
*                                  ! time data
      end if
c
c...  See if normal data
*                                  ! Normal data -- read directly
! MOD TAH 131015: Added type 2: JD/MJD converted to years
      if( field(1).eq.1 .or. field(1).eq.2 ) then
 
          indx = 1
 
          call pos_indx( char_data, field(2), indx )
 
          call read_line( char_data, indx, 'R8', ierr, values(1),
     .             icdum )
          if( field(1).eq.2 ) then ! convert to Dec Years
               call JD_to_Decyrs( values(1), decyrs )
               values(1) = decyrs
          end if
 
*                                                    ! if wanted get sigma
          if( field(3).gt.0 ) then
 
              indx = 1
              call pos_indx( char_data, field(3), indx)
              call read_line( char_data, indx, 'R8', ierr, values(2),
     .                icdum )
! MOD TAH 131015: see if MJD->Years conversion
              if( field(1).eq.2 ) values(2)=values(2)/365.25d0

          end if
      end if

* MOD TAH 200331: Added field 3 which converts JD/MJD to calendar
*         type entry (i.e, field 0)
      if( field(1).eq.3 ) then
         indx = 1
 
          call pos_indx( char_data, field(2), indx )
*         Values(1) is the MJD/JD as it is for field 0. 
          call read_line( char_data, indx, 'R8', ierr, values(1),
     .             icdum )

      end if
 
c
c.... See if we should extract a point type
*                                 !   read a point type
      if( p_field(1).ne.0 ) then
 
          indx = 1
          call pos_indx( char_data, p_field(1), indx )
          call read_line( char_data, indx, 'I4', ierr, point, icdum)
 
      end if
c
c.... See if any unweight flag (point set lower case if unw_flag>0)
*                                  ! get unweight flag
      if( p_field(2).ne.0 ) then
          indx = 1
          call pos_indx( char_data, p_field(2), indx)
          call read_line( char_data, indx, 'I4', ierr, unw_flag, icdum)
          point(2) = unw_flag
 
c
c....     See if we need to change plot point
*                                    ! convert point to lower case
          if( unw_flag.ne.0 ) then
              point(1) = point(1) + 32
          end if
      end if
c
c.... See if point valid ASCII
      if( point(1).gt.64 .or. point(1).lt.0  ) then
*                               ! use '?' for non_valid point type
          point(1) = ichar('?') - 32
      end if
c
c.... Process any errors
 1000 continue
*                             ! error ocurred
      if( ierr.ne.0 ) then
          call report_error('IOSTAT', ierr, 'decod', char_data,
     .        0, 'DECODE_DATA')
      end if
c
c.... Thata all.
      return
      end
 
 
c.......................................................................
 
      subroutine read_file(ema_data, header_only)
c
c     Routine to read data from file into ema.  It will also set
c     up the dynamic mapping of ema for this data set.
c
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                         ! the common block
      include 'plot_com.h'
c
c Variables
c ---------
c ema_data -- an integer array which will used to store the contents
c     of the file and the data to be plotted.  The data to be plotted
c     will be extracted by read_data.
* header_only - Read only header. Not implemented in this routine
c
      integer*4 ema_data(max_plt_space)

      logical header_only
 
c
c
c Local variables
c ---------------
c ierr -- general error variable
c
      integer*4 ierr, trimlen
 
c
c
c Variables used in reading file
c ------------------------------
c idata_array -- an integer array equivalenced to cdata_array
c cdata_array -- a character array equivalenced to idata_array
c eof  -- indicate EOF reached
c i    -- loop counter
c ifbrk -- logical function to determine if program has a break (not used)
c
 
      integer*4 idata_array(int_recl), i
 
      character*(char_recl) cdata_array
 
c
 
      logical eof
 
      equivalence (idata_array,cdata_array)
c
c.... Set file read false
      file_read = .false.
      eof       = .false.
c
c.... open the input file (if one has been specified)
*                                        ! no file name
      if( input_file(1:1).eq.' ' ) then
         write(termlu,'(a)')
     .             " READ_FILE Error: No data file has been given"
*                              ! open and read file
      else
*                                    ! convert name to upper case
         close(200)
         open(200,file=input_file, iostat=ierr, status='old')
 
         call report_error('IOSTAT',ierr,'open',input_file,
     .        0,'READ_FILE')
c
*                                ! no use continuing
         if( ierr.ne.0 ) return
c
c....     Now read the data, Firstly save the header records
          do i = 1, actual_num_headers
              read(200,'(a)', iostat=ierr, err=1000 )
     .            headers(i)
          end do
c
c....     Finish read headers which can be be saved
          do i = actual_num_headers+1, input_num_headers
              read(200,'(a)')
          end do
c
c....     Now start reading data
          num_epochs = 0
          do while ( ierr.eq.0 )
 
              read(200,'(a)', iostat=ierr, err=100)
     .            cdata_array
 
c
c....         See if break issued
C             if( ifbrk() ) then
C                 ierr = -1
C                 call report(' User break: ending read of data file')
C             end if
c
 100          continue
*                                    ! EOF encounter
              if( ierr.eq.-1 ) then
                  eof = .true.
              end if
c
c....         If no error then process the record
              if( ierr.eq.0 .and. 
     .            (cdata_array(1:1).eq.' '.or. ignore_col1) .and.
     .             trimlen(cdata_array).gt.0  ) then
                  num_epochs = num_epochs + 1
c
c....             Save a typical line of the file
                  if( num_epochs.eq.1 ) then
                      typical_record = cdata_array
                  end if
c
c....             Now see if we have too much data
                  if( (num_epochs+1)*(int_recl+9).ge.
*                                          ! too much data, report
     .                max_plt_space ) then
c                                            and set error flag
                      call report_error('Too many data',num_epochs,
     .                    'read',input_file,0,'READ_FILE')
*                                   ! copy this data, but then terminate
                      ierr = -1
*                                   ! read (treat as end of file)
                      eof  = .true.
                  end if
c
c....             Copy to ema
                  do i = 1, int_recl
*                                           ! use ema data because ema not
*                                           ! mapped yet.
                    ema_data((num_epochs-1)*int_recl + i ) =
     .                     idata_array(i)
                  end do
              end if
          end do
c
c....     Finish up the file reading
          close(200)
          file_read = .true.
c
c....     Set up dynamic mapping of ema arrays
          call plt_ema_map
c
c....     Summarise the file
          call sum_file
c
c....     See if any errors
 1000     Continue
*                                               ! error reading file
          if( ierr.ne.0 .and. .not.eof ) then
              call report_error('IOSTAT', ierr,'read',cdata_array,
     .            0,'READ_FILE')
          end if
c
*                  ! file name given
      end if
c
      return
      end
 
 
c......................................................................
 
      subroutine plt_ema_map
c
c     This routine will step up the plotting package mapping of the
c     ema area.  The first (int_recl*num_epochs) is used for the bak_file,
c     the remaining area is used for the x, y, and pt data
c
c Include files
c -------------
*                           ! the plot parameter file
      include 'plot_param.h'
c
*                           ! the plot control common block
      include 'plot_com.h'
c
c All the information we need for this routine is stored in plot_com
c
c
c.... Start setting up the address space
      ibak_array = 1
c
      if( .not.xdiff .and. .not.ydiff ) then
*                                                        ! add bak_array space
      ix_array   = ibak_array + int_recl*num_epochs
c
*                                                    ! 2 real*4 elements for
      iy_array   = ix_array   + 4*num_epochs
c                                                       each epoch
c
*                                                    ! same amount as y_array
      ipt_array   = iy_array   + 4*num_epochs
      end if
c
c.... Make sure we have not gone outside ema space
*                                                         ! not enough room
      if( ipt_array+2*num_epochs.gt.max_plt_space ) then
         write(*,100) num_epochs
 100     format(' There are ',i8,' epochs of data')
         write(*,110) ipt_array+2*num_epochs, max_plt_space
 110     format(/" Out of ema space, space needed is ",i12," words",/,
     .      " Only ",i12," words of ema space available.",
     .      " Increase max_plt_space in &pltpa (the parameter file)")
c
         stop ' Out of ema space in pltsl'
      end if
c
c.... Thats all
      return
      end
 
c......................................................................
 
      subroutine sum_file
c
c     This routine summarizes the contents of the bak_file.  It lists
c     the sites sources and the markov elements which are available
c
c Include files
c -------------
*                          ! the parameter file
      include 'plot_param.h'
c
*                          ! the common file
      include 'plot_com.h'
c
c All the information we need is in the common file
c
c Functions
c ---------
c trimlen -- HP string length utility
c
      integer*4 trimlen
 
c
c
c.... See if file name given
*                                             ! no file yet
      if( trimlen(input_file).eq.0 ) return
c
c.... Start by giving the summary of the solvk solution
      write(termlu,100) input_file(1:trimlen(input_file))
 100  format(/" Summary of data in ",a)
c
      write(termlu,150) num_epochs
 150  format(" There are ",i8," data records")
c
c.... Output the header records
      call out_names(termlu,'header records', headers,
     .    actual_num_headers)
c
      write(termlu,'(" The first data record is:")')
      call report(typical_record)
c
c.... Thats all
      return
      end
 
c........................................................................
 
      subroutine get_field(field, axis_label)
c
c     Subroutine to decode the field information for the plot
c
c Include files
c -------------
*                         ! the parameter file
      include 'plot_param.h'
c
*                         ! the common block
      include 'plot_com.h'
c
c Variables
c ---------
c field  -- the field array:
c     field(1) = type of data (0 is time, 1 is normal)
c     field(2) = the column number for then quanitity to be plotted
c     field(3) = number of colmumns in time field
c axis_label -- the default label to be put on the axis
c
c From the common:
c p_field    -- is the point field array
c     p_field(1) = the field for the point type
c     p_field(2) = the edit flag field
c
      character*(*) axis_label
 
c
      integer*4 field(1)
 
c
c new_label -- dummy label used to construct the axis label
c
      character*(char_recl) new_label
 
c
c Functions
c ---------
c trimlen -- HP utility to return length of string
c
      integer*4 trimlen
 
c
c.... Clear the copy of the axis label
      new_label  = ' '
c
c.... Get the field information
      call get_int(buffer,field,3)
c
c.... Check that valid field type has been entered
*                                                   ! not valid type
! MOD TAH 131015: Allow MJD->DecYear type (2)
! MOD TAH 200331: Allow MJD->Y M D H M sec type (3)
      if( field(1).lt.0 .or. field(1).gt.3  ) then
          call report_error('FIELD',3002,'inputt',buffer,0,'GET_FIELD')
*                          ! will be caught as an error if we try to
          field(1) = -1
*                          ! use this field
      end if
c
c.... Get axis label
      call read_label(buffer, new_label )
c
c.... Add the last '"' to the axis label
      if( trimlen(new_label).gt.0 ) then
          axis_label = '"' // new_label(1:trimlen(new_label)) // '"'
      else
          axis_label = '""'
      end if
c
      return
      end
 
c......................................................................
 
      subroutine help
c
c
c     This routine will list available commands, and give help information
c     from the help_file on specific commands
c
c Include files
c -------------
*                         ! the pltsl parameter file
      include 'plot_param.h'
c
*                         ! the pltsl common block
      include 'plot_com.h'
c
c Variables
c ---------
c iel -- the command number from get_cmd.
c ierr -- the error return from reading files
c eof  -- end of file indicator
c found -- set true when command found
c finished -- set true when when help has finished output from file
c line -- the line read from the help file
c i    -- loop counter
c len  -- length of 'line' or 1 if line is null
c
      integer*4 iel, ierr, i, len
 
c
      logical eof, found, finished
 
c
      character*79 line

* FullHelpName - Full name of help file

      character*128 FullHelpName
c
 
c
c Functions
c ---------
c trimlen -- HP string length utility
c
      integer*4 trimlen
 
c
c.... See if help is requested for a specific command
      call get_cmd(buffer(9:), commands, num_commands, iel)
c
c.... If iel is greater then 0 then give specific help information
      if( iel.gt.0 ) then
c
          call gen_help_name(help_file, fullhelpname)
          open(300,file=fullhelpname,iostat=ierr, status='old')
          if( ierr.ne.0 ) then
              call report_error('IOSTAT',ierr,'open',fullhelpname,
     .            0,'HELP')
c
c....         Set iel = 0 so that standard help massage will be output
              iel = 0
c
*                 ! file openened OK, find help message
          else
c
              eof = .false.
              finished = .false.
              found = .false.
c
c....         Read file until we reach EOF or we have finished outputting
c             the help meassge
              do while ( .not.eof .and. .not.finished )
c
                  read(300,'(a)',iostat=ierr, err=100) line
  100             continue
*                                       ! EOF
                  if( ierr.eq.-1 ) then
                      eof = .true.
*                                       ! see if other error
                  else
                      if( ierr.ne.0 ) then
                          call report_error('IOSTAT',ierr,'read',
     .                        help_file,0,'HELP')
*                                   ! force standard massage
                          iel = 0
                      end if
                  end if
c
c....             See if help command found
*                                       ! check
                  if( ierr.eq.0 ) then
*                                           ! see if command line
                      if( .not.found ) then
                          if( line(3:10).eq.commands(iel) ) then
*                                              ! we found the command
                              found = .true.
                          end if
*                                          ! see if next command found ( end
                      else
*                                          ! of help for current command)
*                                                       ! end
                          if( line(1:2).eq.'**' ) then
                              finished = .true.
*                                          ! write help line
                          else
                              len = max(1,trimlen(line))
                              write(termlu,'(a)') line(1:len)
                         end if
                      end if
*                                  ! no file reading error
                   end if
*                                  ! looping over file
              end do
              close (300)
*                                  ! file opened OK
          end if
*                                  ! iel > 0 (may have been set no to zero
      end if
*                                  ! if file error
c
c.... List the availabke commands to the termlu if iel<=0
      if( iel.le.0 ) then
          write(termlu,300)
 300      format(/" HELP for CPLOTX:",/,
     .            " The runstring for CPLOTX is:",//,
     .            " % CPLOTX <control file or lu> <display>",
     .            " <plot file> <# header records> ",
     .            " <ignore 1> <list of upto 9 more files to plot>",//,
     .            " All runstring parameters are optional")
 
          write(termlu,310)
 310      format(/" All commands must be preceeded by at least one",
     .            " blank character.",/,
     .            " Commands may be abbreviated")
 
          write(termlu,'(" The commands available in plot are:")')
          write(termlu,320) (commands(i),i=1,num_commands)
 320      format(5(1x,a,3x))
 
          write(termlu,330)
 330      format(" Further help may be obtained using",/,
     .           "?  HELP <command> or HELP CPLOTX"/)
c
         call sum_file
      end if
c
      return
      end
 
c......................................................................
 
CTITLE PLOT_FINAL_SETUP
      subroutine plot_final_setup
 
 
*     Routine to do the final setup for the plot programs
 
      include 'plot_param.h'
      include 'plot_com.h'
 
*     Save the record length for the data file
      bak_recl = int_recl
 
c.... Tell the user we are here
      call report('PLOT: (preceed all commands by at least one blank)')
      call report('       (use HELP for instructions, END to quit)')
      if( .not. ignore_col1 ) then
          write(*,*) ' PLOT: (Non-blank in col 1 of data files treated',
     .               ' as comment)'
      end if
 
      return
      end
 
CTITLE POS_INDX
 
      subroutine pos_indx( char_data, field, indx )
 
 
*     Routine to move indx upto the character imediately before the
*     the position we want to use with READLINE
 
*   field   - Field position desired
*   id      - Dummy place holder
*   ierr    - IOSTAT error
*   indx    - Pointer to character position in line
 
 
      integer*4 field, id, ierr, indx
 
*   icdum(30)   - Dummy character
 
      character*1 icdum(30)
 
 
 
      character*(*) char_data
 
***** See if we want to move
 
      if( field.gt.1 ) then
          call multiread( char_data, indx, 'CH', ierr, id, icdum,
     .                     field-1)
      end if
 
***** Thats all
      return
      end
