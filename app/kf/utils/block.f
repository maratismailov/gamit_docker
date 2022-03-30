      program block
 
**************************************************************************
*     Program to read a series of soucres files and to output a file
*     which shows the relation between the subroutine calls
*
*     The runstring of the program is:
*     CI> BLOCK [Input file list] [Output plot file] <max depth> <max length>
*
*     where
*     Input file list is the name of a file containing the list of input
*         files (REQUIRED)
*     Output plot file is the name of a file which can be input to the
*         plotting program PLOT to draw the relational diagram between the
*         subroutines. (REQUIRED)
*     max depth is the maximum depth to which subroutine calls are made.
*         (OPTIONAL, Default is the actual depth of subroutine calls.)
*          RESTRICTION: the max depth value should NOT be be set less than
*         the actual depth of subroutine calls.
*     max length is the number of lines to output on the graphics device.
*         (OPTIONAL: Default is 64).
*
*     T.Herring                    2:50 PM  TUE.,  6  MAY , 1986
*                                 last modified <870421.1330>
*
***************************************************************************
 
 
      include 'block.h'
 
*   depth       - the depth to which the subroutine calls go
*   depth_out   - the depth to be used on the output file
*   entry       - entry in current module
*   i,j,k       - loop counters
*   iel         - a index to routine names in the modules list
*   ierr        - IOSTAT file error
*   decimaltoint    - Character string to integer conversion
*               - utility from HP.
 
*   length      - the length of subroutine call chart
*   length_out  - the length of the output chart (defaults to 50)
*   len_runstring   - Length of the runstring
 
*   lu_user     - users lu number
 
*   max_depth   - the maximum depth of the soubroutine calls
 
*   mod_num     - module number being processed
*   rcpar       - HP runstring reading utility
 
*   trimlen     - HP string length utility
 
      integer*4 depth, depth_out, entry, i,j, iel, ierr,
     .    decimaltoint, length, length_out, len_runstring, lu_user,
     .    max_depth, mod_num, rcpar, trimlen
 
*   eof         - indicates eof
*   fortran     - Indicates if current code is fortran 77.  If not
*               - then code is treated as comments.  We intially
*               - assume that it is fortran.
 
      logical eof, fortran
 
*   name        - name returned with subroutine name in it
 
      character*25 name
 
*   input_file  - file containing list of soucre files
*   output_file - output file cotaining linkage data
*   runstring    - the runstring parameters
 
      character*128 input_file, output_file, runstring
 
*    buffer         - line read from source file
 
 
      character*72 buffer
 
***** START, get the runstring parameters
 
      lu_user = 6
      length_out = 64
      depth_out  = 0
      len_runstring = rcpar(1, input_file )
      if( len_runstring.le.0 ) then
         call proper_runstring('block.hlp','BLOCK',1)
      end if
      len_runstring = rcpar(2, output_file)
      if( len_runstring.le.0 ) then
         call proper_runstring('block.hlp','BLOCK',1)
      end if
 
***** Get the depth to be used on the plots
      len_runstring = rcpar(3, runstring )
 
      if( len_runstring.ne.0 ) then
          depth_out = decimaltoint( runstring, ierr )
          call report_error('decimaltoint', ierr, 'decod',runstring,
     .                       1,'BLOCK' )
      end if
 
***** Get the length of the plots
      len_runstring = rcpar(4, runstring )
 
      if( len_runstring.ne.0 ) then
          length_out = decimaltoint( runstring, ierr )
          call report_error('decimaltoint', ierr, 'decod',runstring,
     .                       1,'BLOCK' )
      end if
 
***** open input and output files
      open(100,file=input_file, iostat=ierr, status='old')
 
      call report_error('IOSTAT',ierr,'open',input_file,1,'BLOCK')
 
      open(200, iostat=ierr, status='SCRATCH' )
 
      call report_error('IOSTAT',ierr,'open','SCRATCH',1,'BLOCK')
 
****  Loop over input file and get the names of the source files
      i = 0
      do  while ( ierr.eq.0 )
          i = i + 1
          if( i.gt.max_sources ) then
              stop ' Too many source files'
          end if
*
          read(100,'(a)',iostat=ierr, end = 100 ) source_names(i)
 
          call trim_lead( source_names(i), ierr )
 
 100      continue
      end do
 
      num_sources = i - 1
 
****  Now we start major loop reading all of the source files and
*     compiling the module and called subroutine names
                                                                 
      write(lu_user,'(/,a,i4,a,a,/,50(2a39,/) )')
     .       'There are ',num_sources,' source files in ',input_file
     .      , (source_names(i),i=1,num_sources)
c     above replaces below to avoid breaking Hollerith rwk 970920
c      write(lu_user,'(/" There are ",i4," source files in ",a,/,
c     .                50(2a39,/) )')
c     .    num_sources, input_file, (source_names(i),i=1,num_sources)
 
 
      num_modules = 0
      num_subr = 0
      routine_numbers(1,1) = 1
 
      do i = 1 ,num_sources
 
*****     open the source file
          open(100,file=source_names(i), iostat=ierr, status='old')
 
          call report_error('IOSTAT',ierr,'open',source_names(i),
     .        0,'BLOCK')
 
          eof = .false.
          if( ierr.ne.0 ) then
              eof = .true.
              write(lu_user,'(" Skipping this file")')
          end if
 
*****     Assume intially that code is fortran
 
          fortran = .true.
 
*****     Now read file
          do while ( .not.eof )
 
              read(100,'(a)', iostat=ierr, err=200 )
     .            buffer
  200         continue
              if( ierr.ne. -1 .and. ierr.ne.0 ) then
                  call report_error('IOSTAT',ierr,'read',
     .                source_names(i), 0,'BLOCK')
              end if
              if( ierr.ne.0 ) then
                  eof = .true.
*                                       ! process line
              else
 
                  call casefold(buffer)
 
*****             Check if fortran.  Routine will comment code if
*                 not fortran
 
                  call check_fortran( buffer, fortran)
 
*****             Process if not comment or zero length
                  if( buffer(1:1).ne.'C' .and. buffer(1:1).ne.'*'
     .            .and. trimlen(buffer).gt.1 ) then

*                     Now clear out the first six characters of the
*                     buffer
                      buffer(1:6) = ' '
 
*****                 See if program
                      call get_name(' PROGRAM',buffer,name)
*                                                    ! PROGRAM
                      if( trimlen(name).gt.0 ) then
                          if( num_modules.ne.0 ) then 
                              write(lu_user,'(a)')
     .                       ' **** WARNING **** >1 Program module '

*****                         Convert program to subroutine (for segmented
*                             programs
 
                              call replace_program( lu_user, buffer )
                          else
                              num_modules = num_modules + 1
                              modules(num_modules) = name
                          end if
                      end if
 
*****                 Replace SEGLD call with direct call to module
 
                      call replace_segld( lu_user, buffer )
 
*****                 Replace LINCA call with direct call to module
 
                      call replace_linca( lu_user, buffer )
 
*****                 See if SUBROUTINE or FUNCTION
 
                      call get_name('SUBROUTINE',buffer,name)
*                                                    ! see if function
                      if( trimlen(name).eq.0 ) then
                          call get_name('FUNCTION',buffer,name)
                      end if
 
*****                 If name present save the results
                      if( trimlen(name).gt.0 ) then
 
*****                     Save number of subroutines
                          routine_numbers(2,num_modules) =
     .                        num_subr - routine_numbers(1,num_modules)
     .                                 + 1
 
*****                     Go to next module number
                          num_modules = num_modules + 1
                          if( num_modules.gt.max_modules ) then
                              stop ' Too many subroutine modules'
                          end if
 
                          modules(num_modules) = name
                          routine_numbers(1,num_modules) =
     .                        num_subr + 1
 
                      end if
 
*****                 See if CALL
                      call get_name(' CALL',buffer,name)
*                                                    !
                      if( trimlen(name).gt.0 ) then
 
*****                     See if we already have a call to
*                         subrouitine
                          do j = routine_numbers(1,num_modules),
     .                           num_subr
                              if( name.eq.subr_names(j) ) then
                                  name = ' '
                              end if
                          end do
                      end if
 
*****                 If still have name save
                      if( trimlen(name).gt.0 ) then
 
                          num_subr = num_subr + 1
                          if( num_subr.gt.max_routines ) then
                              stop ' Too many subroutine calls'
                          end if
 
*****                     Save name
                          subr_names(num_subr) = name
                      end if
 
*                             ! NOT comment
                  end if
*                             ! No error while reading
              end if
 
*                             ! reading over this source file
          end do
 
*****     Save the final routine_numbers entry
          routine_numbers(2,num_modules) = num_subr -
     .        routine_numbers(1,num_modules) + 1
 
*****     OUTPUT a little measage           
          write(lu_user,'(a,i3,a,i4,a,i4,a)') 
     .        'At the end of file ',i,'.  There are ',num_modules
     .       ,' modules and ',num_subr,' calls' 
c         above replaces below -- rwk 970920
c          write(lu_user,'(" At the end of file ",i3,". There are ",
c     .        i4," modules and ",i4," calles")')
c     .        i,num_modules, num_subr
 
*                             ! looping over the files
      end do
 
****  DEBUG write out what we have
                          
      write(lu_user,'(/,a,/,a,/,50(3a25,/))')
     .   'Modules found in source files'
     .  ,'-----------------------------'
     .  ,(modules(i), i=1,num_modules)   
c     above replaces below
c      write(lu_user,'(/"Modules found in source files",/,
c     .                 "-----------------------------",/,
c     .                 50(3a25,/))')
c     .                (modules(i), i=1,num_modules)  
                          
      write(lu_user,'(/,a,/,a,/,50(3a25,/))')
     .   'Calls found in source files'
     .  ,'---------------------------'
     .  ,(subr_names(i), i=1,num_subr)   
c     above replaces below
c      write(lu_user,'(/"Calls found in sources file",/,
c     .                 "---------------------------",/,
c     .                 50(3a25,/))')
c     .                (subr_names(i),i=1,num_subr)
       write(lu_user,'(/,a,/,a,/,50(3a25,/))') 
     .    "Calls found in sources file","---------------------------"
     .   , (subr_names(i),i=1,num_subr)
 
C     write(lu_user,'(" LINKS ",/,50( 5(2i4,7x),/) )')
C    .   ((routine_numbers(i,j),i=1,2), j=1,num_modules)
 
 
*************************************************************************
*
*     NOW BUILD CROSS REFERENCE
*
*************************************************************************
 
      do i = 1, num_subr
          call find_name( subr_names(i), iel )
 
*                                ! subroutine not found in modules
          if( iel.eq.0 ) then
*                                ! list, so add to list
 
              num_modules = num_modules + 1
              if( num_modules.gt.max_modules ) then
                  stop ' Too many modules while building cross ref'
              end if
 
              modules(num_modules) = subr_names(i)
*                                                 ! no calls
              routine_numbers(2,num_modules) = 0
              iel = num_modules
          end if
 
****      Save cross reference
          cross_ref(i) = iel
 
      end do
 
C     write(lu_user,'(" Cross reference ",/,100( 5(2i4,7x),/) )')
C    .    (i,cross_ref(i),i=1,num_subr)
 
*************************************************************************
*
*     NOW WRITE OUT AND DETERMINE THE LINKAGE CHART
*
*************************************************************************
 
*                                               ! Main program
      write(200,'("   1   1  ",a)') modules(1)
      length = 1
      depth  = 2
      max_ depth = 1
 
*                     ! start at first entry of main program
      mod_num = 1
      entry   = 1
 
      stack_length = 0
*                                  ! push main program onto stack
      call push( mod_num,entry )
 
*                              ! process until we get back to main
      do while ( depth.gt.1 )
 
c....     Get module pointed to by current entry in current module
          mod_num = cross_ref(routine_numbers(1,mod_num) + entry -1)
 
c....     write out this module name
          write(200,300) depth, length, modules(mod_num)
  300     format(2i4,2x,a)
 
          max_depth = max(depth, max_depth)
 
c.....    Get where we are pointing to
          depth = depth + 1
          entry = 1
 
c.....    See if we point anywhere
          if( entry.le.routine_numbers(2,mod_num) ) then
 
              call push( mod_num,entry )
 
*                     ! we have reached the end of the call sequence,
          else
*                     ! go back to previous routines
              length = length + 1
              do while ( entry.gt. routine_numbers(2,mod_num) .and.
*                                         ! get out if we are finished
     .                   depth.gt.1 )
                  depth = depth - 1
*                                         ! the stack should be empty
                  if( depth.gt.1 ) then
*                                         ! if depth gets back to 1
                      call pop( mod_num, entry )
                      entry = entry + 1
                  end if
 
              end do
 
*                                          ! push for later use
              call push( mod_num, entry )
 
          end if
      end do
 
***** Increment length to allow for one blank line at the top of the plot
      length = length + 1
 
***** Finsihed first part of the process.  Now set up the .plt file
 
      write(lu_user,'(/" Maximun depth and length are ",2i4)')
     .    max_depth, length
 
      if( depth_out.eq.0 ) then
          depth_out = max_depth
      end if
 
      if( length_out.eq.0 ) then
          length_out = length
      end if
 
      call block_program( output_file, depth_out, length_out )
 
***** Thats all (Close all of the units)
      close(100)
      close(200)
      close(300)
 
      end
 
c--------------------------------------------------------------------------
 
      subroutine get_name(search, buffer, name )
c
c**********************************************************************
c     Routine to get the 'name' following the string 'search' from
c     buffer or from the next records in the input file
c
c     T.Herring                    9:50 AM  WED.,  7  MAY , 1986
c                                 last modified <870421.1330>
c
c**********************************************************************
c
 
*   end     - end of character position of name following 'search'
*   i       - Loop counter
*   ib      - counter used to find first non blank character
*   iel     - start of 'search' in buffer
*   ierr    - error flags
*   len_search  - length of the search string
*   num_open    - NUmber of open ('s and )'s before SEARCH
 
*   pos_comma   - Position of first comma (used to handle cases such
*           - as PROGRAM name,5,99
*   trimlen - HP string length utility
 
      integer*4 end, i, ib, iel, ierr, len_search, num_open, pos_comma,
     .    trimlen
 
*   found   - indicates when  'name' is found
 
      logical found
 
*   buffer  - buffer read from file
*   name    - name following 'search' (null string if search not
*           - found)
*   search  - the string to be searched for
 
      character*(*) buffer, name, search
 
*   new_buffer  - copy of buffer used for manipluations
 
      character*72 new_buffer
 
***** START, see if we can find search
*                 ! clear name
      name = ' '
 
***** Repace first comma
      pos_comma = index( buffer,',' )
      if( pos_comma.ge.7 ) then
          buffer(pos_comma:pos_comma) = '!'
      end if
 
      iel = index(buffer,search)
 
***** see where any '!' is.  If it before the 'searrch' then this is a comment
      end = index(buffer,'!')
*                                          ! 'search' was a comment
      if( end.lt.iel .and. end.gt.0 ) then
          iel = 0
      end if
 
***** Now search for pairs of quoutes and make sure that the 'SEARCH'
*     label does not appear inside quotes (both double and single)
 
      call find_quotes( buffer, iel, '''' )
      call find_quotes( buffer, iel, '"' )
 
***** Find the start of the fortran line and make sure this matches
*     iel
      if( trimlen(buffer).gt.7 ) then
          ib = 7
          do while ( buffer(ib:ib).eq.' ' )
              ib = ib + 1
          end do
 
*....     ib points to first charcter after column 7, match sure
*         this value matches iel
C COMMENTED OUT THIS CHECK TO HANDLE CASE STATEMENT.  REPLACED WITH
C CODE BELOW.
C         if( iel.ne.ib ) THEN   ! We will still except this line
C                                ! if the first part of string is IF
C             if( buffer(ib:ib+1).ne.'IF' ) then
C                 iel = 0        ! Not a valid string
C             end if
C         end if
 
*         Now if the first character does not point to search then
*         check for closing ().  Should handle IF and CASE statememts
 
          if( iel.ne.ib .and. iel.ge.7 ) then
              num_open = 0
              do i = 1,iel-1
                  if( buffer(i:i).eq.'(' ) num_open = num_open + 1
                  if( buffer(i:i).eq.')' ) num_open = num_open - 1
              end do
 
*             Reset iel if () do not close
              if( num_open.ne.0 ) iel = 0
          end if
 
*         Now if there is
 
 
      end if
 
****  If valid 'search' then process
*                         ! we found the search string, get the name
      if( iel.gt.0 ) then
 
          len_search = trimlen(search)
 
*****     Look for a valid terminator ( '(', '!' or blanks)
 
          end = index(buffer(iel+1:),'(' ) + iel -1
*                                 ! look for '!'
          if( end.eq.iel-1 ) then
              end = index( buffer(iel+1:),'!' ) + iel - 1
          end if
 
*                                 ! get length of string
          if( end.eq.iel-1 ) then
              end = trimlen(buffer)
          end if
 
****      Now move string between end of 'search' and end
          new_buffer = ' '
          if( end.gt. iel+len_search+1) then
              new_buffer = buffer(iel+len_search+1:end)
          end if
 
****      See if new buffer contains name
*                                               ! name must be on next line
          if( trimlen(new_buffer).eq.0 )  then
 
              found = .false.
 
              do while ( .not.found )
 
*****             Get next line of file
                  read(100,'(a)', iostat = ierr, err=1000)
     .                buffer
                  call casefold(buffer)
 
*****             Skip over comments and null strings
                  if( buffer(1:1).ne.'C' .and. buffer(1:1).ne.'*'
     .                .and. trimlen(buffer).gt.7 ) then
 
*****                 Get rid of leading blanks
                      iel = 7
                      do while ( buffer(iel:iel).eq.' ' )
                          iel = iel + 1
                      end do
 
*****                 Found start of name, now find end of name
                      end = index(buffer,'(') - 1
                      if( end.eq.-1 ) end = index(buffer,'!') - 1
                      if( end.eq.-1 ) end = trimlen(buffer)
 
*****                 Copy name
                      new_buffer = buffer(iel:end)
                      found = .true.
*                             ! not a comment
                  end if
*                             ! while not found
              end do
*                             ! name from buffer was null
          end if
 
*****     Save the name (Trim any leading blanks as we save the save
*         the name)
 
*                                               ! trim leading blanks
          IF( trimlen(new_buffer).gt.0 ) THEN
*                      ! ib is the character number of the first no blank
              ib = 1
*                      ! character.
              do while ( new_buffer(ib:ib).eq.' ' )
                  ib = ib + 1
              end do
 
              name = new_buffer(ib:trimlen(new_buffer))
 
          ELSE
              name = ' '
          END IF
 
*                             ! search not found
      end if
 
***** Thats all
*                             ! Process any error
1000  continue
 
      call report_error('IOSTAT',ierr,'read',buffer,1,'GET_NAME')
 
      return
      end
 
c-------------------------------------------------------------------------
 
      subroutine find_name( name,iel )
 
c************************************************************************
c     routine to find 'name' in the list of modules.  iel is returned
c     zero if name not found
c
c     T.Herring                   10:26 AM  WED.,  7  MAY , 1986
c                                 last modified <870421.1330>
c*************************************************************************
 
 
      include 'block.h'
 
*   i       - a loop counter
*   iel     - the entry in the modules list of 'name'
*   match   - the index returned from index
 
      integer*4 i, iel, match
 
*   found   - indicates match found
 
      logical found
 
*   name    - the name to be found in modules
 
 
      character*(*) name
 
***** START, loop over modules to see if we can name in modules
      found = .false.
      i     = 0
      iel   = 0
 
      do while ( i.lt.num_modules .and. .not. found )
 
          i = i + 1
          match = index(modules(i),name)
*                                 ! found
          if( match.gt.0 ) then
              found = .true.
              iel = i
          end if
      end do
 
***** Thats all
      return
      end
 
c------------------------------------------------------------------------
 
      subroutine push( mod_num, entry )
 
c************************************************************************
c     Routine to push mod_num and entry onto the stack
c
c     T.Herring                   10:35 AM  WED.,  7  MAY , 1986
c                                 last modified <870421.1330>
c***********************************************************************
 
 
      include 'block.h'
 
*   entry       - the entry in the list of routines called by
*               - this module
*   mod_num     - the number of module being processed
 
 
      integer*4 entry, mod_num
 
***** START, increment stack length and push
      stack_length = stack_length + 1
 
*                                             ! stack overflow
      if( stack_length.gt. max_stack ) then
 
          stop ' Stack overflow'
 
      end if
 
      stack(1,stack_length) = mod_num
      stack(2,stack_length) = entry
 
****  Thats all
      return
      end
 
c------------------------------------------------------------------------
 
      subroutine pop( mod_num, entry )
 
c**********************************************************************
c     Routine to pop the stack
c
c     T.Herring                   11:47 AM  WED.,  7  MAY , 1986
c                                 last modified <870421.1330>
c**********************************************************************
 
 
      include 'block.h'
 
*   entry       - entry popped from stack
*   mod_num     - number of module
 
      integer*4 entry, mod_num
 
***** START, take value off top of stack
 
      if( stack_length .le.0 ) then
          stop ' Stack underflow'
      end if
 
      mod_num = stack(1,stack_length)
      entry   = stack(2,stack_length)
 
      stack_length = stack_length - 1
 
***** Thats all
      return
      end
 
c----------------------------------------------------------------------
 
      subroutine report_error(type,ierr,operation,file,terminate,prog)
c
c     Routine to report a file manipulation error and possibly stop
c     program from running.  The error meassage is printed to the
c     log lu.
c
c Variables
c ---------
c type -- the type of error (usually FMGR or IOSTAT)
c ierr -- the error which ocurred
c operation -- what we were doing at the time
c file -- the file being manipulated
c terminate -- if non-zero causes the program to stop.
c prog -- the program which was running.
c
      character*(*) type, operation, file, prog
 
c
      integer*4 ierr, terminate
 
c
c Local variables
c ---------------
c trimlen -- HP utility for length of string
c error -- logical which indicates if there is an error
c
      integer*4 trimlen
c
      logical error
c
c.... see if error  (for FMGR less than 0 is error, for anything
c     else a non zero value is error
      error = .false.
      if( type(1:2).eq.'FM' ) then
         if( ierr.lt.0 ) error = .true.
      else
         if( ierr.ne.0 ) error = .true.
      end if
c
*                               ! no error
      if( .not.error ) return
c
c.... Report error
      write( *  ,100) type, ierr, operation, file(1:trimlen(file)),
     .   prog
 100  format(/,1x,a," error ",i4," occurred ",a,"ing file ",a,
     .   " in routine ",a)
c
c.... see if we should stop
      if( terminate.gt.0 ) then
         write( * ,110) prog
 110     format(1x,a," terminating")
         stop ' File error reported by REPORT_ERROR'
      end if
c
      return
      end
 
c--------------------------------------------------------------------------
 
      subroutine block_program( output_file, max_width, max_height )
*
*************************************************************************
*
*     Program to read the block_program.data file and output a command
c     file for plot which will show the reslationship between the
c     elements of the solve program.
c
c     T.Herring                   10:47 AM  TUE., 22  APR., 1986
c                                 Last modified <870421.1330>
c************************************************************************
 
 
*   decimaltoint - HP function for ASCII to int value
*   default_height - the height of the plot as orginally specified
*               - Initially equal to max_height.  Maxheight will
*               - change if the plot overflows its boundaries.
 
*   ierr        - error flag
*   len         - length of string
*   pos(2)      - the position of the subroutine name
 
*   psuedo_height - the psuedo height of the paper to make the
*               - characters fit in their boxes
*   psuedo_width - the psuedo width of the paper to make the
*               - charcaters fit in their boxes
 
*   max_height  - the maximum number of subroutine calls in program
*   max_width   - the maximum level of subroutine nesting in prog.
*   trimlen     - HP string length utility
 
      integer*4 default_height, ierr, pos(2),
     .    psuedo_height, psuedo_width, max_height, max_width, trimlen
 
*   x_offset    - view offset in the x direction
*   y_offset    - view offset in the y direction
 
*   x_scale     - scaling factor from x position to view
*   y_scale     - scaling factor from y position to view
 
*   x_length    - the view length in x
*   y_length    - the view length in y
 
*   view(4)     - the view size to be used for each subroutine name
 
      real*4 x_offset, y_offset, x_scale, y_scale, x_length, y_length,
     .    view(4)
 
*   sub_name    - subroutine name
 
      character*24 sub_name
 
*   temp_buff   - temporary storage of subroutine name
 
      character*22 temp_buff
 
*   output_file    - name of output file
 
      character*64 output_file
 
*   buffer      - line read from file
 
 
      character*80 buffer
 
***** Open the output file
 
      open(300,file=output_file, iostat=ierr, status='unknown')
      write(*,*) 'opened ', output_file
 
      call report_error('IOSTAT',ierr,'open',output_file, 1,
     .    'block_program')
 
****  write header to block_program.plt
      write(300,'("  charsz  1.2")')
      write(300,'("  Font  1 ")       ')
 
      psuedo_width = 2.9*25.* max_width
      psuedo_height = 5.0*max_height/0.77
 
      write(300,'("  ps_size  ",2i8)') psuedo_width, psuedo_height
 
      write(300,'("  scale  0 63.8  0 4.2 ")')
 
      x_scale  =  (0.98/max_width)
      x_length =   0.95 *x_scale
      x_offset =   0.1  *x_scale
 
      y_scale  =  (0.99/max_height)
      y_length =   0.80 *y_scale
      y_offset =  -0.1  *y_scale
 
      default_height = max_height
 
***** Now read the input file
 
*                   ! the scratch file
      rewind(200)
 
      do  while ( ierr.eq.0 )
 
          read(200,'(a)', iostat=ierr, err=1000, end=900 ) buffer
 
*****     Get the subroutine postion
          read(buffer,*) pos
 
*****     Get label
          temp_buff = buffer(10:)
 
          call trim_lead(temp_buff,ierr)
 
****      If no error continue
          if( ierr.eq.0 ) then
 
              sub_name = '"' // temp_buff(1:trimlen(temp_buff))
     .                       // '"'
 
****          Get the view based on pos
              view(1) = (pos(1)-1)*x_scale + x_offset
              view(2) = view(1) +  x_length
 
*****         See if we will go outsite the plot, If we do erase plot
*             (which will eject the paper on the 7150) and reset
*             max_height
              if( max_height - pos(2) - 1.le.0 ) then
                  write(300,'("  ERASE")')
                  max_height = max_height + default_height - 1
              end if
 
              view(3) = (max_height - pos(2) - 1)* y_scale + y_offset
              view(4) = view(3) +  y_length
 
*****         write out the line to plot file
              write(300,'("  view ",4f10.4)') view
              write(300,'("  label  0.2 0.5 1 0 ",a)')
     .                sub_name(1:trimlen(sub_name))
 
*****         Draw box  
                                  
              write(300,'("  xmn 0 0 ",3H" ",/,"  ymn 0 0 ",3H" ")')
              write(300,'("  xmx 0 0 ",3H" ",/,"  ymx 0 0 ",3H" " )')
c   
c             above replaces below (rwk 970920)
c              write(300,'("  xmn 0 0 ",3H" ",/,"  ymn 0 0 ",3H" ",/,
c     .                    "  xmx 0 0 ",3H" ",/,"  ymx 0 0 ",3H" " ) ')
c             the following commented out earlier
c             write(300,'("  xaxis_mn 0 0 ",3H" ")')
 
          end if
      end do
c
 900  continue
c
      write(300,'("  end")')
      close(100)
      close(200)
c
c**** END
 1000 continue
      if( ierr.ne.-1) then
          call report_error('IOSTAT',ierr,'read','file',0,
     .        'block_program')
 
      end if
 
      end
 
c.........................................................................
 
      subroutine trim_lead(buffer,ierr)
c
c     Routine to remove leading blanks from string.
c     See description of ierr returns
c
c Variables
c ---------
c buffer -- buffer to have its leading blanks removed
c ierr   -- a returned error number:
c           ierr = 0 if OK
c           ierr = -2 if string is all blanks
c           ierr = -3 if string is too long to be copied
c
      character*(*) buffer
 
c
      integer*4 ierr
 
c
c new_buffer -- dummy buffer used to copy buffer needed becuase of bug
c     in fortran complier
c ilen  -- length of buffer
c ib    -- counter used to find first non-blank character
c
      character*80 new_buffer
 
c
      integer*4 ilen, ib
 
c
c Functions used
c --------------
c Trimlen -- HP utility
c
      integer*4 trimlen
 
c
c.... Firstly make sure buffer is not too long to be copied
c     through 'new_buffer'
      ierr = 0
      ilen = trimlen(buffer)
c
      if( ilen.eq.0 )then
*                              ! set string empty error
         ierr = -2
         return
*                              ! see if too long
      else
*                              ! string too long
         if( ilen.gt.80 ) then
            ierr = -3
            return
         end if
      end if
c
c.... String is fine now find leading blanks
      ib = 1
*                                          ! check for blank
      do while ( buffer(ib:ib).eq.' ' )
         ib = ib + 1
      end do
c
c.... Now copy string back into itself. Due to compiler bug we must
c     use the intermediate buffer new_buffer
*                           ! kill off leading blanks
      if( ierr.ge.0 ) then
         new_buffer = buffer(ib:ilen)
         buffer     = new_buffer
      end if
c
      return
      end
 
*---------------------------------------------------------------------
 
      subroutine find_quotes( buffer, iel, quote )
 
**********************************************************************
*     Subroutine to see whether the start of the string we are
*     searching for is between a pair of quotes.  The start of
*     string we are looking for is in position IEL.  The charater
*     string QUOTE contains the single character delimiter.
*
*     If IEL is between quotes, then its value will b set to zero
*     indicating that no valid string was found.
*
*     T.Herring                   11:19 AM  TUE.,  2  SEPT, 1986
*
**********************************************************************
 
 
*   current_pos - the current position in the string
*               - as we search for pairs of the quotes.
*   end_quote   - the position in the string of the ending quote
*   iel         - starting position of current SEARCH string
 
*   start_quote - the position in the string of the starting pair
*               - of quotes.
 
*   trimlen     - HP utility ot return length of string.
 
      integer*4 current_pos, end_quote, iel, start_quote, trimlen
 
*   searching_for_quotes    - Indicates that we are still seaching
*               - for pairs of quotes
 
      logical searching_for_quotes
 
*   buffer      - the line read from the source file
 
      character*(*) buffer
 
*   quote       - the quote character (either ' or " )
 
 
      character*1 quote
 
***** Now search for pairs of quoutes and make sure that the 'SEARCH'
*     label does not appear inside quotes (both double and single)
 
      current_pos = 0
      searching_for_quotes = .true.
 
      do while ( searching_for_quotes )
 
          start_quote = index(buffer(current_pos+1:),quote)  +
     .                   current_pos
 
*                                               ! single quote found, see
          IF( start_quote.gt.current_pos ) THEN
*                                      ! where the end quote is.
 
              end_quote = index(buffer(start_quote+1:),quote) +
     .                     start_quote
 
*                                                  ! set end_quote to
              if( end_quote.eq.start_quote ) then
*                                                    ! end of string.
                  end_quote = trimlen(buffer)
                  searching_for_quotes = .false.
              end if
 
*****         Now see if our 'SEARCH' string appears inside these quotes
 
              if( iel.gt.start_quote .and. iel.lt.end_quote ) then
*                             ! 'SEARCH' inside string, so ignore.
*                             ! Indicates that no valid string was found
                  iel = 0
*                                                 ! We are finished
                  searching_for_quotes = .false.
 
              end if
*****         Update the current position in the string so that we can
*             search for the next pair of quotes.
 
              current_pos = end_quote + 1
*                                                       ! we are at end of
              if( current_pos.ge.trimlen(buffer) ) then
*                                                       ! string
                  searching_for_quotes = .false.
              end if
 
*                         ! No more quotes found, so stop searching
          ELSE
 
              searching_for_quotes = .false.
 
          END IF
 
      END DO
 
***** Thats all
      RETURN
      END
 
 
CTITLE CHECK_FORTRAN
 
      subroutine check_fortran( buffer, fortran )
 
 
*     Routine to check if current code is fortran.  If it is not
*     (based on directive ASMB) then the code is commented out.
*                                  9:34 AM  TUE., 21  APR., 1987
 
*   fortran     - Indicates if fortran code
 
      logical fortran
 
*   buffer      - Line read from source file
 
 
      character*(*) buffer
 
***** See if fortran dirertive
 
*                                       ! fortran again, set true
      if( buffer(1:3).eq.'FTN' ) then
          fortran = .true.
      end if
 
*                                       ! ASMB, set fortran false.
      if( buffer(1:4).eq.'ASMB' ) then
          fortran = .false.
      end if
 
***** If not fortran then commnet
      if( .not. fortran ) then
          buffer(1:1) = '*'
      end if
 
***** Thats all
      return
      end
 
CTITLE REPLACE_PROGRAM
 
      subroutine replace_program( lu_user, buffer )
 
 
*     Routine to replace the PROGRAM statement with a SUBROUTINE statement
*     if more than one program module is found in the set of source codes.
*     This should allow segmented programs to be processed.
*                                  9:47 AM  TUE., 21  APR., 1987
 
*   full_len    - Full length of the string or used length
*   i,j     - Loop counters
*   ierr    - Error from READ_LINE
*   indx    - Index for READ_LINE
*   lu_user - The user's LU number
*   pos_program - Position of the program statement
*   trimlen - HP function to return used length of string
 
      integer*4 full_len, i, ierr, indx, lu_user, 
     .    trimlen
 
 
      character*7 name
 
*   buffer  - the line read from the source file
 
 
      character*(*) buffer
 
***** Get the first thing on the line and see if it program
 
      indx = 1
      call read_line(buffer,indx,'CH',ierr,i,name)
 
*                                     ! We need to replace
      if( name.eq.'PROGRAM' ) then
 
          full_len = min(len(buffer)-3, trimlen(buffer) )
          write(lu_user,100) buffer(1:full_len)
  100     format(' Replacing: ',a)
 
*         Now replace, Move string up
 
          do i = full_len+3, indx, -1
              buffer(i:i) = buffer(i-3:i-3)
          end do
 
*         Now replace PROGRAM with SUBROUTINE
 
          buffer(indx-7:indx+2) = 'SUBROUTINE'
 
*         Tell user what we did
          write(lu_user,120) buffer(1:full_len+3)
  120     format(' With     : ',a)
      end if
 
***** Thats all
      return
      end
 
CTITLE REPLACE_SEGLD
 
      subroutine replace_segld( lu_user,buffer )
 
 
*     Routine to remove SEGLD from a line, so that the call looks
*     like CALL module directly.
*     Routine assumes that the module name is given as 6HNAME, thus
*     the line of code should look like:
 
*     CALL SEGLD(6Hname..,....).
 
*     This line will be replaced by:
 
*     CALL name..!..... NOTE the comma is replace by ! to make a
*     valid line of fortran.
 
*   i,j     - Loop counters
*   len     - Length of the string
*   lu_user - User LU for messages
*   pos_comma   - Position of first comma
*   pos_segld   - Position of SEGLD in string
*   trimlen - Used portion of string
 
      integer*4 i, len, lu_user, pos_comma, pos_segld, trimlen
 
*   buffer  - Line read from source file
 
      character*(*) buffer
 
****  See if SEGLD in string
 
      pos_segld = index( buffer,'SEGLD(6H' )
 
*                                 ! Replace
      if( pos_segld.ge.7 ) then
          len = trimlen(buffer)
          write(lu_user,100) buffer(1:len)
 100      format(' Replacing : ',a)
 
          do i = pos_segld, len-8
              buffer(i:i) = buffer(i+8:i+8)
          end do
 
*         Now replace first comma
          pos_comma = index( buffer,',' )
          if( pos_comma.gt.0 ) then
              buffer(pos_comma:pos_comma) = '!'
          end if
 
*         Tell user new string
          write(lu_user,120) buffer(1:len-8)
  120     format(' With      : ',a)
 
      end if
 
***** Thats all
      return
      end
 
CTITLE REPLACE_LINCA
 
      subroutine replace_linca( lu_user,buffer )
 
 
*     Routine to remove LINCA from a line, so that the call looks
*     like CALL module directly.
*     Routine assumes that the module name is given as 6HNAME, thus
*     the line of code should look like:
 
*     CALL LINCA(6Hsegmnt,6Hsubnam,.....
 
*     This line will be replaced by:
 
*     CALL subnam!..... NOTE the comma is replace by ! to make a
*     valid line of fortran.
 
*   i,j     - Loop counters
*   len     - Length of the string
*   lu_user - User LU for messages
*   pos_comma   - Position of first comma
*   pos_linca   - Position of LINCA in string
*   trimlen - Used portion of string
 
      integer*4 i, len, lu_user, pos_comma, pos_linca, trimlen
 
*   buffer  - Line read from source file
 
      character*(*) buffer
 
****  See if LINCA in string
 
      pos_linca = index( buffer,'LINCA(6H' )
 
*                                 ! Replace
      if( pos_linca.ge.7 ) then
          len = trimlen(buffer)
          write(lu_user,100) buffer(1:len)
 100      format(' Replacing : ',a)
 
*                                        ! Allow space for Segment name
          do i = pos_linca, len-17
              buffer(i:i) = buffer(i+17:i+17)
          end do
 
*         Now replace first comma
          pos_comma = index( buffer,',' )
          if( pos_comma.gt.0 ) then
*                                            
              buffer(pos_comma:pos_comma) = '!'
          end if
 
*         Tell user new string
          write(lu_user,120) buffer(1:len-17)
  120     format(' With      : ',a)
 
      end if
 
***** Thats all
      return
      end
