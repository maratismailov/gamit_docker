CTITLE READ_CTOG_CMDS
 
      subroutine read_ctog_cmds

      implicit none
 
*     This rotuine will read the the ctogobs commands and act accordingly.  It will
*     exit with either the END command or when EOF is reached.  The command file
*     is opened and read here.  The command file uses the standard command rules. 

* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES  - None


* LOCAL PARAMETERS

 
* LOCAL VARIABLES
 
*   iel     - command number
*   indx    - Position in string.
*   jndx    - Temp value for get_cmd
*   ierr    - IOSTAT error on read
*   trimlen - Length of string
 
      integer*4 iel, indx, jndx, ierr, trimlen
 
*   inline  - Line read from command
 
      character*256 inline
 
*   first_word  - Word of line
 
      character*(ctog_cmd_len) first_word
 
*   finished    - Indicates we are finished
 
      logical finished
 
****  Start, open the command file and loop until finised

      open(100, file= ctogobs_com_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open', ctogobs_com_file, 0,
     .                  'read_ctog_cmds')

      if( ierr.ne.0 ) then
          write(*,100) ctogobs_com_file(1:max(1,
     .                           trimlen(ctogobs_com_file)))
 100      format(' Problem with command file ',a,/,
     .           ' **WARNING** Program running with default settings')
          RETURN
      end if

      finished = .false.
      do while ( .not. finished )
 
          read(100,'(a)', iostat=ierr ) inline
          if( ierr.ne.0 ) then
              finished = .true.
          else
 
*****         See if we have command
              if( inline(1:1).eq.' ' .and. trimlen(inline).gt.0) then
                  indx = 1
                  call GetWord( inline, first_word, indx)
                  call casefold(first_word)
                  jndx = 1
                  call get_cmd(first_word, ctog_commands,
     .                num_ctog_cmds, iel, jndx)
                  if( iel.gt.0 ) then 
                      call proc_ctog_cmd(inline, indx, iel, finished,
     .                    first_word )
                  else
                      if( iel.eq.-1 ) then 
                          write(*,200) inline(1:trimlen(inline))
 200                      format(' COMMAND NOT FOUND: ',a)
                      else if( iel.eq.-2 ) then
                          write(*,220) inline(1:trimlen(inline))
 220                      format(' AMBIGUOUS COMMAND:',a)
                      else if( iel.ne.0 ) then
                          call report_error(' ',ierr,'decod',
     .                          inline,0,'read_ctog_commands')
                      end if
                  end if

*                     ! Processing command
              end if
*                     ! Error OK
          end if
*                     ! Looping until EOF or END
      end do
 
****  Thats all
      return
      end
 
