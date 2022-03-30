CTITLE 'PROPER_RUNSTRING'
 
      subroutine proper_runstring( help_file, program, kill)

      implicit none
 
 
*     Routine to ouput the runstring for nutcorr
*     To use this subroutine the user needs to set the enviroment
*     variable HELP_DIR. This directory is then preprended to
*     the file name
 
*           ierr        - IOSTAT error (opening help file)
*   trimlen             - Length of string
*   kill                - indicates iof program should be killed.
*                         0 means do not stop program
*                        >0 kill program
*                        <0 kill program but list help in 24 line
*                           segment (with q to quit)
*   len_dir             - Length of directory name
*   cnt                 - Keeps track of number of lines written
*                        so that user can be quizzed to continue if
*                        kill set < 0
 
      integer*4 ierr, trimlen, kill, len_dir, cnt
 
*            line       - Line read from input.
 
      character*100 line
 
*             help_file - Name of help file
*       program         - Name of program.
 
      character*(*) help_file, program

*     Full_help_file  - Help file name with directory added.

      character*256 Full_help_file

*     ans - Answer to query for more lines
      character*4 ans
 
***** open the help file.  First construct name
      call getenv('HELP_DIR',Full_help_file)
      len_dir = trimlen(Full_help_file)
      if( len_dir.gt.0 ) then
           if( Full_help_file(len_dir:len_dir).eq.'/' ) then
               Full_help_file(len_dir+1:) = help_file
           else
               Full_help_file(len_dir+1:) = '/' // help_file
           endif
      else
           Full_help_file = help_file
      end if
 
      open(999, file=full_help_file, iostat=ierr, status='old')
*     NOD TAH 980813: If file not found tell user that to do
      if( ierr.ne.0 ) then
          write(*,110)
 110      format(/,'**ERROR** Opening help file.  Check that',
     .             ' enviroment variable HELP_DIR points to',/,
     .             '          directory with help files.')
      end if
      call report_error('IOSTAT',ierr,'open',full_help_file, 0,
     .                  program)
 
*                             ! Loop until EOF
      cnt = 0
      do while ( ierr.eq.0 )
          read(999,'(a)', iostat=ierr) line
          cnt = cnt + 1
          if( ierr.eq.0 ) write(*,'(a)') line(1:max(trimlen(line),1))
          if( mod(cnt,24).eq.0 .and. kill.lt.0 ) then

*             see if we should continue list
              write(*,120)
 120          format('More? (q to quit) ',$)
              read(*,'(a)') ans
              if( ans(1:1).eq.'q' .or. ans(1:1).eq.'Q' ) then
                 ierr = -1
              end if
          end if
      end do
 
***** close the input
      close(999)
 
****  Now see if we should stop
      if( kill.ne.0 ) then
          write(*,200) program(1:max(trimlen(program),1))
  200     format(' PROGRAM ',a,' terminating')
          stop ' Incomplete runstring'
      end if
 
****  Thats all
      return
      end
 
 
