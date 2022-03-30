CTITLE 'REPORT_ERROR'
 
      subroutine report_error(type,ierr,action,file, kill, prog)
  
      implicit none

*     Routine to Report that an error has occurred and possibly
*     stop the program (depending on the value of Kill: no kill if
*     0, otherwize kill)
 
*         ierr      - IOSTAT error
*         kill      - Kill program if value is not zero.
*         len_type  - Used length of type
*         len_action    - Used length of action
*         len_file  - Used length of file
*         len_prog  - Used length of prog
*         trimlen   - Get length of string
 
      integer*4 ierr, kill, len_type, len_action, len_file, len_prog,
     .    trimlen
 
*             type  - Type of error
*   action          - Thing being done when error occurred
*   file            - Name of file being acted on
*   prog            - Program/subroutine name of calling
*                   - routine.
 
      character*(*) type, action, file, prog
      
*   Additional variables needed for use of report_stat with this
*   routine

*   rcpar   - Get the name of the program calling this routine
*   len_mod - Length of the module name
*   jerr    - IOSTAT error creating message line
*   message - The concatinated error message
*   module  - name of program

      integer*4 rcpar , len_mod, jerr
      character*128 module
      character*128 message
      
 
****  If error is not zero, print message
      if( ierr.ne.0 ) then          
          call check_ascii( type )
          call check_ascii( action )
          call check_ascii( file )
          call check_ascii( prog )

          len_type = max(1,trimlen(type))
          len_action = max(1,trimlen(action))
          len_file = max(1,trimlen(file))
          len_prog = max(1,trimlen(prog))
 
          if( len_file+len_action+len_type+len_prog.lt.55 ) then
              write(*,100) type(1:len_type), ierr,
     .                     action(1:len_action),
     .                     file(1:len_file), prog(1:len_prog)
  100         format(1x,a,' error ',i6,' occurred ',a,'ing ',a,
     .                    ' in ',a)
          else
              write(*,120) type(1:len_type), ierr,
     .                     action(1:len_action),
     .                     file(1:len_file), prog(1:len_prog)
  120         format(1x,a,' error ',i6,' occurred ',a,'ing',/,
     .                  a,/,' in ',a)
          end if

****      Set up the report_stat string
          len_mod = rcpar(0, module)
          if( len_mod.le.0 ) then
              module = 'MAIN'
              len_mod = 4
          end if
          write(message,140, iostat=jerr) type(1:len_type), 
     .                                action(1:len_action)
  140     format(a,' error ',a,'ing file')
          if( kill.ne.0 ) then
              call report_stat('fatal',module, prog, file, 
     .                         message,ierr)
          else
              call report_stat('warning',module, prog,file, 
     .                          message,ierr)
          end if
 
****      Now see if we should kill
          if( kill.ne.0 ) then
              write(*,200) prog(1:len_prog)
  200         format(' Program terminating in routine ',a)
              stop ' Stopped from report_error'
          end if
      end if
 
***** Thats all
      return
      end
 
