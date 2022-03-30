CTITLE RUN_GLOUT
      subroutine run_glout ( out_file, opts )

      implicit none 
 
 
*     This routine sets up the GLOUT run string and runs the program
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   idev        - Output device
*   ierr        - Error flag
*   ifive(5)    - 5 rmpar parameters passed back from GLOUT
*               - (NOT USED)
*   len         - length of program name
*   opts        - Options for output
*   trimlen     - HP function for length of string
*   offset_com - Command file offset (Unit 100)
 
      integer*4 ierr, ifive(5), len, opts, trimlen, offset_com,
     .          len_out

*   out_file    - Name of the output file

      character*(*) out_file
 
*   prog_root   - Just the program name
 
      character*9 prog_root
 
*   prog_name   - Program name and runstring
 
      character*128 prog_name
 
      data prog_root / 'glout' /
 
***** Add in the runstring for GLOUT
 
      len = trimlen( prog_root )
      len_out = trimlen( out_file )
      write(prog_name,100, iostat=ierr )
     .    prog_root(1:len), out_file(1:len_out), opts,
     .    glb_com_file(1:max(1,trimlen(glb_com_file)))
 
  100 format(a,' ',a,' ',(i6,' '),a)
 
      call report_error('IOSTAT',ierr,'writ','GLOUT RUNSTRING',1,
*                                     ! kill, if error
     .                  'RUN_GLOUT')
 
***** Now run the program
 
      call execute( prog_name, ifive, 1,6, offset_com )
 
***** Thats all
      return
      end
 
