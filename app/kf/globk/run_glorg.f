CTITLE RUN_GLORG
      subroutine run_glorg ( cmd_file, out_file, opts )

      implicit none 
 
 
*     This routine sets up the glorg run string and runs the program
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   idev        - Output device
*   ierr        - Error flag
*   ifive(5)    - 5 rmpar parameters passed back from glorg
*               - (NOT USED)
*   len         - length of program name
*   opts        - Options for output
*   trimlen     - HP function for length of string
*   offset_com - Command file offset (Unit 100)
 
      integer*4 ierr, ifive(5), len, opts, trimlen, offset_com,
     .          len_out

*   out_file    - Name of the output file
*   cmd_file    - Name of the command file

      character*(*) out_file, cmd_file
 
*   prog_root   - Just the program name
 
      character*9 prog_root
 
*   prog_name   - Program name and runstring
 
      character*128 prog_name
 
      data prog_root / 'glorg' /
 
***** Add in the runstring for glorg
 
      len = trimlen( prog_root )
      len_out = trimlen( out_file )
      write(prog_name,100, iostat=ierr )
     .    prog_root(1:len), out_file(1:len_out), opts,
     .    glr_cmd_file(1:max(1,trimlen( glr_cmd_file))),
     .    glb_com_file(1:max(1,trimlen(glb_com_file)))
 
  100 format(a,' ',a,' ',(i6,' '),a,1x,a)
 
      call report_error('IOSTAT',ierr,'writ','GLORG RUNSTRING',1,
*                                     ! kill, if error
     .                  'RUN_GLORG')
 
***** Now run the program
 
      call execute( prog_name, ifive, 1,6, offset_com )
 
***** Thats all
      return
      end
 
