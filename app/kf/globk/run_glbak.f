CTITLE RUN_GLBAK
      subroutine run_glbak

      implicit none 
 
 
*     This routine sets up the GLBAK run string and runs the program
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   ierr        - Error flag
*   ifive(5)    - 5 rmpar parameters passed back from GLBAK
*               - (NOT USED)
*   len         - length of program name
*   trimlen     - HP function for length of string
*   offset_com - Command file offset (Unit 100)
 
      integer*4 ierr, ifive(5), len, trimlen, offset_com
 
*   prog_root   - Just the program name
 
      character*9 prog_root
 
*   prog_name   - Program name and runstring
 
      character*128 prog_name
 
      data prog_root / 'glbak' /
 
***** Add in the runstring for GLBAK
 
      len = trimlen( prog_root )
      write(prog_name,100, iostat=ierr )
     .    prog_root(1:len),
     .    glb_com_file(1:max(1,trimlen(glb_com_file)))
 
  100 format(2(a,','))
 
      call report_error('IOSTAT',ierr,'writ','GLBAK RUNSTRING',1,
*                                     ! kill, if error
     .                  'RUN_GLBAK')
 
***** Now run the program
 
      call execute( prog_name, ifive, 1, 6, offset_com )
 
***** Thats all
      return
      end
 
