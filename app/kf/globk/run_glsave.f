CTITLE RUN_GLSAVE
      subroutine run_glsave

      implicit none 
 
*     This routine sets up the GLSAVE run string and runs the program
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   ierr        - Error flag
*   ifive(5)    - 5 rmpar parameters passed back from GLFOR
*               - (NOT USED)
*   len         - length of program name
*   trimlen     - HP function for length of string
*   offset_com - Command file offset (Unit 100)
 
      integer*4 ierr, ifive(5), len, trimlen, offset_com
 
*   prog_root   - Just the program name
 
      character*9 prog_root
 
*   prog_name   - Program name and runstring
 
      character*128 prog_name
 
      data prog_root / 'glsave' /
 
***** Add in the runstring for GLSAVE
 
      len = trimlen( prog_root )
      write(prog_name,100, iostat=ierr )
     .    prog_root(1:len),
     .    glb_com_file(1:max(1,trimlen(glb_com_file)))
 
  100 format(2(a,','))
 
      call report_error('IOSTAT',ierr,'writ','GLSAVE RUNSTRING',1,
*                                      ! kill, if error
     .                  'RUN_GLSAVE')
 
***** Now run the program
 
      call execute( prog_name, ifive, 1, 6, offset_com )
 
***** Thats all
      return
      end
 
