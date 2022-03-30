CTITLE RUN_GLINIT
      subroutine run_glinit

      implicit none 
 
 
*     This routine sets up the GLINIT run string and runs the program
*
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   copy_crt    - copy of crt_unit
*   copy_prt    - copy of prt_unit
*   copy_log    - copy of log_unit
*   ierr        - Error flag
*   ifive(5)    - 5 rmpar parameters passed back from GLINIT
*               - (NOT USED)
*   len         - length of program name
*   trimlen     - HP function for length of string
*   offset_com - Command file offset (Unit 100)
 
      integer*4 copy_crt, copy_prt, copy_log, ierr, ifive(5), len,
     .    trimlen, offset_com
 
c     integer*4 system
 
*   prog_root   - Just the program name
 
      character*16 prog_root
 
*   copy_mar_file   - Copy of the markov file name
*   copy_prt_file   - Copy print file name
*   copy_log_file   - Copy log file name
 
      character*256 copy_mar_file, copy_prt_file, copy_log_file

*   copy_comopt -- Copy of command option
      character*256 copy_comopt
 
*   prog_name   - Program name and runstring
 
      character*256 prog_name
 
c     logical opened
c     character*128 file
 
      data prog_root / 'glinit' /
 
***** Print message
      write(crt_unit, 100)
  100 format(/' Initalising Global solution'/)
 
***** save unit numbers since these are not passed to GLINIT which
*     creates the common
 
      copy_crt = crt_unit
      copy_prt = prt_unit
      copy_log = log_unit
      copy_mar_file = glb_mar_file
      copy_prt_file = glb_prt_file
      copy_log_file = glb_log_file
      copy_comopt = comopt

*     See if common file name is passed.  If it is not then use default
      if( trimlen(glb_com_file).eq.0 ) then
          glb_com_file = glb_com_default
      end if
      if( trimlen(sort_file).eq.0 ) then
          sort_file = ''' '''
      end if
      if( trimlen(eq_inp_file(1)).eq.0 ) then
          eq_inp_file(1) = ''' '''
      end if
      if( trimlen(glb_svs_file).eq.0 ) then
          glb_svs_file = ''' '''
      end if
 
***** Add in the runstring for GLINIT

      len = trimlen( prog_root )
      write(prog_name,200, iostat=ierr )
     .    prog_root(1:len), list_file(1:max(1,trimlen(list_file))),
     .    glb_com_file(1:max(1,trimlen(glb_com_file))),
     .    sort_file(1:max(1,trimlen(sort_file))),
     .    sort_direction,
     .    eq_inp_file(1)(1:max(1,trimlen(eq_inp_file(1)))),
     .    glb_svs_file(1:max(1,trimlen(glb_svs_file)))
 
  200 format(4(a,' '),i2,' ',a,1x,a)
 
      call report_error('IOSTAT',ierr,'writ','GLINIT RUNSTRING',1,
*                                      ! kill, if error
     .                  'RUN_GLINIT')
 
***** Now run the program
 
      call execute( prog_name, ifive, 1,100, offset_com )
 
***** Check return values to make sure that GLINIT run OK
 
*                                     ! Run had problems, so program
      if( ifive(2).ne.0 ) then
          stop ' GLOBK Stop: Problem running GLINIT'
      end if
 
***** Now open and read the globk common
 
      call rw_globk_common( 'R' )
 
***** Now restore the unit numbers
      crt_unit = copy_crt
      prt_unit = copy_prt
      log_unit = copy_log
      glb_mar_file = copy_mar_file
      glb_prt_file = copy_prt_file
      glb_log_file = copy_log_file
      comopt = copy_comopt 

***** Thats all
      return
      end
 
