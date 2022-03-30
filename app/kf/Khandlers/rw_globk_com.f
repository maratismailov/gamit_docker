CTITLE RW_GLOBK_COMMON
 
      subroutine rw_globk_common(option)
 

      implicit none
c
c     Routine to read, write and close the GLOBK common block
c     so that it can be passed between the programs which form
c     the KALMAN filter suite of programs
 
*                                 ! the globk parameter file
      include '../includes/kalman_param.h'
c
*                                 ! the globk control common block
      include '../includes/globk_common.h'
c
c
c Variables
c ---------
c option -- this variable determines which function the subroutine
c     will perform.  This options are:
c     'R' -- read the common.  The common file will be opened
c     first if this option is invoked and the file has not yet been
c     opened.
c     'W' -- this option will write the file and then close it.
c     If the file does not yet exist, it will be created by this
c     routine.
* MOD TAH 051211: Added new option to allow glorg to write a new common
*     file
*     'N' -- Write new file and ignores any open status
c
c Local variables
c ---------------
c ierr -- an error variable to check that the file manipulations are
c     okay.
c common_open -- local variable which indicates if the common file
c     is open
c icrt -- the terminal being used -- we must use this alias because
c     if the common is being read and an error ocurrs we do not know
c     the user terminal
c
 
      character*(*) option
 
*   opt     - Upper case version of option
      character*1 opt
 
c
*   ierr    - FMGR file error flag
*   trimlen - HP trimlen function.
 
      integer*4 ierr, trimlen
 
c
 
      logical common_open
 
      save common_open
c
      data common_open / .false. /
c
c.... Set the default name for the common (ie if it was not passed)
      if( trimlen(glb_com_file).eq.0 ) then
          glb_com_file = glb_com_default
      end if
 
***** Check option
      opt = option
      call casefold ( opt )
 
c.... Firstly see if the common file is open at the moment
*                                                   ! create the file
      if( .not.common_open .and. opt.eq.'W' ) then
 
*         Try to create the common
          call create_GLOBK( glb_com_file )
          common_open = .true.
      end if
* MOD TAH 051211: See if N status is set
      if( opt.eq.'N' ) then
*         Try to create the common
          call create_GLOBK( glb_com_file )
          common_open = .true.
          opt = 'W'
      end if
c
*                                                  ! open the file
      if( .not.common_open .and. opt.eq.'R' ) then
 
          call open_GLOBK( glb_com_file, ierr )
          if ( ierr.ge. 0 ) common_open = .true.
 
      end if
 
*     Read each of the blocks
      call rw_globk_block( opt,'CONTROL', glb_control, ierr)
      call report_error('FMP', ierr, opt, glb_com_file,
     .                   1, 'rw_globk_common')
 
      call rw_globk_block( opt,'MARKOV' , glb_markov , ierr)
      call report_error('FMP', ierr, opt, glb_com_file,
     .                   1, 'rw_globk_common')
 
 
      call rw_ema_block ( opt,'EMA'     , glb_ema, ierr)
      call report_error('VREAD/WRITE', ierr, opt, glb_com_file,
     .                   1, 'RW_GLOBK_COMMON')
 
***** If we write the file, then close as well
 
      opt = option
      call casefold( opt )
      if( opt.eq.'W' .or. opt.ne.'N' ) then
          call close_globk( ierr )
          call report_error('FmpClose',ierr,'clos',glb_com_file,0,
     .                      'rw_globk_common')
 
          common_open = .false.
 
      end if
c
      return
      end
 
