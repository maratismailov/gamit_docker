CTITLE RW_COMMON
 
      subroutine rw_common(option)

      implicit none
 
c
c     Routine to read, write and close the SOLVK common block
c     so that it can be passed between the programs which form
c     the KALMAN filter suite of programs
 
*                                 ! the solvk parameter file
      include '../includes/kalman_param.h'
c
*                                 ! the solvk control common block
      include '../includes/solvk_common.h'
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
      if( trimlen(com_file).eq.0 ) then
          com_file = Default_com_file
      end if
 
c.... Firstly see if the common file is open at the moment
*                                  ! open the file
      if( .not.common_open ) then
 
*         Try to create the common
          call create_SOLVK( com_file )
          common_open = .true.
 
      end if
c
c.... Common file is now opened or created, update the open flag
 
      common_open = .true.
 
      call rw_solvk_block( option,'CONTROL', control, ierr)
      call report_error('FMP',ierr,option,com_file,1,'RW_COMMON')
 
      call rw_solvk_block( option,'MARKOV' , markov , ierr)
      call report_error('FMP',ierr,option,com_file,1,'RW_COMMON')
 
***** If we write the file, then close as well
 
      opt = option
      call casefold( opt )
      if( opt.eq.'W' ) then
          call close_solvk( ierr )
          call report_error('FmpClose',ierr,'clos',com_file,0,
     .                      'RW_COMMON')
 
      common_open = .false.
 
      end if
c
      return
      end
 
 
