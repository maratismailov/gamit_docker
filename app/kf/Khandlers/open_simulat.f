CTITLE OPEN_SIMULATION
 
      subroutine open_simulation(SO_name, ierr)
 

      implicit none
 
*     This subroutine is used in the simulation software to ensure
*     that only KalObs files specifically setup for simulations can
*     be opened.  To create such a simulation file, use the program
*     COSI to COpy the SImulation file from a standard KalObs file.
*
*     This routine simplies calls OPEN_KALOBS, reads the values
*     block to get the expt_title and checks that !SIMULATION appears
*     appears in the title.  If it does not the user's program is
*     stoped becuase of the risk of overwritting a valid KalObs file.
*
*                                  9:47 AM  FRI., 19  JUNE, 1987
 
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
 
*   ierr    - Error returned from OPEN_KALOBS.  This error
*           - is simply passed bak to the user's program.
*   trimlen - HP function for length of string
 
      integer*4  ierr, trimlen
 
*   SO_name - Name of the simulation file to be opened.
 
 
      character*(*) SO_name
 
***** Try to open the KalObs file
 
      call open_Kalobs(SO_name, ierr)
 
*     Now read the values block, and check that expt_title is OK
 
      if( ierr.ge.0 ) then
          call rw_KalObs_block('R','Values',values, ierr, 0)
 
*         Check experiment title
*                                                       ! Not a simulation
          if( expt_title(13:23).ne.'!SIMULATION' ) then
*                                           ! file so stop program
              write( * ,100) SO_name(1:max(1,trimlen(SO_name)))
  100         format(//' **** DISASTER **** ',a,' is not a simulation',
     .                 ' file',/,' Use COSI to make a copy for',
     .                 ' simulations'//)
              stop ' Not a simulation file. Stop in OPEN_SIMULATION'
          end if
      end if
 
***** Thats all
      return
      end
