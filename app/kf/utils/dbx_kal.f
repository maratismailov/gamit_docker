
      program dbx_kal
 
*     This program simply reads a KalObs file header and re-writes
*     to disk.  The program is meant to run under dbx so that 
*     sundry quantities in the KalObs header can be changed. 

*     The first use was to change unalias source names.
*
*     Runstring
*     % dbx_kal KalObs_file 
*
*     where KalObs_file is the name of the KalObs file to be
*         updated
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'

* ierr     - IOSTAT error
* rcpar    - Gets runstring
* len_run  - Length of string

      integer*4 ierr, rcpar, len_run
 
* Update   -  Logical to say to write the header.  Should be
*             set true in dbx 
 
      logical update
 
*   KO_name     - Name of KalObs file
*   runstring   - Runstring read from command line
 
 
      character*64 KO_name, runstring
 
****  Decode the runstring
 
      Update = .false. 
 
      len_run = rcpar(1,runstring)
      if( len_run.gt.0 ) then
          KO_name = runstring
      else
          call proper_runstring('dbx_kal',6,1)
      end if
 
***** Try to open the KalObs file
      call open_KalObs( KO_name, ierr)
      if( ierr.ne.0 ) stop ' dbx_kal: Error opening KalObs file'
 
      call rw_KalObs_header('R', ierr)
      call report_error('FmpRead',ierr,'read','KalObs header',1,
     .                  'dbx_kal')
 
***** Write out KalObs header, is Update has been changed.
*     Use:     
*     Assign update = 1   to set to true.
      if( update ) then
          call rw_KalObs_header('W',ierr)
      else
          call proper_runstring('dbx_kal',6,0)
      end if

      call close_KalObs(ierr)
 
      end
 
 
