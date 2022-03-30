CTITLE RW_KALOBS_HEADER
 
      subroutine rw_KalObs_header(option,ierr)

      implicit none
c
c     Routine to read or write the KalObs header records.  OPTION
c     is a character variable.  If OPTION(1:1) = 'W', the record
c     will be written, if 'R', read.  Assumes KalObs is open.
c
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
 
 
      integer*4 ierr
 
 
      character*(*) option
 
 
      character*1 opt
 
c
c.... Get first character of option
      opt = option(1:1)
      call CaseFold(opt)
c
***** Now call utility to read/write each of the blocks in the
*     header
 
      call rw_KalObs_block( opt,'VALUES',values  , ierr, 0)
      call rw_KalObs_block( opt,'NAMES', names   , ierr, 0)
      call rw_KalObs_block( opt,'APR'  , aprioris, ierr, 0)
 
***** Thats all
      return
      end
 
 
