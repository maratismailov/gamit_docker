
CTITLE COMPRESS_site
 
      subroutine compress_site( num_obs_site )
 
      implicit none

*     This routine will remove from the list of sites those
*     which have no observations.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
  
* PASSED VARIABLES

*  num_obs_site - Number of obsverations at this site

       integer*4 num_obs_site(qnum_sites)
 
* LOCAL VARIBALES
 
*   i,j     - Loop counters
 
      integer*4 i,j

***** First get the total number of observations.  If it zero do
*     nothing since this is an old format h-file.

      j = 0

      do i = 1, qnum_sites
         j = j + num_obs_site(i)
      end do

*     If there is no observations, then just return.

      if( j.eq.0 ) RETURN
 
****  Scan over list of satelites finding those with no observations
 
      i = 0
      do while (i.lt.qnum_sites )
          i = i + 1
*                                 ! Remove from list
          if( num_obs_site(i).eq.0 ) then
              do j = i, qnum_sites-1
                  num_obs_site(j) = num_obs_site(j+1)
                  qsite_names(j) = qsite_names(j+1)
                  qfull_names(j) = qfull_names(j+1)
                  qrecv_ty(j) = qrecv_ty(j+1)
                  qrecv_sn(j) = qrecv_sn(j+1)
                  qrecv_fw(j) = qrecv_fw(j+1)
                  qante_ty(j) = qante_ty(j+1)
                  qradome_ty(j) = qradome_ty(j+1)
                  qante_sn(j) = qante_sn(j+1)
                  qant_mod(j) = qant_mod(j+1)
              end do
 
*             Decrement number and pointer
              qnum_sites = qnum_sites - 1
              i = i - 1
          end if
      end do
 
****  Thats all
      return
      end
 
 
 
