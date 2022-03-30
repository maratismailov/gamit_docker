CTITLE INIT_QAPR_COV

      subroutine init_qapr_cov
      
      implicit none

*     Routine to clear the apriori constraints arrays

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* LOCAL VARIABLES
*  i,j,k   - Loop counters

      integer*4 i, j, k
 
*     Loop over the 3x3's for sites
     
      do i = 1, qnum_sites
         do j = 1, 3
            do k = 1, 3
               qapr_cov(k,j,i) = 0.d0
            end do
         end do
      end do
      
      qnum_apr_cov = 0

*     Loops of the 6x6 for satellite IC's.  Radiation aprioris
*     goto into the diagonal elements      
      do i = 1, qnum_svs
         do j = 1, 6
            do k = 1,6
               qapr_svs(k,j,i) = 0.d0
            end do
         end do
      end do
      
      qnum_apr_svs = 0

*     Loop over those things are just diagonals               
      do i = 1,qnum_parn
         qapr_diag(i) = 0
      end do
      
      qnum_apr_diag = 0
      
****  Thats all
      return
      end
            
