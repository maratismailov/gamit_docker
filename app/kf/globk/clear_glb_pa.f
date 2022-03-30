CTITLE CLEAR_GLB_PARN
 
      subroutine clear_glb_parn
 
      implicit none 

 
*     routine to clear all of the parameter counting arrays
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   i,j,k   - loop counters
 
 
      integer*4 i,j
 
***** Loop over all parn variables
 
      do i = 1, gnum_sites
          parn_axo(1,i) = 0   ! axis offset
          parn_axo(2,i) = 0   ! axis offset rate
          parn_atm(i)   = 0   ! Atmospheric delay

          do j = 1,3
*                                         ! position
              parn_site(j,1,i) = 0
*                                         ! position rate
              parn_site(j,2,i) = 0
*                                         ! h/l/lag
              parn_tid (j,i)   = 0
*
              parn_log(j,i)  = 0          ! Earthquake log terms
          end do
      end do
 
      do i = 1, gnum_sources
*                                 ! Ra and dec
          do j = 1,2
              parn_source(j,1,i) = 0
              parn_source(j,2,i) = 0
          end do
      end do

      do i = 1, gnum_svs
          do j = 1,max_svs_elem
             parn_svs(j,i) = 0
          end do
      end do
 
*                                 ! Site origin translation
      do i = 1,3
          parn_tran(i,1) = 0
          parn_tran(i,2) = 0
      end do

*                                 ! Scale
      parn_scale(1) = 0
      parn_scale(2) = 0
 
*                                 ! RA origin
      parn_rao = 0
 
*                                 ! Nutation angles
      do i = 1,8
          parn_nut_ang(i) = 0
      end do

      do i = 1,4
          parn_eor_ut1(i) = 0
      end do
      do i = 1,6
          parn_eor_xy(i)  = 0
      end do
      do i = 1, max_glb_sites
          do j = 1,12
              parn_eor_etd(j,i) = 0
          end do
      end do
 
*                                 ! nutation coefficients
      do i = 1, max_nut_coeff
          parn_nut_coeff(1,i) = 0
          parn_nut_coeff(2,i) = 0
      end do
 
*                                 ! Extended earth tide coefficients
      do i = 1, max_etd_sites
          do j = 1, max_etd_coeff
              parn_etd_coeff(1,j,i) = 0
              parn_etd_coeff(2,j,i) = 0
          end do
      end do

      do i = 1, max_ut1_coeff
         parn_ut1_coeff(1,i) = 0
         parn_ut1_coeff(2,i) = 0
      end do

      do i = 1, max_xy_coeff
         parn_xy_coeff(1,i) = 0
         parn_xy_coeff(2,i) = 0
      end do
 
      parn_gamma = 0
 
****  Thats all
      return
      end
 
 
