CTITLE GET_MAX_GLB_DERIV
 
      subroutine GET_MAX_GLB_DERIV

      implicit none
 
 
*     Routine to compute maximum number of derivatie per observation
*     we may have
 
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
      include '../includes/globk_markov.h'
 
*   i,j     - Loop counters
*   tmax    - Temporary max (while we scan over all coefficeints)
 
 
      integer*4 i,j,k, tmax
      integer*4 log_num  ! number of earthquakes that could effect one
                         ! site
      integer*4 log_max  ! Mamximum value of log_num
 
****  Start with default and starting adding up
 
      max_glb_deriv = 5
 
*     Add up nutation series partials (since these will exceed the
*     value above
 
      do i = 1, max_nut_coeff
          do j = 1,2
              if( parn_nut_coeff(j,i).ne.0 ) then
                  max_glb_deriv = max_glb_deriv + 1
              end if
          end do
      end do

****  Check ut1 coeff partials

      tmax = 5
      do i = 1, max_ut1_coeff
         do j = 1,2
            if( parn_ut1_coeff(j,i).ne.0 ) then
               tmax = tmax + 1
            end if
         end do
      end do

      max_glb_deriv = max(max_glb_deriv, tmax )

****  Check xy coeff partials

      tmax = 5
      do i = 1, max_xy_coeff
         do j = 1,2
            if( parn_xy_coeff(j,i).ne.0 ) then
               tmax = tmax + 1
            end if
         end do
      end do

      max_glb_deriv = max(max_glb_deriv, tmax )

****  Check etd coeff partials

      do k = 1, min(max_etd_sites,gnum_sites)
          tmax = 5
          do i = 1, max_etd_coeff
             do j = 1,2
                if( parn_etd_coeff(j,i,k).ne.0 ) then
                   tmax = tmax + 1
                end if
             end do
          end do

          max_glb_deriv = max(max_glb_deriv, tmax )
      end do

****  MOD TAH 090713: Add in apr_rot contributions (there are 2 per station 
*     coordinate, coordinate does not depend of rotation about its axis).

      tmax = 1  ! Allow 1 for coordinate 
      if( parn_wob(1).gt.0 .or. parn_ut1(1).gt. 0 ) tmax = 3  ! EOP adds 2

      if( parn_rot(k,1).gt.0 .or. parn_rot(k,2).gt.0 ) tmax = tmax+2

*     Now do translation 
      if( parn_tran(1,1).gt.0 ) tmax = tmax + 1
      if( parn_scale(1).gt.0  ) tmax = tmax + 1



c      max_glb_deriv = max_glb_deriv + tmax

****  MOD TAH 090708: Add in contributions from log coefficients: 
*     
      log_max = 0
      log_num = 0
      do k = 1, gnum_sites
          if( k.gt.1 ) then
*             Here we try to determine number of earthquakes that could 
*             effect a single site. 
* MOD TAH 160324: Due to inter mixxing of site names, go back to until
*             4-character site code changes
              if( parn_log(1,k).gt.0 ) then
                  log_num = 1
              else
                  log_num = 0
              endif
              i = 1
              do while ( gsite_names(k-i)(1:4).eq.gsite_names(k)(1:4)
     .                   .and. k-i.gt.1 .and. parn_log(1,k).gt.0 )
                  if( parn_log(1,k-i).gt.0 ) log_num = log_num + 1
                  log_max = max(log_max,log_num)
                  i = i + 1
              end do
          end if
          log_max = max(log_max,log_num)
      end do

      max_glb_deriv = max(max_glb_deriv,tmax) + log_max*3
      write(*,120) max_glb_deriv, tmax, log_max 
 120  format(' Maximum derivatives found ',I3,
     .       ' with ',I3,' coordinate and ',i3,' EQ log terms')


***** That should do us
      return
      end
 
 
