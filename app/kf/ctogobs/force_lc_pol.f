CTITLE FORCE_LC_POLY
 
      subroutine force_lc_poly(pod, apl,bpl, apr, bpr, 
     .                         dobs, sobs)

      implicit none
 
*     This routine will estimate the offset between observables
*     while forcing the cofficients of the polynomial to be 
*     same.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES

*  pod  - Order of the polynomial being estimated for 
*         dimensioning (2=linear)

      integer*4 pod

*  apl(max_patch_poly,max_patch_poly), bpl(max_patch_poly) -- 
*          Polynomial estimates from left segment.
*  apr(max_patch_poly,max_patch_poly), bpr(max_patch_poly) -- 
*          Polynomial estimates from right segment.
*  dobs, sobs  - Estimates and sigma of the jump accross the
*          bias flag
      
      real*8 apl(max_patch_poly,max_patch_poly), bpl(max_patch_poly),
     .       apr(max_patch_poly,max_patch_poly), bpr(max_patch_poly),
     .       dobs, sobs

* LOCAL VARIABLES

*  i,j,k    - Loop counters
*  ipivot   - Pivot elements for inversion.

      integer*4 i,j,k, ipivot(max_patch_poly+1)

*  apc(max_patch_poly+1,max_patch_poly+1), bpc(max_patch_poly+1) -- 
*          Polynomial estimates from combined data.  The zero order
*          element is the difference Pleft - Pright
*  vpl(max_patch_poly), vpr(max_patch_poly) - Estimates of polynomial
*          coefficients.  (Will be updated in invert_vis)
*  chi  -- chi**2 once polynomials are forced to be the same
*  prechi  - Prefit chi**2
*  dchi    - Change in chi**2 for estimate
*  scale(max_patch_poly)  - scale factor
*  apart(max_patch_poly)  - Partial derivatives
*  avi(max_patch_poly+1,max_patch_poly) - Partials matrix by inverse of
*          covariance matrix
      
      real*8 apc(max_patch_poly+1,max_patch_poly+1), 
     .       bpc(max_patch_poly+1),
     .       vpl(max_patch_poly), vpr(max_patch_poly),
     .       vpc(max_patch_poly+1),
     .       chi, pre_chi, dchi, scale(max_patch_poly+1), 
     .       apart(max_patch_poly,max_patch_poly+1),
     .       avi(max_patch_poly+1,max_patch_poly),
     .       debug(2)


****  Start, initialize the estimator
      debug(1) = bpl(1) - bpr(1)
      debug(2) = sqrt(apl(1,1)+apr(1,1))

      do i = 1, pod+1
         bpc(i) = 0.d0
         do j = 1, pod+1
            apc(i,j) = 0.d0
         end do
      end do

      pre_chi = 0.d0

****  Now form partials for Left segment
      do i = 1, pod
         do j = 1, pod+1 
            apart(i,j) = 0.d0
         end do
      end do
      apart(1,1) = 1.d0
      do i = 1, pod-1
         do j = 1, pod-1
            apart(i+1,j+2) = 1.d0
         end do
      end do

***** Copy solution over before inverting the covariance matrix
      do i = 1, pod
         vpl(i) = bpl(i)
      end do

      call invert_vis(apl,vpl, scale, ipivot, pod, max_patch_poly,1)

****  Now compute prefit chisquared.
      do i = 1, pod
         pre_chi = pre_chi + bpl(i)*vpl(i)
      end do

****  Now form the increment to the normal equations. First do ATV-1
      do i = 1, pod+1
         do j = 1, pod
            avi(i,j) = 0.d0
            do k = 1, pod 
               avi(i,j) = avi(i,j) + apart(k,i)*apl(k,j)
            end do
         end do
      end do

****  Now complete the computation by post-multipling by AT and bpl
      do i = 1,pod+1
         do j = 1, pod+1
            do k = 1, pod
               apc(i,j) = apc(i,j) + avi(i,k)*apart(k,j)
            end do
         end do
      end do

*     Increment soltion vector
      do i = 1, pod+1
         do j = 1, pod
            bpc(i) = bpc(i) + avi(i,j)*bpl(j)
         end do
      end do

****************************************
****  Now form partials for Right segment
      do i = 1, pod
         do j = 1, pod+1 
            apart(i,j) = 0.d0
         end do
      end do
      apart(1,2) = 1.d0
      do i = 1, pod-1
         do j = 1, pod-1
            apart(i+1,j+2) = 1.d0
         end do
      end do

***** Copy solution over before inverting the covariance matrix
      do i = 1, pod
         vpr(i) = bpr(i)
      end do

      call invert_vis(apr,vpr, scale, ipivot, pod, max_patch_poly,1)

****  Now compute prefit chisquared.
      do i = 1, pod
         pre_chi = pre_chi + bpr(i)*vpr(i)
      end do

****  Now form the increment to the normal equations. First do ATV-1
      do i = 1, pod+1
         do j = 1, pod
            avi(i,j) = 0.d0
            do k = 1, pod 
               avi(i,j) = avi(i,j) + apart(k,i)*apr(k,j)
            end do
         end do
      end do

****  Now complete the computation by post-multipling by AT and bpr
      do i = 1,pod+1
         do j = 1, pod+1
            do k = 1, pod
               apc(i,j) = apc(i,j) + avi(i,k)*apart(k,j)
            end do
         end do
      end do

*     Increment soltion vector
      do i = 1, pod+1
         do j = 1, pod
            bpc(i) = bpc(i) + avi(i,j)*bpr(j)
         end do
      end do

**** *******************************
*     Now complete the solution
      do i = 1, pod+1
         vpc(i) = bpc(i)
      end do

      call invert_vis(apc,vpc, scale, ipivot, pod+1, max_patch_poly+1,1)

****  Compute the new chi squared for the solution
      dchi = 0
      do i = 1, pod+1
         dchi = dchi + vpc(i)*bpc(i)
      end do

      chi = (pre_chi - dchi)/(2*pod)
      if( chi.lt.1.d0 ) chi = 1.d0

*     Now get the estimate of the jump and its sigma
      dobs = vpc(1) - vpc(2)
      sobs = sqrt((apc(1,1)+apc(2,2)-2.d0*apc(1,2))*chi)

****  Thats all
      return
      end
