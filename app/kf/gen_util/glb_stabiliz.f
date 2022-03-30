CTITLE GLB_STABILIZE
 
      subroutine glb_stabilize( parn, num, ndimp, cov_parm, ndimc,
     .                          nused, sol_parm, origin_gain,
     .                          used, translation)
 
      implicit none 

 
*     Routine to stabilize the translation orgin of either the
*     station positions or of the RA origin.
*
*     This subroutine is a mini kalman filter with one observation,
*     namely, the sum of the parameters given in PARN should be
*     zero.
*     The used array indicates (by bit being set) which sites are to
*     used in the origin.  Translation gives the value that the
*     origin should be fixed to.
 
*   i,j,k   - Loop counters
*   iel,jel - Matrix row numbers
*   ndimc   - Dimension of the covariance matrix
*   ndimp   - the dimension of the first index in parn
*   num     - Number of parameters to be stabilized
*   nused   - Number of rows used in cov_parm
*   parn(ndimp,1)   - Parameter numbers of the parameters
*   used(1) - Bit mapped array giving which sites or sources
*           - should be used in the origin
 
      integer*4 i,j,iel, ndimc, ndimp, num, nused,
     .    parn(ndimp,1), used(1)
 
*   cov_parm(ndimc,ndimc)   - Covariance matrix to be stabilized
*   origin_gain(1)          - Kalman gain vector for the
*                           - origin of the system
*   sol_parm(ndimc)         - solution vector
*   sol_total               - sum of the sol_parm vector
*   total                   - intermediate result, sum of the
*                           - elements in a column
*   translation             - Value to which the origin should
*                           - be fixed (i.e., average correction
*                           - to the values)
      real*8 cov_parm(ndimc,ndimc), origin_gain(1), sol_parm(ndimc),
     .    sol_total, total, translation
 
*   kbit                    - Bit checking function
 
      logical kbit
 
 
 
***** First check to see if all parameters are estimated
 
      do i = 1, num
          if( parn(1,i).eq.0 .and.
*                                                ! Not all estimated
     .        kbit( used,i )        ) RETURN
      end do
 
***** Form the kalman gain for stabilization
 
      do i = 1, nused
 
          origin_gain(i) = 0.d0
          do j = 1, num
*                                           ! Include
              if( kbit( used,j ) ) then
                  iel = parn(1,j)
                  origin_gain(i) = origin_gain(i) + cov_parm(iel,i)
              end if
          end do
      end do
 
****  Now sum the orgin gain
      total     = 0
      sol_total = 0
      do j = 1, num
*                                            ! Include
          if( kbit( used,j ) ) then
              iel = parn(1,j)
              total     = total     + origin_gain(iel)
              sol_total = sol_total + (sol_parm(iel)+translation)
          end if
      end do
 
***** Now get the adjustment to the parmeters
 
      do i = 1,nused
          sol_parm(i) = sol_parm(i) - origin_gain(i)*sol_total/total
      end do
 
***** Now adjust the covariance matrix
      do i = 1, nused
          do j = i, nused
              cov_parm(j,i) = cov_parm(j,i) -
     .                        origin_gain(j)*origin_gain(i)/total
              cov_parm(i,j) = cov_parm(j,i)
          end do
      end do
 
***** Thats all
      return
      end
 
      subroutine stabilize(parn, num, ndimp, cov_parm, ndimc, sol_parm)
 
 
c
c     routine to stabilize the kalman filter solution if all
c     clock and/or site parameters have been estimated.  This
c     subroutine removes the singularity in the solution for
c     these parameters where we need to assume that some parameters
c     are 'fixed' i.e., one clock and/or one site.
c     Note: the apriori constrains placed on these parameters ensures
c     that the adjusts to them have zero mean.  Therefore we only
c     need to manipulate the covariance matrix to get the correct
c     sigmas.
c
c Variables
c ---------
c parn -- the parameter numbers for the paramaters being stabilized
c num  -- the number of parameters to be stabilized
c ndimp -- the dimension of the first index in parn
c cov_parm -- the covariance matrix from the solution.
c ndimc -- the dimension of cov_parm
c sol_parm -- the parameter adjustment estimates
c
      integer*4 ndimp, ndimc, parn(ndimp,1), num
 
c
      real*8 cov_parm(ndimc,1), sol_parm(1)
 
c
c
c Local variables
c ---------------
c iel -- index for matrix access
c jel -- another index for matrix access
c ind_cov -- a column index for cov_parm
c total -- intermediate computation value
c sol_total -- intermediare computation value
c
 
      integer*4 iel
 
*   i,j,k   - Loop counters
      integer*4 i,j
 
c
      real*8 total, sol_total
 
c
c Scratch common
c --------------
c ident -- identification for use of common
c origin_gain -- the kalman gain vector for the aplliction of the
c     zero mean shift.
c
      integer*4 ident, dumout
 
c
      real*8 origin_gain(1)
 
c
      common ident, dumout, origin_gain
 
c
c
c.... First check if we should stabilize .i.e, are all parameters on
      do i = 1, num
*                                       ! at least one value not estimated
         if( parn(1,i).eq.0 ) return
c                                         so donot stabilize
      end do
c
c.... OK all are turned on
c
c.... Form the Kalmam gain for the stabilization
      iel = parn(1,1)
      do j = 1, ndimc
         call dwsum(origin_gain(j),cov_parm(iel,j),ndimp, num)
      end do
c
c.... Now sum orgin_gain
      call dwsum(total,origin_gain(iel),ndimp, num)
c
c.... Now adjust the parameter estimates
      call dwsum(sol_total, sol_parm(iel),ndimp, num)
      do i = 1, ndimc
         sol_parm(i) = sol_parm(i) - origin_gain(i)*sol_total/total
      end do
c
c.... Now adjust the covariance matrix
      do i = 1, ndimc
*                        ! loop over all elements for these parameters.
         do j = 1, ndimc
            cov_parm(j,i) = cov_parm(j,i) - origin_gain(j)*
     .         origin_gain(i)/total
         end do
      end do
c
c.... Thats all
      return
      end
 
