CTITLE APPLY_COND
 
      subroutine apply_cond(parn, ndim, pmu_parts, num, cov_parm,
     .            numc, cond_part, origin_gain, cat, used )
 
      implicit none 

*     Routine to apply conditions so that the net rotation and
*     translation of the coordinate system will be be zero.
 
*   ndim                - Second dimension for parn array (
*                       - 2 is postion and velocity, 1 for
*                       - postion only (solvk)
*   num                 - Number of sites in the parn array
*   numc                - NUmber of rows and columns in the
*                       - covariance matrix.
*   parn(3,ndim,num)    - Parameter numbers for the sites
*   used(*)             - Bit mapped array saying which sites
*                       - to use.
 
      integer*4 ndim, num, numc, parn(3,ndim,num), used(*)
 
*   cov_parm(numc,numc) - Covariance matrix.
*   pmu_parts(3,3,num)  - Polar motion/UT1 partials
*   cond_part(6,num)    - Condition partials.
*   origin_gain(numc,6) - Computed Kalman Gain for applying the
*                       - conditions.
*   cat(numc,6)         - Temporary storage of:
*                       -  m   T
*                       - C   A
*                       -  m+1
 
      real*8 cov_parm(numc,numc), pmu_parts(3,3,num), cond_part(6,num),
     .    origin_gain(numc,6), cat(numc,6)
 
* LOCAL VARIABLES
 
*   i,j,k,l         - Loop counters
*   ir,ic           - Row and col counters
*   iel             - Generic element in matrix
*   lel             - last elment used in cond_part
*   ipivot(6)       - Used by invert_vis
 
      integer*4 i,j,k,l, ir,ic, iel, lel, ipivot(6)
 
*   kbit            - Check to see if bit is set.
 
      logical kbit
 
*   ptp(6,6)            - Normal equations for the pmu_part inversion
*   scale(6)            - Used by invert_vis
*   tmp                 - Temporary summation of covariance decrament.
 
 
      real*8 ptp(6,6), scale(6), tmp 

      lel = 0
 
****  Start, first computed the normal eqautions for the pmu_part
*     i.e., we want to compute the generalize inverse of pmu_part
      do i = 1,6
          do j = 1,6
              ptp(i,j) = 0.d0
          end do
      end do
 
*     Now start forming the normal equations
      do i = 1, num
*                                         ! We are using this site
          if( kbit(used,i) ) then
 
*             Translations first
              ptp(1,1) = ptp(1,1) + 1
              ptp(2,2) = ptp(2,2) + 1
              ptp(3,3) = ptp(3,3) + 1
 
*             Translation/rotation cross terms
              do j = 1,3
                  do k = 1,3
                      ptp(j,k+3) = ptp(j,k+3) +
     .                    pmu_parts(j,k,i)
                      ptp(k+3,j) = ptp(j,k+3)
                  end do
              end do
 
*             Rotatations
              do j = 1,3
                  do k = 1,3
                      do l = 1,3
                          ptp(j+3,k+3) = ptp(j+3,k+3) +
     .                        pmu_parts(j,l,i)*pmu_parts(k,l,i)
                      end do
                  end do
              end do
          end if
      end do
 
****  Now invert the normal equations
      write(*,150) 'PTP Before inversion', ptp
 150  format(/,a,6(/,6F15.6))
      call invert_vis(ptp, origin_gain,scale, ipivot, 6,6,0)
      write(*,150) 'PTP After inversion', ptp
 
****  Now form the condition partials
      iel = parn(3,1,num)

      if( iel.eq.0 ) RETURN

      do i = 1,iel
        do j = 1,6
              cond_part(j,i) = 0.d0
          end do
      end do
 
***   Now fill in the non-zero values
      do i = 1, num
          if( kbit(used,i) ) then
              do j = 1,6
                  do k = 1,3
                      iel = parn(k,1,i)
                      if( iel.gt.0 ) then
****                      Save the last element used in cond_part
                          lel = iel
*                         Do the translation part
                          cond_part(j,iel) = ptp(j,k)
 
*                         Now do the rotation
                          do l = 1,3
                              cond_part(j,iel) = cond_part(j,iel)+
     .                                    ptp(j,l+3)*pmu_parts(l,k,i)
                          end do
                      end if
                  end do
              end do
          end if
      end do
 
****  Now we have partials of conditions.  Now set the variances
*     of these linear combinations to zero.
*     First form Cov_parm*cond_part(R)
      do ir = 1,numc
          do ic = 1,6
              cat(ir,ic) = 0.d0
 
              do i = 1,lel
                  cat(ir,ic) = cat(ir,ic) +
     .                cov_parm(ir,i)*cond_part(ic,i)
              end do
          end do
      end do
 
*                m   T
****  Now form AC   A
*                m+1
      do ir = 1,6
          do ic = 1,6
              ptp(ir,ic) = 0.d0
              do i = 1,lel
                  ptp(ir,ic) = ptp(ir,ic) +
     .                cond_part(ir,i)*cat(i,ic)
              end do
          end do
      end do
 
****  Now invert
      write(*,150) 'Data Covariance before inversion', ptp
      call invert_vis(ptp, origin_gain,scale, ipivot, 6,6,0)
      write(*,150) 'Data Covariance after  inversion', ptp
 
****  Now finishup the origin gain
      do ir = 1, numc
          do ic = 1,6
              origin_gain(ir,ic) = 0
              do i = 1, 6
                  origin_gain(ir,ic) = origin_gain(ir,ic) +
     .                cat(ir,i)*ptp(i,ic)
              end do
          end do
      end do
 
****  Now finally update covarirance matrix
      do ir = 1,numc
          do ic = 1,numc
              tmp= 0.d0
              do i = 1,6
                  tmp = tmp + origin_gain(ir,i)*cat(ic,i)
              end do
              cov_parm(ir,ic) = cov_parm(ir,ic) - tmp
              if( ir.eq.ic .and. cov_parm(ir,ic).lt.0.d0 ) then
                  write(*,250) ir, cov_parm(ir,ic), tmp
 250              format('Negative diagonal Row ',i4,' After&before',
     .                  2d20.8)
              end if
          end do
      end do
 
****  Thats all.
      return
      end
 
 
 
 
