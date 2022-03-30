CTITLE SOURCE_PARTIAL
 
      subroutine source_partial( comp, pn, part_pnt, a_part, source)

      implicit none 
 
*     Routine to compute the partial derivatives with respect to the
*     parameters in the global solution.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
 
*   comp    - Component of position (1=ra, 2=dec)
*   i,j     - loop counter
*   pn      - Local parameter number
*   part_pnt(2,1)   - pointers to the partial derivatives
 
*   source  - Global source number
 
      integer*4 comp, i,j, pn, part_pnt(2,1), source
 
*   a_part(1)   - Array for partial derivatives
*   coeff_deriv(2,max_nut_coeff)  - Nutation series coefficient
*               - partials
*   amp_deriv(2,4,max_nut_coeff)  - Nutation ampltides partials
*   deriv(2)    - Derivatives for nutation angles
*   epoch       - Epoch of experiment for computing series partials
*   fcn_period  - FCN period in main memory
 
      real*8 a_part(1), coeff_deriv(2,max_nut_coeff), deriv(2), epoch,
     .    fcn_period
     
C, amp_deriv(2,4,max_nut_coeff) - Not used
 
*   kbit        - Bit checking function
 
      logical kbit
 
 
 
***** Start adding in partials, global position first
 
*     See if source to be used
      if( .not.kbit(guse_source,source) ) RETURN
 
      call add_partial(indx_pnt(pn), part_pnt,
     .                 parn_source(comp,1,source), a_part, 1.d0 )
 
*     See if we should add translation (only for RA)
 
      if( comp.eq.1 ) then
          call add_partial(indx_pnt(pn), part_pnt, parn_RAO,
     .                     a_part, 1.d0 )
      end if
 
*     Now check for nutation angle coefficients ( Only needed if not
*     estimated)
 
*                                 ! Add nutation angle partials
      if( .not.nut_ang_est ) then
 
          call nut_ang_part( source, deriv, apr_val_source(1,1,source),
     .                       comp, gepoch_expt )
 
*                             ! Nutation DEPS, and DPSI
          do i = 1, 2
              call add_partial(indx_pnt(pn), part_pnt,
     .                         parn_nut_ang(i), a_part, deriv(i) )
          end do
 
*         Now check on the nutation series coefficient parts
 
*                                 ! Gepoch_epxt is EMA
          epoch = gepoch_expt
          fcn_period = nut_period
          call nut_coeff_part(deriv, epoch, fcn_period, coeff_deriv)
c         call gamp_parts( epoch, fcn_period, amp_deriv )
 
          do i = 1, max_nut_coeff
*                             ! In and Out of phase
              do j = 1,2
                  call add_partial( indx_pnt(pn), part_pnt,
     .                 parn_nut_coeff(j,i), a_part, coeff_deriv(j,i))
              end do
c             do j = 1,4
c                 coeff_deriv(j) = deriv(1)*amp_deriv(1,j,i) +
c    .                             deriv(2)*amp_deriv(2,j,i)
c             end do
c             do j = 1,4
c                 call add_partial( indx_pnt(pn), part_pnt,
c    .                 parn_nut_coeff(j,i), a_part, 
c    .                 coeff_deriv(j))
c             end do
          end do
 
      end if
 
***** Thats all, later should add series partials
      return
      end
 
