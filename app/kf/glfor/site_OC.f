CTITLE SITE_OC
 
      subroutine site_OC ( comp, site, sol_obs )

      implicit none  
 
*     Routine to compute site o minus c.  Initially only accounts
*     for position and velocity.  Later should include any affects
*     of changes in PM/UT1 if these values are not estimated.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

      include '../includes/glb_hdr_def.h'
 
*   comp        - Componeent (1=X,2=Y,3=Z)
*   site        - global site number
 
      integer*4 comp, site
 
*   sol_obs     - Parameter to be updated
 
      real*8 sol_obs

* LOCAL PARAMETERS

*   deriv(3) - partials of comp with respect to pole pos and ut1
*   dpmu(3)  - Changes in the pole position estimates
*   dsol(3)  - Changes in position due to non-secular terms

      real*8 deriv(3), dpmu(3), dsol(3)
*   reppmu(3) - Reported PMU error (one set, value needs to change
*               for next error).
      real*8 reppmu(3)

*   kbit     - Test bits:
      logical kbit

      save reppmu
 
***** Remove the C part of the estimates
      if( comp.le.3 ) then
          sol_obs = sol_obs - apr_val_site(comp,1,site) -
     .              apr_val_site(comp,2,site)*(gepoch_expt-
     .                                    site_epoch(site))/365.25d0

          call eval_nonsec(site, gepoch_expt, num_nonsec,
     .         param_nonsec, apr_val_nonsec, dsol, tran_est)
          sol_obs = sol_obs - dsol(comp)

*         See if we need to update polar motion UT1.  If pole components
*         have not been estimated, then force rotation.
          call pmu_part(site, deriv, apr_val_site(1,1,site),
     .                  comp, gut1_apr )
	  dpmu(1) = 0.d0
	  dpmu(2) = 0.d0
	  dpmu(3) = 0.d0
	  if( cwob_apr(1).ne.0 .and. .not.kbit(pmu_est,1)) 
     .                dpmu(1) = gwob_apr(1) - cwob_apr(1)
          if( cwob_apr(2).ne.0 .and. .not.kbit(pmu_est,2) ) 
     .                dpmu(2) = gwob_apr(2) - cwob_apr(2)
          if( cut1_apr(1).ne.0 .and. .not.kbit(pmu_est,3) ) 
     .                dpmu(3) = gut1_apr(1) - cut1_apr(1)
          if( abs(dpmu(3)).gt.10 ) then
              if( dpmu(3).ne.reppmu(3))
     .        write(*,120) pmu_est, cwob_apr(1:2), gwob_apr(1:2),
     .                      cut1_apr(1), gut1_apr(1)
 120          format('* LARGE dUT1 : pmu_est ',o3,' CWOB ',2F10.3,
     .               ' GWOB ',2F10.3,' CUT1 ',F15.3, ' GUT1 ',F15.3)
              reppmu(3) = dpmu(3)
              dpmu(3) = 0    ! Gross error so don't apply
          endif
          sol_obs = sol_obs - dpmu(1)*deriv(1) - dpmu(2)*deriv(2)
     .                      - dpmu(3)*deriv(3) 
      else
          sol_obs = sol_obs - apr_val_site(comp-3,2,site)
      end if
 
***** Thats all
      return
      end
 
