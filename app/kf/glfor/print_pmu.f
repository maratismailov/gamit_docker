CTITLE PRINT_PMU
 
      subroutine print_pmu( unit , cov_parm, sol_parm, ngp)

      implicit none 
 
*     This routine will print the running estimate of the
*     polat motion/UT1 adjustments. (if they are being estimated)
*
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
* PASSED VARIABLES
 
*   unit        - Output unit
*   ngp         - Number of global parameters
 
      integer*4 unit, ngp
 
*   cov_parm(ngp,ngp)   - Covariance matrix
*   sol_parm(ngp)       - Solution vector
 
      real*8 cov_parm(ngp,ngp), sol_parm(ngp)
 
* LOCAL VARIABLES
 
*   i,j         - Loop counters
*   np          - Parameter number (temp)
*   date(5)     - Date of experiment
 
      integer*4 i, np, date(5)
 
*   pmu(2,3)    - Value and +- on current PMU Values.
*   sectag      - Seconds tag on determination
 
      real*8 pmu(2,3), sectag
 
*   pmu_estimated   - True if some component of PMU hase been
*               - estimated
 
      logical pmu_estimated
 
****  Loop over the pmu parameters and collect values
      pmu_estimated = .false.
      np = parn_wob(1)
      if( np.gt.0 ) then
          pmu_estimated = .true.
          pmu(1,1) = sol_parm(np)
          pmu(2,1) = sqrt(cov_parm(np,np))
      end if
      np = parn_wob(2)
      if( np.gt.0 ) then
          pmu_estimated = .true.
          pmu(1,2) = sol_parm(np)
          pmu(2,2) = sqrt(cov_parm(np,np))
      end if
      np = parn_ut1(1)
      if( np.gt.0 ) then
          pmu_estimated = .true.
          pmu(1,3) = sol_parm(np)
          pmu(2,3) = sqrt(cov_parm(np,np))
      end if
 
****  Now see if we should output
      if( pmu_estimated ) then
          call jd_to_ymdhms( gepoch_expt, date, sectag)
          write(unit,100) (date(i), i=1,3), pmu
 100      format('Orient_adj. (mas) ',I4,2i3, ' dX/Y +- ',
     .           2(F7.2,1x,F7.2,1x),' dUT1 ',
     .            (F7.2,1x,f7.2))
      end if
 
****  Thats all
      return
      end
 
 
