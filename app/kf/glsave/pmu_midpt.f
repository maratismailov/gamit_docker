CTITLE PMU_MIDPT

      subroutine pmu_midpt

      implicit none 

*     Routine to convert single day polar motion/UT1 into single 
*     entry multiday polar motion/UT1 when the output epoch of the
*     soltion does not match the PMU epoch (happens when -M is set
*     in glsave or midp option is selected in org or prt_opts.

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'

* LOCAL VARIABLES
      integer*4 i, j
      logical kbit

     
***** Polar motion/UT1 values
* MOD TAH 070823: If we are outputing at mid-point convert the
*     dailty polar motion/UT1 to multiday with the correct 
*     epoch (this needs to be the last value)
      print *,'PMU ',gepoch_expt,gepoch_out, num_mul_pmu 
      stop
      if( abs(gepoch_expt-gepoch_out).gt.1.d-3 .and. 
     .    num_mul_pmu.eq. 0) then
          do i = 1,2
             parn_mul_pmu(1,i,1) = parn_wob(i)
             apr_val_mul_pmu(1,i,1) = apr_val_wob(i) + gwob_apr(i)

             parn_mul_pmu(2,i,1) = parn_wob(i+2)
             apr_val_mul_pmu(2,i,1) = apr_val_wob(i+2) + gwob_apr(i+2)
          end do
          parn_mul_pmu(1,3,1) = parn_ut1(1)
          apr_val_mul_pmu(1,3,1) = apr_val_ut1(1) + gut1_apr(1)
          parn_mul_pmu(2,3,1) = parn_ut1(2)
          apr_val_mul_pmu(2,3,1) = apr_val_ut1(2) + gut1_apr(2)
*         Set the number
          num_mul_pmu = 1
          gmul_pmu_ep(1) = gepoch_expt  ! Must be referenced to original time
          do i = 1,4
             parn_wob(i) = 0  ! Stops regular output
          end do 
          do i = 1,2
             parn_ut1(i) = 0
          end do
          do j = 1,gnum_sites
             if( kbit(guse_site,j) ) call sbit(mul_pmu_site(1,1),j)
          enddo
          do j = 1, gnum_svs
             call sbit(mul_pmu_svs(1),j)
          enddo

      end if

****  Thats all
      return
      end


