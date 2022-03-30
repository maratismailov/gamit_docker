       subroutine fit_igs_clk( ins, params_cse, par_flag_cse,
     .                        L1_cyc_cse, L2_cyc_cse, 
     .                        ctol_cse,  type)


*     This subroutine will fit linear models to the clock
*     estimates or apply the fits to the range clocks to
*     phase clock values.  The notion here is to improve
*     the overall allignment of phase and range clocks by
*     alligning the offsets from the range clocks with the
*     phase clocks. 

* INCLUDES

      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include '../includes/mfile_def.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
*   par_flag_cse(num_param,num_ep)         - parameter estimate flag
*   ins  -- Summary file unit number
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number

      integer*4 par_flag_cse(num_param,num_ep), ins,
     .          ctol_cse(num_chan, num_cfiles, num_ep)

*   params_cse(num_param,num_ep)    - Clock parameter estimates (L1 cycles)
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement


      real*8 params_cse(num_param,num_ep),
     .    L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep)

* type  -- P or R for range or phase clocks.  When R clocks are used
*     polynomials are fitted, when P clocks are passed, polynomials
*     are fit and the offsets between the ranges and phase clocks
*     applied to the phase clocks and the number of cycles (the latter
*     should ensure that the total clock+phase remains unchanged)

      character*1 type

* LOCAL VARIABLES
* i,j,k  -- Loop counters
* ep     -- Epoch counter
* num_stat -- Number of values in clock estimate
* ipivot(2) -- Pivot elements

      integer*4 i,j,k, ep, num_stat, ipivot(2) 

* norm_clk(2,2) -- Clock polynomial normal equations
* sol_clk(2), bvec_clk(2) -- Solution vectors for clock solution
* dtx -- Epoch offset to get uncorrelated time
* dt  -- Time offset (sec)
* dclk -- Average clock offset
* dphs(2)  -- Phase offsets (not used in initial version)
* sum_stat -- RMS statistics summation
* rms      -- RMS scatter 
* scale(2) -- Scaling for inversion 

      real*8 norm_clk(2,2), sol_clk(2), bvec_clk(2),
     .       dtx, dt, dclk, dphs(2), sum_stat, rms, scale(2)

* dclk_save(max_cfiles+max_gprn) -- Save clock values for output
      real*4 dclk_save(max_cfiles+max_gprn)

* name  -- Name of site or satellite
* name_save(max_cfiles+max_gprn) -- Saved names for output

      character*4 name, name_save(max_cfiles+max_gprn)

* Variables used in fitting diffeences (results saved in save_drms(i)

      integer*4 num_var
      real*8 sum_var, sum_dmean, dmean, rms_dmean, prev_clk
              
****  Fit the linear polynomials to all the clocks (and compute the
*     RMS about the fits).
      write(ins,110) type
 110  format('CLOCK FITS FOR ',a,' Data type',/,
     .       'C T Name      #       Offset (ns)     dRate (D-14)',
     .       '     RMS (ns)    Rate (D-14)')
      do i = 1, num_param

         sum_stat = 0
         num_stat = 0
         do j = 1,2
            sol_clk(j) = 0.d0
            do k = 1, 2
               norm_clk(j,k) = 0.d0
            end do
         end do

*        OK, now loop over all of the data
         do ep = 1, num_ep
            if( par_flag_cse(i,ep).eq.0 ) then
*               Estimate is good, accumulate into normal equations
                num_stat = num_stat + 1
                dt = (ep-1)*sampling + 0.001d-3
                sum_stat = sum_stat + params_cse(i,ep)**2
                do j = 1, 2
                   sol_clk(j) = sol_clk(j) + 
     .                          params_cse(i,ep)*dt**(j-1)
                   do k = 1, 2
                      norm_clk(j,k) = norm_clk(j,k) + 
     .                                dt**(j-1)*dt**(k-1)
                   end do
                end do
C               if( i.le.num_cfiles ) then
C                  write(*,120) ep, cf_codes(i),
C    .                       params_cse(i,ep)/fClk*vel_light 
C120               format('APClk ',i5,1x,a4,1x,F15.2, ' m')
C               endif
            end if
         end do

****     OK: Now solve the system of equations
         do j = 1, 2
            bvec_clk(j) = sol_clk(j)
         end do

*        Generate name of parameter
         if( i.le.num_cfiles ) then
            name = cf_codes(i)
         else
* MOD TAH 200628: Added GNSS code.
            write(name,210) sv_gnss, prn_list(i-num_cfiles)
 210        format(a1,I2.2)
         end if
         name_save(i) = name
         if( num_stat.gt.2 ) then
             call invert_vis(norm_clk, sol_clk, scale, ipivot, 2,2,1)

*            Now compute postfit RMS
             sum_stat = sum_stat - sol_clk(1)*bvec_clk(1) -
     .                             sol_clk(2)*bvec_clk(2)
             rms = sqrt(abs(sum_stat)/num_stat)

             write(ins,220) type, name, num_stat, (sol_clk(1)/fClk)*1.d9
     .                  , (sol_clk(2)/fClk)*1.d14,
     .                    (rms/fClk)*1.d9, 
     .                    (apr_clk_poly(2,i)+sol_clk(2)/fClk)*1.d14
 220         format('C ',a1,1x,a4,i6,F20.3,F15.3, F10.2,3x,F15.3)
         else
             write(ins,230) type, name
 230         format('C ',a1,1x,a4,' No Data')
             rms = 999.99d0
         end if

*****    Now: If this is range, save the offset (at the uncorrelated
*        epoch) and rate.  If it is phase, compute the offset needed
*        and apply this value
* MOD TAH 200504: Only adjust clock if RMS to linear fit is small.
         if( type.eq.'R' ) then
             if( num_stat.gt.2 .and. rms/fClk.lt.100.0d-9 ) then 
                 dtx = -norm_clk(1,2)/norm_clk(2,2)
                 sol_clk(1) = sol_clk(1) + dtx*sol_clk(2)
             else
                 dtx = 0.0d0
                 sol_clk(1) = 0
                 sol_clk(2) = 0
             end if
             do j = 1, 2
                save_clk(j,i) = sol_clk(j)
             end do
             save_epc(i) = dtx
             save_rms(i) = rms
         else

*            Now compute and apply offset.
             if( num_stat.gt.2 ) then
                 dtx =  -norm_clk(1,2)/norm_clk(2,2)
                 sol_clk(1) = sol_clk(1) + dtx*sol_clk(2)
             else
                 dtx = 0.0d0
                 sol_clk(1) = 0
                 sol_clk(2) = 0
             end if

*            Compute the value of the offset
             dclk = sol_clk(1) - save_clk(1,i)
             dclk_save(i) = dclk

*            Save the final clock values to be used.
             save_clk(2,i) = sol_clk(2)
             save_epc(i) = dtx
             save_rms(i) = rms

*            Now apply the clk offset and the phase adjustment
*            (For the moment, skip the cycle update until we see
*            how well this works
             do ep = 1, num_ep
                params_cse(i,ep) = params_cse(i,ep) - dclk
             end do
         end if
      end do

*     For the P phase case, output the offsets applied
      if( type.eq.'P' ) then
          write(ins,350) 
 350      format('A P S/S      Dclk (ns)      dL1 (cyc)  dL2 (cyc)') 
          do i = 1, num_param

*            Get the offset 
             name = name_save(i)
             dclk = dclk_save(i)
             if( i.gt. num_cfiles ) then   ! Reference to satellite freq
                dphs(1) = dclk*fL1(i-num_cfiles)/fClk
                dphs(2) = dclk*fL2(i-num_cfiles)/fClk
             else
                dphs(1) = dclk
                dphs(2) = dclk*fL2(1)/fClk
             endif 

             write(ins,360) type, name, (dclk/fClk)*1.d9, 
     .                    dphs(1), dphs(2)
 360         format('A ',a1,1x,a4,1x,F10.3,1x,
     .              2F10.3)
          end do
      end if


***   Now compute the rms of the differences.  This is a better measure
*     of the stability of the clock for older data.
      write(ins,410)
 410  format('C R Site     #      DiffRMS   TotalRMS   ',
     .       '   DiffRMS   TotalRMS ',/,
     .       'C                         (cycles)       ',
     .       '        (ns)')
      do i = 1, num_cfiles
         sum_var = 0
         num_var = 0
         sum_dmean = 0
         prev_clk = 0.d0
*        Find first clock value
         k = 1
         prev_clk =  params_cse(i,1)        
         do while ( par_flag_cse(i,k).ne.0 .and. k.lt.num_ep)
             k = k + 1
         enddo 
      
         do ep = k+1, num_ep
            if( par_flag_cse(i,ep).eq.0 ) then
                num_var = num_var + 1

* MOD TAH 0001718: Remove the average slope while computing the
*               mean
                sum_dmean = sum_dmean + (prev_clk-
     .                               params_cse(i,ep))
                sum_var = sum_var + (prev_clk-
     .                               params_cse(i,ep))**2
                prev_clk = params_cse(i,ep)
c               write(*,998) i,ep, params_cse(i,ep), par_flag_cse(i,ep)
c998            format('CLK: ',i3,i5,1x,F20.2,1x,i6)
           end if
         end do
         if( num_var.gt.2 ) then
            dmean = sum_dmean/num_var
            rms_dmean = sqrt(abs(sum_var/num_var-dmean**2))
            save_drms(i) = rms_dmean

         else
            save_drms(i) = 999.99
         end if
         write(ins,420) type, cf_codes(i), num_var, 
     .                     save_drms(i), save_rms(i),
     .                    (save_drms(i)/fClK)*1.d9,
     .                    (save_rms(i)/fClk)*1.d9
 420     format('C ',a1,1x,a4,i6,1x,2f11.2,1x,2f11.2)

      end do

****  Thats all
      return
      end

CTITLE APP_PREFIT_CLK
 
      subroutine app_prefit_clk(L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, params_cse, par_flag_cse  )
 
      implicit none
 
*     Routine to take updated linear fits to the station clocks
*     and apply them to clock estimates and range and phase
*     residuals to make it look like the polynomials were used
*     in model.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/cfile_def.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param, num_ep)
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse(num_chan, num_cfiles, num_ep)  - L1 range residuals
*                   - cycles at L1
*   L2r_rng_cse(num_chan, num_cfiles, num_ep)  - L2 range residuals
*                   - cycles at L2
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)
 
* LOCAL VARIABLES
      integer*4 i,j   ! Loop counters
      integer*4 lv    ! Satellite list number
      integer*4 ep    ! Epoch number
      logical   kbit  ! Test bit status

      real*8 dt       ! Time offset from epoch 1 (plus small offset
                      ! to dt**n
      real*8 dClk     ! Clock rate offset (cycles at fClk)


****  Update the fitted polynonomial values (save_clk has the estimated
*     adjustments to the apr_clk_poly values so we apply and then apply 
*     the adjustment to the phase and range data.
      do i = 1, num_cfiles
*        Update the apr_clk_poly values for the new estimates (only
*        rate is updated)
         apr_clk_poly(2,i) = apr_clk_poly(2,i) + save_clk(2,i)/fClk  

****     Update the clock_estimates (params_cse) and the phase and
*        range values
         do ep = 1, num_ep 
* MOD TAH 200507: Change to apply to the mid-point epoch of the
*           data (no mean shift).
*           dt = (ep-1)*sampling 
            dt = ep*sampling-save_epc(i) 
            dClk = save_clk(2,i)*dt
            params_cse(i,ep) = params_cse(i,ep) - dClk
*           Now remove from range and phase residuals as well
*           Station clock so apply to all channels
            do j = 1, actual_max_chan
*               If bit 30 is on then no-data in this channel
                lv = ctol_cse(j,i,ep)
* MOD TAH 200507: Changed to use fL1/2(lv) rather than fL1/2(j). 
*               Probably no effect because as coded all frequencies
*               are the same.  Test on zero values is needed for
*               single frequency data.
                if( .not.kbit(data_flag_cse(j,i,ep),30) ) then
                    L1r_phs_cse(j,i,ep) = L1r_phs_cse(j,i,ep) - 
     .                                    dClk*fL1(lv)/fClk
                    if(  L2r_phs_cse(j,i,ep).ne.0 )
     .              L2r_phs_cse(j,i,ep) = L2r_phs_cse(j,i,ep) - 
     .                                    dClk*fL2(lv)/fClk
                    L1r_rng_cse(j,i,ep) = L1r_rng_cse(j,i,ep) - 
     .                                    dClk*fL1(lv)/fClk
                    if(  L2r_rng_cse(j,i,ep).ne.0 )
     .              L2r_rng_cse(j,i,ep) = L2r_rng_cse(j,i,ep) - 
     .                                    dClk*fL2(lv)/fClk
                endif
            end do
         end do
         save_clk(2,i) = 0.d0  ! Reset value
      end do
*
*     Thats all
      return
      end

 
