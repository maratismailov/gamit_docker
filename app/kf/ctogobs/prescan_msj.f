CTITLE PRESCAN_MSJ

      subroutine prescan_msj(L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, params_cse, par_flag_cse  )
 
      implicit none
 
*     Routine to look for milli-second jumps in the clocks
*     used in the model and to remove msec jumps from clocks
*     and range and phase data.  Routine only called if prefit_clk
*     is being used.
*     Routine added: TAH 200507.
 
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
      integer*4 lv    ! Satellite index number.
      integer*4 ep    ! Epoch number
      logical   kbit  ! Test bit status

      real*8 ms_offsets  ! accumulated time offset due to 
                      ! milli-second jumps (msec)
      real*8 jump_ms  ! Estimate jump in data between current
                      ! and previous epoch. 
      real*8 prev_clk ! Normally params_cse(i,ep-1) but if there
                      ! are missing epochs might be earlier.(cyc)
      real*8 dms      ! Difference between jump and 1 msec (cyc) 


* MOD TAH 200507: Pre-scan clocks to see if there are millisecond
*     jumps in the clocks.  Correct RINEX will not have this problem
*     but it does occur.  When found correct both range and phase
*     assuming RINEX conversion did not account for jumps.

      do i = 1, num_cfiles
*        Loop over epochs (start from ep 2, take epoch 1 as reference
         ms_offsets = 0
         prev_clk = params_cse(i,1)
         do ep = 2, num_ep 
             if( .not.kbit(par_flag_cse(i,ep),1) ) then
                jump_ms  = nint((params_cse(i,ep) - 
     .                           prev_clk)/(fCLk*1.d-3))
                dms = (params_cse(i,ep)-prev_clk)-
     .                 jump_ms*(fCLk*1.d-3)
                if( jump_ms.ne.0.0d0 ) then
                   ms_offsets = ms_offsets + jump_ms
                   write(*,50) ep,  cf_codes(i), ms_offsets, 
     .                         dms/fClk*vel_light
 50                format('At Epoch ',i5,' MSec jump at ',a4,
     .               ' Total ', f6.1,' Residual ',f10.2,
     .               ' m PARAMS_CSE')
                endif
                prev_clk = params_cse(i,ep)
             endif
*            If there are millisecond jumps accumulating, 
*            correct the clock, range and phase data
             if( ms_offsets.ne.0 ) then 
*               Correct clock to make continous (cycles)
                params_cse(i,ep) = params_cse(i,ep) -
     .                             ms_offsets*(fClk*1.d-3)
*               Now correct the range and phase values for this change
*               in the clock
                do j = 1, actual_max_chan               
                   lv = ctol_cse(j,i,ep)
                   if( .not.kbit(data_flag_cse(j,i,ep),30) ) then
                      L1r_rng_cse(j,i,ep) = L1r_rng_cse(j,i,ep) - 
     .                                ms_offsets*(fL1(lv)*1.d-3)
                      if( nol1only ) 
     .                L2r_rng_cse(j,i,ep) = L2r_rng_cse(j,i,ep) -
     .                                ms_offsets*(fL2(lv)*1.d-3)
C MOD TAH 200508: Don't apply to phase?  Offset will be caught later
C                     if they are present.
C                     L1r_phs_cse(j,i,ep) = L1r_rng_cse(j,i,ep) - 
C    .                                ms_offsets*(fL1(lv)*1.d-3)
C                     if( nol1only )
C    .                L2r_phs_cse(j,i,ep) = L2r_rng_cse(j,i,ep) -
C    .                                ms_offsets*(fL2(lv)*1.d-3)
                   end if
                end do
             end if
         end do
      end do

****  That all
      return
      end

