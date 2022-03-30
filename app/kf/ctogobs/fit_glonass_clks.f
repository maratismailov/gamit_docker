CTITLE FIT_GLONASS_CLKS
 
      subroutine fit_glonass_clks( L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, params_cse, par_flag_cse  )

      implicit none
 
*     Routine to estimate GLONASS clocks baecuase the values
*     coming from model and saved in the cfiles seem to be very
*     bad )satellites possible have wrong sign).
 
 
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
 
*   L1r_rng_cse(num_chan, num_cfiles, num_ep)  - L1 range residuals
*                   - cycles at L1
*   L2r_rng_cse(num_chan, num_cfiles, num_ep)  - L2 range residuals
*                   - cycles at L2
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)
 
* LOCAL VARIABLES
      integer*4 ep    ! Epoch counter
     .,         i,j   ! Loop counter
     .,         lv    ! Satellite list number
     .,         av_num   ! Number of values on station clock sum
     .,         ns    ! Satellite position in parameter array

      logical data_OK  ! Function to test if data flags are OK.

      real*8  res     ! Range residual (cycles, ion_free)
     .,       av_sum, av_var   ! Sum of residuals and res^2
     .,       av_clk, av_std   ! Station clock estimate and std dev.
     .,       rng_omc     ! Function to return range residual.

 

****  Code assumes that satellite clock hae be be changed in sign.
*     We then estimate the receiver clock

      call set_rng_mask( rng_mask )

      write(*,'("RANGE MASK ",o10)' ) rng_mask
      do ep = 1, num_ep

*         Bug fixed in makej so clock now has correct sign.
          do j = 1, num_sat
             ns = num_cfiles + j
             apr_clk_val(ns) = params_cse(ns,ep)
             params_cse(ns,ep) = apr_clk_val(ns)
          end do

****      Now loop over stattions to get estimate of station 
*         clock
          do i = 1, num_cfiles
*            Compute range O-C
             av_sum = 0.d0
             av_var = 0.d0
             av_num = 0
             do j = 1, actual_max_chan
                if( L1r_rng_cse(j,i,ep).eq.0 .or. 
     .              L2r_rng_cse(j,i,ep).eq.0 ) then
                   call sbit(data_flag_cse(j,i,ep),30,1)
                endif

                if( data_ok(data_flag_cse(j,i,ep),0, rng_mask) ) then
                   lv = ctol_cse(j,i,ep)    ! Satellite number
                   if( lv.gt.0 ) then 
                      ns = num_cfiles + lv 
*                     Compute residual with zero stattion clock
                      res = rng_omc(L1r_rng_cse(j,i,ep),
     .                              L2r_rng_cse(j,i,ep),
     .                      0.d0,apr_clk_val(ns), fL1(lv), fL2(lv) )
                      av_sum = av_sum + res
                      av_var = av_var + res**2
                      av_num = av_num + 1
C                     if( ep.lt.10 )
C    .                print *,'EP RES ', ep, i, lv, res,  
C    .                    apr_clk_val(ns), ' RCV Clk ', params_cse(i,ep)
                    end if
                end if
             end do
             if( av_num.gt.0 ) then
                av_clk = av_sum/av_num
                av_std = sqrt((av_var - av_clk**2*av_num)/(av_num-1))
                if( abs((params_cse(i,ep)-av_clk)/fClk*1.d6).gt.1.0 )
     .          print *,'EP AVCLK ',ep, i, cf_codes(i), av_clk, av_std,          
     .               ' m, RCV ', params_cse(i,ep), ' m ',
     .               (params_cse(i,ep)-av_clk)/fClk*1.d6,' dC us'
                params_cse(i,ep) = av_clk
                par_flag_cse(i,ep) = 0
             end if
          end do
      end do

****  Thats all
      return
      end
 

                   
  

