CTITLE PRE_SCAN_CLK
 
      subroutine prescan_clk(ep,  L1r_rng_cs, L2r_rng_cs,
     .                ctol_cs, data_flag_cs, dclk_next, par_flag_cs )

      implicit none
 
*     This routine will pre-scan the range residuals at this epoch
*     and see how well they match the current observations.  The idea
*     is to catch "large" jumps in the receiver clocks and or
*     bad range data. We do this by first averaging over all data at
*     each site to get the satellite clocks (some clocks may be flagged
*     because non of the data passed the quality check.  This is true
*     for the first observation.) Once the satellite clocks are computed
*     we average over satellites at aeach station to get station clock.
*     The residuals are then checked to see if we have a large jump in
*     station clocks or some bad data.  Bad ranges are flagged here but
*     never unflagged.  (We could unflag once the final clocks are
*     computed if they appear OK at that point.)
 
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   ctol_cs(num_chan, num_cfiles) - Conversion from channel
*                 - number to local list of PRN's
*   data_flag_cs(num_chan, num_cfiles) - Data flags (to see if any
*                 - good.)
*   par_flag_cs(num_param)  - Parameter quality flag so that we
*                 - set the jump bit if needed.
*   ep      - Epoch number being processed.
 
 
 
      integer*4 ctol_cs(num_chan, num_cfiles),
     .    data_flag_cs(num_chan, num_cfiles),
     .    par_flag_cs(num_param), ep
 
*   L1r_rng_cs(num_chan, num_cfiles)    - L1 range measurements (L1 cycles)
*   L2r_rng_cs(num_chan, num_cfiles)    - L2 range measurements (L2 cycles)
*   dclk_next(num_cfiles)   - Expected change from the apriori for the
*                             next clock epoch

 
      real*8 L1r_rng_cs(num_chan, num_cfiles),
     .    L2r_rng_cs(num_chan, num_cfiles), 
     .    dclk_next(num_cfiles)
 
* LOCAL VARIABLES
 
*    i,j,k         - Loop counters
*    lv          - satellite number
*    ns          - parameter number of satellite
*    num_gmean_svs(max_gsvs)      - Number of values in mean by
*                  satellite for data that pass the tolerance test
*    num_bmean_svs(max_gsvs)      - Number of values in mean by
*                  satellite for data that do not pass the tolerance
*                  test
*    num_gmean(max_cfiles) - Number of values in means by receiver for
*                  - "good" data (see notes on gmean_rcv)
*    num_bmean(max_cfiles) - Number of values in means by receiver for
*                  - "bad" data (see notes on gmean_rcv)
 
*    num_avail(max_cfiles)- number of values at site which did
*                         - or did not trip the jump condition
*    mean_sites(max_cfiles) - List of corresponding site numbers
*                         - after the mean list have been sorted.
*    ir   - Reference site for the smallest clock jump (used when seeing
*           if a clock jumped or not).
*    ms   - Site with the samllest jump (used in checking range quality).
*    ch   - Channel number
*    ltoc - Integer*4 function to return the channel number of a given
*           satellite list entry
*    nc   - short variale name for num_cfiles (keeps lines shorter)
 
      integer*4 i,j, lv, ns, num_gmean_svs(max_gsvs),
     .    num_bmean_svs(max_gsvs),
     .    num_gmean(max_cfiles), num_bmean(max_cfiles),
     .    num_avail(max_cfiles), mean_sites(max_cfiles),
     .    ms, ir, ch, ltoc, nc
 
*   data_ok     - Logical function which is true if data OK
*   one_site_set    - Indicates at least one site tripped the
*               - jump detector
*   site_set(max_cfiles) - Indicates that this speific site failed the
*               - the clock jump detector
*   kbit        - Check if bit set.
 
      logical data_ok, one_site_set, site_set(max_cfiles), kbit
 
*   prefit(max_gsvs,max_cfiles)     - prefit residuals boxed for
*               - easy analysis
*   site_res(max_gsvs, max_cfiles)  - Residuals after removing the
*               - mean satellite clocks
*   res         - Generic residual values
*   gmean_svs(max_gsvs)    - Mean value of satellite clock (cycles)
*   bmean_svs(max_gsvs)    - Mean value of satellite clock (cycles)
*   grms_svs(max_gsvs)     - RMS (unweighted) fit to the clock value for
*                  data that pass tolerance.
*   brms_svs(max_gsvs)     - RMS to clock fit (for data that do not
*                  pass tolerance.
*   gwght_svs    - Sum of the weights for doing mean_svs calcuation
*                 (used instead of simple mean becuase of the
*                  different quality ground stations)
*   bwght_svs    - Sum of the weights for data that do not pass
*                  tolerance test on satellite
*   gmean_rcv(max_cfiles)    - Mean receiver clock are updating
*               - satellite clocks using "good" data (i.e., data which
*               - which did not exceed the range tolerance jump.
*   grms_rcv(max_cfiles)     - RMS scatter of "good" ranges to a station
*   bmean_rcv(max_cfiles)    - Mean receiver clock are updating
*               - satellite clocks using "bad" data (i.e., data which
*               - which did exceed the range tolerance jump.  This does
*               - not mean that the data bad per se just that all of the
*               - tripped the rnage jump dectector (by comparing good
*               - bad we will determine if this is bad data or a jump in
*               - the clock
*   brms_rcv(max_cfiles)     - RMS scatter of "bad" ranges to a station
 
*   rng_omc     - Function to return range residual
*   jump_mag     - Magntiude of jump (set based on clock noise
*                 and a mimimum value (0.1 usec))

*   ms_offsets(max_cfiles) -- Millisecond offsets in the clock.  Saved
*    in milliseconds
*   jump_ms  -- Jump in clock to see if clock to millisecond jump
*   dms      -- Difference between millisec and jump.
 
      real*8 prefit(max_gsvs,max_cfiles),
     .    site_res(max_gsvs, max_cfiles), res, gmean_svs(max_gsvs),
     .    bmean_svs(max_gsvs), grms_svs(max_gsvs), brms_svs(max_gsvs),
     .    gwght_svs, bwght_svs,
     .    gmean_rcv(max_cfiles), bmean_rcv(max_cfiles),
     .    grms_rcv(max_cfiles), brms_rcv(max_cfiles),
     .    rng_omc, jump_mag, res_ms, ms_offsets(max_cfiles),
     .    jump_ms, dms

      real*8 summ, sumv, meanv, rmsv
      integer*4 numm

      save ms_offsets

****  If this is the first epoch, reset the millisecond offsets
      if( ep.eq.1 ) then
         do i = 1, num_cfiles
            ms_offsets(i) = 0.d0
         end do
      end if

****  Correct for any millioffsets accumulated so far.  (We remove
*     from dclk_next and the range data)
C      reset_tol = 0     ! Commented out 030127 TAH
C     Un-commented code ! TAH 040917.  Not clear why reset to zero?
      do i = 1, num_cfiles
****     See if there are any new millisecond offsets
         jump_ms  = nint(dclk_next(i)/(fCLk*1.d-3)-ms_offsets(i))
         dms = dclk_next(i)-(ms_offsets(i)+jump_ms)*(fClk*1.d-3)
         if( abs(dms).lt.reset_tol .and. jump_ms.ne.0.d0 ) then
*            OK, we have a new millisecond jump
             call sbit(par_flag_cs(i),3,1)
             ms_offsets(i) = ms_offsets(i) + jump_ms
             write(*,90) ep,  cf_codes(i), jump_ms, dms/fClk*1.d6
 90          format(' At Epoch ',i5,' Milli-second jump at ',a4,
     .           ' Size ', f11.1,' Residual ',f10.6,
     .           ' usec (dclk_next)')
         end if
                                                          

         if( ms_offsets(i).ne.0 ) then 
             dclk_next(i) = dclk_next(i) - ms_offsets(i)*(fClk*1.d-3)
*            Now correct the range values for this change in the clock
             do j = 1, actual_max_chan               
                lv = ctol_cs(j,i)
                if( .not.kbit(data_flag_cs(j,i),30) ) then
                    L1r_rng_cs(j,i) = L1r_rng_cs(j,i) -
     .                                ms_offsets(i)*(fL1(lv)*1.d-3)
                    if( L2r_rng_cs(j,i).ne.0.d0 )
     .              L2r_rng_cs(j,i) = L2r_rng_cs(j,i) -
     .                                ms_offsets(i)*(fL2(lv)*1.d-3)
                end if
             end do
         end if
      end do
 
***** Scan mean and rms of pre-scan range residuals to see if we
*     might need to update the dclk_next value

      do i = 1, num_cfiles
         summ = 0.0
         sumv = 0.0
         numm = 0
         do j = 1, actual_max_chan
            if( data_ok(data_flag_cs(j,i),0, rng_mask) ) then
                lv = ctol_cs(j,i)
*               compute range O-C
                ns = num_cfiles + ctol_cs(j,i)
                if( kbit(status_rep,16) )
     .          write(*,105) ep, cf_codes(i), prn_list(ctol_cs(j,i)),
     .               L1r_rng_cs(j,i),L2r_rng_cs(j,i),
     .                apr_clk_val(i),dclk_next(i),apr_clk_val(ns),
     .                fL1(lv)*1e-6, fL2(lv)*1.e-6
 105            format('In epoch ',i5,1x,a4,1x,' PRN ',i2.2, 
     .                 ' Rs ',2F15.3,' Clk ',3F15.3,' Fs ',2F10.3)
                res = rng_omc(L1r_rng_cs(j,i),L2r_rng_cs(j,i),
     .                apr_clk_val(i)+dclk_next(i),apr_clk_val(ns),
     .                fL1(lv), fL2(lv) )
                if( kbit(status_rep,16) )
     .          write(*,110) ep, cf_codes(i),
     .                prn_list(ctol_cs(j,i)),
     .                res/fL1(lv)*vel_light, 
     .                apr_clk_val(ns)/fL1(lv)*vel_light
 110            format('At epoch ',i5,' Pre-Range  ',
     .                 a4,' PRN',i2.2,' Residual ',F15.2,' m,',
     .                 ' SV Clock ',f10.2,' m')
                summ = summ + res
                sumv = sumv + res**2
                numm = numm + 1
             end if
          end do
          if ( numm.gt.1 ) then
              meanv = summ/numm
              rmsv = sqrt(sumv/numm - meanv**2)
                if( kbit(status_rep,16) )
     .          write(*,115) cf_codes(i), numm,  
     .                dclk_next(i)/fL1(lv)*vel_light,
     ,                meanv/fL1(lv)*vel_light,
     .                rmsv/fL1(lv)*vel_light
 115          format('Clock Adj ',a,1x,I3,' dC ',3F15.2,' m')
              if( rmsv .lt. 100.d0 .and.
     .            abs(meanv).gt.rng_res_min ) then
                  dclk_next(i) = dclk_next(i) + meanv
              end if
          else
             meanv = 0.0
             rmsv = 0.0
          endif
      end do

 
***** Scan all the range residuals first to remove any large
*     outliers.  The initial clock estimates are from model for
*     stations and the broadcast ephemeris for the satellites
*     Each of these should be accurate within a few hundred 
*     meters.
      do i = 1, num_cfiles
         do j = 1, actual_max_chan
            if( data_ok(data_flag_cs(j,i),0, rng_mask) ) then
                lv = ctol_cs(j,i)
*               compute range O-C
                ns = num_cfiles + ctol_cs(j,i)       
                res = rng_omc(L1r_rng_cs(j,i),L2r_rng_cs(j,i),
     .                apr_clk_val(i)+dclk_next(i), apr_clk_val(ns),
     .                fL1(lv), fL2(lv) )
                if( ep.le.0 ) then
                   print *,'EP: I,J ',ep, i,j, lv, ns
                   print *,'RNG VAL ', L1r_rng_cs(j,i),L2r_rng_cs(j,i)
                   print *,'CLK VAL ',apr_clk_val(i), dclk_next(i), 
     .                  apr_clk_val(ns), res
                end if
                res_ms = nint(res/(fL1(lv)*1.d-3))
                jump_mag = max(rng_res_min,rng_res_tol*rng_noise(i))
                jump_mag = min(jump_mag, rng_res_max)     
                if( abs(res).gt.jump_mag ) then                                
                     call sbit(data_flag_cs(j,i),13,1)
                     write(*,120) ep, cf_codes(i),
     .                     prn_list(ctol_cs(j,i)),
     .                     res/fL1(lv)*vel_light, 
     .                     apr_clk_val(ns)/fL1(lv)*1.d6
 120                 format(' At epoch ',i5,' BAD Pre-RNG data for ',
     .                      a4,' PRN',i2.2,' Residual ',F15.2,' m,',
     .                      ' SV Clock ',f10.3,' us')
                     call sbit(data_flag_cs(j,i),13,1)
                 endif
             end if       
          end do
      end do     

***** Clear the prefit box before we start
      nc = num_cfiles
      do i = 1, num_cfiles
          do j = 1, num_sat
              prefit(j,i) = 0.d0
          end do
      end do
 
*     Initialize the parameter quality flag for this epoch.  By
*     setting bit 1, we make the parameter estimate as no-data.
      do i = 1, num_param
          par_flag_cs(i) = 1
      end do
 
*
****  Construct a uniform matrix of range residuals index by
*     list and site.  We can then easily where the contributions are.
*     If L2_range is available use an ionosphere free combination.
 
      do i = 1, num_cfiles
          do j = 1, actual_max_chan
              if( data_ok(data_flag_cs(j,i),0, rng_mask) ) then
                  lv = ctol_cs(j,i)
*                 compute range O-C
                  ns = num_cfiles + ctol_cs(j,i)
C                 res = rng_omc(L1r_rng_cs(j,i),L2r_rng_cs(j,i),
C    .                apr_clk_val(i),apr_clk_val(ns),
C    .                fL1(lv), fL2(lv)  )
* MOD TAH 931222: Use the value of the clock offset from the c-file
*                 as the first guess.
                  res = rng_omc(L1r_rng_cs(j,i),L2r_rng_cs(j,i),
     .                apr_clk_val(i)+dclk_next(i),apr_clk_val(ns),
     .                fL1(lv), fL2(lv) )
 
*                 Becuase we use zero residual to indicate that there
*                 no data, change any exactly zero residual to 0.001
*                 cyles
                  if( res.eq.0 ) res = 0.001d0
 
*                 Save the residual
                  prefit(ctol_cs(j,i),i) = res
              end if
          end do
      end do
 
****  Now scan down each satellite.  If all is Ok then all the values
*     from all the sites should be about the same if not we may have
*     a site jump.
 
      one_site_set = .false.
      do i = 1, num_cfiles
          site_set(i) = .false.

*         Check to see if there is a jump in clock between the apr_clk_val
*         based on the previous epoch and on the value that was read from
*         the cfile
          jump_mag = max(rng_jump_min,rng_jump_tol*apr_clk_sd(i))
          if( abs(dclk_next(i)).gt.jump_mag ) then
              call sbit(par_flag_cs(i),2,1)
              write(*,100) ep, cf_codes(i), dclk_next(i)/fClk*1.d6
 100          format(' At Epoch ',i4,' Clock jump at ',a4,' Size ',
     .            f11.3,' usec (dclk_next)')
              apr_clk_val(i) = apr_clk_val(i) + dclk_next(i)
          end if
      end do
 
      do i = 1, num_sat
 
*         Loop over all the sites
          gmean_svs(i) = 0.d0
          bmean_svs(i) = 0.d0
          grms_svs(i)  = 0.d0
          brms_svs(i)  = 0.d0
          gwght_svs = 0.d0
          bwght_svs = 0.d0
          num_gmean_svs(i) = 0
          num_bmean_svs(i) = 0
          do j = 1, num_cfiles
 
*             Increment mean if value is not zero (no data) and
*             this is not one of values which exceeds to tolerance
C             jump_mag = max(rng_jump_min,
C    .                   rng_jump_tol*sqrt(apr_clk_sd(j)**2 +
C    .                                     apr_clk_sd(nc+i)**2))
              jump_mag = rng_jump_min
 
*             if we have data and it passed the tolerance then
*             accumulate with the good data, otherize put in the
*             bad data (may be OK, we will use the rms to tell).
*             NOTE: The rms is unweighted
              if( prefit(i,j).ne.0 .and.
     .            abs(prefit(i,j)).lt.jump_mag) then
                  gmean_svs(i) = gmean_svs(i) +
     .                           prefit(i,j)/apr_clk_var(j)
                  grms_svs(i) = grms_svs(i) + prefit(i,j)**2/
     .                           apr_clk_var(j)
                  gwght_svs = gwght_svs + 1.d0/apr_clk_var(j)
                  num_gmean_svs(i) = num_gmean_svs(i) + 1
              end if
*             If it does not pass tolerance then accumulate under
*             bad data
              if( prefit(i,j).ne.0 .and.
     .            abs(prefit(i,j)).ge.jump_mag) then
                  bmean_svs(i) = bmean_svs(i) +
     .                           prefit(i,j)/apr_clk_var(j)
                  brms_svs(i) = brms_svs(i) + prefit(i,j)**2/
     .                           apr_clk_var(j)
                  bwght_svs = bwght_svs + 1.d0/apr_clk_var(j)
                  num_bmean_svs(i) = num_bmean_svs(i) + 1
              end if
          end do
 
*         Compute the mean satellite clock
          if( num_gmean_svs(i).gt.0 ) then
              gmean_svs(i) = gmean_svs(i) / gwght_svs
              grms_svs(i) = sqrt(grms_svs(i)/gwght_svs -
     .                           gmean_svs(i)**2+1.d0)
             
C             write(*,920) i, gmean_svs(i)/fClk*1.d6, num_gmean_svs(i)
C920          format(i3,f10.2,1x,i3)
          end if
 
****      See if we need to use the "bad" satellite clock data (i.e.,
*         data which did not pass requirements
          if( num_gmean_svs(i).eq.0 .and. num_bmean_svs(i).gt.0 ) then
              bmean_svs(i) = bmean_svs(i) / bwght_svs
              gmean_svs(i) = bmean_svs(i)
              brms_svs(i) = sqrt(brms_svs(i)/bwght_svs -
     .                          bmean_svs(i)**2+1.d0)
 
*             Print a warning if no good data and brms is large
              if( brms_svs(i).gt. rng_jump_min ) then
                  write(*,940) ep, prn_list(i),
     .                   brms_svs(i)/fClk*1.d6,
     .                   bmean_svs(i)/fClk*1.d6, num_bmean_svs(i)
 940              format('**WARNING** Epoch ',i5,' PRN ',i2.2,
     .                   ' RMS on clock fit ', F10.2, ' usec',
     .                   ' Mean ',F10.2, ' usec, number ',I3)
              end if
          end if
 
*****     Now compute the residuals at each site (These will be
*         used to check for clock jumps and bad range data)
          if( num_gmean_svs(i).gt.0 .or. num_bmean_svs(i).gt.0 ) then
              do j =  1, num_cfiles
                  jump_mag = max(rng_jump_min,
     .                           rng_jump_tol*sqrt(apr_clk_sd(j)**2 +
     .                                             apr_clk_sd(nc+i)**2))
                  if( prefit(i,j).ne.0 ) then
                      res = prefit(i,j) - gmean_svs(i)
                      site_res(i,j) = res
                      if( abs(res).gt. jump_mag) then
                          site_set(j) = .true.
                          one_site_set = .true.
                      end if
                  end if
              end do
          else
*             Clear all of the site residuals for this satellite
              do j = 1, num_cfiles
                 site_res(i,j) = 0
              end do
          end if
      end do
 
****  now scan accross the site_res line (by site to see who is bad)
*     Do this if there are appears to be jump.  If not continue
*     processing
 
      if( .not. one_site_set ) RETURN
 
      if( kbit(status_rep,1) ) write(*,190)
 190  format('** Processing clock jump value')
 
      do j = 1, num_cfiles
          gmean_rcv(j) = 0.d0
          grms_rcv(j)  = 0.d0
          num_gmean(j) = 0
          bmean_rcv(j) = 0.d0
          brms_rcv(j)  = 0.d0
          num_bmean(j) = 0
          num_avail(j) = 0
 
*         Accumulate the mean range and rms for the data passing
*         and failing the rnage jump detector
 
          do i = 1, num_sat
*                                 ! We have data, see if station failed
              if( prefit(i,j).ne.0 ) then
 
*                 Get the jump magnitude for this site.  Accumulate the
*                 good and "bad" data for this site.  (Bad just may mean
*                 that there has been a jump in the clock.)
                  jump_mag = max(rng_jump_min,
     .                           rng_jump_tol*sqrt(apr_clk_sd(j)**2+
     .                                             apr_clk_sd(nc+i)**2))
                  res = site_res(i,j)
                  if( abs(res).gt. jump_mag) then
*                                       ! accumulate jump data
                      bmean_rcv(j) = bmean_rcv(j) + res
                      brms_rcv(j) = brms_rcv(j) + res**2
                      num_bmean(j) = num_bmean(j) + 1
*                                       ! accumulate good non-jump data
                  else
                      gmean_rcv(j) = gmean_rcv(j) + res
                      grms_rcv(j) = grms_rcv(j) + res**2
                      num_gmean(j) = num_gmean(j) + 1
                  end if
                  num_avail(j) = num_avail(j) + 1
              end if
          end do
 
*         Compute mean if we have data
          if ( num_gmean(j).gt.0 ) then
              gmean_rcv(j) = gmean_rcv(j)/num_gmean(j)
              grms_rcv(j) = sqrt( grms_rcv(j)**2/num_gmean(j) -
     .                            gmean_rcv(j)**2+1.d0 )
          end if
          if ( num_gmean(j).eq.0 .and. num_bmean(j).gt.0 ) then
              bmean_rcv(j) = bmean_rcv(j)/num_bmean(j)
              gmean_rcv(j) = bmean_rcv(j)
              brms_rcv(j) = sqrt( brms_rcv(j)**2/num_bmean(j) -
     .                            bmean_rcv(j)**2+1.d0 )
          end if
      end do
 
*     Find the two largest jumps.  If the largest is much larger
*     than the others, then just one jump probably.
      call jmp_sort( num_cfiles, gmean_rcv, mean_sites )
 
*     Find the site with the samllest jump (but one which is non-zero
*     since zero means no data)
      ir = 0
      do i = 1,num_cfiles
         if( gmean_rcv(i).ne.0 ) ir = i
      end do
      if( ir.eq.0 ) ir = num_cfiles
 
*     Assign jump to first site and see if there are any others
      do i = 1, num_cfiles-1
          if( abs(gmean_rcv(i)-gmean_rcv(ir)).gt.
     .        max(rng_jump_tol*sqrt(apr_clk_var(mean_sites(i))+
     .                          apr_clk_var(mean_sites(ir))),
     .            rng_jump_min)   ) then
              if( kbit(status_rep,2) ) 
     .                     write(*,200) ep, cf_codes(mean_sites(i)),
     .                     gmean_rcv(i)/fClk*1.d6
 200          format(' At Epoch ',i4,' Clock jump at ',a4,' Size ',
     .            f11.3,' usec')
              apr_clk_val(mean_sites(i)) = apr_clk_val(mean_sites(i)) +
     .                                     gmean_rcv(i)
 
*             Set the jump indicator in the parameter flag
              call sbit(par_flag_cs(mean_sites(i)),2,1)
          end if
 
****      Now check for "bad" range data.  Here we look at the bad data
*         if there was good data.
          if( num_bmean(mean_sites(i)).gt.0 ) then
              ms = mean_sites(i)
 
*             Loop on the editing routine for this data
              j = 99
              do while ( j.gt.0 )   
                  lv = ctol_cs(j,i)
                  jump_mag = max(rng_res_tol*rng_noise(ms),
     .                           rng_res_min ) 
                  call preedit_range(prefit(1,ms), site_res(1,ms),
     .                   gmean_rcv(ms), num_sat, jump_mag, j )
                  if( j.gt.0 ) then
*                                             ! This stops the data
                     prefit(j,ms) = 0.d0
*                                                 ! data being further
*                                                 ! considereed.
                     ch = ltoc(ctol_cs(1,ms), j, actual_max_chan)
                     call sbit(data_flag_cs(ch,ms),13,1)
                     write(*,400) ep, cf_codes(ms),
     .                  prn_list(j),
     .                 (site_res(j,ms)-gmean_rcv(ms))/fL1(lv)*vel_light
 400                 format(' At epoch ',i5,' BAD Range data for ',A4,
     .                      ' PRN',i2.2,' Residual ',F10.2,' m')
                  end if
              end do
          end if
 
      end do
 
****  Thats all
      return
      end
 
CTITLE PREEDIT_RANGE
 
      subroutine preedit_range(prefit, site_res, gmean_rcv, 
     .                         num_sat, tol, edit )

      implicit none
 
*     Routine to iteratively check the  range residuals and edit out the
*     "bad" ones.  If no bad data is found then edit is returned as -1.
 
* PASSED VARIABLES
 
*   num_sat     - Number of satellites
*   edit        - Satellite number to edit.  If all the data looks
*               - OK then edit is returned -1
 
      integer*4 num_sat, edit
 
*   prefit(num_sat)     - Prefit residuals before svs clocks are
*               - are corrected (a zero value here mean no data)
*   site_res(num_sat)   - Range residuals after correcting for
*               - for satellite clocks.
*   gmean_rcv   - Mean offset for the reciever.
*   tol         - Number of cycles of error allowed for this station
*               - this is the product of the rng_res_tol and the
*               - rng_noise.
 
      real*8 prefit(num_sat), site_res(num_sat), gmean_rcv, tol
 
* LOCAL VARIABLES
 
*   num         - Number of ranges used in checking the quality
*   i           - Loop counter
*   is          - Channel with biggest residual in the site_res list
*   ip          - Channel with biggest residual in the prefit list (the
*                 reason we check both is becuase the sattelit clocks
*                 be be distorted wehn one of the sites is left out.
 
      integer*4 num, i, is, ip
 
*   means       - Mean of the residuals for site_res
*   meanp       - Mean of the residuals for prefit
*   biggests    - Biggest residual from the mean for the site_res
*   biggestp    - Biggest residual from the mean for the prefits.
 
      real*8 means, meanp, biggests, biggestp
 
****  Compute the mean and see who has the biggest resiudak
      edit = -1
      means = 0
      meanp = 0
      num  = 0
      is = 1
      ip = 1
      do i = 1, num_sat
          if( prefit(i).ne.0 ) then
              means = means + (site_res(i)-gmean_rcv)
              meanp = meanp + prefit(i)
              num = num + 1
          end if
      end do
 
***   See if we have enough data to check.
      if ( num.ge.2 ) then
          means = means/ num + gmean_rcv
          meanp = meanp/ num
          biggests = 0
          biggestp = 0
          do i = 1, num_sat
              if( prefit(i).ne.0 .and.
     .            abs(site_res(i)-means).gt.biggests ) then
                      biggests = abs(site_res(i)-means)
                      is = i
              end if
              if( prefit(i).ne.0 .and.
     .            abs(prefit(i)-meanp).gt.biggestp ) then
                      biggestp = abs(prefit(i)-meanp)
                      ip = i
              end if
          end do
 
****      See if biggest is too big
          if( biggests.gt.tol .and. biggestp.gt.tol ) then
              edit = is
              if( is.ne.ip ) then
                   write(*,100) is, ip
 100               format('**WARNING** The following edit for channel ',
     .                     i2,' in site_res different from channel ',
     .                     i2,' in prefit')
              end if
          end if
      end if
 
****  Thats al
      return
      end
 
