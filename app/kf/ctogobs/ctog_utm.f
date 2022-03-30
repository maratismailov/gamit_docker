 
CTITLE READ_CFDATA
 
      subroutine read_cfdata(L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse , data_flag_cse, bf_type_cse, params_cse, 
     .    par_flag_cse, azel_cse , pf_dphs_cse, svcL1_ce,  serr )
 
      implicit none

*     Routine to read all of the cfile data into memory
*     so that the receiver and ground clocks can be computed
*     and approximate values for the one-way cycles at L1
*     and L2 computed.
 
* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include '../includes/mfile_def.h'
      include '../includes/gobs_header.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)
*   par_flag_cse(num_param)  - parameter estimate flag (zero means OK)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number

*   serr            - Maximum IOSTAT error to occurr
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    par_flag_cse(num_param,num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep), serr
 
*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement.  Units are
*                   - L1 cycles (may be fractional for half wavelength
*                   - units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)   - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)   - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse(num_chan, num_cfiles, num_ep)   - L1 range residuals
*                   - cycles at L1
*   L2r_rng_cse(num_chan, num_cfiles, num_ep)   - L2 range residuals
*                   - cycles at L2
*   params_cse(num_param,num_ep)                - Clock correction for
*                     sites read from cfile (L1 cycles)
*   pf_dphs_cse(num_chan, num_cfiles, num_ep)   - L1 adjustments to the
*                     omc to account for postfit residuals (L1 cycles)
*   svcL1_ce(num_sat,num_ep) - Part of satellite clock removed by model
*                     (L1 cycles)
  
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep), 
     .    pf_dphs_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param,num_ep)
      real*4
     .    svcL1_ce(num_sat,num_ep)
     
*  azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
*                      angles to be written to phs_res_root files.
*                      (Real*4, radians)

      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)  
 
* LOCAL VARIABLES
 
*   i,j,k   - Loop counters
*   ierr    - IOSTAT error.
*   num_good_rng - Counts number of good observations at each
*              epoch (used to see if we have a reliable clock
*              estimate)
*   rng_mask_loc - Local copy of the range mask
*   mfs          - Site number in mfile
*   trimlen      - Length of string

      integer*4 i,j,k, ierr, num_good_rng, rng_mask_loc, 
     .          mfs, trimlen
 
*   clk_not_found - Indicates that we have not found the
*        station clock initial value yet.
*   kbit          - Checks to see if bit is set
*   data_OK       - Function to see if data is OK.
 
      logical clk_not_found, kbit, data_OK

*   dt  - Time since start of experiment (epoch number -1)*
*         sampling interval
*   dphs_L1 - Postfit correction to phase at L1 (for postfit residuals)

      real*8 dt, dphs_L1

*   last_good_ep -- Epoch of the last good epoch of data
*   min_data_spacing -- Minimum spacing between epochs of good data
*                       (used to check sampling interval)
*   num_nolcaut -- Number delete due to LC_AUTCLN requirements.

      integer*4 last_good_ep, min_data_spacing, num_nolcaut

*   report_line -- Character string for report_stat messages

      character*128 report_line
 
***** Initialize the apriori values of the station clock
      
      do i = 1, num_cfiles
          rclock_s(i) = 0.d0
      end do
 
      do i = 1, num_param
          init_clk_val(i) = 0.d0
      end do
 
      actual_max_chan = 0
      num_obs = 0
      min_cfile_elev = 1.57

* MOD TAH 180312: Iniitialize the clock parameter flag (par_flag_cse)
*     to show no data  (changed as values are added).
      do i = 1, num_cfiles+num_sat
         do j = 1, num_ep
            par_flag_cse(i,j) = 1   ! Bit 1 set showing no data originally
         end do
      end do 

 
***** Loop over all cfiles reading all the data we need

      call set_rng_mask( rng_mask_loc )
      serr = 0

      do j = 1, num_cfiles
          write(*,100) cf_codes(j), rcv_types(j)
 100      format(' Reading data for ',a4,' Reciever ',a4)

* MOD TAH 091228: Makes sure that data types are correct when
*         L1 only data is being processed.

          do i = 1, num_gdata_types(j)

* MOD TAH 091228: Check that L2 data types are correct if 
*            we are processing L1 only data.  Save and then 
*            see if needs to be updated.
             if( .not. nol1only ) then  ! L1 only so sets 
                                        ! type 2 and 4 to zero
                 if( data_types_s(i,j).eq.2 .or. 
     .               data_types_s(i,j).eq.4 ) then
                     data_types_s(i,j) = 0
                 endif
                 do k = 1,num_sat
                     if( data_types_s(i,j).eq.0 .and. 
     .                   lambda(k,i,j).ne.0 ) then
                        lambda(k,i,j) = 0
                     endif
                 end do 
             end if
          end do
 
*         See if we are to compute postfit residuals
          if( use_postfit ) then
              call get_mfs( cf_codes(j), mfs )
              if ( mfs.gt.0 ) then
                   call set_np( mfs )
               else
                   write(report_line,115) cf_codes(j) 
                   write(*,'(a)') report_line(1:trimlen(report_line))
 115               format('Site ',a4,' NOT found in m-file')        
                   call report_stat('warning','autcln','read_cfdata',
     .                              ' ',report_line,0)
               end if
          end if
          
*         Set to say we have not found the clock values at this
*         site yet.
          clk_not_found = .true.
          last_good_ep = 0 
          min_data_spacing = num_ep
          num_nolcaut = 0

          do i = 1, num_ep
 
*             Get the time block
              dt = (i-1)*sampling

              call read_cf4(100+j,'ALL', ierr)
              serr = max(serr, ierr )
*             If error reading cfile; Kill the whole run.
              call report_error('IOSTAT',ierr,'read',cfiles(j),
     .                          1,'READ_CFDATA/read_cf4')

*             Update actual number of channels used
              actual_max_chan = max(actual_max_chan,cf_msat)

* MOD TAH 990920: Check on the miniumum data spacing
              if( cf_msat.gt.0 ) then
*                 OK, we have data, set the minimum spacing
                  if( last_good_ep.gt.0 ) then
                      min_data_spacing = min(min_data_spacing,
     .                                      (i-last_good_ep))
                  end if
                  last_good_ep = i
              end if

 
****          Now read the indivual data records
              num_good_rng = 0
              do k = 1, cf_msat
                  call read_cf5(100+j, 'ALL', ierr)

*                 If error reading cffile; Kill the whole run.
                  call report_error('IOSTAT',ierr,'read',cfiles(j),
     .                              1,'READ_CFDATA/read_cf5')

* RWK 150204: Settng ctol_cse needed to index fL1 and fL2 so move from lower down                                      
*                 save the channel number.
                  ctol_cse(k,j,i) = prntol(cf_iprn)
cd                  print *
cd     .              ,'READ_CFDATA k j i cf_iprn ctol_cse '
cd     .              ,k,j,i,cf_iprn,ctol_cse(k,j,i) 

                  num_obs = num_obs + 1 
*                 Save values we want.
*** WARNING **** This code assumes dual frequency data and
*                at least one range measurement
* MOD TAH 990512: Mod to allow single frequency L1 only data.
                  L1r_phs_cse(k,j,i) = cf_omcs(1)
                  if( data_types_s(2,j).ne.0 ) then
                     L2r_phs_cse(k,j,i) = cf_omcs(2)
                  else
                     L2r_phs_cse(k,j,i) = 0.d0     
                  end if
                  L1r_rng_cse(k,j,i) = cf_omcs(3)
*                 See if we have L2 range measurement
                  if( data_types_s(4,j).eq.4 ) then
                      L2r_rng_cse(k,j,i) = cf_omcs(4)
                  else
                      L2r_rng_cse(k,j,i) = 0.d0
                  end if
                 
****              Now see if we compute postfit residuals
                  if( use_postfit ) then
                      call comp_dphs(i, mfs, j, dphs_L1)
                      L1r_phs_cse(k,j,i) = L1r_phs_cse(k,j,i) -
     .                                     dphs_L1
                      if( L2r_phs_cse(k,j,i).ne. 0 ) then
                         L2r_phs_cse(k,j,i) = L2r_phs_cse(k,j,i) - 
     .                 dphs_L1*fL2(ctol_cse(k,j,i))/fL1(ctol_cse(k,j,i))
                      end if
                      L1r_rng_cse(k,j,i) = L1r_rng_cse(k,j,i) - 
     .                                     dphs_L1
                      if( L2r_rng_cse(k,j,i).ne.0 ) then
                          L2r_rng_cse(k,j,i) = L2r_rng_cse(k,j,i)
     .             - dphs_L1* fL2(ctol_cse(k,j,i))/fL1(ctol_cse(k,j,i))
                      end if
                      
* MOD TAH 970828: If we are going to normal point then save this adjustment
                      if( np_size.gt.1 ) then
                          pf_dphs_cse(k,j,i) = dphs_L1
                      end if
                  end if

* MOD TAH 970828: If phs_res_root given, then save the azimuth and elevation
*                 as well
                  if( (trimlen(phs_res_root).gt.0) .or.
     .                (apply_phs_clk .and. pc_max_iter.gt.1)  ) then
                      azel_cse(1,k,j,i) = cf_azimuth(1)
                      azel_cse(2,k,j,i) = cf_elev(1)
                  end if 
                  
*                 Clear the number of cycles
                  L1_cyc_cse(k,j,i) = 0.d0
                  L2_cyc_cse(k,j,i) = 0.d0
 
*                 Save the error flag and channel to list value
                  call cerrfl_to_df(cf_ierfl, cf_data_flag, 
     .                 data_flag_cse(k,j,i))

* MOD TAH 000718: See if we are to ignore original bias flags
                  if( .not.use_orig_bf ) 
     .                        call sbit(data_flag_cse(k,j,i),32,0)
                  if( kbit(data_flag_cse(k,j,i),32) ) then
                      call sbit(bf_type_cse(k,j,i), 1,1)
                  end if

* MOD TAH 950818: Check to see if the L2 and L1 range are the
*                 same.  If so kill data.
* MOD TAH 960122: Check to see if L2 range is more than 100 meters
*                 for L1 range. (This is second chekc below).
                  if( data_types_s(4,j).eq.4 ) then
                      if( cf_obsv(3).eq.cf_obsv(4) ) then
                          call sbit(data_flag_cse(k,j,i),4,1)
                      end if
                      if( abs(cf_obsv(3)-cf_obsv(4)).gt.100.d0 ) then
                          call sbit(data_flag_cse(k,j,i),4,1)
                      end if
                   endif
* MOD TAH 090124:  Also use this flag if no L2 data and in lc_autcln
*                  mode (ie., resolve_wl is true).  Test is on
*                  non-full cycle L2 or no P2 data.
                   if( resolve_wl .and. 
     .                (lambda(prntol(cf_iprn),2,j).ne.-1 .or.
     ,                 lambda(prntol(cf_iprn),4,j).ne. 1) ) then
                       call sbit(data_flag_cse(k,j,i),4,1)
                       num_nolcaut = num_nolcaut + 1
                   end if


* MOD TAH 080512: See if L2 phase is zero and we are not allowing
*                 L1 only data to be processed.  If so flag the data
                  if( nol1only .and.
     .                L2r_phs_cse(k,j,i).eq. 0.d0 ) then
                      call sbit(data_flag_cse(k,j,i),22,1)
                  end if

*                 See if this observation is less than the cleaning
*                 minimum elevation angle
C                 if( cf_elev(1).lt.site_celev(j) ) then
C                     call sbit(data_flag_cse(k,j,i), 24,1)
C                 end if

*                 See if observation will be below the final elevation
*                 cutoff. Set here since this will be considered good 
*                 data until cfiles are written
C                 if( cf_elev(1).lt.site_oelev(j) ) then
C                     call sbit(data_flag_cse(k,j,i), 27,1)
C                 end if

* MOD TAH 990519: Replaced the above code with a call to check_elev
*                 which will check azimuth mask as well as standard
*                 cutoffs.
                  call check_elev(j, cf_elev(1), cf_azimuth(1),
     .                            data_flag_cse(k,j,i))

*                 Keep track of solve minimum angle incase we don't
*                 set one (i.e., use_gamit_elev yes command)
                  if( cf_elev(1).lt. min_cfile_elev .and.
     .               (cf_ierfl.eq.0 .or. cf_ierfl.eq.10) )
     .                 min_cfile_elev = cf_elev(1)

*                 If this is minimac data, set the ignore minimac
*                 range flag (bit 26)
                  if( rcv_types(j).eq.'MIN ' ) then
                      call sbit(data_flag_cse(k,j,i), 26,1)
                  end if

*                 See if we are pre-editing this data.
                  call set_pre_edit( data_flag_cse(k,j,i),cf_iprn,j,i)

*                 See if we should edit based on low SNR
*                 If the cfile values are zero, then we have no information
*                 so set value to 9 (so wont be deleted)
                  if( cf_isnr(1).eq.0 ) cf_isnr(1) = 9
                  if( cf_isnr(2).eq.0 ) cf_isnr(2) = 9
                  if( cf_isnr(1).lt. site_snr(1,j) ) then
                      call sbit(data_flag_cse(k,j,i), 1,1)
                  end if
                  if( cf_isnr(2).lt. site_snr(2,j) ) then
                      call sbit(data_flag_cse(k,j,i), 1,1)
                  end if

*                 See if L2 is bad, but L1 may be OK (set bit 22 for
*                 this case.
                  if( cf_isnr(1).ge. site_snr(1,j) .and.
     .                cf_isnr(2).lt. site_snr(2,j) ) then
                      call sbit(data_flag_cse(k,j,i), 22,1)
                  end if

*                 See if we have good data (may move this lower later)
                  if( data_OK(data_flag_cse(k,j,i),0,rng_mask_loc ) )
     .                num_good_rng = num_good_rng + 1

* MOD TAH 990917: Save the satellite clock estimate.  We only need
*                 to save the offset part since the non-offset part has
*                 already been removed.
* MOD TAH 051116: Only save the satellite clock if it non-zero. Certain
*                 data errors in model result in no satellite clock value
*                 and we want to aviod these
* MOD TAH 180311: Clock should be saved at the reference clock frequency
*                 set in set_freqs.
C                 if( cf_svcepc.ne.0 )
C    .            params_cse(num_cfiles+prntol(cf_iprn),i) = 
C    .                           cf_svcepc*fL1(ctol_cse(k,j,i)) 
                  if( cf_svcepc.ne.0 ) then
                     params_cse(num_cfiles+prntol(cf_iprn),i) = 
     .                           cf_svcepc*fClk 
* MOD TAH 180312: Set parameter flag to zero to show an estimate
                     par_flag_cse(num_cfiles+prntol(cf_iprn),i) = 0
                  end if
  
cd        print *,'READ_CFDATA cf_svcepc cf_iprn k j i ctol_cse f params '
cd     .    , cf_svcepc,cf_iprn,k,j,i,ctol_cse(k,j,i),fl1(ctol_cse(k,j,i))
cd     .   ,  params_cse(num_cfiles+prntol(cf_iprn),i) 


* MOD TAH 000222: Save the modeled portion of the satellie clock
*                 if we are going to write the igs clock file.
* MOD TAH 061210: Only save value if range data is good.  (Calculation
*                 in model can be wrong when data not good).
                  if( write_igs_clk .and. 
     .              data_OK(data_flag_cse(k,j,i),0,rng_mask_loc)) then
                      svcL1_ce(prntol(cf_iprn),i) = cf_svcL1
                  end if
 
* RWK 150204: Move this code to the top of the loop so that it's available for fL1, fL2             
*                 save the channel number.
*                 ctol_cse(k,j,i) = prntol(cf_iprn)
*                         ! Looping over channeles
              end do
 
*             Set the data flag for the rest of data slots    at this
*             site and epoch to show no data
              do k = cf_msat+1, num_chan
                  call sbit(data_flag_cse(k,j,i),30,1)
                  ctol_cse(k,j,i) = 0
              end do
 
*             Save needed information from time block
*             See if we should set the apriori clock value
              if( num_good_rng.gt.0 .and. clk_not_found ) then
                  rclock_s(j) = cf_rclock*fClk - 
     .                (apr_clk_poly(2,j)* dt +
     .                 apr_clk_poly(3,j)*(dt)**2 +
     .                 apr_clk_poly(4,j)*(dt)**3)*fClk
                  init_clk_val(j) = rclock_s(j)
                  apr_clk_epoch(j) = i
                  clk_not_found = .false.
              end if

*             Save the site clock in the clock parameter array.
              if( num_good_rng.gt.0 ) then
                  params_cse(j,i) = cf_rclock*fClk -
     .                (apr_clk_poly(2,j)* dt +
     .                 apr_clk_poly(3,j)*(dt)**2 +
     .                 apr_clk_poly(4,j)*(dt)**3)*fClk
* MOD TAH 180312: Set parameter flag to zero to show an estimate
                  par_flag_cse(j,i) = 0
              else

*                 If this is just a gap in the data then extrapolate
*                 form the last value
                  if( i.gt.1 ) then
                      params_cse(j,i) = params_cse(j,i-1) 
                  else
                      params_cse(j,i) = 0.d0
                  end if
              end if
          
          end do

*         Check the data spacing that we finally measured.
          if( orig_sampling(j).ne.min_data_spacing*cf_inter .or.
     .        orig_sampling(j).eq.0 ) then
*             OK, they do not match. 
              write(*,290) cf_codes(j), orig_sampling(j),
     .                     min_data_spacing*cf_inter
 290          format('**WARNING** Sampling interval does not match ',
     .               'at ',a4,' Original cf_inter of ',i4,
     .               ' s reset to ',i4, ' s')
              orig_sampling(j) = min_data_spacing*cf_inter
          end if

          close(100+j, iostat=ierr) 
          serr = max(serr,ierr)
          call report_error('IOSTAT',ierr,'clos',cfiles(j),0,
     .                      'read_cfdata')
          if( num_nolcaut.gt.0 ) then
              write(report_line,295) num_nolcaut, cf_codes(j)
 295          format('LC_AUTCLN: ',i6,' data deleted at ',a,
     .               ' due to 1/2 L2 phase or no P2 range')
              write(*,'(a)') report_line(1:trimlen(report_line))
              call report_stat('warning','autcln','read_cfdata',
     .                              ' ',report_line,0)
          end if  
      end do
 
****  Thats all
      write(*,300) actual_max_chan, num_chan, num_obs
 300  format(' Maximum number of channels used was ',i3,' Using ',i3,
     .       ' in current run.',/,
     .       ' There are ',i8,' one-way observations in these files')
      if( actual_max_chan.gt.num_chan ) then
          call report_stat('fatal','autcln','read_cfdata',' ',
     .                      'Too many channels in receiver',0)
      end if
      
      return
      end
 
CTITLE PROCESS_RNG
 
      subroutine process_rng(L1r_phs_cse, L2r_phs_cse,
     .    L1r_rng_cse, L2r_rng_cse, L1_cyc_cse, L2_cyc_cse,
     .    ctol_cse, data_flag_cse, params_cse, par_flag_cse  )
 
      implicit none
 
*     Routine to extract from the range and phase residuals the
*     station and satellite clocks and to compute the numbers
*     of carrier phase cycles by epoch.
 
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
 
*   iter       - Number of iterations to get convergence on clock
*               - statistics and values
*   i           - Epoch loop counter
*   trimlen     - Length of string
 
 
      integer*4 iter, i, j,  trimlen

*   dclk_next(max_cfiles+max_gsvs) - Estimate of the change in the clock for
*                 next epoch (initally based on the cfile estimates
*                 of rclock)

      real*8 dclk_next(max_cfiles+max_gsvs)
 
*   converged   - Logical indicating that clock solution converged
 
      logical converged
 
****  set the data mask to be used for edting the range data
      call set_rng_mask( rng_mask )
 
****  Keep the clock estimation until we seem converged on
*     statistics and values
      converged = .false.
      iter = 0
 
      do while ( .not.converged .and. iter.lt.max_rclk_iter )
          iter = iter + 1
*****     Start by initializing the reciever clocks with the values
*         obtained from the Cfile or from the previous iteration.
 
          call init_clk( iter )
          call init_rms( sum_rng_var, sum_rng_num, num_cfiles)
          
          do j = 1, num_param
             dclk_next(j) = 0.d0
          end do
 
*         Start looping over all epochs of data

          do i = 1, num_ep
          
* MOD TAH 990917: Get the apriori values for the clock from the
*            parameter estimates.  Now pre-loaded during read
             do j = 1, num_param
                apr_clk_val(j) = params_cse(j,i) - dclk_next(j)
             end do          
 
*             Set up the normal equations
* MOD TAH 990917: Run init_norm with clock known since we now get
*             values from rclock and svcoff in the cfiles.
              call init_norm(i, .true.)
 
*             Pre-scan the residuals to see if we have any jumps in the
*             clocks
C             write(*,7999) i,(params_cse(j,i), j=1,num_param)
C7999         format('PARAM ',i5,50(F11.1,1x))
C             write(*,8000) i, (dclk_next(j),j=1,num_cfiles)
C8000         format('RCK cyc',i5,20(F11.1,1x))
              call prescan_clk(i, L1r_rng_cse(1,1,i),
     .                        L2r_rng_cse(1,1,i),
     .                        ctol_cse(1,1,i), data_flag_cse(1,1,i),
     .                        dclk_next, par_flag_cse(1,i) )

*             Before updating the clock estimates get the change to 
*             next epoch
              if( i.lt.num_ep ) then
                  do j = 1, num_cfiles

*                    Only set the jump to the next epoch if we are
*                    past the first valid obervation at this site.
                     if( i.ge. apr_clk_epoch(j) ) then
                         dclk_next(j) = params_cse(j,i+1) - 
     .                                  params_cse(j,i)
                     else
                        dclk_next(j) = 0.d0
                     end if
                  end do
              end if
 
*             Now get the new estimates at this epoch and sum the
*             range residual rms estimates.
              call est_clk_rng(i, L1r_rng_cse(1,1,i),
     .                        L2r_rng_cse(1,1,i),
     .                        ctol_cse(1,1,i), data_flag_cse(1,1,i),
     .                        params_cse(1,i), par_flag_cse(1,i))
C             write(*,8010) i, 
C    .               (sqrt(sum_rng_var(j)/(max(1,sum_rng_num(j))))/
C    .                fClk*vel_light,
C    .                j=1,num_cfiles)
C8010         format('RMS',i5,20(F7.1,1x))
 
          end do
 
*****     Now scan the clock parameter estimates and see how many
*         jumps of one millisec there are.  Theses millisecond jumps
*         be removed from the range residuals here.  When the Gobs file
*         written they will be recorded as range cycle slips.  (By using
*         as a range slip there will be no phase shift associated with
*         them
         
          call get_rng_resets(L1r_rng_cse, L2r_rng_cse, data_flag_cse,
     .                     params_cse, par_flag_cse, ctol_cse )
 
*****     Update the initial clock estimates for the next iteration
          call update_init_clk( params_cse, par_flag_cse)
 
*****     Compute the new statistics and see if we have converged
          call update_stats(params_cse, par_flag_cse, converged)
 
*****     Tell user what is happening
          call report_clk_rng( iter, converged, 6 )
*                 ! Until converged (or outof iterations)
      end do
 
      if( trimlen(rng_clk_root).gt.0 ) then
          call write_clk(rng_clk_root, params_cse, par_flag_cse)
      end if
      
****  Thats all
      return
      end
 
CTITLE CHECK_ELEV 

      subroutine check_elev( ns, elev, azimuth, data_flag)
 
      implicit none

*     Routine to check if a data point should be edited due
*     to elevation angle.

* INCLUDES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
* elev, azimuth -- Elevation and azimuth (rads)

      real*4 elev, azimuth

* ns            -- site number
* data_flag     -- Data flag.  Bits 24 and 27 might be set

      integer*4 ns, data_flag

* LOCAL VARIABLES

* i             -- Loop counter
      integer*4 i, na

* min_elev      -- Estimate of minimum elevation to use (rad)
      real*4 min_elev

****  See if the az_mask has been given for this site
      min_elev = -1
      na = num_azmask(ns)
      if( num_azmask(ns).gt.0 ) then

*         Scan along Azimuth values given and see if we fall
*         inside range
          do i = 1, num_azmask(ns)-1
             
             if( azimuth.ge.az_mask(1,i,ns) .and.
     .           azimuth.le.az_mask(1,i+1,ns) ) then
                 min_elev = az_mask(2,i,ns)
             end if
          end do
      end if

*     If min_elev is still -1 then the azimuth is not included
*     range, so use the station default value
      if( min_elev.eq.-1 ) then
          min_elev = site_celev(ns)
      end if

*     See if the observation is below this value
      if( elev.lt. min_elev ) call sbit(data_flag, 24,1)

*     Now see if we should set the output elevation cutoff.  This is
*     not needed anymore becuase solve can change the cutoff
      if( num_azmask(ns).eq.0 .and. min_elev.lt.site_oelev(ns) ) then
         if( elev.lt. site_oelev(ns) ) call sbit(data_flag, 27,1)
      end if
      if( num_azmask(ns).gt.0 ) then
         if( elev.lt. min_elev ) call sbit(data_flag, 27,1)
      end if

****  Thats all
      return
      end

