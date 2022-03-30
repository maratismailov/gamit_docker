CTITLE UPDATE_CFILES
 
      subroutine update_cfiles(L1r_phs_cse, L2r_phs_cse, 
     .              L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse,
     .              ctol_cse, data_flag_cse,
     .              params_cse, par_flag_cse, pf_dphs_cse )

      implicit none
 
*     This routine will open and read the cfiles in the ctogobs run ad
*     write new cfiles with (initally) the cfile series incremented and
*     the number of cycles updated. Any new bias flags will also be
*     written to the new files.
 
*     In this version, the number of partials to set to zero to speed
*     up the reading and writting of files (and to save space).
 
* INCLUDES FILES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
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
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
*   pf_dphs_cse(num_chan, num_cfiles, num_ep)   - L1 adjustments to the
*                     omc to account for postfit residuals (L1 cycles)
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep),
     .    pf_dphs_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)

* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   trimlen - Length of string
*   i,j,k       - Loop counters
 
      integer*4 ierr, trimlen, i
c     integer*4 lnblnk --this not currently used

      integer*4 fmprename 

*   new_name    - Name of the new-cfile
*   old_name    - Old name or renamed cfile if we are replacing
*                 cfile
*   gen_name    - Name of cfile while it is being generated (used
*                 when overwritting the cfile.
 
      character*256 new_name, old_name, gen_name

*   replace_cfile - Logical which indicates that we are to replace
*                   the cfile i.e., cfile renamed, new one made and
*                   then is all is OK old one deleted.

      logical replace_cfile

***** See if we need new cfiles
      if( trimlen(new_cfile_series).eq.0 ) RETURN

***** Set the phase masks

      call set_phs_mask(phs_mask, phs_bias_mask ) 

***** See if user wants the first bias flag removed.
      if( remove_first_bias ) then
          call remove_fb(ctol_cse, data_flag_cse)
      end if

****  Loop over the cfiles.  Start by opening and creating new cfiles
      do i = 1, num_cfiles

*         Get the code and series number
          call gen_new_cfname(cfiles(i),new_cfile_series,new_name,
     .                        replace_cfile)

*         If we are replacing, rename the old cfile
          if( replace_cfile ) then
              gen_name = new_name(1:trimlen(new_name)) // '.new'
          else       
*             Just assign the cfile name
              gen_name = new_name
          end if
          old_name = cfiles(i)
          call open_cf(100+i,old_name , ierr)
 
*         If there has not been an error try to create new cfile
          write(*,120) cfiles(i)(1:trimlen(cfiles(i))),
     .                new_name(1:trimlen(new_name))
 120      format(' Updating ',a,' to ',a)
          if( ierr.eq.0 ) then
              call create_cf(200,gen_name,ierr)
          end if
 
****      If still no error then copy header information over to new
*         cfile
 
          if( ierr.eq.0 ) then
              call copy_cfheads( 100+i, 200, i, -1, ierr)
          end if
 
****      If still no error then copy over data
          if( ierr.eq.0 ) then
              call copy_cfdata( 100+i,200, i, L1r_phs_cse, L2r_phs_cse, 
     .              L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .              ctol_cse, data_flag_cse,
     .              params_cse, par_flag_cse, pf_dphs_cse, ierr )
          end if

****      If we replaced the cfile and there is no error, then
*         remobe the old file 
          close(200)
          if( ierr.eq.0 .and. replace_cfile ) then
              close(100+i,status='delete',iostat=ierr) 
              call report_error('IOSTAT',ierr,'delet',old_name,0,
     .                          'update_cfiles')
              ierr = fmprename(gen_name,new_name,'K')
c This doesn't work because of blanks in the filename...
c             ierr = rename(gen_name,new_name) 
c This did work until g77 0.5.20 made rename a subroutine...
c              ierr = rename(gen_name(:lnblnk(gen_name))
c     .                  ,new_name(:lnblnk(new_name)))
              call report_error('IOSTAT',ierr,'renam',gen_name,0,
     .                          'update_cfiles')
          else
              close(100+i,iostat=ierr)
          end if
      end do
 
****  Thats all
      return
      end
 
CTITLE GEN_NEW_CFNAME
 
      subroutine gen_new_cfname(name,series,new_name, replace)

      implicit none
 
*     Routine to generate the name of the new cfiles bases on the name
*     of the old one and the new series letter.  A series of '+' will
*     increment the cfile file series
 
* PASSED VARIABLES
 
*   name        - Name of the old cfile
*   series      - Series for new cfiles (+ mean increment the
*               - series)
*   new_name    - Name of the new cfile
 
      character*(*) name, series, new_name

*   replace     - True if the old series equals the new series.
*                 (returned by this routine)

      logical replace
 
* LOCAL VARAIBLES
 
*   len_name    - Name of the old cfile
*   trimlen     - Length of string
*   i           - Loop counter
 
      integer*4 len_name, trimlen, i
 
*   old_series  - Old cfile seris
*   new_series  - new cfile series based on the series type
*               - passed.
 
      character*4 old_series, new_series
 
****  Start by finding the last '/'in the name
      len_name = trimlen(name)
      i = len_name
      do while (i.gt.1 .and. name(i:i).ne.'/')
          i = i - 1
      end do
 
*     Adjust i so that it points to the 'c'
      if( name(i:i).eq.'/' ) i = i+1
 
 
*     Now generate the new name
      old_series = name(i+5:i+5)
      if( series(1:1).eq.'+' ) then
          if( ichar(old_series(1:1)).ge.ichar('a') ) then
              new_series = char(ichar(old_series(1:1))+1)
          else
              new_series = 'a'
          end if
      else
          new_series = series
      end if

*     See if we have explicitlty said to use old name.
      if( series(1:1).eq.'.' ) then
          new_series = old_series
      end if
 
      new_name = name(1:i+4) // new_series(1:1) // name(i+6:)

*     See if we are overwritting old c-files
      if( new_series(1:1).eq.old_series(1:1) ) then
          replace = .true.
      else
          replace = .false.
      end if
 
****  Thats all
      return
      end
 
CTITLE  COPY_CFHEADS
 
      subroutine copy_cfheads( inu, outu, cf, nparts, ierr )

      implicit none
 
*     This routine will copy the header records of the cfile on
*     inu to the one on outu whle optionally changing the number
*     partials written to the output file.
 
 
* INCLUDES FILES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   inu, outu   - Input and output unit numbers
*   cf          - Cfile number
*   nparts      - variable to specfy number of partials to
*               - to be written.  If it is zero then none are
*               - written.  Any other value will cause them all
*               - all be written
*   ierr        - IOSTAT error either on reading or writing.
 
      integer*4 inu, outu, cf, nparts, ierr
 
* LOCAL VARIABLES
 
*   date(5)     - Date that file written
*   k           - temporay storage for number of text records
*   j           - Loop counter
 
      integer*4 date(5), k, j, trimlen
 
*   sectag      - seconds tage of time
 
      real*8 sectag

      character*12 full_ver, hsver
 
*     Read the first cfile block
      call read_cf1( inu, 'ALL', ierr )
 
*     Update the history entries on the file
      if( cf_ntext.lt.cf_maxtxt) then
          k = cf_ntext + 1
          cf_ntext = k
          call systime( date, sectag)
          full_ver = hsver(ctogobs_version)
          write(cf_text(k),100) caprog_name, 
     .             full_ver(1:trimlen(full_ver)), date, sectag
 100      format(a8,': vers. ',a,' run ',i4,'/',i2,'/',i2,
     .            1x,i2,':',i2,1x,F5.2)
      end if

      if( cf_ntext.lt.cf_maxtxt) then
          k = cf_ntext + 1
          cf_ntext = k
          write(cf_text(k),120) caprog_name, (ion_rw_tol(j,cf),j=1,3), 
     .                           dd_wl_tol, dd_lc_tol
 120      format(a8,': Tol: LG ',3f6.2,' WL ',3f6.2,' LC ',3F6.2) 
      end if

      if( cf_ntext.lt.cf_maxtxt) then
          k = cf_ntext + 1
          cf_ntext = k
          write(cf_text(k),140) caprog_name, max_wl_ret, max_dd_ret,
     .           max_lg_use, tol_one_way_fix
 140      format(a8,': Returns: WL ',I4,' LC ',I4,' LG ',I4,
     .           ' OW ',I4,' sec')
      end if

      if( cf_ntext.lt.cf_maxtxt) then
          k = cf_ntext + 1
          cf_ntext = k
          write(cf_text(k),160) caprog_name, dchi2_ratio, dchi2_min_val,
     .           dchi2_max_sep
 160      format(a8,': ChiTol: Ratio ',f5.2,' Min ',f5.2,
     .           ' Max ', f7.1,' cyc')
      end if
      if( cf_ntext.lt.cf_maxtxt .and. np_size.gt.0 .and.
     .    np_start.gt.0 ) then
          k = cf_ntext + 1
          cf_ntext = k
          write(cf_text(k),180) caprog_name, np_size, np_start
 180      format(a8,': Data Normal pointed with ',i3,' epochs,',
     .           ' Starting epoch ',i4)
      end if

      if( cf_ntext.lt.cf_maxtxt .and. apply_phs_clk ) then 
          k = cf_ntext + 1
          cf_ntext = k
          write(cf_text(k),190) caprog_name
 190      format(a8,': Phase clocks applied to residuals')
      end if

 
      call write_cf1(outu, 'ALL',ierr)
 
****  Now read the second type record
      call read_cf2(inu, 'ALL', ierr)
 
*     If nparts was passed as zero, update the record here sso
*     that partals will not be written
      if( nparts.eq.0 ) cf_npart = 0
      
*     Check to see if we are outputing normal points.  If we
*     are then say the original sampling interval was the
*     normal point spacing
      if( np_size.gt.0 .and. np_start.gt.0 ) then
          cf_ircint = sampling*np_size 
      end if           

***** Update the elevation cutoff to be consistent with the
*     output cutoff.
      cf_elvcut = site_oelev(cf)*180.d0/pi

* MOD TAH 190722: Check the cf_lambda values before writing.  
*     Unless L1 only is being processed, update any L1 only
*     factors to L1/L2.  Including L1 only in dual frequency
*     results, even when all the data is deleted causes large 
*     (10^2) nrms values in solve.
      if( nol1only ) then  ! No L1 only data being processed
          do k = 1, cf_nsat
             if( k.eq.0 ) then
                write(*,220)  cf_lambda(1,:), cf_dattyp
 220            format('Updating  cf_lambda. cf_dattyp from ',
     .                 5I4,' TYPE ',5I3)
             endif           
             if( cf_lambda(k,2).ne.-1 ) then
                 cf_lambda(k,2) = -1
                 cf_dattyp(2) = 2
                 cf_lambda(k,4) = +1
                 cf_dattyp(4) = 4
             endif
          enddo
      end if

* MOD TAH 200511: Update the frequencies so that network values
*     are saved.  Avoids problems with first or last c-file 
*     zero frequency due to no data (possible with GNSS processing
*     rate for GPS).
      do k = 1, cf_nsat
         cf_fL1(k) = fL1(k)
         cf_fL2(k) = fL2(k)
      enddo

*     Get the correct number of data types for this data.
      call write_cf2(outu, 'ALL', ierr)
 
****  Now do the third record
      call read_cf3(inu,'ALL',ierr)
      call write_cf3(outu,'ALL' ,ierr)
 
****  Thats al
      return
      end
 
CTITLE COPY_CFDATA
 
      subroutine copy_cfdata( inu, outu, cf, L1r_phs_cse, L2r_phs_cse, 
     .              L1_cyc_cse, L2_cyc_cse,  L1r_rng_cse, L2r_rng_cse,
     .              ctol_cse, data_flag_cse,
     .              params_cse, par_flag_cse, pf_dphs_cse, ierr )

      implicit none
 
*     Routine to loop through the cfile reading the data and
*     writing it with the number of cycles updated.
 
 
* INCLUDES FILES
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   inu, outu       - Input and output units
*   cf              - cfile number
*   ierr            - IOSTAT error
 
      integer*4 inu, outu, cf, ierr
 
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
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
*   pf_dphs_cse(num_chan, num_cfiles, num_ep)   - L1 adjustments to the
*                     omc to account for postfit residuals (L1 cycles) 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep),
     .    pf_dphs_cse(num_chan, num_cfiles, num_ep),
     .    params_cse(num_param, num_ep)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters                              
*   lv      - satellite number
*   elev_dist(20) - distribution of data with
*             elevation angle (5 deg bins)
*   bin     - Bin number for particular elevation.
*   trimlen - Length of string.
*   start_ep, tot_npep  - Start and total number of normal
*             point epochs.

*   act_msat  - Actual number of satellites observed.  Used
*               when outputing normal points and cf_msat is
*               set to zero for all non-normal point epochs.

*   cs       - Parameter number for satellite clock estimate.
 
      integer*4 i,j, lv, elev_dist(20), bin, prn, trimlen,
     .          start_ep, act_msat, cs

*   push_bias(max_gsvs) - Logical to indicate to df_to_cerrfl that
*               we have shown as marginal a data point with a bias
*               flag and that the bias flag should be pushed to next
*               good data point
*   kbit    - Test to see if a bit is set.
*   data_OK - Logical function returns true if data does not
*             have any bits from the mask set.
*   good_bf - Logical function that returns true if data point
*             is good and has a bias flag. (mask used should not
*             have bias flag bits (31 and 32) set.
*   write_np - Set true if we are writting out normal points.

      logical push_bias(max_gsvs) , kbit, data_ok, good_bf,
     .        write_np

*   drng    - Difference between L1 range used and value in cfile
*             If they differ then we have added millisecond clock
*             jumps.
*   ms_dclk_sec - Change in clock time in seconds
*   ms_dclk_L1cyc - Change in clock time in L1 cycles
*   ms_dclk_L2cyc - Change in clock time in L2 cycles
*   ms_dclk_meters - Change in clock time in meters

*   dnp_omc(4)  - Difference in the observed omc and normal pointed
*                 omc.  (used to make the new obervable values)
*   npep_est    - Estimate of normal point epoch.  (If this is integer
*                 then we are on normal point).

      real*8 drng, ms_dclk_sec, ms_dclk_L1cyc, ms_dclk_L2cyc, 
     .       ms_dclk_meters, cf_minel, dnp_omc(4), npep_est 

* MOD TAH 180307: Added ion_delay to applied to the omc values for
*     use in solve (for non-Glonass systems, should not affect results)
      real*8 ion_delay(4)  ! Ion delays for L1 L2 phase and P1 P2 range.

****  Initialize the push biases and elev distribution
      do i = 1, num_sat
         push_bias(i) = .false.
      end do
      do i = 1,20
         elev_dist(i) = 0
      end do

****  Set up for outputing normal points
      if( np_size.gt.0 .and. np_start.gt.0 ) then
          write_np = .true.
          start_ep = abs(np_start)
      else
          write_np = .false.
          start_ep = 1
      end if

****  Loop over all the epochs in the data
      cf_minel = 90.d0

      i = 0
*     Do as a while loop so that we can exit if there is an error. 
      do while( i.lt. num_ep )
          i = i + 1
          call read_cf4( inu, 'ALL', ierr)

*         See if this is a normal point epoch if we are forming
*         normal points
          act_msat = cf_msat
          if( write_np ) then
              npep_est = float(i-start_ep-int(np_size/2))/np_size 
              if( npep_est.ge.0.d0 .and.
     .            abs(npep_est-int(npep_est)).lt.1.d-6 ) then
*                 This is an epoch at which we should output a normal
*                 point.
                  cf_msat = act_msat
              else
*                 This is not an epoch for output of normal points
                  cf_msat = 0
              end if
          end if
 
*         Now loop over the channels
          do j = 1, act_msat
              call read_cf5( inu, 'ALL', ierr)
                                                    
              lv = ctol_cse(j,cf,i)

*             Now update the number of cycles
*             If we are writing normal points then compute the
*             new observable
              if( write_np .and. cf_msat.ne.0 ) then
                  dnp_omc(1) = L1r_phs_cse(j,cf,i) - cf_omcs(1)
                  dnp_omc(2) = L2r_phs_cse(j,cf,i) - cf_omcs(2)
                  dnp_omc(3) = L1r_rng_cse(j,cf,i) - cf_omcs(3)
                  dnp_omc(4) = L2r_rng_cse(j,cf,i) - cf_omcs(4)

                  cf_omcs(1) = L1r_phs_cse(j,cf,i)
                  cf_omcs(2) = L2r_phs_cse(j,cf,i)
                  cf_omcs(3) = L1r_rng_cse(j,cf,i)
                  cf_omcs(4) = L2r_rng_cse(j,cf,i)
                  
*                 If we used postfit residuals; add back the adjustment
*                 so that normal points referr to original aprioris.
*                 NOTE: We do not change dnp_omc because this is added
*                 back to the raw measurements.
                  if( use_postfit ) then
                      cf_omcs(1) = cf_omcs(1) + pf_dphs_cse(j,cf,i)
                      cf_omcs(2) = cf_omcs(2) + pf_dphs_cse(j,cf,i)*
     .                                          fL2(lv)/fL1(lv)
                      cf_omcs(3) = cf_omcs(3) + pf_dphs_cse(j,cf,i)
                      cf_omcs(4) = cf_omcs(4) + pf_dphs_cse(j,cf,i)*
     .                                          fL2(lv)/fL1(lv)
                  end if                

*                 For the range make sure that we don't have any
*                 millisecond jumps.
                  if( remove_ms_jump.and.
     .                abs(dnp_omc(3)).gt.100.d0 ) then
*                     We have some millisecond jumps.  Since these
*                     have been removed in range residuals, add back
*                     now so that the epoch will be computed correctly.
                      ms_dclk_sec  = nint(dnp_omc(3)/fL1(lv)*1.d3)*1.d-3
                      ms_dclk_meters = ms_dclk_sec*vel_light
                      ms_dclk_L1cyc = ms_dclk_sec * fL1(lv)
                      ms_dclk_L2cyc = ms_dclk_sec * fL2(lv)
                      cf_omcs(3) = cf_omcs(3) - ms_dclk_L1cyc 
                      cf_omcs(4) = cf_omcs(4) - ms_dclk_L2cyc 
                   else
                      ms_dclk_L1cyc = 0.d0
                      ms_dclk_L2cyc = 0.d0
                      ms_dclk_meters = 0.d0
                   end if

*****              Now replace the obervables themselves with the normal
*                  point value.  (This will only be checked by running 
*                  ctox -> model.)
                   cf_obsv(1) = cf_obsv(1)+ dnp_omc(1)
                   cf_obsv(2) = cf_obsv(2)+ dnp_omc(2)
                   cf_obsv(3) = cf_obsv(3)+ dnp_omc(3)/fL1(lv)*vel_light
     .                                    - ms_dclk_meters
                   cf_obsv(4) = cf_obsv(4)+ dnp_omc(4)/fL2(lv)*vel_light 
     .                                    - ms_dclk_meters
              end if


*             Now complete the normal processing of the data
C MOD TAH 190704: use re-scaled versions
              cf_omcs(1) = cf_omcs(1) + L1_cyc_cse(j,cf,i)
              cf_omcs(2) = cf_omcs(2) + L2_cyc_cse(j,cf,i)
C MOD TAH 190710: Removed change.  Fixed the GLONASS frequencies instead
C             cf_omcs(1) = L1r_phs_cse(j,cf,i) + L1_cyc_cse(j,cf,i)
C             cf_omcs(2) = L2r_phs_cse(j,cf,i) + L2_cyc_cse(j,cf,i)
C             cf_omcs(3) = L1r_rng_cse(j,cf,i)
C             cf_omcs(4) = L2r_rng_cse(j,cf,i)
 
              cf_obsv(1) = cf_obsv(1) + L1_cyc_cse(j,cf,i)
              cf_obsv(2) = cf_obsv(2) + L2_cyc_cse(j,cf,i)

*             Now check to see if we have changed the range
*             to account for millisecond jumps.  Do this only for the
*             first satellite (need to write the type 4 record out
*             with the update clock value.
              if( abs(cf_omcs(3)-L1r_rng_cse(j,cf,i)).gt.1.d0 .and.
     .            remove_ms_jump ) then

*                 There appears to be millisecond jumps.  Compute
*                 milliseconds
                  drng = L1r_rng_cse(j,cf,i) - cf_omcs(3)
                  ms_dclk_sec  = nint(drng/fL1(lv)*1.d3) * 1.d-3
                  ms_dclk_meters = ms_dclk_sec*vel_light
                  ms_dclk_L1cyc = ms_dclk_sec * fL1(lv)
                  ms_dclk_L2cyc = ms_dclk_sec * fL2(lv)

*                 correct the range 
                  cf_omcs(3) = cf_omcs(3) + ms_dclk_L1cyc
                  cf_obsv(3) = cf_obsv(3) + ms_dclk_meters

*                 See if L2 range needs to be done.
                  if( L2r_rng_cse(j,cf,i).ne.0.d0 ) then
                      cf_omcs(4) = cf_omcs(4) + ms_dclk_L2cyc
                      cf_obsv(4) = cf_obsv(4) + ms_dclk_meters
                  end if

*                 Now if we have not changed the epoch, then
*                 change now
                  cf_rclock = cf_rclock + ms_dclk_sec
                  cf_sod    = cf_sod + ms_dclk_sec
              end if
              
* MOD TAH 970107: Check to see if we are applying the phase clock
*             to the data.  There are some problems with gamit here.
*             we should update rclock, but this value is not used
*             to correct the phase, so it will no good.  Instead just
*             apply clock corrections to omcs values
              if( apply_phs_clk ) then
                  cs = num_cfiles + ctol_cse(j,cf,i)
                  cf_omcs(1) = cf_omcs(1) - (params_cse(cf,i) -
     .                                       params_cse(cs,i))           
                  cf_omcs(2) = cf_omcs(2) - (params_cse(cf,i) -
     .                               params_cse(cs,i))*(fL2(lv)/fL1(lv))
                  cf_omcs(3) = cf_omcs(3) - (params_cse(cf,i) -
     .                                       params_cse(cs,i))           
                  cf_omcs(4) = cf_omcs(4) - (params_cse(cf,i) -
     .                               params_cse(cs,i))*(fL2(lv)/fL1(lv)) 
              end if

* MOD TAH 180307: Use the phase data to compute the ion-delay for phase
*             and range and apply the values to cf_omcs and cf_obsv
*             The channel numner lv is need for GLONASS processing.
              if( app_ion ) then 
                  call comp_ion( ion_delay, cf_omcs, cf_obsv, lv )
              endif

*             If this is first channel write type 4 record
              if( j.eq.1 ) then
                  call write_cf4( outu, 'ALL',ierr )
              end if

****          Now see if we are changing the elevation cutoff
*             for the output.  The "mask" tells if we are checking
*             the cfile min_elev and only when we are using the
*             min_ctog_elev do we allow this change on output.
              if( .not.kbit(phs_mask,15) ) then 
*                 GAMIT elevation cutoff not being used, so see
*                 current measurement is below the output cutoff
*                 NOTE: df_to_cerrfl will handle the bias flag pushing.
                  if( cf_elev(1).lt.site_oelev(cf) ) then

****                  See if data currently marked "good".  If it is and
*                     it has a bias flag then set that we must push_bias
*                     to the next good point.
                      
                      if( good_bf(data_flag_cse(j,cf,i),0,
     .                            phs_mask) ) then
                            push_bias(prntol(cf_iprn)) = .true.
                      end if
                      call sbit(data_flag_cse(j,cf,i),24,1)
                  end if
              end if

****          Convert the gobs data flag back to a solve errfl
              call df_to_cerrfl(data_flag_cse(j,cf,i), cf_ierfl,
     .              push_bias(prntol(cf_iprn)), phs_mask )

* MOD TAH 050410: Save the autcln data flag
              cf_data_flag = data_flag_cse(j,cf,i)

*             accumulate elevation distriubution.  Only if we will
*             really be outputing these values.
              if( data_OK(data_flag_cse(j,cf,i),0,phs_mask) .and.
     .            cf_msat.gt.0 ) then
                  bin = cf_elev(1)*180/pi/5 + 1
                  if( cf_elev(1)*180./pi.lt.cf_minel ) then
                      cf_minel = cf_elev(1)*180./pi
                  end if
                  elev_dist(bin) = elev_dist(bin) + 1
              end if

*             Only output if we have satellites
              if( cf_msat.gt.0 ) then
                   call write_cf5( outu,'ALL',ierr )
              end if

*             See if we should write the single difference files
*             Generated i root to file name has been given and
*             this is one of of the selected sites (and greater
*             than site 1 since we difference against this site.)
              if( trimlen(sng_diff_root).gt.0 .and.
     .            kbit(phs_res_sites,cf) .and. cf.gt.1 .and. 
     .            cf_msat.gt.0 ) then
                   prn = cf_iprn
                   call write_sd_file(i, cf, prn, 
     .                L1r_phs_cse, L2r_phs_cse,
     .                L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse,
     .                ctol_cse, data_flag_cse,
     .                params_cse, par_flag_cse )
              end if

*             See if we should write the single difference files

          end do

*         Make sure we write the type 4 record when no satellites
          if( act_msat.eq.0 ) then 
              call write_cf4( outu, 'ALL',ierr )
          end if

****      See if an error has ocurred
          if( ierr.ne.0 ) i = num_ep
      end do

****  Now write out the eleveation angle distribution
      if( cf.eq.1 ) then
*         write the header line
          write(uns,200) ((j-1)*5,j*5, j =1,18)
 200      format('ELEVATION ANGLE HISTOGRAM',/,
     .           'SITE',18(1x,i2,'-',i2),' Min (dg)')
      end if
      write(uns,220) cf_codes(cf),(elev_dist(i),i=1,18), cf_minel
 220  format(A4,18(1x,I4,1x),1x,F5.2)
 
****  Thats all
      return
      end

CTITLE DF_TO_CERRFL

      subroutine df_to_cerrfl(data_flag, cf_ierfl, push_bias, mask)

      implicit none

*     Routine to convert the gobs data flag back to a cfile errflag
*     Becuase we may have data flaged bad with a bias flag, we
*     need to wait for first good data after flag to set the bias
*     flag in GAMIT

* PASSED VARIABLES

* data_flag - GOBS data flag bit mapped.
* cf_ierfl  - cfile error flag (may be changed)
* mask      - Mask to see if gobs data is good (ignores the bias
*             flags)

      integer*4 data_flag, mask 

      integer*2 cf_ierfl 

* push_bias - True if we are waiting for a good data point to
*      show a bias flag

      logical push_bias

* LOCAL VARIABLES	

* data_OK - Gobs functon to see if data is OK.
* kbit    - Tests if bit is set

      logical data_OK, kbit


*     Convert the gobs data flag back into the closest GAMIT
*     error flag.
      cf_ierfl = 0
      if( .not.data_OK(data_flag, 0,mask) ) then
          cf_ierfl = -1 

*         If we were using GAMIT elevation cutoff then
*         check this bit
          if( kbit(mask,15) .and.
     .        kbit(data_flag,15)) cf_ierfl = 4
*         See if we were using cview edits
          if( kbit(mask,14) .and.
     .        kbit(data_flag,14)) cf_ierfl = -1	
          if( kbit(data_flag,24)) cf_ierfl = 4
          if( kbit(data_flag,2) ) cf_ierfl = 3
          if( kbit(data_flag,4) ) cf_ierfl = 2
      end if
        
*     See if should output a bias flag. 
      if( data_OK(data_flag, 0,mask) .and.
     .   (kbit(data_flag,31) .or. kbit(data_flag,32)) ) then
         cf_ierfl = +10
      end if

*     See if we really can tell cview that there is a bias flag
      if( cf_ierfl.eq.0 .and. push_bias ) then
          push_bias = .false.
          cf_ierfl = +10
      end if

      if( cf_ierfl.eq.10 ) call sbit(data_flag,31,1)

***** Thats all
      return
      end
