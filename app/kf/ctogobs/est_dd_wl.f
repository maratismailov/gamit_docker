CTITLE est_dd_wl
 
      subroutine est_dd_wl(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, bf_type_cse, azel_cse) 

      implicit none
 
*     This routine will use the bf_table from flat_dd to estimate
*     number of cycles in the last bias flags for each site and
*     satellite.  MW-WL biases are estimated by station and satellite with the
*     satellite offsests depending on dcb type for the receiver.  
*     Double difference mabuguities are formed based on baseline length and 
*     sigma of double diffs.  

*     Sequence is:
*     Estimate one-way ambuguities from LC, EX-WL (L2-L1), and MW-WL (different
*        data sigmas and baseline weighting are used for the three estimates).
*     For MW-WL estimate corrections to the DCBs (iterated with cycles removed).
*     Approximately resolve L1 and L2 cycles using LC and MW-WL and apply to
*        phase data.  Remove integers from estimates.
*     Form the DD biases to be passed to solve and assess which ones can be reliably
*        set to integer values.


* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep)
     .,    ctol_cse(num_chan, num_cfiles, num_ep)
     .,    bf_type_cse(num_chan, num_cfiles, num_ep)

*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation
       real*4 azel_cse(2,num_chan, num_cfiles, num_ep)
 
* LOCAL VARIABLES
      integer*4 shortest_pair(2)  ! Two stations with the shortest 
                                  ! baseline length
     .,         llbf_sv     ! Satellite for long last bias flag at shortest
                            ! baseline pair of statons
     .,         llbf_dur(2) ! Duration of long last bias flag for pair
     .,         llbf_bf(2)  ! BF_TABLE entry of long last bias for pair

      logical all_resolved  ! Set true if all bias flags are resolved. Mainly
                            ! use to indicate too many bias flags
     .,       single_ref    ! Set true if single reference satellite can be found (visible
                            ! at all sites)

      character*256 stat    ! Status message saying which site and satellite to start
                            ! with for wide-lane bias fixing.

****  Start: Compute the baseline lengths betweeen all the sites
*     and return the pair of stations with the shortest length.  
      call get_blens( shortest_pair)

***** Remove data that has no double differences
      call scan_nodd( ctol_cse, data_flag_cse, bf_type_cse) 

****  Rescan the bias flags to make the bf_table (this is in case some
*     bias flags were removed in pf_remove_bias operation.  This operation may
*     move to before trim_oneways in ctogobs main.
      call get_biases( data_flag_cse, ctol_cse, all_resolved, 'REP' )

*     See which station and satellite we will choose as base.  The
*     measure here is the longest last bias flag.  Satellite and 
*     length is returned
      single_ref = .true.
      call find_long_last_bf( shortest_pair, llbf_sv, 
     .                        llbf_dur,llbf_bf, single_ref)

***** See if we found a satellite
      if( llbf_sv.eq.0 ) then
          single_ref = .false.
      endif
      if( .not.single_ref ) then
         call find_long_last_bf( shortest_pair, llbf_sv, 
     .                        llbf_dur,llbf_bf, single_ref)
      end if

*     See which is better (more quality tests could be add here
*     at the moment use just the length)
      wl_ref_svs  = llbf_sv
      if( llbf_dur(1).gt.llbf_dur(2) ) then
         wl_ref_site = shortest_pair(1)
         wl_ref_dur  = llbf_dur(1)
         wl_ref_bf   = llbf_bf(1)
      else
         wl_ref_site = shortest_pair(2)
         wl_ref_dur  = llbf_dur(2)
         wl_ref_bf   = llbf_bf(2)
      endif
      write(stat,120) cf_codes(wl_ref_site), prn_list(wl_ref_svs),
     .                wl_ref_bf, wl_ref_dur
 120  format('WL DD Reference site ',a4,' PRN_',i2.2,' BF ',i4,
     .       ' Duration ',I4,' Epochs')
      call report_stat('status','autcln','Est_dd_wl',' ',stat,0)

****  Now generate estimates for LC, EX-WL (L2-L1) and MW-WL.
      call est_comb(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, azel_cse,'LC')
     
      call est_comb(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, azel_cse,'EXWL')

*     Do the MWWL case last because we will use the covariance matrix
*     from this solution when we form double differences
      call est_comb(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, azel_cse,'MWWL')


****  OK With these biases, compute the number of ambiquities to add
*     to the data.
      call update_dd_wl( L1_cyc_cse, L2_cyc_cse,  
     .    ctol_cse, data_flag_cse)

****  OK now set up final DD ambiguties
      call fin_gamit_dd(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse)


****  Thats all
      return
      end
 

CTILE UPDATE_DD_WL

      subroutine update_dd_wl( L1_cyc_cse, L2_cyc_cse, 
     .    ctol_cse, data_flag_cse) 

      implicit none

*     Routine to form update the number of cycles based on the 
*     direct estimates of the last bias flags.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES

*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)

*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep)


* LOCAL VARIABLES
      integer*4 ns, lv  ! Site and SVS
     .,         p  ! Parameter number
     .,         IMW  ! Integer MW-WL computed with corrected MW
     .,         dMW  ! Values arounf IMW (-1,0,1)
     .,         dernum  ! dMW value associated with best fit.
     .,         ltoc ! function converting SV list to channel
     .,         ch   ! channel number

      real*8 MW, MWC  ! MW-WL and its corrected value for site and SVS
     .,      dermin   ! Minimum deviation of N1 from integer
     .,      derraw(3)   ! Raw values of differences for dMW values
     .,      dN1, dN2 ! Offset of L1 and L2 cycles
     .,      LC       ! LC estimate
     .,      fN1      ! Estimate of number of L1 cycles needed

* MOD TAH 200617: Added frequnency ratios.
      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases. 
                 ! Estimated MW-WL with nominal frequencies will be 
                 ! integer (MW-WL)*fL1u/fL1

      character*4 mode  ! Use of CAL or RAW MW-WL values



***** Loop over each site and satellite and update the number of
*     cycles needed.  We have MW-WL estimate of L2-L1 cycles and 
*     we have the LC estimate:
*     LC = lcf1*N1 + lcf2*N2 and
*     MW = N2-N1  -> N2 = MW+N1
*     LC = lcf1*N1 + lcf2*(MW+N1)  -> N1 = (LC - lcf2*MW)/(lcf1+lcf2)

*     MWC is MW with the site and satellite corrections applied
*     Loop over station and satellites
      do ns = 1, num_cfiles
         do lv = 1, num_sat   

* MOD TAH 200617: Compute frequency ratios (fr1=fr2 for Glonass) 
            fr1 = fL1u(lv)/fL1(lv)
            fr2 = fL2u(lv)/fL2(lv)
                   
*           Make sure we have data
            p = (ns-1)*num_sat + lv
            if( ambnum_mw(p).gt.0 ) then 
*              Compute estimates of N1 and N2
               MW = ambest_mw(p)
               MWC = MW - (mwbias_site(ns)+
     .                     mwbias_svs(lv,site_dcb(ns)+1))
               LC = ambest_lc(p)

*              Try with correct value first the straigth estimates
               mode = 'CAL'
*              If this close to integer, use this value
               dermin = 1.e6
               dernum = 0.d0
* MOD TAH 200617: Convert so integer
               IMW = nint(MWC*fr1)
               do dMW = -1,1
                   fN1 = (LC-lcf2(lv)*(IMW+dMW))/(lcf1(lv)+lcf2(lv))
                   derraw(dMW+2) = fN1-nint(fN1*fr1)/fr1 
                   if( abs(derraw(dMW+2)).lt.dermin ) then
                       dermin = abs(derraw(dMW+2))
                       dernum = dMW
                   endif
               enddo

****           See which was best
               dMW = IMW+dernum
               dN1 = nint((LC-lcf2(lv)*dMW)/(lcf1(lv)+lcf2(lv))*fr1)/fr1
               dN2 = dN1 + dMW/fr1
               if( dN1.ne.0 .or. dN2.ne.0 ) then
*                  Update our saved values
                   ambest_lc(p) = ambest_lc(p) - 
     .                (lcf1(lv)*dN1 + lcf2(lv)*dN2)
                   ambest_mw(p) = ambest_mw(p) - dMW
                   ambest_ex(p) = ambest_ex(p) - dMW
                   
*                  Apply the corrections to the data
                   call update_cyc(-dN1,-dN2, ns, lv, 1, num_ep, 0,
     .                   L1_cyc_cse, L2_cyc_cse, 
     .                   ctol_cse, data_flag_cse, 'N' )
               endif
               write(*,220) cf_codes(ns),site_dcb(ns), prn_list(lv), 
     .                    dN1, dN2, dermin, ambest_mw(p), ambest_ex(p),
     .                    ambest_lc(p), MWC-dMW
 220           format('OWRES ',a4,1x,i1,1x,I2.2,' dN12 ',2F6.2,1x,
     .                ' DER ',F7.3,1x,' MW ',F7.3,' EX ',F7.3,
     .                ' LC ',F7.3,' MWC ',F7.3,' cyc')
             endif
          enddo
      end do

****  Thanks all
      return
      end



CTILE EST_COMB

      subroutine est_comb(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse, azel_cse,
     .    comb_type)

      implicit none

*     Routine to form direct estimate of for comb_type of
*     LC  -- LC combination
*     EXWL -- L2-L1 difference
*     MWWL -- MW wideline (L2-L1) cycles.  DCB Biases are also estimated
*        in the MWWL case
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES

*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)

*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
 
 
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)

*   azel_cse(2,num_chan, num_cfiles, num_ep)     - Azimuth and elevation

      real*4 azel_cse(2,num_chan, num_cfiles, num_ep)

      character*(*) comb_type ! combination LC, EXWL, MWWL      
* NOTE to g77 users: This construct does not seem to work in g77 when
* used in the concatenation at line 401 below though it represents
* a standard and safer declaration of the variable.  Since g77 (gcc <4.0)
* is no longer supported, we recommend that you upgrade to the currently
* supported compiler gfortran (gcc 4.2).  If g77 is used,  replace the 
* above line with 
*      character*4 comb_type    

* LOCAL VARIABLES

*   obs(5)  - MW-WL, L2-L1, L1, L2, LC residuals in cycles
*   wgh(5)     - Weight of measurement (inversely proporational to
*             1/sin(el)**2 and phase noise on satellite

      real*8 obs(5,4), wgh(5,4), tot_wgh

*   OK      - Logical: Set OK if data is OK
      logical OK(4)

      integer*4 ep    ! Epoch being extracted

      integer*4 ch    ! Channel number at this particular epoch
     .,      ltoc     ! Function to return channel number of specific
                      ! satellite list number
      integer*4 i,j,k,l, in,jn , ent, bl_ent
      integer*4 np(4) 
     .,   nump          ! Number of bias parameters

      real*8 dobs      ! Observation combination 
      real*8 blen      ! Baseline length 
      real*8 scale(max_neq)
      real*8 ap(4)     ! DD operator [1 -1 -1 1]
     .,      maxdiag  ! Largest elements on diagonal

      integer*4 ipivot(max_neq),
     .       sum_num(max_cfiles*max_gprn)

      real*8 neqb(max_cfiles+3*max_gprn, max_cfiles+3*max_gprn),
     .       bvcb(max_cfiles+3*max_gprn),
     .       apbias(max_cfiles+3*max_gprn)
      real*8 res, rec, obc
      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases.

      integer*4 off, it, p

      logical kbit, done


      logical data_OK ! Function which returns true to data passes
                      ! the phs_mask on the data flag
     .,       good_bf ! Function to see if good bias flags
      logical ibfd(max_cfiles*max_gprn) ! Logical set true when an implicit bias
                          ! have been removed (to avoid doing it multiple
                          ! times when a new bias flag have been found).

***** OK: All done, Start formal estimation
      call report_stat('status','autcln','EST_COMB',' ',
     .                 'Start '//comb_type,0)
      nump = num_cfiles*num_sat
      do  i = 1, nump
          sol_eq(i) = 0.d0
          sum_num(i) = 0.d0
          do j = 1,nump
             norm_eq(i,j) = 0.d0
          end do
      end do

      do i = 1,num_cfiles
         site_dcb(i) = -1
      end do

      ap(1) =  1.0d0
      ap(2) = -1.0d0
      ap(3) = -1.0d0
      ap(4) =  1.0d0

****  Now build up the normal equations
      do ep = 1, num_ep
         do i = 1, nump
            ibfd(i) = .false.
         end do

         do i = 1, num_cfiles-1
            do j = 1, num_sat-1

*              Loop over all epochs and satellite combinations 
*              solving biases implictly when needed.
               ch =  ltoc( ctol_cse(1,i,ep),j, 
     .                     actual_max_chan)

*              See if bias flag on good data that we need to remove
*              bias flag.
               if( ch.gt.0 ) then
                 if( good_bf(data_flag_cse(ch,i,ep),0,phs_mask) ) then
                     np(1) = (i-1)*num_sat+j
                     call implicit_bf(np(1),num_cfiles*num_sat, 
     .                    max_neq, norm_eq, sol_eq,1)
                     ibfd(np(1)) = .true.
                 endif
               endif 

*              Data at this time is in the interval of the last
*              bias flag so see which DD ambs it contributes to
               call get_obs(ep, i,j,  
     .            L1r_phs_cse, L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .            L1r_rng_cse, L2r_rng_cse, data_flag_cse, ctol_cse,
     .            azel_cse, obs(1,1), wgh(1,1), OK(1))

*              If OK:then this oneway matches and has good data
               if( OK(1) ) then 

*                  Set dcb type for site
                   if( site_dcb(i).eq.-1 ) then
                      site_dcb(i)  = 0
                      if( kbit(data_flag_cse(ch,i,ep),28) ) then
                         site_dcb(i)  = 1
                       endif
                       if( kbit(data_flag_cse(ch,i,ep),29) ) then
                         site_dcb(i)  = 2
                       end if
                   end if

*                  Find the next station
                   do k = i+1, num_cfiles
                      ch =  ltoc( ctol_cse(1,k,ep),j, 
     .                            actual_max_chan)

                      if( ch.gt.0 ) then

*                        See if implicit bias solution needed for this
*                        site/sv combination
                         np(2) = (k-1)*num_sat+j
                         if( good_bf(data_flag_cse(ch,k,ep),
     .                                          0,phs_mask) .and.
     .                        .not. ibfd(np(2)) ) then
                            call implicit_bf(np(2),num_cfiles*num_sat, 
     .                           max_neq, norm_eq, sol_eq,1)
                            ibfd(np(2)) = .true.
                         endif
                         OK(2) = data_OK(data_flag_cse(ch,k,ep),
     .                                 0,phs_mask)
                      else
                          OK(2) = .false.
                      endif

                      if ( OK(2) ) then
*                         Set dcb type for site
                          if( site_dcb(k).eq.-1 ) then
                             site_dcb(k)  = 0
                             if( kbit(data_flag_cse(ch,k,ep),28) ) then
                                site_dcb(k)  = 1
                             endif
                             if( kbit(data_flag_cse(ch,k,ep),29) ) then
                                site_dcb(k)  = 2
                             end if
                         end if
                         call get_obs(ep, k, j, L1r_phs_cse, 
     .                      L2r_phs_cse, L1_cyc_cse, L2_cyc_cse,
     .                      L1r_rng_cse, L2r_rng_cse, 
     .                      data_flag_cse, ctol_cse,azel_cse, 
     .                      obs(1,2), wgh(1,2), OK(2))

*                        Now find next satellites
                         do l = j+1,num_sat
*                           Check station 1 to see of OK
                            ch =  ltoc( ctol_cse(1,i,ep),l, 
     .                               actual_max_chan)
                            OK(3) = .false.
                            OK(4) = .false.
                            if( ch.gt.0 ) then
**                             See if implicit bias solution
*                              needed here
                               np(3) = (i-1)*num_sat+l
                               if( good_bf(data_flag_cse(ch,i,ep),
     .                                          0,phs_mask) .and.
     .                              .not. ibfd(np(3)) ) then
                                  call implicit_bf(np(3),
     .                                 num_cfiles*num_sat,  
     .                                 max_neq, norm_eq, sol_eq,1)
                                  ibfd(np(3)) = .true.
                               endif

                               OK(3)=data_OK(data_flag_cse(ch,i,ep),
     .                                 0,phs_mask)
                            endif

                            ch =  ltoc( ctol_cse(1,k,ep),l, 
     .                               actual_max_chan)
                            if( ch.gt.0 ) then
                               np(4) = (k-1)*num_sat+l
                               if( good_bf(data_flag_cse(ch,k,ep),
     .                                          0,phs_mask) .and.
     .                             .not. ibfd(np(4)) ) then 
                                  call implicit_bf(np(4),
     .                                 num_cfiles*num_sat,  
     .                                 max_neq, norm_eq, sol_eq,1)
                                  ibfd(np(4)) = .true.
                               endif

                               OK(4)=data_OK(data_flag_cse(ch,k,ep),
     .                                 0,phs_mask)

                            endif
*                           If both OK, continue
                            if( OK(3) .and. OK(4) ) then

                               call get_obs(ep, i, l, L1r_phs_cse, 
     .                           L2r_phs_cse, L1_cyc_cse,
     .                           L2_cyc_cse, L1r_rng_cse,  
     .                           L2r_rng_cse,data_flag_cse,ctol_cse, 
     .                           azel_cse, obs(1,3), wgh(1,3), OK(3))
                               call get_obs(ep, k, l, L1r_phs_cse, 
     .                           L2r_phs_cse, L1_cyc_cse,
     .                           L2_cyc_cse, L1r_rng_cse,  
     .                           L2r_rng_cse,data_flag_cse,ctol_cse, 
     .                           azel_cse, obs(1,4), wgh(1,4), OK(4))

*****                          We have data, so increment normal equations
*                              Based on comb_type get the wgh and dobs
                               ent = bl_ent(i,k)
                               blen = baselens(ent)
                               if( comb_type.eq.'LC' ) then 
                                  tot_wgh = 1.d0/(1.d0/wgh(5,1)+
     .                                            1.d0/wgh(5,2)+
     .                                            1.d0/wgh(5,3)+
     .                                            1.d0/wgh(5,4))
*                                 Down-weight the longer baselines.  Longest
*                                 baseline is 0.01 in weight (10% sigma)
                                  tot_wgh = tot_wgh
C     .                                *sqrt(1.d-4/max(blen/12.d6,1.d-4))
                                  dobs = (obs(5,1)-obs(5,2))-
     .                                   (obs(5,3)-obs(5,4))
                               elseif( comb_type.eq.'EXWL' ) then
                                  tot_wgh = 1.d0/(1.d0/wgh(2,1)+
     .                                            1.d0/wgh(2,2)+
     .                                            1.d0/wgh(2,3)+
     .                                            1.d0/wgh(2,4))
*                                 For EXWL downweight proportional to length
*                                 Longest baseline is 1.d-4 (1% sigma).
                                  tot_wgh = tot_wgh*
     .                                 (msig_wl_tol/
     .                                  max(blen/12.d6,msig_wl_tol))
*    .                                 (exp(-(blen/100.d3)**2))
                                  dobs = (obs(2,1)-obs(2,2))-
     .                                   (obs(2,3)-obs(2,4))

                               else     ! Default is MWWL combination.
                                  tot_wgh = 1.d0/(1.d0/wgh(1,1)+
     .                                            1.d0/wgh(1,2)+
     .                                            1.d0/wgh(1,3)+
     .                                            1.d0/wgh(1,4))
                                  dobs = (obs(1,1)-obs(1,2))-
     .                                   (obs(1,3)-obs(1,4))
                               endif

                               np(1) = (i-1)*num_sat+j
                               np(2) = (k-1)*num_sat+j
                               np(3) = (i-1)*num_sat+l
                               np(4) = (k-1)*num_sat+l


*                              Increment equations
                               if( ep.lt.0 ) then
                                  write(*,405) ep, comb_type, i,j, 
     .                               k, l, tot_wgh,
     .                               (jn, (wgh(in,jn), in=1,5),jn=1,4)
 405                              format('WGH ',i5,1x,a4,' SS ',4i3,1x,
     .                               F8.3,4(' SW ',I2,5(1x,F8.3)))
                               end if
                                
                               do in = 1,4
                                   sol_eq(np(in)) = sol_eq(np(in))+
     .                                ap(in)*dobs*tot_wgh
                                   sum_num(np(in)) = sum_num(np(in))+1
                                   do jn = 1,4
                                       norm_eq(np(in),np(jn)) = 
     .                                    norm_eq(np(in),np(jn)) +
     .                                    ap(in)*ap(jn)*tot_wgh
                                   enddo
                               enddo
                            endif   ! OK(3) and OK(4)
                         enddo      ! Loop on l satellite
                      endif         ! Second station OK
                   enddo            ! Loop on k station
               endif                ! First station/satellite OK
            end do                      ! First satellite
         end do                         ! First station
      enddo                            ! Looping over epochs

*     Make sure all dcb's have been set (deleted stations will not
*     have values (see to zero so index will be OK).
      do i = 1,num_cfiles
         if( site_dcb(i).eq. -1 ) site_dcb(i) = 0
      end do


C     write(*,410) num_cfiles, num_sat, num_cfiles*num_sat
C410  format('AMB: ',3i6)
      if( 1.eq.2 ) then 
         do i = 1,num_cfiles*num_sat
             if( norm_eq(i,i).eq.0 ) norm_eq(i,i) = 1.d0
             write(*,420) comb_type,(norm_eq(i,j),
     .                               j=1,num_cfiles*num_sat)
 420         format('NEQ',a,1x,1600E18.10)
         enddo
         do i = 1, num_cfiles*num_sat
             write(*,440) comb_type, sol_eq(i)
 440         format('SOL',a,1x,1600E18.10)
         enddo
      endif


***** Now apply basic constraints; One reference site and satellite
      write(*,480) wl_ref_site, cf_codes(wl_ref_site)
 480  format('Constraining site ',i4,2x,a4)
      maxdiag = 0.d0
      do i = 1, nump
         if( norm_eq(i,i).gt.maxdiag ) maxdiag = norm_eq(i,i)
      end do
      do i = 1,num_sat
         in = (wl_ref_site-1)*num_sat+i
         norm_eq(in,in) = norm_eq(in,in) + 2*maxdiag
      enddo

      do i = 1,num_cfiles
        if( i.ne.wl_ref_site) then
           in = (i-1)*num_sat + wl_ref_svs
           norm_eq(in,in) = norm_eq(in,in) + 2*maxdiag  !  1.d10
        endif
      enddo

      write(*,485) wl_ref_svs, prn_list(wl_ref_svs), maxdiag
 485  format('Constraining sv ',i4,' PRN ',i2.2,' MaxDiag ',E12.4)

*     Now apply a loose constraint for dependent biases
      do i = 1, nump 
C        norm_eq(i,i) = norm_eq(i,i) + 0.1d0   !  *maxdiag
         norm_eq(i,i) = norm_eq(i,i) + 1.d-6*maxdiag
      end do
      
      call invert_vis(norm_eq,sol_eq,scale,ipivot,nump, max_neq,1)

*
      if( 1.eq.2 ) then
         do i = 1,num_cfiles*num_sat
             write(*,520) comb_type, (norm_eq(i,j),
     .                    j=1,num_cfiles*num_sat)
 520         format('INV',a,1x,1600E18.10)
         enddo
         do i = 1, num_cfiles*num_sat
             write(*,540) comb_type, sol_eq(i)
 540         format('ANS',a,1x,1600E18.10)
         enddo
      endif

****  Now print results and list statostics
      do i = 1, num_cfiles
         do j = 1, num_sat
            in = (i-1)*num_sat + j
            write(*,600) comb_type, in, cf_codes(i), site_dcb(i), 
     .                   prn_list(j), sol_eq(in),
     .                   sqrt(norm_eq(in,in)), sum_num(in)
 600        format('OWBIAS Type ',a4,1x,i5,2x,a4,1x,i1,1x,I2.2,
     .              F9.3,' +- ',F7.3,1x, I6)

*           Based on comb_type; save the values
            if( comb_type.eq.'LC' ) then
                ambest_lc(in) = sol_eq(in)
                ambsig_lc(in) = sqrt(norm_eq(in,in))
            elseif( comb_type.eq.'EXWL' ) then
                ambest_ex(in) = sol_eq(in)
                ambsig_ex(in) = sqrt(norm_eq(in,in))
            else
                ambest_mw(in) = sol_eq(in)
                ambsig_mw(in) = sqrt(norm_eq(in,in))
                ambnum_mw(in) = sum_num(in)
            endif

         enddo
      enddo 

****  If this is EX-WL save the covariance matrix so that we can
*     compute DD sigmas (weights need full covariance matrxi)
      if( comb_type.eq.'EXWL' ) then
         do i = 1, num_cfiles*num_sat
            do j = 1, num_cfiles*num_sat
               Cov_OWEx(i,j) = norm_eq(i,j)
            end do
         end do
      end if

****  Now if we doing the MWWL estimates then compute the site and
*     svs bias by DCB type.

      if( comb_type.eq.'MWWL' ) then

****      For the MW-WL case now add estimation of biases on each satellite by
*         receiver type
          do i = 1, num_cfiles+3*num_sat
              apbias(i) = 0.d0
          enddo 

*****     Iterate the solution
          done = .false.
          it = 0
          do while ( .not. done )
             it = it + 1
             if( it.ge.10 ) done = .true.  ! Check below on adjustments 
             do i = 1, num_cfiles+3*num_sat
                 bvcb(i) = 0.d0
                 do j = 1, num_cfiles+3*num_sat
                    neqb(i,j) = 0.d0
                 end do
             end do
*            Loop over one-ways
             do i = 1, num_cfiles
                off = site_dcb(i)*num_sat
                do j = 1, num_sat
                   fr1 = fL1u(j)/fL1(j)
                   p = (i-1)*num_sat + j
                   np(1) = i
                   np(2) = num_cfiles+off+j
*                  If we have data
                   if( sum_num(p).gt.0 ) then
                       obc = sol_eq(p) - (apbias(np(1))+apbias(np(2)))
                       rec = obc-nint(obc*fr1)/fr1
                       tot_wgh = 1.d0/(norm_eq(p,p)+1.d-4)
                       do k = 1,2
                          in = np(k)
                          bvcb(in) = bvcb(in) + rec*tot_wgh
                          do l = 1,2
                             neqb(np(k),np(l)) = neqb(np(k),np(l))+
     .                                           tot_wgh
                          enddo
                       end do
                   end if
                end do
             end do

****         Formed the normal equations, now solve system with weak constraint
             nump = num_cfiles+3*num_sat
             do i = 1,nump
                neqb(i,i) = neqb(i,i) + 1
             enddo

             call invert_vis(neqb,bvcb,scale,ipivot,nump,
     .                   max_cfiles+3*max_gprn,1)

****         Update the aprioris, output the results and iterate
             res = 0.d0
             do i = 1, nump
                apbias(i) = apbias(i) + bvcb(i)
                if( abs(bvcb(i)).gt.res ) res = abs(bvcb(i))
             enddo
             write(*,620) it, res
 620         format('OWB',I2.2,' Max Adj ',F8.3)
             if( res.lt.1.d-4 ) it = 10

          end do

****      Now save the results for later use
          do i = 1, num_cfiles
             mwbias_site(i) = apbias(i)
          enddo
          do i = 0,2
             do j = 1,num_sat
                p = num_cfiles + i*num_sat + j
                mwbias_svs(j,i+1) = apbias(p)
            enddo
          enddo

****      Now save the results for later use
          write(*,630) it
 630      format('DCB Bias updates: Iter ',i3)
          do i = 1, num_cfiles
             write(*,640) cf_codes(i), site_dcb(i), 
     .                 apbias(i), bvcb(i), sqrt(neqb(i,i))
 640         format('DCB_BIAS: Site ',a4,1x,i1,1x,
     .         ' WLB ',2F8.3, ' +- ',F8.3,' cyc')
          enddo
          do i = 0,2
             do j = 1,num_sat
                p = num_cfiles + i*num_sat + j
                write(*,650)prn_list(j), i, apbias(p), 
     .             bvcb(p), sqrt(neqb(p,p))
 650            format('DCB_BIAS:    PRN ',I2.2,1x,i1,4x,
     .               ' WLB ',2F8.3, ' +- ',F8.3,' cyc')
             enddo
          enddo 

      endif

****  Return thats all
      return
      end

CTITLE fin_gamit_dd
 
      subroutine fin_gamit_dd(L1r_phs_cse, L2r_phs_cse,
     .    L1_cyc_cse, L2_cyc_cse, L1r_rng_cse, L2r_rng_cse, 
     .    ctol_cse, data_flag_cse)

      implicit none
 
*     Routine to write out the double difference ambiquities to 
*     be used by solve.  Scheme here is to generate all the baselines
*     from a single station and then to replace longer baselines with
*     shorter ones.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'ctogobs_com.h'

* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
*   par_flag_cse(num_param, num_ep)     - Parameter estimate quality
*                   - flags.
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)

*   L1_cyc_cse(num_chan, num_cfiles, num_ep)    - Number of cycles
*                   - needed for each L1 phase measurement .  May be
*                   - fracttional for half cycle units.
*   L2_cyc_cse(num_chan, num_cfiles, num_ep)    - number of cycles
*                   - needed for each L2 phase measurement
*   L1r_phs_cse(num_chan, num_cfiles, num_ep)  - L! phase residuals
*                   - cylces at L1
*   L2r_phs_cse(num_chan, num_cfiles, num_ep)  - L2 phase residuals
*                   - cycles at L2
*   L1r_rng_cse, L23_rng_cse -- L1 and L2 range residuals.
*   params_cse(num_param, num_ep)       - Clock parameter estimates
*                   - by epoch.
      real*8 L1_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L2_cyc_cse(num_chan, num_cfiles, num_ep),
     .    L1r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L2r_phs_cse(num_chan, num_cfiles, num_ep),
     .    L1r_rng_cse(num_chan, num_cfiles, num_ep),
     .    L2r_rng_cse(num_chan, num_cfiles, num_ep)
 

* LOCAL VARIABLES

      integer*4 nblen   ! number of baselines
     .,   bfi((max_cfiles*(max_cfiles+1))/2) ! Index of baselines, 
                        ! set 0 when baseline included in sort list
     .,   bfj((max_cfiles*(max_cfiles+1))/2) ! Sort list of baseline lengths
     .,   i,j,k         ! Loop countres
     .,   pair(2)       ! Sites in current baseline
     .,   last          ! BF table entries of the last bias on each satellite
                        ! in entries being tested
     .,   numout, numfxd  ! Number of DD biases writen and of these the number
                        ! fixed.
     .,   last_bfs(max_gprn,max_cfiles)  ! List of last bias flags in the bf_table each
                        ! site and satellite
     .,   lbf_status(max_gprn,max_cfiles)  ! Status of last bias flags, Set -1 if no 
                        ! data, 0 initially and then 1 once it hs been used
     .,   nnew, nuse    ! Number of satellites new for baseline (not used at either
                        ! site) and used (used at both sites)
     .,   svnew(max_gprn),svuse(max_gprn) ! List of satellites new for baseline and 
                        ! used for baselined.  Only nnew-1 ambiguities are added at
                        ! each baseline, but differences can be formed with nnew+nuse-1
                        ! satellites
     .,   js, ks, jt    ! Satellite numbers
     .,   num_last_bfs  ! Number of last bias flags
     .,   bl_ent        ! Function to return entry in baseline table for sites i,j
     .,   longcom_dur, longcom_bf ! Longest duration and bf_table number for satellite
                        ! visible at all sites
     .,   longsng_dur, longsng_bf ! Longest duration and bf_table number for satellite
                        ! that is visible only at some stations
     .,   longone_dur, longone_bf ! Longest duration and bf_table number for satellite
                        ! that is visible only at some stations
     .,   bf            ! Generic bf_table entry number
*    .,   namb          ! Entry in ambiquity tables (computed from site and satellite 
*                       ! number)
*    .,   amb_tab(5,max_cfiles*max_gprn)  ! Ambiquity table.  Entries in table are:
*                       ! 1 -- Site 1 (starts as ref_site)
*                       ! 2 -- Site 2 (Second site in baseline)
*                       ! 3 -- SVS  1 (starts as Ref_svs)
*                       ! 4 -- SVS  2 (second satellite in DDiff)
*                       ! 5 -- Status: 1 if amb OK, -1 if reference site/svs
*                       !      -2 No data on satellite
     .,   kbf(4)        ! The four bf_table entries that make up double differnce
                        ! ambiquity being output
     .,   cbl           ! Baseline number
     .,   chk_overlap   ! Function to check overlap of double diffs
     .,   overlap, max_ovl, max_sv  ! Overlap, max value and satellite associated
                        ! maximum value
     .,   bovl          ! overlap value for best sigma combination
     .,   eps, epe      ! Epoch start and stop of overlap
     .,   DD(2,max_cfiles-1)  ! Baselines that go into ambiquities for
                        ! each satellite.  Site entries set negative if not needed
     .,   DDBL(max_cfiles-1,max_gprn)  ! Baseline numbers associated with DD and
                        ! saved by satellite
     .,   usedDD(max_cfiles-1)  ! Flag to show baseline has been used.
     .,   svs_list(max_gprn)    ! Sorted list of amount of satellite data for each
                        ! baseline.
     .,   num_list(max_gprn)    ! Number of data per satellite
     .,   expected_namb ! Expected number of ambiguities
     .,   nzero         ! Number of zeroes
     .,   usedSite(max_cfiles)  ! Status of site in forming baselines
                        ! 0 - Site not used yet, -1 no data on PRN, used
     .,   minblindx     ! Minimum baseline index in sorted list (proxy for length)
     .,   ib, kb, nb    ! Sites with shortest baselines
 

      real*8  WLAmb         ! Value for widelane ambiquitiy
      real*8  DDLC, DDEX, DDMW, DDMWC, DDMWT ! Double diff LC, EX-WL, 
                        ! MW-WL and MW-WL with dcb correction
      real*8  SLC, SEX, SMW   ! Sigmas for above values.
      real*8  Exwgh      ! EX-WL weight based on length
      real*8  derraw(3)  ! Sigma scale errors for -1,0,1 ambiguity



      real*8 maxlen     ! Maximum length of currently shortest baseline
     .,   confid        ! Sum of confidence for dd wl.


      integer*4 np(4)   ! Parameter numbers in ambiguity estimates.
     .,         jc, kc  ! Counters
     .,         bks  ! Best DD covariance for svs js and ks
     .,         IMW, dernum ! Integer MW and offset value
     .,         derfn1  ! Integer estimates of number of L1 cycles
     .,         dMW     ! Values arounf DMW to be tested.

      real*8 ap(4)      ! DD partial [1 -1 -1 1]
     .,      DDCov      ! Variance of DD 
     .,      bDDCov     ! Best DDcov value
     .,      fN1, sN1   ! Float estimate fo L1 cycles and sigma
     .,      blen       ! Baseline length
     .,      dermin     ! Mimum value of error metric to see if bias fixable

* MOD TAH 200617 Added frequency ratio (only fr1 needed)
      real*8 fr1,fr2  ! Frequency ratio for Glonass: fr1 = fL1u/fL1 and we 
                 ! multiply by the factor to get integer estimate;
                 ! once resolved to an integer; we divide to get 
                 ! fractional cycle to be applied to remapped phases.
 

      logical OK_svs        ! Set true to indicate satellite present at all sites.
     .,   OK_site       ! Set true is site has last bias on all satellites or
                        ! all stations are missing that same satellite
     .,   bldone        ! Set true when baselines all filled

      character*1 type  ! Set X for fixed, R for Free bias

      write(uns,120)
 120  format('ACBIAS',4x,'double_difference biases:')

*     First sort the baselines into increasing length
      nblen = (num_cfiles)*(num_cfiles-1)/2
*     Form an index array of all baselines, when one is included
*     in the sorted list, this entry is set zero so that it 
*     won't be used again.  This covers the cases when some baselines
*     have exactly the same length.
      do i = 1, nblen
          bfi(i) = i
      end do

      ap(1) =  1.d0
      ap(2) = -1.d0
      ap(3) = -1.d0
      ap(4) =  1.d0

****  Check for sites that have been deleted and move these to the
*     end if of the length list
      call get_blens( np )  ! computes baseline lengths and
*               ! returns shorest WL-weighted length 
C     do i = 1, num_cfiles-1
C        do j = i+1, num_cfiles
C           k = bl_ent(i,j)
C           orderlen(k) = baselens(k)
C           if( lc_num(i).eq.0 .or. lc_num(j).eq.0 ) then
C               orderlen(k) = baselens(k) + 200.d6
C           endif
C           if( WL_RMS(i).gt.1.d0 .or. WL_RMS(j).gt.1.d0 ) then
C               orderlen(k) = baselens(k) + 15.d6*
C    .             max(WL_RMS(i),WL_RMS(j))
C           endif
C        end do
C     end do

      do i = 1, nblen
          k = 0
          maxlen = 1.d30
          do j = 1,nblen
             if( orderlen(j).le.maxlen .and. bfi(j).ne.0 ) then
                 maxlen = orderlen(j)
                 k = j
             endif
          end do
          bfj(i) = k
          bfi(k) = 0
      end do

***** Generate the list of last bias flags for each site and satellites
      num_last_bfs = 0
      do i = 1, num_cfiles
         do j = 1, num_sat
             last_bfs(j,i) = 0  ! Clear the point for each site/sv combination
                                ! (If value remains zero then no data on sv)
          end do
      end do
      do i = 1, tot_bf
         call get_last_bias(i,last)
         if( last.eq.i ) then
             last_bfs(bf_table(3,i),bf_table(2,i)) = i
             num_last_bfs = num_last_bfs + 1
         end if
      end do

      write(*,140) tot_bf, num_last_bfs
 140  format('Of ',i5,' Total biases, there are ',i5,' last biases')

***** Find the satellite that is present at all stations and is the longest 
*     length.  We will use this station and satellite as the base
      longcom_dur = 0
      longsng_dur = 0
      longone_dur = 0
      do i  = 1, num_cfiles
         write(*,150) cf_codes(i),(last_bfs(j,i),j=1,num_sat)
* MOD MAF 210701: Updated 32(I5) to 50(I5) to allow for 45 Beidou satellites
c150     format('LBF ',a4,32(I5))
 150     format('LBF ',a4,50(I5))

* MOD TAH 050410: Only consider sites that have data on all biases
        OK_site = .true.
        do j = 1, num_sat
*           Station is missing a satellite, exclude from consideration
            if( last_bfs(j,i).eq.0 ) then
               OK_svs = .false.
*              Check to see if any other stations have the same satellite
               do k = 1, num_cfiles
                  if( last_bfs(j,k).gt.0 ) OK_svs = .true.
               end do
*              If some other sites have this satellite then exclude this
*              station.
               if( OK_svs ) OK_site = .false.
            end if
        end do
*       If site excluded, report now, otherwise use it for check
        if( .not. OK_site ) then
            write(*,160) cf_codes(i)
 160        format('Removing ',a4,' as possible reference site due ',
     .             'to missing satellite')
* MOD TAH 050618: Keep track of longest here in case no station can 
*       see all satellites
            do j = 1, num_sat
               bf = last_bfs(j,i)
               if( bf.gt.0 ) then
                   if( bf_table(4,bf).gt.longone_dur ) then
                       longone_dur = bf_table(4,bf)
                       longone_bf  = bf
                   end if
                endif
             end do

        ELSE
*          Check for longest sequence
           do j = 1, num_sat
               bf = last_bfs(j,i)

               OK_svs = .true.
               if( bf.gt.0 ) then
                  do k = 1, num_cfiles
*                    Only check sites that have at least some data
                     if( last_bfs(j,k).eq.0 .and. LC_num(k).gt. 0 ) 
     .                                              OK_svs = .false.
                  end do
                  
*                 If this satellite is everywhere, check the length
                  if( OK_svs ) then
                      if( bf_table(4,bf).gt.longcom_dur ) then
                          longcom_dur = bf_table(4,bf)
                          longcom_bf  = bf
                      end if
                  else     ! Save the non-all present longest one
                      if( bf_table(4,bf).gt.longsng_dur ) then
                          longsng_dur = bf_table(4,bf)
                          longsng_bf  = bf
                      end if
                  endif
               endif
           end do
        endif
      end do
*
*     Report the longest last bias flag and save the referernces 
      if( longcom_dur.gt.0 ) then
         write(*,200) cf_codes(bf_table(2,longcom_bf)), 
     .                prn_list(bf_table(3,longcom_bf)),
     .                bf_table(1,longcom_bf), bf_table(4,longcom_bf)
 200     format('Longest com BF ',a4,' PRN ',i2.2,' Start ',i5,
     .          ' Dur ',i4)
         ref_site = bf_table(2,longcom_bf)
         ref_svs  = bf_table(3,longcom_bf)
      else if( longsng_dur.gt.0 ) then
         write(*,210) cf_codes(bf_table(2,longsng_bf)), 
     .                prn_list(bf_table(3,longsng_bf)),
     .                bf_table(1,longsng_bf), bf_table(4,longsng_bf)
 210     format('Longest sng BF ',a4,' PRN ',i2.2,' Start ',i5,
     .          ' Dur ',i4)
         ref_site = bf_table(2,longsng_bf)
         ref_svs  = bf_table(3,longsng_bf)
      else if( longone_dur.gt.0 ) then
         write(*,220) cf_codes(bf_table(2,longone_bf)), 
     .                prn_list(bf_table(3,longone_bf)),
     .                bf_table(1,longone_bf), bf_table(4,longone_bf)
 220     format('Longest one BF ',a4,' PRN ',i2.2,' Start ',i5,
     .          ' Dur ',i4)
         ref_site = bf_table(2,longone_bf)
         ref_svs  = bf_table(3,longone_bf)
      else
         write(*,230)
 230     format('NO REF SITE and SVS FOUND: Stop')
         stop 'AUTCLN: No ref site and satellite found'
      endif

****  Initialize the last bias flag status (lbf_status)
      do i = 1, num_cfiles
         do j = 1, num_sat
             if( last_bfs(j,i).eq.0 ) then
*                There is not last bias flag for this 
*                combination
                 lbf_status(j,i) = -1
             else
                 lbf_status(j,i) = 0
                 if( LC_SVS_num(j,i).eq.0 ) then
                    write(*,320) cf_codes(i), prn_list(j),
     .                 lbf_status(j,i) 
 320                format('GAMITDD Warning no data for ',a4,
     .                     ' PRN',i2.2,' Last BF ',i5)
                 endif
             endif
         end do
      end do

      namb = 0
      do j = 1, num_cfiles
         write(*,360) cf_codes(j),(lbf_status(k,j),k=1,num_sat)
         do k = 1, num_sat
            if( lbf_status(k,j).eq.-1 ) then
                namb = namb + 1
            endif
         end do
      end do

*     Now account for missing stations and satellites
      ib = 0
      do j = 1, num_cfiles
*         See if all satellites missing
          ib = 0
          do k = 1, num_sat
             if( lbf_status(k,j).ne.-1 ) ib = ib + 1
          end do
*         Decrement namb if no last bias flags at all
          if( ib.eq.0 ) namb = namb - 1
      enddo
*     Repeat for satellites
      do k = 1, num_sat
         ib = 0
         do j = 1,num_cfiles
            if( lbf_status(k,j).ne.-1 ) ib = ib + 1
          end do
*         Decrement namb if no last bias flags at all
          if( ib.eq.0 ) namb = namb - 1
      enddo

      expected_namb = (num_cfiles-1)*(num_sat-1)-namb
      write(*,325) (num_cfiles-1)*(num_sat-1), namb, 
     .             expected_namb
 325  format('GAMITDD: All biases ',i4,' Missing SVS ',I4,
     .       ' Expected ambiguities ',i4)

C     do i = 1, nblen
C        call ent_bf(bfj(i),num_cfiles, pair)
C        write(*,365) i, bfj(i), pair, cf_codes(pair(1)), 
C    .                cf_codes(pair(2))
C 365    format('BL ',i4,i5, 2i4, 1x, a4,1x,a4)
C     end do

****  Set up baselines to be used for each satellite
      nzero = 0
      do j = 1, num_sat

*        Now generate the list of independent baselines in ascending order
*        Start with the shortest baseline and then get closest stations
*        to end and work out from there finding next closest station.
*
         do i = 1, num_cfiles-1
            DD(1,i) = 0
            DD(2,i) = 0
            usedDD(i) = 0
         end do
         do i = 1, num_cfiles
            usedSite(i) = 0
            if( lbf_status(j,i).eq.-1 ) usedSite(i) = -1
         end do

*        Now fill the DD array with the shortest baselines that can see
*        this satellite
         bldone = .false.
*        Find the first baseline with data on this satellite
         k = 0
         do while ( k.lt. nblen .and. .not. bldone )
             k = k + 1 
             call ent_bf(bfj(k),num_cfiles, pair)
             if( lbf_status(j,pair(1)).eq.0 .and.
     .           lbf_status(j,pair(2)).eq.0 ) then
*                Baseline has good data at both sites so start
*                here.
                 DD(1,1) = pair(1)
                 DD(2,1) = pair(2)
                 usedDD(1) = 1
                 bldone = .true.
                 usedSite(pair(1)) = 1
                 usedSite(pair(2)) = 1
             endif
         end do
         if( .not.bldone) then
             bldone = .true.  ! No use searching for more baselines
         else
             bldone = .false.
         endif

*****    OK Now fill out all the baseline
         nb = 1
         do while ( .not. bldone )
             minblindx = nblen+1
             ib = 0
             do i = 1, num_cfiles
                if( usedSite(i).gt.0 .and.
     .              lbf_status(j,i).eq. 0 ) then
*                  Site has been used and has data on satellite.
*                  Check lengths to unused used site that have 
*                  good data
                   do k = 1, num_cfiles                      
                      if ( usedSite(k).eq.0 .and.
     .                     lbf_status(j,k).eq. 0 ) then
                         cbl = bl_ent(i,k)
*                        Find where this baseline ranks
                         do ks = 1, nblen 
                            if( bfj(ks).eq.cbl ) js = ks
                         end do
                         if( js.lt.minblindx ) then
*                            This is a shorter baseline so save
                             minblindx = js
                             kb = k
                             ib = i
                         end if
                      end if
                   end do
                end if
             end do
*            OK: If ib and kb have been set then we have
*            found a good baseline so add
             if( ib.gt. 0 ) then
                 nb = nb + 1
*                See if we need to reverse the baseline.
                 if( kb.lt.ib ) then
                    k = kb 
                    kb = ib
                    ib = k
                 endif
                 DD(1,nb) = ib
                 DD(2,nb) = kb
                 usedDD(nb) = 1
                 usedSite(ib) = usedSite(ib)+1
                 usedSite(kb) = usedSite(kb)+1
             else
*                Nothing left so we are done
                 bldone = .true.
             endif
         end do
         if( 1.eq.2 )
     .   write(*,328) prn_list(j),(DD(1,k),DD(2,k),usedDD(k), 
     .                             k=1,num_cfiles-1)
 328     format('DD PRN',i2.2,' DDTrip ',50(2I3,2x,i2,1x))

*        Now assign baseline numbers to the DD entries and save
*        by satellite
         do i = 1,num_cfiles-1
            if( DD(1,i).gt.0 ) then
                cbl = bl_ent(DD(1,i),DD(2,i))
                DDBL(i,j) = cbl
            else
                DDBL(i,j) = 0   ! No baseline with this site
            end if
         end do

         if( 1.eq.2 )
     .   write(*,330) prn_list(j),(DDBL(i,j),i=1,num_cfiles-1)
 330     format('DDBL PRN',i2.2,' BL ',50I5)
         do i = 1,num_cfiles-1
            if ( DDBL(i,j).eq.0 ) then
              nzero = nzero + 1
            endif
         end do

      end do
C     write(*,335) nzero
C335  format('GAMITDD: Nzeros in DDBL ',i4)                     

            
****  Now scan down the baselines in increasing length to
*     set up double difference biases
      namb = 0          
      do i = 1, nblen
*        Get the sites in this baseline
         call ent_bf(bfj(i),num_cfiles, pair)

*        OK: For this baseline, generate sorted list of satellites
         do j = 1,num_sat
             kbf(1) = last_bfs(j,pair(1))
             kbf(2) = last_bfs(j,pair(2))
             if( kbf(1).gt.0 .and. kbf(2).gt.0 ) then
                overlap = chk_overlap(kbf(1), kbf(2),0,0,eps, epe)
             else
                overlap = -9999
             end if
*            See where this fits in list
             OK_svs = .false.

             do js = 1, j-1
                 if( overlap.lt.num_list(js) .and. .not.OK_svs ) then
*                    Move list up so we can insert
                     do ks = j-1,js, -1
                        num_list(ks+1) = num_list(ks)
                        svs_list(ks+1) = svs_list(ks)
                     end do
                     num_list(js) = overlap
                     svs_list(js) = j
                     OK_svs = .true.
                  endif
              end do
              if( .not. OK_svs ) then   ! Add to end
                  num_list(j) = overlap
                  svs_list(j) = j
                  OK_svs = .true.
              endif
         end do
         if( 1.eq.2 )
     .   write(*,337) pair,cf_codes(pair(1)),
     .                cf_codes(pair(2)),
     .                (svs_list(j),num_list(j),j=1,num_sat)
 337     format('GAMITDD: BSRT ',2i3,1x,a4,1x,a4,1x,50(I2.2,1x,i5,2x)) 
        
*        Now find all new and used satellite pairs for this
*        baseline
         nnew = 0
         nuse = 0
         do jt = 1, num_sat
*           Get the sorted satellite number.
            j = svs_list(jt)
*           See if we use this baseline for this satellite
            do k = 1, num_cfiles - 1
               if( DDBL(k,j).eq.bfj(i) ) then
*                  This is baseline for this site so add
                   nnew = nnew + 1
                   svnew(nnew) = j
               end if
            end do 

*           See if satellites used before but not one that we just
*           added.
            if( lbf_status(j,pair(1)).gt.0 .and.
     .          lbf_status(j,pair(2)).gt.0  .and.
     .          svnew(nnew).ne. j   ) then
*               Satellite has been used and so can be
*               for double differences
                nuse = nuse + 1
                svuse(nuse) = j
            endif
         end do


****     OK: Now see if we have any new satellites
         if( nnew.gt.0 ) then
*            For each new satellite, find longest overlap
*            others.
* NOTE: Could replace here with code that generates smallest sigma
*       for the double differences!.
             max_ovl = 0
             max_sv = 0
             ks = 0
             do j = 1, nnew
                max_ovl = 0
                js = svnew(j)
                kbf(1) = last_bfs(js,pair(1))
                kbf(2) = last_bfs(js,pair(2))
                bDDCov  = 1.d10

*               Find best choice for second satellite
                do k = j+1,nnew
                    ks = svnew(k)
                    kbf(3) = last_bfs(ks,pair(1))
                    kbf(4) = last_bfs(ks,pair(2))

                    if( kbf(1).gt.0 .and. kbf(2).gt.0 .and. 
     .                  kbf(3).gt.0 .and. kbf(4).gt.0 ) then

                        overlap = chk_overlap(kbf(1), kbf(2), 
     .                        kbf(3), kbf(4),eps, epe)

*                       Compute sigma of this combination
                        np(1) = (pair(1)-1)*num_sat + js
                        np(2) = (pair(2)-1)*num_sat + js
                        np(3) = (pair(1)-1)*num_sat + ks
                        np(4) = (pair(2)-1)*num_sat + ks

                        DDCov = 0.d0
                        do jc = 1,4
                           do kc = 1,4
                              DDCov = DDCov + 
     .                            ap(jc)*norm_eq(np(jc),np(kc))*ap(kc)
                           enddo
                        enddo
                        if( DDCov.lt.bDDCov ) then 
                           bDDCov =  DDCov
                           bks = ks
                           bovl = overlap
                        endif

                        if( overlap.gt. max_ovl ) then
                            max_ovl = overlap
                            max_sv  = ks
                        end if
                    end if
                end do 

*               Now test against reference list to see if smaller
*               sigma is possible

                do k = 1,nuse
                    ks = svuse(k)
                    kbf(3) = last_bfs(ks,pair(1))
                    kbf(4) = last_bfs(ks,pair(2))
                    if( kbf(1).gt.0 .and. kbf(2).gt.0 .and. 
     .                  kbf(3).gt.0 .and. kbf(4).gt.0 ) then

                        overlap = chk_overlap(kbf(1), kbf(2), 
     .                        kbf(3), kbf(4),eps, epe)
*                       Compute sigma of this combination
                        np(1) = (pair(1)-1)*num_sat + js
                        np(2) = (pair(2)-1)*num_sat + js
                        np(3) = (pair(1)-1)*num_sat + ks
                        np(4) = (pair(2)-1)*num_sat + ks

                        DDCov = 0.d0
                        do jc = 1,4
                           do kc = 1,4
                              DDCov = DDCov + 
     .                           ap(jc)*norm_eq(np(jc),np(kc))*ap(kc)
                           enddo
                        enddo
                        if( DDCov.lt.bDDCov ) then 
                           bDDCov =  DDCov
                           bks = ks
                           bovl = overlap
                        endif

                        if( overlap.gt. max_ovl ) then
                            max_ovl = overlap
                            max_sv  = ks
                        end if
                    end if
                end do

*****           It is possible that there is no overlap on
*               this baseline for this new satellite.  This can
*               happen on early segements of data
                if( (max_ovl.eq. 0) .and. js.ne.ks) then
                    max_ovl = -1
                    if( ks.gt.0 ) then
                       max_sv  = ks  ! Use last satellite?
                    else
                       max_sv = svs_list(num_sat)
* MOD TAH 070410: Make sure not same SVS
                       if( max_sv.eq.js ) then
                           max_sv = svs_list(num_sat-1)
                       endif
                    end if
                    if( bks.eq.0 ) bks = max_sv
                endif
 
                if( max_ovl.ne.0 ) then
                    namb = namb + 1
                    amb_tab(1,namb) = pair(1)
                    amb_tab(2,namb) = pair(2)
                    amb_tab(3,namb) = js
                    amb_tab(4,namb) = bks
                    amb_tab(5,namb) = max_ovl    !  bovl
                    lbf_status(js,pair(1))     = 
     .                                   lbf_status(js,pair(1))+1
                    lbf_status(bks,pair(1)) = 
     .                                   lbf_status(bks,pair(1))+1
                    lbf_status(js,pair(2))     = 
     .                                   lbf_status(js,pair(2))+1
                    lbf_status(bks,pair(2)) = 
     .                                   lbf_status(bks,pair(2))+1
                    write(*,340) namb, cf_codes(pair(1)),
     .                           cf_codes(pair(2)),
     .                           prn_list(js),prn_list(bks), bovl,
     .                           sqrt(bDDCov)
 340                format('GAMITDD Amb ',i4,1x,a4,1x,a4,1x,I2.2,
     .                     1x,I2.2,1x,' MAX_OVL ',i5,' Sigma ',F8.3)
                endif
             enddo 
         endif      
      enddo 

      do i = 1, num_cfiles
         write(*,360) cf_codes(i),(lbf_status(j,i),j=1,num_sat)
* MOD MAF 210701: Updated 32(I5) to 50(I5) to allow for 45 Beidou satellites
c360     format('GAMITDD ',a4,1x,32(I5))
 360     format('GAMITDD ',a4,1x,50(I5))
      end do          

*
***** OK: All done, write out table
      numout = 0
      numfxd = 0
      do i = 1,namb

         if( amb_tab(5,i).ne.0 ) then
   
****         Do final calc on wide lane
*            Compute sigma of this combination
             np(1) = (amb_tab(1,i)-1)*num_sat + amb_tab(3,i)
             np(2) = (amb_tab(2,i)-1)*num_sat + amb_tab(3,i)
             np(3) = (amb_tab(1,i)-1)*num_sat + amb_tab(4,i)
             np(4) = (amb_tab(2,i)-1)*num_sat + amb_tab(4,i)

             DDCov = 0.d0
             DDLC = 0.d0
             SLC  = 0.d0
             DDEX = 0.d0
             SEX  = 0.d0
             DDMW = 0.d0
             SMW  = 0.d0
             SEX  = 0.d0
             do jc = 1,4
                DDLC = DDLC + ap(jc)*ambest_lc(np(jc))
                SLC  = SLC  + ambsig_lc(np(jc))**2
                DDEX = DDEX + ap(jc)*ambest_ex(np(jc))
                DDMW = DDMW + ap(jc)*ambest_mw(np(jc))
                SMW  = SMW  + ambsig_mw(np(jc))**2
             	do kc = 1,4
             	   DDCov = DDCov + 
     .       		   ap(jc)*norm_eq(np(jc),np(kc))*ap(kc)
                   SEX   = SEX + ap(jc)*Cov_owex(np(jc),np(kc))*ap(kc)
             	enddo
             enddo

             SLC = sqrt(SLC * DDCov/SMW)
             SEX = sqrt(SEX)
             SMW = sqrt(SMW * DDCov/SMW)

             DDMWC = DDMW - 
     .          ( mwbias_svs(amb_tab(3,i),site_dcb(amb_tab(1,i))+1)-
     .            mwbias_svs(amb_tab(3,i),site_dcb(amb_tab(2,i))+1)-
     .            mwbias_svs(amb_tab(4,i),site_dcb(amb_tab(1,i))+1)+
     .            mwbias_svs(amb_tab(4,i),site_dcb(amb_tab(2,i))+1)  )
             DDMWT = DDMWC
             DDMWC = DDMW


****         OK Now check resolution
             cbl = bl_ent(amb_tab(1,i),amb_tab(2,i))
             blen = baselens(cbl)
             IMW = nint(DDMWC)
             dermin = 1.d8
             do dMW = -1,1   
* RWK 150203: This may need to be reformulated to account for differences
*             in frequencies among Glonass SVs
                 fN1 = (DDLC-lcf2(1)*(IMW+dMW))/(lcf1(1)+lcf2(1))
                 sN1 = sqrt(SLC**2+(lcf2(1)*SMW)**2)/(lcf1(1)+lcf2(1))
*                EX-WL also has a baseline dependent downweight
*                so we do need an extra one here
                 Exwgh = msig_wl_tol/max(blen/12.d6,msig_wl_tol)
*                Exwgh = 1.d0
        

                 derraw(dMW+2) = abs(fN1-nint(fN1))/sN1 + 
     .                           abs(DDMWC-(IMW+dMW))/max(SMW,0.05d0) +
     .                           abs(DDEX-(IMW+dMW))/SEX*Exwgh

                 
                 write(*,380) 'dMW ',dMW,abs(fN1-nint(fN1))/sN1, sN1,
     .                   abs(DDMWC-(IMW+dMW))/SMW, SMW,
     .                   abs(DDMWT-(IMW+dMW))/SMW, SMW,
     .                   abs(DDEX-(IMW+dMW))/SEX*Exwgh, SEX,  Exwgh
 380             format(a,1x,I2,4(2x,F8.4,1x,F6.3),1x,F6.3)
 
                 if( abs(derraw(dMW+2)).lt.dermin ) then

*                    Best fit so far
                     dermin = abs(derraw(dMW+2))
                     dernum = dMW
                     derfn1 = nint(fN1)
                 endif
             enddo

             WLAmb = IMW   + dernum
             DDEX  = DDEX  - WLAmb
             DDMW  = DDMW  - WLAmb
             DDMWC = DDMWC - WLAmb
             DDMWT = DDMWT - WLAmb 
* RWK 150203: This may need to be reformulated to account for differences
*             in frequencies among Glonass SVs
             DDLC  = DDLC  - (lcf1(1)*derfn1+lcf2(1)*(derfn1+WLAmb))

             if ( dernum.eq.-1 ) confid = min(derraw(2)/derraw(1),
     .                                        derraw(3)/derraw(1))
             if ( dernum.eq. 0 ) confid = min(derraw(1)/derraw(2),
     .                                        derraw(3)/derraw(2))
             if ( dernum.eq. 1 ) confid = min(derraw(1)/derraw(3),
     .                                        derraw(2)/derraw(3))
             if ( confid.gt.dchi_wl_tol ) then
                 type = 'X'
                 numfxd = numfxd + 1
             else
                 type = 'R'
             endif

* MOD TAH 070405: Flip the sign of the WLAmb value written out
* Seems to be long standing bug. (WLAmb replaced with -WLAmb).
* Flip the signs of the means as well to avoid confusion.

             write(uns,420) cf_codes(amb_tab(1,i)),
     .                      cf_codes(amb_tab(2,i)),
     .                      prn_list(amb_tab(3,i)),
     .                      prn_list(amb_tab(4,i)), -WLAmb,
     .                      confid, type, 
     .                      site_dcb(amb_tab(1,i)), 
     .                      site_dcb(amb_tab(2,i)),
     .                      -DDMWC, SMW, -DDMWT, SMW, -DDEX, SEX,
     .                      -DDLC, SLC, -derfn1, 
     .                      blen/1000.d0, Exwgh

 420         format('ACBIAS',4x,'B1L2 ',a4,1x,a4,1x,
     .                i2.2,1x,i2.2,1x,F6.1,1x,
     .                F10.2,1x,a1,4x,' ! ',
     .              ' DCB ',2i2,
     .              ' MWC ',F7.3,' +- ',F7.3,
     .              ' MW  ',F7.3,' +- ',F7.3,
     .              ' EX  ',F7.3,' +- ',F7.3,
     .              ' LC  ',F7.3,' +- ',F7.3,
     .              ' dN1 ',I5,' LEN ',F8.2,' km, EWWgh ',F6.3)
             print *,'ACMWC ', cf_codes(amb_tab(1,i)),
     .                      cf_codes(amb_tab(2,i)),
     .                      prn_list(amb_tab(3,i)),
     .                      prn_list(amb_tab(4,i)), real(-WLAmb),
     .                      real(confid), type 

             numout = numout + 1
         else 
             write(*,480) i,cf_codes(amb_tab(1,i)),
     .                      cf_codes(amb_tab(2,i)),
     .                      prn_list(amb_tab(3,i)),
     .                      prn_list(amb_tab(4,i)),
     .                      amb_tab(5,i)
 480         format('GAMITDD Warning: No data on ambiguity ',
     .               i4,1x,a4,1x,a4,1x,' PRN ',2I3.2,1x,i5)
         endif
      end do
*
***** Tell user how many fixed
      write(uns,'(a)') 'ACBIAS    '
      write(uns,520) numout, num_last_bfs, numfxd,
     .               expected_namb
 520  format('NBIAS out ',2i6,' Number fixed ',i6,' Expected ',i6)
      if( lca_type.eq.'DIR' ) then
         write(uns,540) lca_type, dchi_wl_tol, msig_wl_tol 
 540     format('TOLBIAS Type ',a,' Ratio ',F6.2,' IonScale ',e10.3)
      else
         write(uns,550) lca_type, min_wl_tol, dchi_wl_tol, msig_wl_tol,
     .                  mdev_wl_tol 
 550     format('TOLBIAS Type ',a,' Min Number ',i4,' Min_WL_tol ',F6.2, 
     .          'Min WL Sigma and Deviation ',2F6.2,' cyc')
      endif


****  Thats all
      return
      end

