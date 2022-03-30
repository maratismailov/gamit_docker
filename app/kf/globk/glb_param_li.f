CTITLE GLB_PARAM_LIST
 
      subroutine glb_param_list
 
      implicit none 

 
*     This routine forms the list of global parameters to be estimated,
*     for forms the pointers for the state transission elements in this
*     solution.

* MOD TAH 000808: Reset the value which initializes the apriori sigmas
*     and markov variance to 1.d-16 (from 1.d-6 which is too large for 
*     markov process).
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
*   i,j,k       - Loop counters
*   np          - number of global parameters
*   nm          - number of markov parameters
*   nt          - number of transsions rows found
*   nc          - total number of transission columns
*   kt          - Number of sites needed for tides
*   ne          - Earthquake number
 
      integer*4 i,j,k,  np, nm, nt, nc, kt, in, ne
 
*   need_tran   - Indicates that we need a translation parameters
*               - for either the sites or RA.
*   kbit        - Bit checking routine
 
      logical need_tran, kbit

*   dist  - Distance between site and earthquake
      real*8 dist
 
*   apr_temp, mar_temp  - Temp apririo sigma for site position so that all
*                 will be turned on when NEU are specified.
*   apr_temp_save, mar_temp_save - Apriori and Markov to be saved 
*                 in the cov_mar array.
 
      real*4 apr_temp, mar_temp, 
     .       apr_temp_save, mar_temp_save
 
***** State loop over all parameters and seeing which ones are estimated
*     and which ones are markov.  We also save the apriori variances
*     and markov statistics for each parameter.
 
      np = 0
      nm = 0
      nt = 0
      nc = 0
 
*                                 ! Clear all parameter numbers, and
      call clear_glb_parn
*                                 ! state transsions

****  SITE positions and rates
*                                 ! Loop over all sites
      do i = 1, gnum_sites
        IF( KBIT(GUSE_SITE,i) ) then
*                                 ! Loop over XYZ
*         Check the NEU components.  If are are non-zero then
*         set so that all components will be estimated.
          apr_temp = 0.d0
          mar_temp = 0.d0
          do j = 1,3
             if( apr_neu(j,1,i).ne.0 ) then
*                We set a small value here since this value will be
*                saved in the apriori variance array and we do want
*                to effect the results.
                 apr_temp = 1.d-16          
             end if
             if( mar_neu(j,1,i).ne.0 ) then
                 mar_temp = 1.d-16          
             end if
             if( apr_neu(j,1,i).eq. -1.d0 ) apr_neu(j,1,i) = 0.d0
             if( mar_neu(j,1,i).eq. -1.d0 ) mar_neu(j,1,i) = 0.d0
             
          end do

          do j = 1,3
               apr_temp_save = apr_temp + apr_site(j,1,i)
               mar_temp_save = mar_temp + mar_site(j,1,i)
 
               call glb_par_count( apr_temp_save, mar_temp_save,
     .                             np, nm, parn_site(j,1,i))
          end do
*                                 ! Loop over XYZ dot
*         Check the NEU components.  If are are non-zero then
*         set so that all components will be estimated.
          apr_temp = 0.d0
          mar_temp = 0.d0
          do j = 1,3
             if( apr_neu(j,2,i).ne.0 ) then
                 apr_temp = 1.d-16 
             end if
             if( mar_neu(j,2,i).ne.0 ) then
                 mar_temp = 1.d-16 
             end if
             if( apr_neu(j,2,i).eq. -1.d0 ) apr_neu(j,2,i) = 0.d0
             if( mar_neu(j,2,i).eq. -1.d0 ) mar_neu(j,2,i) = 0.d0
          end do


          do j = 1,3
*             Use save values to stop site accumulating.
              apr_temp_save = apr_temp + apr_site(j,2,i)
              mar_temp_save = mar_temp + mar_site(j,2,i)

              call glb_par_count( apr_temp_save, mar_temp_save,
     .                            np, nm, parn_site(j,2,i))
 
*             Now check for state transission entries
              call check_trans(parn_site(j,1,i), parn_site(j,2,i),
     .                          nc, nt, 1.0 )
          end do
        end if
      end do

* MOD TAH 040703: See if atmospheric delay parameters are estimated
      do i = 1, gnum_sites
         if( kbit(guse_site,i) ) then
*           Count in the atmospheric delays
            call glb_par_count( apr_atm(i), mar_atm(i),
     .                        np, nm, parn_atm(i) )
         end if
      end do
        

* MOD TAH 030607: See if we have parameter estimates for log terms.  Only 
* estimate the log terms if there is more than a few days of data
*     after earthquakes
      do i = 1, gnum_sites
        IF( KBIT(GUSE_SITE,i) .and. 
     .      gepoch_end-gepoch_start.gt.2 ) then
*           See if this site is associated with an earthquake
            do ne = 1, num_eq
                call eval_dist(eq_pos(1,ne), apr_val_site(1,1,i), dist)
                if( dist.le.eq_rad(ne) .and. 
     .              gsite_names(i)(7:8).eq.eq_codes(ne)(1:2) ) then

*                   OK, This is efffected by earthquake ne.  See if
*                   log is to be estimated.

                    do j = 1,3 
                       apr_temp_save = sqrt(eq_log_sig(j,ne)**2 +
     .                                      (eq_log_sig(j+3,ne)*
     .                                      (eq_depth(ne)/dist)**2)**2)
                       mar_temp_save = 0.0
                       call glb_par_count( apr_temp_save, mar_temp_save,
     .                        np, nm, parn_log(j,i) )
                    end do
* MOD TAH 050510: Output to stdout and only if value is non-zero
                    if( apr_temp_save.gt.0 ) 
     .              write(*,220) gsite_names(i), dist/1000, 
     .                           apr_temp_save*1000
 220                format('LOGSig: Site ',a8,' Distance ',F8.2,' km, ',
     .                     'LogSig ',F10.2,' mm')
                endif
            end do
        end if
      end do

***** See if we need to have a translation origin.  Needed when all
*     site positions are estimated
 
*                                 ! Loop over XYZ
      do j = 1,3
*         check to see if the RW process is zero, while the
*         white noise process is non-zero.  This way the parameter
*         process noise will be counted.
          mar_temp = mar_tran(j,1)
          if( mar_tran(j,1).eq.0 .and. mar_tran(j,3).ne.0 ) 
     .        mar_temp = 0.001
          call glb_par_count( apr_tran(j,1), mar_temp,
     .                        np, nm, parn_tran(j,1) )
      end do

*     Now count the rates and account for the transitions.
      do j = 1,3
          call glb_par_count( apr_tran(j,2), mar_tran(j,2),
     .                        np, nm, parn_tran(j,2) )
          call check_trans(parn_tran(j,1), parn_tran(j,2),
     .                        nc, nt, 1.0 )
      end do  

* MOD TAH 981215: Count the parameters for rotations
*                                 ! Loop over XYZ rotatation
      do j = 1,3
*         check to see if the RW process is zero, while the
*         white noise process is non-zero.  This way the parameter
*         process noise will be counted.
          mar_temp = mar_rot(j,1)
          if( mar_rot(j,1).eq.0 .and. mar_rot(j,3).ne.0 ) 
     .        mar_temp = 0.001
          call glb_par_count( apr_rot(j,1), mar_temp,
     .                        np, nm, parn_rot(j,1) )
      end do

*     Now count the rates and account for the transitions.
      do j = 1,3
          call glb_par_count( apr_rot(j,2), mar_rot(j,2),
     .                        np, nm, parn_rot(j,2) )
          call check_trans(parn_rot(j,1), parn_rot(j,2),
     .                        nc, nt, 1.0 )
      end do 

****  Now do the scale and rate of change
* MOD TAH 130716: See if white noise scale invoked
      mar_temp =  mar_scale(1)
      if( mar_scale(3).gt.0 .and. mar_scale(1).eq.0 ) then
          mar_temp = 1.d-6   !  Use small value to set parameter
                             ! Actual value will be used in glb_predict)
      end if 
      call glb_par_count( apr_scale(1), mar_temp,
     .                    np, nm, parn_scale(1) )
      call glb_par_count( apr_scale(2), mar_scale(2),
     .                    np, nm, parn_scale(2) )
      call check_trans(parn_scale(1), parn_scale(2),
     .                    nc, nt, 1.0 )

 
***** AXIS offsets and rates
 
      do i = 1, gnum_sites
        if( kbit(guse_site,i) ) then
*                                                           ! Offset
          call glb_par_count( apr_axo(1,i), mar_axo(1,i),
     .                        np, nm, parn_axo(1,i) )
 
*                                                           ! Rate
          call glb_par_count( apr_axo(2,i), mar_axo(2,i),
     .                        np, nm, parn_axo(2,i) )
 
          call check_trans( parn_axo(1,i), parn_axo(2,i),
     .                      nc, nt, 1.0 )
        end if
      end do
 
***** SOURCE positions
 
*                                 ! Loop over all sources
      do i = 1, gnum_sources
        if( kbit(guse_source,i) ) then
*                                 ! Loop over RA and Dec
          do j = 1,2
 
              call glb_par_count( apr_source(j,1,i),mar_source(j,1,i),
     .                            np, nm, parn_source(j,1,i))
 
          end do
*                                 ! Loop over RA Dec dot
          do j = 1,2
              call glb_par_count( apr_source(j,2,i), mar_source(j,2,i),
     .                            np, nm, parn_source(j,2,i))
 
*             Now check for state transission entries
              call check_trans(parn_source(j,1,i), parn_source(j,2,i),
     .                          nc, nt, 1.0 )
          end do
        end if
      end do
 
***** Check for RA origin shift parameter
      need_tran = .true.
      do i = 1, gnum_sources
          if( parn_source(1,1,i).eq.0 ) need_tran = .false.
      end do
 
      if( need_tran ) then
          call glb_par_count( apr_rao, 0.,
     .                        np, nm, parn_rao )
          call check_trans( parn_rao, parn_rao,
*                                              ! Causes zero transission
     .                      nc, nt, -1.0 )
      end if

****  Do Satellite position

      do i = 1, gnum_svs
         do j = 1,max_svs_elem-3  ! Treat PCO differently
            call glb_par_count( apr_svs(j,i), mar_svs(j,i), np, nm,
     .                          parn_svs(j,i) )
         end do
* MOD TAH 1906228: Treat the special case of negative mar_svs elements 
*        for PCO values

         do j = max_svs_elem-2,max_svs_elem  ! Treat PCO differently
            if( mar_svs(j,i).gt.0 ) then     ! Treat normally
               call glb_par_count( apr_svs(j,i), mar_svs(j,i), np, nm,
     .                             parn_svs(j,i) )
            else       ! Case of white noise PCO, set mar_svs to 0
               mar_temp = 0.0 
               call glb_par_count( apr_svs(j,i), mar_temp, np, nm,
     .                             parn_svs(j,i) )
            endif
         end do 
      end do
 
****  Wobble and UT1 parameters
*     Wobble first
*                         ! Loop over the eight components in the model,
* MOD TAH 960506: Last two values are now used for white-noise, sigma
*     for increment.
* MOD TAH 981020: Only use single PMU if we are not estimating multi-pmu
*     values
      if( num_mul_pmu.eq.0 ) then
          do i = 1,4
              call glb_par_count( apr_wob(i), mar_wob(i),
     .                            np, nm, parn_wob(i) )
          end do
 
* NEW CODE ADDED
*         X-pole (only link to x rate term)
* MOD TAH 950417: Changed 1 to days in years because process
*         is mas/day.
          call check_trans(parn_wob(1), parn_wob(3), nc,nt, 365.25)

*         Y-pole (only link to y rate term)
          call check_trans(parn_wob(2), parn_wob(4), nc,nt, 365.25)
 
         
*****     Now do UT1  ( Loop over the component model)
* MOD TAH 960506: Last value used for white-noise increment sigma 
          do i = 1, 2
              call glb_par_count( apr_ut1(i), mar_ut1(i),
     .                            np, nm, parn_ut1(i))
          end do 
*         Add the UT1 Transition
* MOD TAH 190528: If mar_ut1(2) (lod) is negative, de-couple UT1 from
*         integaated LOD (allow zero becasue LOD maybe deterministic.
          if( mar_ut1(2).ge.0 ) then
              call check_trans( parn_ut1(1), parn_ut1(2), nc,nt, 365.25)
          endif
      else

*         Multiple PMU values
          do i = 1,3                   ! Loop X, Y, UT1
             do j = 1, 2               !  Loop offset and Rate
                do k = 1, num_mul_pmu  ! Loop over epoch values
                   if( i.le.2 ) then 
                       call map_wob_indx('MTOS',in,j,i)
                       call glb_par_count( apr_wob(in), 0.0,
     .                               np, nm, parn_mul_pmu(j,i,k))
                   else
                       call glb_par_count( apr_ut1(j), 0.0,
     .                                np, nm, parn_mul_pmu(j,i,k))
                   end if
                end do
             end do
          end do
      end if
 
***** Nutation angles (Should be done similar to wobble, but for the
*     moment just implement random walk)
 
*                     ! Loop over deps, and dpsi, and seasonal model
      do i = 1,2
          call glb_par_count( apr_nut_ang(i), mar_nut_ang(i),
     .                        np, nm, parn_nut_ang(i) )
      end do

***** UT1, xy and extended Earth tides

      do i = 1, 4
         call glb_par_count( apr_eor_ut1(i), mar_eor_ut1(i), np, nm,
     .                       parn_eor_ut1(i))
      end do 
 
      do i = 1, 6
         call glb_par_count( apr_eor_xy(i), mar_eor_xy(i), np, nm,
     .                       parn_eor_xy(i))
      end do 

***** Set the Earth tide number of stations (1 for global tides, and 
*     number of stations for site dependent tides)
      if( glb_glb_tides ) then
          kt = 1
      else
          kt = gnum_sites
      end if

***** Now do extended Earth tide values.

      do k = 1, kt
         do i = 1, 12
             call glb_par_count( apr_eor_etd(i,k), mar_eor_etd(i,k),
     .                       np, nm, parn_eor_etd(i,k))
         end do
      end do


***** nutation series coefficients
 
*                                 ! Loop over each term
      do i = 1, max_nut_coeff
*                                 ! In and out of phase
          do j = 1, 2
              call glb_par_count( apr_nut_coeff(j,i), 0.0,
     .                            np, nm, parn_nut_coeff(j,i) )
          end do
      end do
 
****  Now do tides ( h,l,lag)
      do i = 1, kt
*                                     ! h,l, lag
          do j = 1,3
              call glb_par_count( apr_tid(j,i), 0.0,
     .                            np, nm, parn_tid(j,i) )
          end do
      end do

***** Now UT1, xy and ETD coefficients

      do i = 1, max_ut1_coeff
         do j = 1, 2
              call glb_par_count( apr_ut1_coeff(j,i), 0.0,
     .                       np, nm, parn_ut1_coeff(j,i) )
         end do
      end do

      do i = 1, max_xy_coeff
         do j = 1, 2
              call glb_par_count( apr_xy_coeff(j,i), 0.0,
     .                       np, nm, parn_xy_coeff(j,i) )
         end do
      end do

 
***** Extended earth tide coefficients

      do k = 1, min(kt,max_etd_sites)
*                             ! Loop over each term
          do i = 1, max_etd_coeff
*                             ! In and Out of phase
              do j = 1, 2
                  call glb_par_count( apr_etd_coeff(j,i,k), 0.0,
     .                           np, nm, parn_etd_coeff(j,i,k) )
              end do
         end do 
      end do
 
***** Gamma (last one)
 
      call glb_par_count( apr_gamma, 0.0, np, nm, parn_gamma )
 
***** Now save the totals which we have
 
      num_glb_parn = np
      num_glb_mar  = nm
      num_trans_row  = nt
      num_trans_col  = nc
 
****  Thats all
      return
      end
 
