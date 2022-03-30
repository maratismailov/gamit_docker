CTITLE EST_POS

      subroutine est_pos( ep, data_type, options, an, iter )
* MOD AZ 190305: use built Module to store VMF grid values
      use vmf_type_build

      implicit none

*     This is the main routine for estimating positions of kinematic
*     sites.  Depending on the options passed it will also compute
*     adjustments to the bias parameters.

      include '../includes/const_param.h' 
      include 'track_com.h'
* MOD AZ 190305: use header to include VMF grid path
      include 'vmf_com.h'
* MOD AZ 190305: new variables that merely for function integrity
* dry_map_dump -- dry mapping function from VMF, which remains unused
* wet_map_dump -- wet mapping function from VMF, which remains unused
      double precision :: dry_map_dump, wet_map_dump

* PASSED VARIABLES
* ep -- Epoch number at which to compute position
* iter -- Iteration number (used to change the rank procedure for higher
*         iteration of searching.  Not used for normal position estimation.
* an   -- Number of ambiquity we are searching on or if we are doing a
*         float analysis the sample interval

      integer*4 ep, iter, an

* data_type -- Data types to be used (L1,L2,LC, P1, P2)
* options   -- Options for getting the position estimates.
*              RESBF  -- Resolve bias flags
*              ELEV   -- Compute elevation angles
*              KALF   -- Run as a kalman filter in that the covariance matrix
*                        is propagated forward in time.
*              ION    -- If set, ionospheric delay model also estmated
*              SMOOTH -- If covariance matrix is to be written for the
*                        smoothing run,

      character*(*) data_type, options

* LOCAL VARIABLES
* converged -- Set true when the kinematic position estimate is converged

      logical converged

* geod_coord(3) -- Geodetic coordinates
* geod_ref(3)   -- Geodetic coordinates of reference site
* dNEU(3)       -- NEU difference in position (m)
* cov_xyz(3,3)  -- Position XYZ covarinance matrix 
* cov_neu(3,3)  -- Position NEU covariance
* cov_tmp(3,3)  -- Temporary matrix
* sectag        -- Seconds tag
* fday          -- Fractional day
* curr_mjd      -- MJD of current epoch
* curr_out      -- COrrect output time with clock error included
* sec_in_day    -- Number of seconds since start of day
* av_static(3,num_site) -- Average GEOD lat/long/and height during static times
* av_square(3,num_site) -- Squares of values to compute RMS
* av_sig(3,num_site) -- Sigam of each component based on RMS
* av_init(3,num_site) -- Initial coordinates of points (improves accuracy in RMS
*            calculation)
* adj_xyz(3), adj_neu(3) -- Adjustment to XYZ and NEU to apriori final
*            output.
* dif_xyz(3), dif_neu(3) -- difference between final and initial site coords
*            in XYZ and NEU directions.
* loc_coord  -- co-lat, long height

      real*8 geod_coord(3), geod_ref(3), dNEU(3), cov_xyz(4,4),
     .       cov_neu(4,4), cov_tmp(4,4), 
     .       sectag, fday, curr_mjd, curr_out,
     .       sec_in_day, av_static(3,max_site),av_square(3,max_site),
     .       av_sig(3,max_site), av_sectag, av_mjd, av_init(3,max_site),
     .       adj_xyz(3), adj_neu(3), dif_xyz(3), dif_neu(3), 
     .       loc_coord(3)

      real*8 sig_neu(3) ! Sigams NEU (m)

* i  -- Loop counter
* na -- Parameter number for atmospheric delay
* date(5) -- Calender date of measurements
* pj,pk   -- Pointers to parameter numbers
* doy     -- Day of year
* pnt     -- Pointer to one-way used in double difference
* num_total -- Total one-ways used in forming double diffs
* num_fixd  -- Number of fixed double differences
* ks        -- Full site number of kinematic site
* av_num(num_site) -- Number of values in average
* av_date(5) -- Date of average time of static part
      integer*4 i, j,k, na, date(5), pj, pk, doy, pnt, num_total,
     .          num_fixd, ks, av_num(max_site), av_date(5)

* checked(max_obs) -- Logical array to record which oneways in the
*     double differences have been checked.
* kbit  -- Checks if bit set
* head_out -- Set true once position headers are written.

      logical checked(max_obs), kbit, add_avg, head_out

* atm_delay  -- Estimate atmospheric delay
* atm_total  -- Total estimate of atmospheric delay 
* atm_sig    -- Atmospheric delay sigma
* rot_mat(3,3) -- Rotation matrix from XYZ to NEU
* tran_mat(4,4) -- Rotation + atm delay
* dry_zen, wet_zen -- Total dry and wet delay.
* rho_ha     -- Correlation between height and zenith delay
* zend        -- Zenith distance (rads) 
* undu  -- Geoid height from GPT model.

      real*8 atm_delay, atm_total, atm_sig, rot_mat(3,3), 
     .       dry_zen, wet_zen, tran_mat(4,4), rho_ha, zend, undu 

* lat, ht  -- Latitude and height of site (used for atmospheric delay
*     calculations)
* temp_c, press, rel_hum -- Temperature, press, and rel_hum at site
*     based on seasonal model (currently)
* tbias    -- Bias in surface temperature due to inversions (not used)

      real*8 lat, long, ht, temp_C, press, rel_hum, tbias 

* sig_out -- Sigma of height output

      real*8 sig_out

* Variables for getting average position of site
* st_stopgo(5), en_stopgo(5) - Start and end time of stop go
* st_sgsec, en_sgsec - Start end and end sec

      real*8 avg_xyz(3,max_kine)
      integer*4 avg_num(max_kine), st_stopgo(5), en_stopgo(5)
      real*8 st_sgsec, en_sgsec, st_mjd

* stopgo_type -- Set to K or S depending on whether receiver is static
*     or kinematic
      character*1 stopgo_type

***** OK, The steps we need here are: initialize matrices; get the data and
*     compute the o-minus-c values based on the current apriori position;
*     update the position estimate.  Do the other steps needed as given
*     in the options. 

* MOD TAH 040316: See if kinematic/static mode changes at this epoch
*     (only important if stopgo_mode is .true.)

      do i = 1, num_site
         if( kbit(data_flag_cse(1,i,ep),6) ) then
*            Kinematic moving mode set.  Set the number of static_obs
*            to zero
             static_obs(i) = 0
         else
*            Increment the number of static obs.  Only increment if
*            data is available.
             if( kbit(data_flag_cse(1,i,ep),1) ) 
     .  	  static_obs(i) = static_obs(i) + 1
             if( static_obs(i).eq.1 ) then
                 stopgo_point(i) = stopgo_point(i) + 1
             endif 
         endif
      enddo
  

*     Initialize the solution matrices, and the apriori coordinates
      call init_mat(ep, options, data_type, an)

*     Clear the averaging values
      if( ep.eq.1 ) then
          do i = 1,num_kine
             do j = 1,3
                avg_xyz(j,i) = 0.d0
             end do
             avg_num(i) = 0
          end do
      end if
      if ( ep.eq.1 ) then
         do i = 1, num_site
            stopgo_point(i) = 0
         end do
      end if

****  See if we should write the covariance matrix and solution
      if( index(options,'SMOOTH').gt.0  ) then
          if( an.gt.0 ) then 
              call wr_cov_parm(ep,'R',data_type)
              call smooth( ep )
          end if
      end if

*     See how many data types are used.
      num_dtype = 0
      if( index(data_type,'L1').gt.0 ) num_dtype = num_dtype + 1
      if( index(data_type,'L2').gt.0 ) num_dtype = num_dtype + 1
      if( index(data_type,'LC').gt.0 ) num_dtype = num_dtype + 1
      if( index(data_type,'P1').gt.0 ) num_dtype = num_dtype + 1
      if( index(data_type,'P2').gt.0 ) num_dtype = num_dtype + 1
      if( index(data_type,'PC').gt.0 ) num_dtype = num_dtype + 1


      if( num_dtype.gt.max_dtype ) then
*         Too many data types selected.  Kill run
          write(*,110) num_dtype, data_type, max_dtype
 110      format('**DISASTER** Too many types ',i4,' from ',a,
     .           ' Maximum allowed is only ',i1)
          stop '**DISASTER** Too many data types'
       end if


****  If kinematic positions have not been initially estimated estimate
*     them now.  If only P1 has been selected this is the same solution.
*     In get_kine_apr we use the apriori position given by the user for
*     first epoch, and the previous epochs positions other epochs unless
*     kine_known flag is set, in which case we use the value based initially
*     P1 Solution
      call get_kine_apr(ep)

      if( ep.ge.debug_start .and. ep.le.debug_end ) then 
         write(*,160) ep, kine_known, kbit(kine_OK(1,i),ep),
     .             (i,(curr_site_xyz(j,kine_to_site(i)),j=1,3),
     .             i = 1,num_kine)
 160     format(/,60('-'),/,
     .          'CURR_SITE_XYZ ',I7,1x,2l1,1x,6(I3,3F15.3))
      end if


****  Now we may need to iterate on the solution estimate if the apriori
*     is too far in error.
      converged = .false.
      do while ( .not. converged )
          call update_pos(ep, data_type, options, converged, an,iter)
 
      end do 

      if( index(options,'NoWrite') .eq.0 ) then

*        At this epoch see how many biases we have not fixed.
         num_fixd = 0
         num_total = 0
         head_out = .false.
         do i = 1, max_obs
            checked(i) = .false.
         end do

         do i = 1, num_dd
*           Check the one way to see if fixed 
            do k = 1,4
               pnt = dd_des(k,i)

*              Find the ambiquity number for this one-way
               do j = 1, num_ambs
                  if( bf_ents(1,j).eq.ow_des(1,pnt) .and.
     .                bf_ents(2,j).eq.ow_des(2,pnt) .and.
     .                ep.ge.bf_ents(3,j) .and. 
     .                ep.le.bf_ents(4,j) .and. .not.checked(j)) then
*                     Found the correct one.  See if fixed
                      checked(j) = .true.
                      num_total = num_total + 1
                      if( kbit(bf_ents(5,j),2) ) then
                          num_fixd = num_fixd + 1
                      end if
                  end if
               end do
            end do
         end do

*        See if we need to write atmospheric delay from FIXED sites
         do i = 1, num_site

*           Compute geodetic position
            call XYZ_to_GEOD( rot_mat, site_apr(1,i), geod_coord )

            lat = pi/2.d0-geod_coord(1)
            long = geod_coord(2)
            ht = geod_coord(3)
            if( site_type(i).eq.0 .and. atm_parn(i).gt.0 ) then
                na = atm_parn(i)
                if( index(options,'SMOOTH').gt.0 ) then
                    atm_delay = sol_vec_sav(na)
                    atm_sig = sqrt(cov_parm_sav(na,na))
                else
                    atm_delay = sol_vec(na)
                    atm_sig = sqrt(cov_parm(na,na))
                end if

* MOD TAH 0703515: Compute apriori delay
                rel_hum = 0.d0

*               Get the time: For output we want GPS time (not time on
*               the receiver).
                curr_mjd = ref_start + ((ep-1)*usr_interval)/86400.d0 
     .                + sec_offset(i,ep)/86400.d0
                call mjd_to_ymdhms(curr_mjd, date, sectag)
                call ymd_to_doy( date, doy)
                sec_in_day = date(4)*3600.d0 + 
     .                       date(5)*60.d0 + sectag
               
                if ( atm_mtt ) then
                   call met_seasonal( temp_c, press, rel_hum, tbias, 
     .                   curr_mjd,lat, ht ) 
                                  
*                  Now get the zenith delay terms and map to the elevation angle
                   if( .not.use_atm_ref ) then 
                       call dry_saas_zen( temp_C, lat, ht, press,  
     .                                    dry_zen)
                       call wet_saas_zen( temp_C, lat, ht, rel_hum,
     .                                    wet_zen)
                   else
                       call get_atm( i, curr_mjd, dry_zen, wet_zen, 
     .                            lat, long, ht )
                   end if
                elseif ( atm_gmf ) then
*                  Use GPT and GMF
                   call gpt(curr_mjd,lat,long,ht, press, temp_C, undu) 
                   if( .not.use_atm_ref ) then 
                       call dry_saas_zen( temp_C, lat, ht, press,   
     .                                    dry_zen)
                       rel_hum = ref_rel_hum
                       call wet_saas_zen( temp_C, lat, ht, rel_hum, 
     .                                    wet_zen)
                   else
                       call get_atm( i, curr_mjd, dry_zen, wet_zen, 
     .                               lat, long, ht )
                   end if
* MOD AZ 190305: Additional VMF option and call of related function
                elseif ( atm_vmf ) then
*                  Use VMF
                  if( .not.use_atm_ref ) then
                   call vmf1_grid( VMF1_grid_file, curr_mjd, lat,
     .                             long, ht, 0, dry_map_dump,
     .                             wet_map_dump, dry_zen, wet_zen )
                  else
                   call get_atm( i, curr_mjd, dry_zen, wet_zen,
     .                           lat, long, ht )
                  end if
                else
* MOD AZ 190305: Sign of problematic implementation
                   write(*,399)
 399               format('ERROR: ATM Model wrongly read',
     .                    'ATM_MODELC command')
                   stop 'TRACK: ATM Commands out of order'

                endif

                atm_total = dry_zen + wet_zen + atm_delay

                fday = doy*1.d0 + date(4)/24.d0 + date(5)/1440.d0 +
     .                 sectag/86400.d0 

c                write(*, 400) date(1), doy, sec_in_day, 
c     .                        atm_total*1000, atm_delay*1000., 
c     .                        atm_sig*1000.0,  fday, ep
 400            format('FIXD ATM ', i5,i5,f11.2,' ZEN (mm) ',3F9.2,
     .                1x,f15.9,1x,i6)
* MOD TAH 100712: Output values to data files 
                if( index(posit_type,'DHU').gt.0 .and.iter.ne.0 ) then
                   write(21+3*num_site+i, 850) date, sectag, 
     .                 (0.00, 0.00,j=1,3), 
     .                  0.00 , 0,
     .                  atm_delay*1000., atm_sig*1000.0, fday, ep,
     .                  0, 0, 'F', 0.00
                endif
                if( index(posit_type,'GEOD').gt.0 .and. iter.ne.0 ) then
                   write(21+i, 800) date(1), doy, sec_in_day, 
     .                 90.d0-geod_coord(1)*180/pi, 
     .                 geod_coord(2)*180/pi, geod_coord(3),
     .                 (0.00,j=1,3),  0.00, 0,atm_total*1000., 
     .                  atm_sig*1000.0, fday, ep, 0, 0, 'F' 
                end if
                if( index(posit_type,'NEU').gt.0 .and.iter.ne.0 ) then
                   write(21+num_site+i, 830) date, sectag, 
     .                (0.00, 0.00 ,j=1,3),  0.00, 0,
     .                atm_delay*1000., atm_sig*1000.0, fday, ep,
     .                 0, 0, 'F', 0.00 
                end if  

                if( index(posit_type,'XYZ').gt.0 .and.iter.ne.0 ) then
                   write(21+2*num_site+i, 870) date, sectag, 
     .                (0.00, 0.00 ,j=1,3),  0.00, 0,
     .                atm_delay*1000., atm_sig*1000.0, fday, ep,
     .                 0, 0, 'F', 0.00 
                end if  

            end if

****        Output initial position is this site is fixed:
            if( site_type(i).eq. 0 .and. ep.eq.num_epochs) then
*               Write header
                if( .not.head_out ) then
                    write(lus,500) 
 500                format(/,'Position Estimates and Adjustments ',
     .                     ' at epoch',/, 
     .                     'SITE      X (m)        dX      +-    ',
     .                     '    Y (m)        dX      +-        ',
     .                     'Z (m)        dX      +-      Lat (deg)',
     .                     '       dN (m)   +-    Long (deg)       ',
     .                     'dE (m)    +-   Height (m)    dH      +- ') 

                    head_out = .true.
                endif
                write(lus,510) site_names(i),
     .                       (site_apr(j,i),0.0,0.0,j=1,3),
     .                        lat*180/pi, 0.0, 0.0,
     .                        long*180/pi,  0.0, 0.0,
     .                        ht, 0.0, 0.0
 510            format(a4,3(F14.4,1x,F8.4,1x,F7.4),1x,
     .                 2(F15.9,1x,F8.4,1x,F7.4),1x,
     .                   F10.4,1x,F8.4,1x,F7.4)
            end if

         end do
*
*        Now do kimenatic sites.      
         do i = 1, num_kine

*           Get the current time of the position determination
*           For output we want actual time
            curr_mjd = ref_start + ((ep-1)*usr_interval)/86400.d0 
     .               + sec_offset(kine_to_site(i),ep)/86400.d0
            call mjd_to_ymdhms(curr_mjd, date, sectag)
  
            call ymd_to_doy( date, doy)
            sec_in_day = date(4)*3600.d0 + 
     .                   date(5)*60.d0 + sectag

            fday = doy*1.d0 + date(4)/24.d0 + date(5)/1440.d0 +
     .             sectag/86400.d0 

            ks = kine_to_site(i)

*           Now compute the output position information.
*           If we are doing a smoothing run, then replace the adjustment
*           in kine_out
            if( index(options,'SMOOTH').gt.0 ) then
                do j = 1, 3
                   kine_out(j,i) = kine_out(j,i) - sol_vec((i-1)*3+j) +
     .                             sol_vec_sav((i-1)*3+j)
                   adj_xyz(j) = sol_vec_sav((i-1)*3+j)
                end do
            else
                do j = 1, 3
                  adj_xyz(j) = sol_vec((i-1)*3+j)
                end do
            end if
  
            call XYZ_to_geod(rot_mat, kine_out(1,i), geod_coord)

****        Compute adj to NEU
            do j = 1,3
               adj_neu(j) = 0.d0
               do k = 1,3
                   adj_neu(j) = adj_neu(j)+rot_mat(j,k)*adj_xyz(k)
               end do
            end do


****        Add value into average.  See if consistent with previous values
            if( avg_num(i).gt. 0 ) then
                if( abs(kine_out(1,i)-avg_xyz(1,i)/avg_num(i))
     .                                             .lt.1.d0 ) then
                   add_avg = .true.
                else
                   add_avg = .false.
                end if
            else
                add_avg = .true.
            end if
            if( add_avg ) then
                do j = 1, 3
                    avg_xyz(j,i) = avg_xyz(j,i)+kine_out(j,i)
                end do
                avg_num(i) = avg_num(i) + 1
            end if

*           Pull of the covariance matrix elements
            do j = 1,3
               pj = pos_parn(j,i)
               do k = 1, 3
                  pk = pos_parn(k,i)
                  if( index(options,'SMOOTH').gt.0 ) then
                      cov_xyz(j,k) = cov_parm_sav(pj,pk)
                  else
                      cov_xyz(j,k) = cov_parm(pj,pk)
                  end if
               end do
               pk = atm_parn(ks)
               if( pk.gt.0 ) then
                  if( index(options,'SMOOTH').gt.0 ) then
                      cov_xyz(j,4) = cov_parm_sav(pj,pk)
                      cov_xyz(4,j) = cov_parm_sav(pj,pk)
                  else
                      cov_xyz(j,4) = cov_parm(pj,pk)
                      cov_xyz(4,j) = cov_parm(pj,pk)
                  end if
               else
                  cov_xyz(j,4) = 0.d0
               endif
            end do
            pk = atm_parn(ks)
            if( pk.gt.0 ) then
               if( index(options,'SMOOTH').gt.0 ) then
                   cov_xyz(4,4) = cov_parm_sav(pk,pk)
               else
                   cov_xyz(4,4) = cov_parm(pk,pk)
               end if
            else
               cov_xyz(4,4) = 1.d-10
            endif

*           Make the trans matrix
            do j = 1,3
               do k = 1,3
                  tran_mat(j,k) = rot_mat(j,k)
               end do
            end do
            do j = 1,4
               tran_mat(j,4) = 0.d0
               tran_mat(4,j) = 0.d0
            enddo
            tran_mat(4,4) = 1.d0
               


*           Get the NEU variance-covariance
C           call var_comp(rot_mat, cov_xyz, cov_neu, cov_tmp, 3,3,1)
            call var_comp(tran_mat, cov_xyz, cov_neu, cov_tmp, 4,4,1)

            rho_ha = cov_neu(3,4)/sqrt(cov_neu(3,3)*cov_neu(4,4))

*           See if atmospheric delay estimated
            if( atm_parn(ks).gt.0 ) then
                na = atm_parn(ks)
                if( index(options,'SMOOTH').gt.0 ) then
                    atm_delay = sol_vec_sav(na)
                    atm_sig = sqrt(cov_parm_sav(na,na))
                else
                    atm_delay = sol_vec(na)
                    atm_sig = sqrt(cov_parm(na,na))
                end if

* MOD TAH 111117: If this scale convert to delay change
                if( atm_scale(ks) ) then
                     atm_delay = atm_delay*
     .                  (geod_coord(3)-site_geod(3,1))*1.d-3
                     atm_sig = atm_sig*
     .                  abs((geod_coord(3)-site_geod(3,1))*1.d-3)
                end if

            else
                atm_delay = 0.d0
                atm_sig = 0.d0
            end if

            sig_out = sqrt(cov_neu(3,3))

* MOD TAH 0401316: See what the statis if
            stopgo_type = 'K'
            if( static_obs(ks).gt.0 ) stopgo_type = 'S'
*           If first epoch save the initial site coordinates
            if( ep.eq.1 ) then
                do j = 1,3
                   av_init(j,ks) = geod_coord(j)
                enddo
            end if

            if( index(posit_type,'GEOD').gt.0 .and.
     .          iter.ne.0 ) then
               if ( sig_out.le.out_sig_limit) then

* MOD TAH 0703515: Compute apriri delay
                 lat = pi/2.d0-geod_coord(1)
                 long = geod_coord(2)
                 ht = geod_coord(3)
                 rel_hum = 0.d0

c                call ymdhms_to_mjd(date,sectag, trn_time)
                 if( atm_mtt ) then 
                    call met_seasonal( temp_c, press, rel_hum, tbias, 
     .                      curr_mjd,lat, ht ) 
                                  
*                   Now get the zenith delay terms and map to the 
*                   elevation angle
                    if( .not.use_atm_ref ) then 
                       call dry_saas_zen( temp_C, lat, ht, press,
     .                                    dry_zen)
                       call wet_saas_zen( temp_C, lat, ht, rel_hum,
     .                                    wet_zen)
                    else
                       call get_atm( ks, curr_mjd, dry_zen, wet_zen, 
     .                            lat, long, ht )
                    end if
                 elseif ( atm_gmf ) then
*                  Use GPT and GMF
                   call gpt(curr_mjd,lat,long,ht, press, temp_C, undu) 
                   if( .not.use_atm_ref ) then 
                       call dry_saas_zen( temp_C, lat, ht, press,   
     .                                    dry_zen)
                       rel_hum = ref_rel_hum
                       call wet_saas_zen( temp_C, lat, ht, rel_hum, 
     .                                    wet_zen)
                   else
                       call get_atm( ks, curr_mjd, dry_zen, wet_zen, 
     .                               lat, long, ht )
                   end if
                 elseif ( atm_vmf ) then
*                  Use VMF
                    if( .not.use_atm_ref ) then
                    call vmf1_grid( VMF1_grid_file, curr_mjd, lat,
     .                             long, ht, 0, dry_map_dump,
     .                             wet_map_dump, dry_zen, wet_zen )
                    else
                    call get_atm( ks, curr_mjd, dry_zen, wet_zen,
     .                               lat, long, ht )
                    end if
                 else
                    write(*,799)
 799                format('ERROR: ATM Model wrongly read',
     .                    'ATM_MODELC command')
                    stop 'TRACK: ahha2 Commands out of order'

                 endif 
 
                 atm_total = dry_zen + wet_zen + atm_delay

                 write(21+ks, 800) date(1), doy, sec_in_day, 
     .             90.d0-geod_coord(1)*180/pi, 
     .             geod_coord(2)*180/pi, geod_coord(3),
     .             (sqrt(cov_neu(j,j))*100.0,j=1,3), 
     .             wrms_dd_site(ks), num_dd_site(ks),
     .             atm_total*1000., atm_sig*1000.0, fday, ep,
     .             num_total, num_total-num_fixd, stopgo_type 
 800           format(i5,i5,f15.6,2f15.9,f12.4,3f6.1,f6.1,i3,
     .                2F9.2,1x,f17.11,1x,i6,2I4,1x,a)
               endif
           end if

*          Accumulate averages and/or output the values.  Also check
*          the last epoch.
           if( (static_obs(ks).eq.0 .or. ep.eq.num_epochs) .and.
     .          .not.static(i)  ) then
*              Changed from static to kinematic.  Output any 
*              old values
               if( av_num(ks).gt.0 ) then
                   do j = 1,3
                      av_static(j,ks) = av_static(j,ks)/av_num(ks) 
                      av_sig(j,ks) = sqrt(
     .                              (av_square(j,ks)/av_num(ks)
     .                              -av_static(j,ks)**2)/av_num(ks))
                      av_static(j,ks) = av_static(j,ks) + 
     .                                  av_init(j,ks)
                   end do
                   av_mjd = curr_mjd - av_num(ks)*usr_interval/
     .                                 (2*86400.d0)
                   call mjd_to_ymdhms(av_mjd,av_date, av_sectag)
* MOD TAH 050123: Change to give start and stop times
                   st_mjd = curr_mjd - av_num(ks)*usr_interval/
     .                                 86400.d0
                   call mjd_to_ymdhms(curr_mjd, en_stopgo, en_sgsec)
                   call mjd_to_ymdhms(st_mjd, st_stopgo, st_sgsec)

                   write(lus, 810) site_names(ks), stopgo_point(ks), 
     .                 st_stopgo, st_sgsec, 
     .                 en_stopgo(4), en_stopgo(5), en_sgsec,
     .                 90.d0-av_static(1,ks)*180/pi, 
     .                 av_static(2,ks)*180/pi, av_static(3,ks),
     .                 (av_sig(j,ks)*180/pi,j=1,2), av_sig(3,ks),
     .                  av_num(ks)
 810               format('STATIC ',a4,' PT',i5.5,1x,i4,4i3,1x,
     .                1x,f5.1,' to ',2i3,1x,f5.1,
     .                ' LLH ',2f15.9,f11.4,1x,
     .               ' +- ',2F12.9,' deg ',F8.4,' (m), Num ',i5)
*                  Now clear the values
                   do j = 1,3
                      av_static(j,ks) = 0.d0
                      av_square(j,ks) = 0.d0
                   end do
                   av_num(ks) = 0
               end if
            else
****           Accumulate values
               do j = 1,3 
                  av_static(j,ks) = av_static(j,ks) + 
     .                   (geod_coord(j)-av_init(j,ks))
                  av_square(j,ks) = av_square(j,ks) + 
     .                   (geod_coord(j)-av_init(j,ks))**2
               enddo
               av_num(ks) = av_num(ks)+1

            endif
 
            if( index(posit_type,'NEU').gt.0 .and.iter.ne.0
     .          .and. sig_out.lt.out_sig_limit  ) then
*              Get the NEU difference in position 
* MOD TAH 0901112: Added use of ref_xyz instead of first site
               if( ref_xyz(1).eq.0 ) then
                   do j = 1,3
                      ref_xyz(j) = site_apr(j,1)
                   end do
               end if

               call XYZ_to_geod(rot_mat, ref_xyz, geod_ref)
               dNEU(1) = -earth_rad*(geod_coord(1)-geod_ref(1))
               dNEU(2) =  earth_rad*sin(geod_ref(1))*
     .                              (geod_coord(2)-geod_ref(2))
               dNEU(3) = (geod_coord(3)-geod_ref(3))
               write(21+num_site+ks, 830) date, sectag, 
     .             (dNEU(j), sqrt(cov_neu(j,j)),j=1,3), 
     .             wrms_dd_site(ks), num_dd_site(ks),
     .             atm_delay*1000., atm_sig*1000.0, fday, ep,
     .             num_total, num_total-num_fixd, stopgo_type, rho_ha  
 830           format(i5,4i3,1x,F10.6, 3(F15.4,1x,F9.4,1x),
     .              F8.2,1x,i3,1x, 2F9.2,1x,f17.11,1x,i6,2i4,1x,a,
     .              1x,F6.3)
            end if
            if( index(posit_type,'DHU').gt.0 .and.iter.ne.0
     .          .and. sig_out.lt.out_sig_limit  ) then
*              Rotate the dXYZ from the initial coordimates into
*              dNEU
               do j = 1,3
                   dif_xyz(j) = kine_out(j,i) - site_apr(j,ks)
                   sig_neu(j) = sqrt(cov_neu(j,j))
               end do
               call rotate_geod(dif_xyz, dif_neu, 'XYZ','NEU',
     .               site_apr(1,ks) , loc_coord, rot_mat)
               write(21+3*num_site+ks, 850) date, sectag, 
     .             (dif_NEU(j)*1000, sig_neu(j)*1000,j=1,3), 
     .              wrms_dd_site(ks), num_dd_site(ks),
     .             atm_delay*1000., atm_sig*1000.0, fday, ep,
     .             num_total, num_total-num_fixd, stopgo_type, rho_ha 
 850           format(i5,4i3,1x,F10.6, 3(F15.2,1x,F9.2,1x),
     .              F8.2,1x,i3,1x, 2F9.2,1x,f17.11,1x,i6,2i4,
     .              1x,a,1x,F6.3)
            endif

            if( index(posit_type,'XYZ').gt.0 .and.iter.ne.0
     .          .and. sig_out.lt.out_sig_limit  ) then
*              Output XYZ coordinates

               write(21+2*num_site+ks, 870) date, sectag, 
     .             (kine_out(j,i), sqrt(cov_xyz(j,j)), j=1,3), 
     .              wrms_dd_site(ks), num_dd_site(ks),
     .             atm_delay*1000., atm_sig*1000.0, fday, ep,
     .             num_total, num_total-num_fixd, stopgo_type, rho_ha 
 870           format(i5,4i3,1x,F10.6, 3(F15.4,1x,F9.4,1x),
     .              F8.2,1x,i3,1x, 2F9.2,1x,f17.11,1x,i6,2i4,
     .              1x,a,1x,F6.3)
            endif

*           Accummulate the phase residual statistics and 
*           if write_res is true, write the residuals out
            if( iter.ne.0 ) then
                call out_resid(ep, i)
            end if

*           If this is the last epoch output average values
            if( ep.eq.num_epochs ) then
                if( avg_num(i).gt.0 ) then
                   write(*,900)  site_names(kine_to_site(i)), 
     .                  avg_num(i), (avg_xyz(j,i)/avg_num(i),j=1,3)
 900               format('SITE ',a8,' Average XYZ from ',i8,' epochs ',
     .                     3F15.4,' m')
                else
                   write(*,910) site_names(kine_to_site(i))
 910               format('SITE ',a8,' No average position ')
                endif

****            OK: Now output the final positions with adjustments
                write(lus,510) site_names(ks),
     .                (kine_out(j,i),adj_xyz(j),sqrt(cov_xyz(j,j)),
     .                 j=1,3),
     .                 90.d0-geod_coord(1)*180/pi,adj_neu(1),
     .                 sqrt(cov_neu(1,1)),
     .                 geod_coord(2)*180/pi, adj_neu(2),
     .                 sqrt(cov_neu(2,2)),
     .                 geod_coord(3),adj_neu(3), 
     .                 sqrt(cov_neu(3,3))
            end if  
C           write(301,915) ep, site_names(ks),
C     .           (kine_out(j,i),adj_xyz(j),sqrt(cov_xyz(j,j)),
C     .            j=1,3)
C 915       format(I5,1x,a4,3(1x,F16.4,1x,F10.4,1x,F10.4,1x))
         end do
      end if

****  Thats all
      return 
      end
     
CTITLE INIT_MAT

      subroutine init_mat( ep, options, data_type, ani )

      implicit none

*     Routine to initialize the estimation matrices needed for getting position
*     esimates
* MOD TAH 100518: Included changing the back solution start so that the aprori
*     position sigams are constistent with the propagated process noise. 
*     (needed for the timedep_procns command)
* MOD TAH 100710: Modified time dependence timedep_procns so that only last one
*     is applied to a site

      include 'track_com.h'

* PASSED VARIABLES
* ep -- Epoch number at which to compute position.  Used to see if this is the
*       first epoch
* ani -- Sampling being used in the solution.  If an is negative then solution
*       is running backwards

      integer*4 ep, ani

* options   -- Options for getting the position estimates.
*              RESBF  -- Resolve bias flags
*              ELEV   -- Compute elevation angles
*              ION    -- Initialize matrices for ion delay estimation
*              FLOAT  -- Floating ambiquities for non-resolved ambiquities
* data_type  -- Type of data to use (L1,L2/LC). If LC used, estimate is 
*               non-integer in L1 slot

      character*(*) options, data_type

* LOCAL VARIABLES
* i,j   -- Loop counters
* ks    -- Site number (from full list) of the current kinematic site
* pn    -- Generic parameter number
* ne    -- Number of ambiguity estimates to output
* as    -- Ambiquity parameter "slot" parameter number to use
* an    -- Absolute value of an

* Lim_amb -- Limit on replacing ambiguities (non-zero when smooth)
* kdep  -- Index for site dependent noise

      integer*4 i,j, k, ks, pn, ne, as, an, lim_amb, kdep

* pos_max_parm -- Possible maximum number of parameters (Added 090221
*     to make sure cov_parm matrix is initialized correctly when row
*     and column are cleared).  Value is num_parm+ne*num_amb
*     (ne = 1 for LC and ne = 2 for L1+L2)
      integer*4 pos_max_parm

* base_len -- real*8 function to return baseline length between two sites
* blen     -- Length between two sites
* curr_mjd -- nominal time

      real*8 base_len, blen, curr_mjd

      real*8 dr, drt  ! Computed from radial change and dt (for mar_atm_hgt)

      real*8 apr_site_prop(3)  ! XYZ propagated apriori sigma
      integer*4 dep   ! Difference in epoch numbers accross time=dependent
                      ! interval. (Needed for variance at end of data)

* kbit   -- Checks status of bit

      logical kbit

****  OK, clear and initialize the matrixes.
      if( ep.eq.0 ) RETURN
      curr_mjd = ref_start + ((ep-1)*usr_interval)/86400.d0 

C     if( index(options, 'KALF').eq.0 .or. ep.eq.1 ) then
      if( index(options, 'KALF').eq.0 .or. 
     .    ((ep.eq.1 .and. ani.gt.0) .or. 
     .     (ep.eq.num_epochs .and. ani.lt.0) ) ) then

*         Clear the matrices and assign the apriori variances
*         to to the covariance matrix.
          do i = 1, num_parm
             do j = 1, num_parm
                cov_parm(i,j) = 0.d0
             end do
* MOD TAH 090219: Put weak constraint on all parameters
             cov_parm(i,i) = 1.d2
          end do

*         Now assign the apriori variances to the diagonal.
          do i = 1, num_kine
             ks = kine_to_site(i)
*            Compute what we should use; start with apriori
*            and then add
             if( .not. (ep.eq.1 .and. ani.gt.0) ) then
*                Propagate the noise to the end of the data
                 do j = 1,3
                    apr_site_prop(j) = apr_site(j,ks) +
     .                  num_epochs*mar_site(j,ks)
                 end do
*                See if time dependent process noise
                 kdep = 0
                 do k = 1, num_tdep
                    if( nst_tdep(k).eq.-1 .or. 
     .                  nst_tdep(k).eq.ks ) then
                        dep = (mjd_tdep(2,k)-mjd_tdep(1,k))*86400.0d0/
     .                         usr_interval
                        kdep = k
                    endif
                 end do
* MOD TAH 100710: Only apply the last element
                 if( kdep.gt.0 ) then 
                    do j = 1,3
                        apr_site_prop(j) = apr_site_prop(j) +
     .                      abs(dep)*mar_tdep(j,kdep)
                    end do
                  end if

*                 Now make sure values are not too large
                  do j = 1,3
                     if( apr_site_prop(j) .gt.10.0 ) 
     .                 apr_site_prop(j) = 10
****                 If user gave aposterori sigmas then use these
*                    values
                     if( pst_site(j,i).ge.0 ) 
     .                 apr_site_prop(j) = pst_site(j,i)

                  end do 
             end if

             do j = 1, 3
                pn = pos_parn(j,i) 
*               If this is a forward run just use the apr_site value
*               but if back solution use the accumlated noise (set
*               upper limit).
                if( (ep.eq.1 .and. ani.gt.0) ) then 
                   cov_parm(pn,pn) = apr_site(j,ks)
                else
                   cov_parm(pn,pn) = apr_site_prop(j)
                end if
             end do
          end do

*         Assign the variances for the atmospheric delay
          do i = 1, num_site
             pn = atm_parn(i)
             if( pn.gt.0 ) then
                 cov_parm(pn,pn) = apr_atm(i)
             end if
          end do

*****     Now initialize the solution vector
          do i = 1, num_parm
             sol_vec(i) = 0.d0
          end do

****      Set the number of static observations to zero 
          do i = 1, num_site
             static_obs(i) = 0
          end do

      else 
*         We are running a Kalman filter solutions; in this case
*         we propagate the covariance matrix forward.
          do i = 1, num_kine
             ks = kine_to_site(i)
             do j = 1, 3
                pn = pos_parn(j,i)
* MOD TAH 040316: Only add process noise if number of static obs is
*               1 or less
                if( static_obs(ks).le.1 .or. 
     .             .not. stopgo_mode ) then 
                   cov_parm(pn,pn) = cov_parm(pn,pn) + mar_site(j,ks)
                else
                   cov_parm(pn,pn) = cov_parm(pn,pn) + 
     .                                     mar_site(j,ks)/stopgo_dvar
                endif
                
             end do
* MOD TAH 100515: See if we need to add additional time dependent
*            process noise
* MOD TAH 100710: Only apply the last entry
             kdep = 0
             do k = 1, num_tdep
                if( nst_tdep(k).eq.-1 .or. nst_tdep(k).eq.ks ) then
                    if( curr_mjd.ge.mjd_tdep(1,k) .and. 
     .                  curr_mjd.le.mjd_tdep(2,k) ) then
****                   OK Correct site and in time range so add noise
                       kdep = k
                    end if
                end if
             end do
* MOD TAH 1007110: See if timedep_procn is to be added
             if( kdep.gt.0 ) then
                do j = 1, 3
                   pn = pos_parn(j,i)
                   if( pn.gt.0 ) 
     .             cov_parm(pn,pn) = cov_parm(pn,pn) + 
     .                               mar_tdep(j,kdep)
                end do
             end if
          end do

*         Assign the variances for the atmospheric delay
          do i = 1, num_site
             pn = atm_parn(i)
             if( pn.gt.0 ) then

                 cov_parm(pn,pn) = cov_parm(pn,pn) + mar_atm(i)
             end if
          end do

****      For kinematic sites see if there is a height rate term
          do i = 1, num_kine
             ks = kine_to_site(i)
             pn = atm_parn(ks)
*            See if we have a height rate tems
             if( mar_atm_hgt(ks).gt.0 .and. pn.gt.0 ) then

*               Get height rate
                if( ep.eq.1 ) then
                    dr = abs(sqrt(kine_xyz(1,i,ep+1)**2+
     .                            kine_xyz(2,i,ep+1)**2+
     .                            kine_xyz(3,i,ep+1)**2)-
     .                       sqrt(kine_xyz(1,i,ep  )**2+
     .                            kine_xyz(2,i,ep  )**2+
     .                            kine_xyz(3,i,ep  )**2))
                    drt = usr_interval
                elseif (ep.eq.num_epochs) then
                    dr = abs(sqrt(kine_xyz(1,i,ep  )**2+
     .                            kine_xyz(2,i,ep  )**2+
     .                            kine_xyz(3,i,ep  )**2)-
     .                       sqrt(kine_xyz(1,i,ep-1)**2+
     .                            kine_xyz(2,i,ep-1)**2+
     .                            kine_xyz(3,i,ep-1)**2))
                    drt = usr_interval
                else
                    dr = abs(sqrt(kine_xyz(1,i,ep+1)**2+
     .                            kine_xyz(2,i,ep+1)**2+
     .                            kine_xyz(3,i,ep+1)**2)-
     .                       sqrt(kine_xyz(1,i,ep-1)**2+
     .                            kine_xyz(2,i,ep-1)**2+
     .                            kine_xyz(3,i,ep-1)**2))
                    drt = 2*usr_interval
                endif
                cov_parm(pn,pn) = cov_parm(pn,pn) + 
     .                            mar_atm_hgt(ks)*(dr/drt)**2
c                print *,'ATM: ',ep, dr/drt, 
c     .               sqrt(mar_atm_hgt(ks)*(dr/drt)**2)*1000
             endif
          enddo

*****     Now initialize the solution vector.  In this case we only
*         zero the elements of the kinematic position if the stochastics
*         are large (in which case the apriori is being updated with
*         the previous solution).  The atmospheric and bias parameters
*         delays are propagated forward in the analysis
          do i = 1, num_kine
             do j = 1, 3
                if( .not.static(i) ) then 
                    sol_vec(3*(i-1)+j) = 0.d0
                end if
             end do
          end do
      end if

****  See if ION is being estimated.  If so set the apriori covariance
*     matrix and solution vector
      if( index(options,'ION').gt.0 ) then

*         Loop over the pairs of stations
          do i = 1, num_site
             ion_vec(i) = 0.d0
             do j = 1, num_site

*               Get the baseline length
                blen = base_len(ep, i,j)
                ion_parm(i,j) = ion_var + 
     .                          ion_corr*exp(-(blen/ion_tau)**2)
             end do
          end do
      end if
                                                    
****  Set up for float solution.  Here the parameter list changes with 
*     time as biases are added and removed. Only execute at start.
      if( index(options,'FLOAT') .gt. 0 .and.
     .    ((ep.eq.1 .and. ani.gt.0) .or. 
     .     (ep.eq.num_epochs .and. ani.lt.0) ) ) then

*         See how many ambiquities we need to estimate
          an = abs(ani)
          ne = 0 
          if( index(data_type,'L1').gt.0 .or.
     .        index(data_type,'LC').gt.0     ) ne = ne + 1
          if( index(data_type,'L2').gt.0     ) ne = ne + 1

* MOD TAH 090221: Compute maxium number of parameters
          pos_max_parm = min(num_parm + ne*num_ambs,max_parm)

*         First check to see if we have any ambiguity parameters that we
*         have finished data with. Check the limits for the two possible
*         directions.
* MOD TAH 071028: A quick fix that needs update later
c         lim_amb = 0
c         if( index(options,'SMOOTH').gt.0 ) then
c             lim_amb = num_ambs
c         endif
* MOD TAH 080923: Something wrong with setting lim_amb so revert to using
*         old code which uses num_ambs
c         lim_amb = num_ambs

c         do i = 1, lim_amb   ! num_ambs
C            if( ep.gt.bf_ents(4,i) .and. amb_parn(1,i).gt.0 ) then
c            if( (ep.gt.bf_ents(4,i) .and. amb_parn(1,i).gt.0 .and.
c    .            ani.gt.0) .or.
c    .           (ep.lt.bf_ents(3,i) .and. amb_parn(1,i).gt.0 .and.
c    .            ani.lt.0) ) then

*                OK, we have passed the end of the ambiguity.  See how
*                many values to output
c                do j = 1, ne
c                   ambiq_float(j,i) = sol_vec(amb_parn(j,i))
c                   ambiq_flsig(j,i) = sqrt(cov_parm(amb_parn(j,i), 
c    .                                               amb_parn(j,i))) 
c                end do

c                if( ep.eq.bf_ents(4,i) .or. 
c    .               ep.eq.bf_ents(3,i) )
c    .           write(*,410) i, ep, (bf_ents(j,i),j=1,5),
c    .                amb_parn(1,i), (sol_vec(amb_parn(j,i)),
c    .                sqrt(cov_parm(amb_parn(j,i), amb_parn(j,i))),
c    .                         j = 1, ne)
c410             format('BF_FLOAT ',I4,i7,i3,' PRN ',I2.2,1x,2I7,1x,
c    .                  I3,' NP ',I3,' Estimated ADJUST. (m) ',
c    .                  2(F10.3,1x,F10.3))

****             Now clear the entry to allow it be reused. NOTE: the
*                covariance and solution are re-initialized when new
*                parameter introduced
c                do j = 1, ne
*                   Now clear the pointer so that it can used again
* MOD TAH 090110: Removed clearing parameter slot
C                   parn_to_amb(amb_parn(j,i)) = 0
*                   Set parameter number to negative value so we know
*                   what it referrs to.  parn_to_amb is used to decide
*                   on open slots.
c                   amb_save(j,i) = amb_parn(j,i)
* MOD TAH 090110: Removed clearing parameter slot
C                   amb_parn(j,i) = 0 ! -amb_parn(j,i)
c                end do
c            end if
c         end do

****      OK, Now that we have output and reinitialized any ambiguities
*         whose data is finished, see if we have new values that we have to
*         add.  
          as = num_parm
          do i = 1, num_ambs
             if( .not.kbit(bf_ents(5,i),2) ) then
*                See if parameter is assigned.  If not then add
                 do j = 1, ne
                    if( amb_parn(j,i).eq.0 ) then
*                       No parameter assigned to assign now
                        as = as + 1
                        amb_parn(j,i) = as
                        amb_save(j,i) = as
                        parn_to_amb(as) = i
* MOD TAH 090221: Make sure we clear row (use pos_max_parm not num_parm
*                       which does not yet include ambiguities).
                        do k = 1, pos_max_parm
                           cov_parm(k,as) = 0.d0
                           cov_parm(as,k) = 0.d0
                        end do
                        cov_parm(as,as) = 1.d0
                        sol_vec(as) = 0.d0
C                       write(*,420) i, (bf_ents(k,i),k=1,5), as, 
C    .                               num_parm
C420                    format('BF_FLOAT Adding new parameter: BF ',i4,
C    .                     i2,' PRN ',i2.2,1x,2i7,1x,i3, ' Parameter ',
C    .                     i4,' of ',i4,' parameters')   
                    endif
                 end do
             endif
          end do
          if( as.gt.num_parm ) then
             write(*,430) as - num_parm
 430         format('BF_FLOAT: Ambiguities added ',I4,' parameters')
             num_parm = as
             if( num_parm.gt.max_parm ) then
                call report_error('MAX EXCEEDED',num_parm,
     .               'Parameter Estimates','est_pos',1,'est_pos')
             endif
          endif


C            if( .not.kbit(bf_ents(5,i),2) .and. 
C    .           ep-bf_ents(3,i).lt.an .and. ep-bf_ents(3,i).ge.0 ) then
C            if( .not.kbit(bf_ents(5,i),2) .and.
C    .          ((ep-bf_ents(3,i).lt.an .and. ep-bf_ents(3,i).ge.0 
C    .           .and. ani.gt.0) .or.
C    .          (ep-bf_ents(4,i).gt.-an .and. ep-bf_ents(4,i).le.0 
C    .           .and. ani.lt.0) )  ) then 
*                OK: This ambiquity is unknown and starts at this epoch
*                Add it to the parameter list.
C                as = 0
C                do j = non_amb_parm+1, num_parm
C                   if( parn_to_amb(j).eq.0 .and. as.eq.0 ) then
*                       OK, this slot is open.  Save the slot number so
*                       that we can use (save value before so that it
*                       is incremented.
C                       as = j - 1 
C                   end if
C                end do

*                see if we found a slot for this ambiquity
C                if( as.eq.0 ) then
*                    No, slot found so add to end of list
C                    as = num_parm
C                end if
*                OK: set up the parameter number and initialization
C                do j = 1, ne
C                   as = as + 1
C                   parn_to_amb(as) = i
C                   amb_parn(j,i) = as
C                   amb_save(j,i) = as
C                   if( as.gt.num_parm ) num_parm = as
C                   if ( as.eq.0 )
C    .              write(*,450) i, (bf_ents(k,i),k=1,5), as, num_parm
C450                format('BF_FLOAT Adding new parameter: BF ',i4,
C    .                     i2,' PRN ',i2.2,1x,2i7,1x,i3, ' Parameter ',
C    .                     i4,' of ',i4,' parameters')   

C                   do k = 1, num_parm
C                      cov_parm(k,as) = 0.d0
C                      cov_parm(as,k) = 0.d0
C                   end do
C                   cov_parm(as,as) = 1.d0
C                   sol_vec(as) = 0.d0
C                end do

*                If we have increases the parmeter number, save new value
c                 if( as.gt.num_parm ) num_parm = as
c            end if
c         end do
      end if

****  Thats all
      return
      end

CTITLE FLOAT_END

      subroutine float_end( data_type )

      implicit none

*     Routine to save and report biases at end of solution
      include 'track_com.h'

* PASSED VARIABLES
* data_type -- Type of data used in solution

      character*(*) data_type

* LOCAL VARIABLES
* i, k -- loop counter
* ne -- Number of entries

      integer*4 i, k, ne


*     See how many ambiquities we need to estimate
      ne = 0 
      if( index(data_type,'L1').gt.0 .or.
     .    index(data_type,'LC').gt.0     ) ne = ne + 1
      if( index(data_type,'L2').gt.0     ) ne = ne + 1

      do i = 1, num_ambs
         if( amb_parn(1,i).gt.0 ) then
*            OK, we have passed the end of the ambiguity.  See how
*            many values to output
             do k = 1, ne
                ambiq_float(k,i) = sol_vec(amb_parn(k,i))
                ambiq_flsig(k,i) = sqrt(cov_parm(amb_parn(k,i), 
     .                                           amb_parn(k,i))) 
             end do

             write(*,410) i, num_epochs, site_names(bf_ents(1,i)),
     .           (bf_ents(k,i),k=2,5),
     .            amb_parn(1,i),(sol_vec(amb_parn(k,i)),
     .            sqrt(cov_parm(amb_parn(k,i), amb_parn(k,i))),
     .                     k = 1, ne)
 410         format('BF_FLOAT ',I6,1x,i7,1x,a4,' PRN ',I3.2,1x,2I7,1x,
     .               o2.2,' NP ',I3,' Estimated Adjust. (m) ',
     .               2(F10.3,1x,F10.3))

         end if
      end do

****  Save the total number of parameters that we have encountered
      tot_parm = num_parm

****  Thats all 
      return
      end

CTITLE SETUP

      subroutine setup( option )

      implicit none

*     Estimation setup routine.  Here we compute the number of parameters
*     that we need to estimate.

      include 'track_com.h'

* PASSED VARIABLES
* option  -- Character string with options.  The choices are
*     EXATM  -- excludes the atmospheric delay from being included
*               in the parameter estimation
*     FLOAT  -- Clears the arrays for float bias parameter estimation. 

      character*(*) option

* LOCAL VARIABLES
* i,j  -- loop counters
* np   -- Parameter number counter

      integer*4 i, j, np

****  Clear the counters for the statistics of the double differences
      rms_dd_avg = 0.d0
      num_dd_avg = 0

****  OK, Find out how many kinematic sites we have.
      num_kine = 0
      do i = 1, num_site
         if( site_type(i).eq.1 ) then
             num_kine = num_kine + 1
*            Save mapping of kinematic site number to all sites numbers
             kine_to_site(num_kine) = i
         end if
      end do

****  Clear all the parameter number arrays.  These are pointers to the
*     parameter number of different types of parameters
* MOD TAH 080626: Fixed limits for kinematic parameter numbers (num_kine
*     instead of num_site)
      do i = 1, num_kine
         do j = 1, 3
            pos_parn(j,i) = 0
         end do
      end do
      do i = 1, num_site
         atm_parn(i) = 0
      end do

*     Start counting up the number of parameters and save their
*     pointers to the parameter numbers
      np = 0
      do i = 1, num_kine
         do j = 1, 3
            np = np + 1
            pos_parn(j,i) = np
         end do
      end do
         

*     Now see if we have atmospheric delays turned on
      if( index(option,'EXATM').eq.0 ) then
         do i = 1, num_site
            if( apr_atm(i).gt.0 ) then
                np = np + 1
                atm_parn(i) = np
            end if
         end do
      end if

*     See if we should clear the parameter arrays for a float solution
      if( index(option,'FLOAT').ne.0 ) then
          do i = 1, num_ambs
             do j = 1,2
                amb_parn(j,i) = 0
                amb_save(j,i) = 0
             end do
          end do
      end if

****  Save the final number of parameters
      num_parm = np
      non_amb_parm = np

      write(*,220) num_kine, num_parm
 220  format('There are ',i2,' kinematic sites in this analysis; ',
     .        i4,' parameters per epoch to be estimated')

*     Thats all
      return
      end

CTITLE GET_KINE_APR

      subroutine get_kine_apr(ep)

      implicit none

*     Routine to load the apriori position estimates for the kinematic
*     sites position.  If this is the first epoch, it uses the position
*     given by the user; for other epochs is uses the last kinematic 
*     position.

      include 'track_com.h'

* PASSED VARIABLES
* ep -- Epoch number at which to compute position

      integer*4 ep

* LOCAL VARIABLES
* i, j  -- Loop counters
* ks    -- Site number corresponding to kinematic site i 

      integer*4 i,j, k, ks

      logical kbit

****  OK, See if first epoch
      if( ep.eq.0 ) RETURN

      if( ep.eq.1 ) then
*         Just copy the user site apriori positions accross
          do i = 1, num_site
             do j = 1, 3

*               Set position to apriri for first epoch
*               Then if we know the kinematic positions
*               save these values, overwriting the apriori values
                curr_site_xyz(j,i) = site_apr(j,i)
                if( kine_known ) then
                    do k = 1, num_kine
                       ks = kine_to_site(k)
                       if( ks.eq.i ) then
                          curr_site_xyz(j,i) = kine_xyz(j,k,ep)
                       end if
                    end do
                endif
             end do
          end do
      else 

*         If the kinematic trajectory is not known then copy the 
*         coordinates from the last epoch.  If kine is known then
*         use the value for this epoch
*        (could add velocity later if needed)
          do i = 1, num_kine
             ks = kine_to_site(i)
c             if( static_obs(i).le.1 ) then
                do j = 1, 3
                   if( kine_known .or. kbit(kine_OK(1,i),ep) ) then 
                       curr_site_xyz(j,ks) = kine_xyz(j,i,ep)
                   else
                       curr_site_xyz(j,ks) = kine_xyz(j,i,ep-1)
                   end if                    
                end do
c             end if
          end do
      end if
C     if( ep.ge.debug_start .and. ep.le.debug_end ) then 
C        write(*,100) ep, kine_known, kbit(kine_OK(1,i),ep),
C    .             (i,(curr_site_xyz(j,kine_to_site(i)),j=1,3),
C    .             i = 1,num_kine)
C100     format(/,60('-'),/,
C    .          'CURR_SITE_XYZ ',i5,1x,2l1,1x,6(I3,3F15.3))
C     end if

****  Thats all
      return
      end

CTITLE UPDATE_POS

      subroutine update_pos(ep, data_type, options, converged, an, iter) 

      implicit none

*     This is the main estimation and final bias flag estimation routine.
*     It is first called with data_type P1 to generate a rough trajectory
*     for the kinematic sites.  

      include '../includes/const_param.h'
      include '../includes/xfile_def.h'
      include 'track_com.h'

* PASSED VARIABLES
* ep -- Epoch number at which to compute position
* iter -- Iteration count when searching biases.
* an   -- Ambiquity number we are searching on or if doing a
*         FLOAT solution, the sampling interval

      integer*4 ep, iter, an

* data_type -- Data types to be used (L1,L2,LC, P1, P2)
* options   -- Options for getting the position estimates.
*              RESBF  -- Resolve bias flags
*              ELEV   -- Compute elevation angles
*              KALF   -- Run as a kalman filter in that the covariance matrix
*                        is propagated forward in time.
*              SEARCH -- Option to Search the bias-parameter ambiquities
*                        to fix the biases to final values
*              FLOAT  -- Non-resolved biases estimated as float parameters
*              SMOOTH -- Indicates that covariance matrices should be written
*                        if solution is being run backwards.

      character*(*) data_type, options

* converged -- Set true when the kinematic position estimate is converged

      logical converged

* LOCAL VARIABLES
* no, ns, nd, nc -- Short name counters for number of one-ways, single differences,
*     double differences, and number of sd at each site.
* pn        -- Current PRN being processed
* i,j        -- Loop counters
* ref_sd     -- Pointer to single difference that is used to form double differences
* ks         -- kinematic site number
* np         -- Parameter number
* mc         -- Number of channels used in computing the station clock
*               epoch correction

      integer*4 no, ns, nd, nc, ni, pn, i, j, k, ref_sd, ks, np, 
     .          lv, mc, nt, is

* rcv_time  -- MJD of time given by the receiver clock
* trn_time  -- MJD of tranmission time as given by satellite clock (needs to be
*              corrected for satellite clock error).
* rcv_sec   -- Seconds from start of SP3 file for time
* trn_sec   -- Seconds from start of SP3 file for tme
* svs_ear(3) -- Satellite coordinates for current PRN and site
* range(2)   -- Range computed to satellite at L1 and L2 (m). (Difference due to
*               antenna offsets at the two frequencies)
* elev       -- Real*8 version of elevation angle (deg)
* az         -- Azimuth (real*8) (deg)

      real*8 rcv_time, trn_time, svs_ear(3), range(2), phase(2),
     .       elev, az, trn_sec

* sol_test(max_parm) -- Test solution vector to make sure that we have converged
* sol_upd(max_obs)   -- Update of the solution vector for the filter update
* scale(max_obs)     -- Scale factors for inversion
* acat(max_obs, max_obs) -- The ACAT matrix (partials by parameter 
*                       covariance by partials transposed)
* AC(max_parm,max_obs)   -- Computation of AC matrix
* ACI(max_site,max_obs) -- Computation of AC matrix for ionosheric delay
* acat_ion(max_obs,max_obs) -- ACAT matrix for ion estimation
* dummy              -- Dummy value passed to invert_vis since we do not  
*                       need solution in inversion
* total_adj          -- Total magntiude of adjust to position.  Needed to 
*                       test convergence
* used_cutoff        -- Used elevation angle cutoff.  (Accounts for small
*                       reduction needed if position changes)
      real*8 sol_test(max_parm), sol_upd(max_obs), scale(max_obs),
     .       acat(max_obs, max_obs), AC(max_parm,max_obs), dummy, 
     .       total_adj, used_cutoff, ACI(max_site, max_obs),
     .       acat_ion(max_obs, max_obs) 

* site_clk  -- Estimate of error in site clock in meters
* avclk     -- Average value of psuedorange minus range with satellite clock
*              correction
* omc(max_chan) -- Range obs-minus-computed
* prev_range(max_chan)  -- Range to satellite from previous iteration.
* dry_map   -- Dry mapping function, passed from theory into partials
* wet_map   -- Wet mapping function (use for partials)
* max_err   -- Max dd error when checking for outliers (mm)

      real*8 site_clk, avclk, omc(max_chan), prev_range(max_chan),
     .       dry_map, wet_map, max_err, var_sum 

      real*8 ssqr(max_site), swgh(max_site) ! Sum of res^2/wgh and 1/wgh for
                                     ! WRMS calcuations
     .,      obvar(max_obs)          ! Apriori variances for data (compare to
                                     ! post-fit residuals).

* ipivot(max_obs)    -- Pivot elements for inversion 
* trimlen  -- Length of string
* ndf, sdf, adf      -- Number of degrees of freedom adjustment, Total, Sites and
*             atmosphere
* msi -- Millisecond offsets on phase (uZ problem)
* max_ind  -- Index of maximum error while checking postfit double differences

      integer*4 ipivot(max_obs), trimlen, ndf, sdf, adf, msi, max_ind

      integer*4 date(5)
      real*8 sectag

      integer*4 PtoL   ! Generate the satellite list number from the PRN
C    .,         lv     ! Satellite list number 

* ow_used(max_obs)  -- Logical to indicate that a one-way measurement
*     has been used in a double difference
* old_ow  -- Set true when ever all one-ways have been used
* data_OK -- Logical function returns true if data OK (based on mask to
*            test bits)
* dd_OK   -- Check on dd size

      logical ow_used(max_obs), old_ow, data_OK, dd_OK

* line   -- Output line 
      character*512 line


****  First need to generate a set of residuals to the current position
*     estimates.
*     See if option to reduced elevation angles cutoff (to account for
*     possible position changes has been passed).  (0.003 deg is about
*     a 1 km station position error).

      if( index(options,'RedE').gt.0 ) then
          used_cutoff = elev_cutoff - 0.003
      else
          used_cutoff = elev_cutoff
      end if

      no = 0
      ni = 0

*     Set the number of degrees of freedom adjustment.  This is only
*     an approximate calcualtion
      sdf = 0
      adf = 0
      do i = 2, num_site
         if( mar_site(1,i).gt.(0.001)**2 ) sdf = 3
         if( mar_atm(i).gt.(0.001)**2 ) adf = 1
      end do
      ndf = sdf + adf 
      nt = 0

      do i = 1, num_site
         
*        Loop over the satellites at this epoch
         rcv_time = ((ep-1)*usr_interval+sec_offset(i,ep))/86400.d0+
     .                ref_start

*        Compute all the satellite clock corrections
         call comp_svs_clk( rcv_time )

         nc = 0

*        Iterate to get the receiver clock correction
         site_clk = 0
         if( ep.ge.debug_start .and. ep.le.debug_end ) then
            write(*,900) i, site_names(i), ep, num_chan_se(i,ep),
     .            ctop_cse(1:num_chan_se(i,ep),i,ep)
 900        format('DEBUG: Site ',I2,1x,a,' EP ',I6,' Chans ',I4,
     .            ' PRNS ',12I4.2)
         end if  

         do j = 1, num_chan_se(i,ep)
            prev_range(j) = 20.d6
            pn = ctop_cse(j,i,ep)
            lv = PtoL(pn)
            trn_time = rcv_time -
     .      	       prev_range(j)/vel_light/86400.d0
            trn_sec  = ((ep-1)*usr_interval+sec_offset(i,ep)) +
     .      	       ref_sec - prev_range(j)/vel_light
            call theory(trn_time, trn_sec, i, j, 
     .      	 curr_site_xyz(1,i), svs_ear, range, phase, elev,
     .      	 az, dry_map, wet_map, 'S', ep)
            prev_range(j) = range(1)
         end do

*        Run the iteration twice, with then the final calculation.
         do k = 1, 2
            avclk = 0.d0
            mc = 0
            do j = 1, num_chan_se(i,ep)
               if( data_OK(data_flag_cse(j,i,ep),data_mask) ) then
                  pn = ctop_cse(j,i,ep)
                  lv = PtoL(pn)
                 mc = mc + 1
                  trn_time = rcv_time -
     .                       prev_range(j)/vel_light/86400.d0
                  trn_sec  = ((ep-1)*usr_interval+sec_offset(i,ep)) +
     .                       ref_sec - prev_range(j)/vel_light
                                                                
                  call theory(trn_time, trn_sec, i, j, 
     .                 curr_site_xyz(1,i), svs_ear, range, phase,  
     .                 elev, az, dry_map, wet_map, 'S', ep)
                                                  
*                 Now get the average difference between psuedorange
*                 and range
                  omc(j) = P1o_all_cse(j,i,ep) - range(1) + 
     .                     svs_clk(lv)*vel_light
                  if( ep.ge.debug_start .and. ep.le.debug_end ) then
                      write(*,110) ep,i,j, pn, trn_sec, omc(j),
     .                      P1o_all_cse(j,i,ep), range(1), 
     .                      svs_clk(lv)*vel_light
 110                  format('DB OMC EP, site CH PRN ',i6,3i4,
     ,                       ' TRN_SEC ',F16.6,' OMC(PN) ',F12.3,
     .                       ' P1 Obs/Theory ',2F14.3,' SV Clk ',F14.3) 
                  end if
*
* MOD TAH 040321: Check that OMC is not in error by millisecs.
                  msi = nint((omc(j)/vel_light)/0.001d0)
                  avclk = avclk + omc(j)
                  prev_range(j) = range(1)
               else
                  omc(j) = P1o_all_cse(j,i,ep) - range(1) + 
     .                     svs_clk(lv)*vel_light
               end if
            end do

*           Check that the clock is OK:
* MOD TAH 040723: Compute avclk with mc (actual number of channels used) not
*           the nominal number of channels.  Calculation can only be done if
*           mc is greater than zero.
C           avclk = avclk/num_chan_se(i,ep)
            if( mc.gt.0 ) then
               avclk = avclk/mc
               call check_clk(ep,i,omc,avclk)
            end if

*           Get the average site clock
            if( mc.gt.0 ) then
               site_clk =  avclk
*              Correct the psuedoranges and the epoch offset for the
*              site clock and re-compute the rcv_time
               sec_offset(i,ep) = sec_offset(i,ep) - 
     .              site_clk/vel_light
               do j = 1, num_chan_se(i,ep)
                   P1o_all_cse(j,i,ep) = P1o_all_cse(j,i,ep) - site_clk
* MOD TAH 030727:  Added all observables for update so that the Widelanes
*                  will be consisent.
                   P2o_all_cse(j,i,ep) = P2o_all_cse(j,i,ep) - site_clk
                   L1o_all_cse(j,i,ep) = L1o_all_cse(j,i,ep) - 
     .                                           site_clk*fR1/vel_light
                   L2o_all_cse(j,i,ep) = L2o_all_cse(j,i,ep) - 
     .                                           site_clk*fR2/vel_light
              end do
     
              rcv_time = ((ep-1)*usr_interval+
     .                    sec_offset(i,ep))/86400.d0+ref_start
C               rcv_time = ((ep-1)*usr_interval+
C     .                    sec_offset(i,ep)-site_clk/vel_light)/86400.d0+
C     .                    ref_start
            end if
            if ( ep.ge.debug_start .and. ep.le.debug_end ) then
               write(*,120)  ep, site_names(i),k, sec_offset(i,ep),
     .             site_clk, (omc(j)-site_clk,
     .                        data_OK(data_flag_cse(j,i,ep),data_mask),
     .                        j=1,num_chan_se(i,ep))
 120           format('CLOCKS: Epoch ',I7,1x,a4,' Iter ',i2,
     .             ' Offset ',d15.8,' sec, ',e14.3,' m, OmC (m) ',
     .             (20(1x,F10.4,L1)))  

C              write(*,*) 'RANGES: ',(P1o_all_cse(j,i,ep)+
C    .                    svs_clk(pn)*vel_light,
C    .                    prev_range(j), j=1, num_chan_se(i,ep))
            end if

         end do
         
****     Now we are ready to start final calculation                 
         do j = 1, num_chan_se(i,ep)

*           Compute the theoretical range (use transmitt time as the
*           receive time minus the psuedorange)
            pn = ctop_cse(j,i,ep)
            lv = PtoL(pn)

c           trn_time = rcv_time - P1o_all_cse(j,i,ep)/vel_light/86400.d0
            trn_time = rcv_time - prev_range(j)/vel_light/86400.d0
            trn_sec  = ((ep-1)*usr_interval+sec_offset(i,ep)) +
     .                 ref_sec - prev_range(j)/vel_light
      
            call theory(trn_time, trn_sec, i, j,
     .                  curr_site_xyz(1,i), 
     .                  svs_ear,  range, phase, elev, az, dry_map, 
     .                  wet_map, 'F', ep) 

*           Save the elevation angle
            elev_cse(j,i,ep) = elev
            az_cse(j,i,ep) = az

*           See if below elevation cuttoff.
            if( elev.lt.used_cutoff ) then
                call sbit(data_flag_cse(j,i,ep),5,1)
            else
                call sbit(data_flag_cse(j,i,ep),5,0)
            end if

*           Now compute the one-way o-minus-c and partials array
            if( ep.ge.debug_start .and. ep.le.debug_end ) then
                if(data_OK(data_flag_cse(j,i,ep),data_mask) ) nt=nt+1
                write(*,130) ep, site_names(i),nt, 
     .                      ctop_cse(j,i,ep), data_flag_cse(j,i,ep),
     .                      data_mask
 130            format('OW: Epoch ',I7,1x,a4,' Cnt ',i3,' PRN',i3.2,
     .                 ' DataFlag ', o6,'o Mask ',o6,'o')

            end if

            if( data_OK(data_flag_cse(j,i,ep),data_mask) ) then
* MOD TAH 090112: Pass wet mapping function as partial.      
                call get_ow_omc(ep, i, j, pn, data_type, no, nc, range, 
     .                  phase, elev, dry_map, wet_map, svs_ear)
                if( index(options,'ION').ne.0 ) then
                    call get_ow_ion(ep, i,j, ni, elev)
                end if
            end if
         end do
         num_ow_by_site(i) = nc
         if( ep.ge.debug_start .and. ep.le.debug_end ) then
             write(*,150) ep, i, (curr_site_xyz(j,i),j=1,3)
 150         format('THEORY XYZ: ',i6,' Site ', i3,' XYZ ',3F20.3)
         endif

      end do
  
*     Save the number of one-way observations
      num_ow = no
      num_ion = ni

*     Set all one-ways as un-used
      do i = 1, num_ow
         ow_used(i) = .false.
      end do
      
****  Now form the double differences.  First difference
*     between the sites using the first site as reference
*     Loop over oneways above the first site.  Do this in groups of the
*     number types that we are using.
      ns = 0
      do i = num_ow_by_site(1)+1, num_ow, num_dtype

*        See if we can find the satellite for this measurment
*        to all sites preceding this one
         do j = 1, i, num_dtype
            if( ow_des(2,j).eq. ow_des(2,i) .and.
     .          ow_des(1,j).ne. ow_des(1,i)   ) then
*               Satellites match, increment the number of
*               single differences and the save the pointers
*               to the one-ways description array
                do k = 1, num_dtype
                   ns = ns + 1
                   sd_des(1,ns) = i+k-1
                   sd_des(2,ns) = j+k-1
                end do
            end if
         end do
      end do

*     Now search over all possible combinations to get the double
*     differences.  Again do the search in groups of the number of
*     data types being used.
*     Save number of single differences
      num_sd = ns
c     write(*,800) num_sd, (sd_des(1,j), sd_des(2,j),j=1,num_sd)
c800  format('Num SD ',i3,20(2i4,3x))


****  Now form the double differences.  Start differencing
*     from first entry
      ref_sd = 1
      nd = 0
      do j = 1, num_sd - num_dtype, num_dtype
          ref_sd = j
          do i = j+num_dtype, num_sd, num_dtype

*           Make sure that the reference station pair has not
*           changed.  Check both station pairs.
            old_ow = ow_used(sd_des(1,ref_sd)) .and.
     .               ow_used(sd_des(1,i))      .and.
     .               ow_used(sd_des(2,ref_sd)) .and.
     .               ow_used(sd_des(2,i))
            if( ow_des(1,sd_des(1,ref_sd)).eq.
     .          ow_des(1,sd_des(1,i))     .and.
     .          ow_des(1,sd_des(2,ref_sd)).eq. 
     .          ow_des(1,sd_des(2,i))     .and. .not.old_ow ) then

*               OK, sites are the same.  Now save the pointers
*               to the four one-ways that make up the double 
*               difference.
                do k = 1, num_dtype
                   nd = nd + 1
                   dd_des(1,nd) = sd_des(1,i+k-1)
                   dd_des(2,nd) = sd_des(2,i+k-1)
                   dd_des(3,nd) = sd_des(1,ref_sd+k-1)
                   dd_des(4,nd) = sd_des(2,ref_sd+k-1)
                   ow_used(sd_des(1,i+k-1)) = .true.
                   ow_used(sd_des(2,i+k-1)) = .true.
                   ow_used(sd_des(1,ref_sd+k-1)) = .true.
                   ow_used(sd_des(2,ref_sd+k-1)) = .true.
                end do

*               From the dd_des array we can now form the values
*               for the ionospheric delay observables
                do k = 1, 4
                   dd_dsi(k,nd/num_dtype) = (dd_des(k,nd)+num_dtype-1)/
     .                                        num_dtype
                end do

            end if
         end do
      end do

****  Save the number of double differences
      num_dd = nd

*     We now have the double differences.  Now form the double differences
*     and the partials array.
      call form_dd(ep)

*     Now do a rough check on omc to make sure there are no
*     gross errors.
      dd_OK = .true.
      do i = 1, num_dd
         if( abs(sol_obs(i)).gt.100.d0 .and. 
     .       index(options,'EDIT').gt.0       ) then

* MOD TAH 020626: Check to see if a mulitple of milliseconds
*            off.
             msi = nint((sol_obs(i)/vel_light)*1000.d0)
             if( abs(sol_obs(i)-msi*vel_light/1000.d0).gt.100.d0)then
             
                dd_OK = .false.
*               Mark the contributing one-way as bad
                ns =  ow_des(1,dd_des(1,i))
                lv =  ow_des(2,dd_des(1,i))
                call mjd_to_ymdhms((ep-1)*usr_interval/86400.d0
     .                      +ref_start,date, sectag)

                write(*,320) ep, date, sectag, site_names(ns), 
     .                       lv, sol_obs(i)
 320            format('+ Edit at Epoch ',i6,1x,I4,4(1x,i2.2),1x,f6.2,
     .                 ' Site ',a4,' PRN ', I3,' Residual ',
     .                  F12.2,' m')
                call mark_slip(ep, ns, lv, 2, sol_obs(i))
             else
*               Seem to have a millisec slip
                ns =  ow_des(1,dd_des(1,i))
                lv =  ow_des(2,dd_des(1,i))
                write(*,330) ep, site_names(ns), lv, sol_obs(i),
     .                       msi
 330            format('+ Millisecond slip at Epoch ',i6,' Site ',a4,
     .                 ' PRN ',I3,' Residual ',F12.2,' m, ',
     .                 i4,' ms')
                sol_obs(i) = sol_obs(i)-msi*vel_light/1000.d0
*               Now correct the original data:
                call fix_ms(ep,i,ns,lv, msi)
             end if
         end if
      end do
      if( .not.dd_OK ) then
         write(*,*) '**ERROR** Bad double differernces at ',ep
      end if
 

* DEBUG: See the ep range to something resonable to get output
      if( (ep.ge.debug_start .and. ep.le.debug_end) ) then
          do i = 1, num_dd
             write(*,340) ep, i, (dd_des(j,i),j=1,4), sol_obs(i)
 340         format('DD EP ',i6,' # DD_DES ',5i4,1x,' OminusC ',F15.3,
     .              ' m')
          end do
      end if

      if( .not.dd_OK ) num_dd = 0

*     Now do the filtering: First form the ACAT matrix.  Save the AC part of the 
*     matrix
      do i = 1, num_parm
         do j = 1, num_dd
            AC(i,j) = 0.d0
            do k = 1, num_parm
               AC(i,j) = AC(i,j) + apart(k,j)*cov_parm(i,k)
            end do
         end do
      end do

*     Now form ACAT matrix
      do i = 1, num_dd
         do j = 1, num_dd
             acat(i,j) = 0.d0
             do k = 1, num_parm
                acat(i,j) = acat(i,j) + AC(k,i)*apart(k,j)
             end do
         end do
      end do

      if( ep.le.0 ) then
         do i = 1,num_dd
            write(*,350) 'CobsM',i,(cov_obs(i,j),j=1,num_dd)
 350        format(a,1x,i3,60(E13.4,1x))
         end do
         do i = 1,num_parm
            write(*,350) 'PARMM',i,(cov_parm(i,j),j=1,num_parm)
        end do
         do i = 1,num_parm
            write(*,350) 'APART',i,(apart(i,j),j=1,num_dd)
         end do
         do i = 1,num_parm
            write(*,350) 'AC   ',i,(AC(i,j),j=1,num_dd)
         end do
         do i = 1,num_dd
            write(*,350) 'ACAT ',i,(ACAT(i,j),j=1,num_dd)
         end do


      end if


*     Now add to data covariance and invert
      do i = 1, num_dd
         obvar(i) = cov_obs(i,i)
         do j = 1, num_dd
            cov_obs(i,j) = cov_obs(i,j) + acat(i,j)
         end do
      end do 
      if( ep.le.0 ) then
         do i = 1,num_dd
            write(*,350) 'CobsP',i,(cov_obs(i,j),j=1,num_dd)
         end do
      end if

      call invert_vis(cov_obs, dummy, scale, ipivot, num_dd, max_obs, 0)

*     Now finish forming the Kalman Gain matrix
      do i = 1, num_parm
         do j = 1, num_dd
            kgain(i,j) = 0.d0
            do k = 1, num_dd
               kgain(i,j) = kgain(i,j) + AC(i,k)*cov_obs(k,j)
            end do
         end do
      end do
      if( ep.le.0 ) then
         do i = 1,num_parm
            write(*,350) 'Kgain',i,(kgain(i,j),j=1,num_dd)
         end do
      end if

****  Now get the difference between the o-minus-c and the current solution
*     estimate (note: the position corrections are always zero.  This is mainly
*     for the atmospheric delay corrections)
      do j = 1, num_dd
         sol_upd(j) = sol_obs(j)
         do i = 1, num_parm
            sol_upd(j) = sol_upd(j) - apart(i,j)*sol_vec(i)
         end do
      end do

*     Update the parameter estimates.  Save this in a temporary vector
*     until we are sure we have converged
      total_adj = 0
      do i = 1, num_parm
          if( cov_parm(i,i).lt.0 ) then
              print *,'BAD COV_PARM Ep ',ep,' Param ',i, 
     .              ' Val ',cov_parm(i,i)
              stop 'Track: Negative covariance element'
          endif
      enddo

      do i = 1, num_parm
         sol_test(i) = sol_vec(i)
         do j = 1, num_dd
            sol_test(i) = sol_test(i) + kgain(i,j)*sol_upd(j)
         end do
                         
*        Only add contribution if station coordinate
         do j = 1, num_kine
            do k = 1,3
               if( i.eq.pos_parn(k,j) ) then
                   total_adj = total_adj + sol_test(i)**2
               end if
            end do
         end do
      end do
c     write(*,999) ep, total_adj, (sol_test(i),i=1,num_parm)
c999  format('TEST:  ',i6,20F10.3)
 
      if( ep.ge.debug_start .and. ep.le.debug_end ) then
         write(*,996) ep, (sol_test(i), i=1,num_parm)
 996     format('SOLVEC: ',i6,(20F10.3))
      endif
 
****  Compute the postfit residuals
      do j = 1, num_dd
         sol_upd(j) = sol_obs(j)
         do i = 1, num_parm
            sol_upd(j) = sol_upd(j) - apart(i,j)*sol_test(i)
         end do
      end do

****  If we estimating the ionospheric delay, start the filtering
*     process.

      if( index(options,'ION').gt.0 ) then 

*         Form the double differences of the one way ionospheric
*         delay values
          num_ddion = num_dd/num_dtype
          call form_iondd

*         Start forming the ACAT matrix first with AC which we
*         use later as well.
          do i = 1, num_site
             do j = 1, num_ddion
                ACI(i,j) = 0.d0
                do k = 1, num_site
                   ACI(i,j) = ACI(i,j) + apion(k,j)*ion_parm(i,k)
                end do
             end do
          end do
          do i = 1, num_ddion
             do j = 1, num_ddion
                acat_ion(i,j) = 0.d0
                do k = 1, num_site
                   acat_ion(i,j) = acat_ion(i,j) + ACI(k,i)*apion(k,j)
                end do
             end do
          end do

*         Now add to the data covariance
          do i = 1, num_ddion
             do j = 1, num_ddion
                ion_cov(i,j) = ion_cov(i,j) + acat_ion(i,j)
             end do
          end do 
          call invert_vis(ion_cov, dummy, scale, ipivot, num_ddion, 
     .                    max_obs, 0 ) 

*         Finish forming the kgain matrix
          do i = 1, num_site
             do j = 1, num_ddion
                kg_ion(i,j) = 0.d0
                do k = 1, num_ddion
                   kg_ion(i,j) = kg_ion(i,j) + ACI(i,k)*ion_cov(k,j)
                end do
             end do
          end do

*         OK, Now form the initial estimate and residuals
          do i = 1, num_site
             ion_vec(i) = 0.d0
             do j = 1, num_ddion
                ion_vec(i) = ion_vec(i) + kg_ion(i,j)*ion_obs(j)
             end do
          end do

*         Now generate the residuals and the mapping from changes
*         in the one-ways to residual one-ways
          do i = 1, num_ddion
             do j = 1, num_site
                ion_obs(i) = ion_obs(i) - apion(j,i)*ion_vec(j)
             end do
          end do
          
          if( ep.ge.debug_start .and. ep.le.debug_end ) then 
             write(line,380) ep, (ion_vec(j),j=1,num_site)
 380         format('ION_EST: ',i6,' Est ',6F12.3)
             write(line(trimlen(line)+2:),385) 
     .                 (ion_obs(j),j=1,num_ddion)
 385         format('Residuals ',50F8.2)
             if( num_ddion.gt.0 )
     .       write(*,'(a)') line(1:trimlen(line))
          end if

*         Now for the mapping from changes in the one-ways to 
*         changes in the residuals
          do i = 1, num_ddion
             do j = 1, num_ddion
                ion_to_res(i,j) = 0.d0
                do k = 1, num_site
                   ion_to_res(i,j) = ion_to_res(i,j) + 
     .                               apion(k,i)*kg_ion(k,j)
                end do
             end do
          end do
      end if

* DEBUG:
      if( ep.ge.debug_start .and. ep.le.debug_end ) then
         write(*,430) ep, num_dd, (sol_upd(j),j=1,num_dd)
 430     format('Post-fit Residuals (m): ',i6,1x,i4,100(F10.4))
      end if

****  Update the station coordinates
      do i = 1, num_kine
         ks = kine_to_site(i)
         if( 1.eq.2 )
     .   write(*,800) ep, data_type(1:2), ks, 
     .       (curr_site_xyz(j,ks), sol_test(pos_parn(j,i)),
     .        curr_site_xyz(j,ks)+sol_test(pos_parn(j,i)),j=1,3)
 800     format('POSUPD ',I7,1x,a2,1x,i2,' Coords ',
     .           3(F14.4,1x,F9.4,F14.4))
 
         do j = 1, 3
            np = pos_parn(j,i)
            if( np.gt.0 ) then
                curr_site_xyz(j,ks) =  curr_site_xyz(j,ks) + 
     .                                 sol_test(np)
            end if
         end do
         if( ep.ge.debug_start .and. ep.le.debug_end ) then
            write(*,997) ep, ks, (curr_site_xyz(j,ks),j=1,3),
     .           (sol_test(pos_parn(j,i)),j=1,3)
 997        format('UPDATED: ',i6,i3,3F13.3,' DP ',3F8.3)
         end if
      end do

****  OK: See if we have converged.  Check the size of the adjusments to 
*     to the positions (converge to 10 m).  Return for next solution 
*     unless we are searching the ambiquiuty space and/or running a 
*     float solution.
      if( total_adj.gt. 10.d0 .and. 
     .    index(options,'SEARCH').eq.0  .and.
     .    index(options,'FLOAT').eq.0) then
          converged = .false.
          RETURN
      end if


****  OK, solution looks good.  Now save the final answer and update the
*     covariance matrix
      do i = 1, num_parm
         sol_vec(i) = sol_test(i)
         do j = 1, num_parm
            do k = 1, num_dd
               cov_parm(i,j) = cov_parm(i,j) - kgain(i,k)*AC(j,k)
            end do
         end do
      end do

      if( ep.ge.debug_start .and. ep.le.debug_end ) then
         write(*,1501) ep, (cov_parm(i,i), i=1,num_parm)
1501     format('COVDIAG: ',i6,(50e13.5))
      endif

****  Accumulate statistics by site
      do is = 1, num_site
         ssqr(is) = 0    ! Sum (res/sig)^2
         swgh(is) = 0    ! Sum (1/sig)^2
         num_dd_site(is) = 0
      end do

*     Compute the final residuals and RMS
      rms_dd = 0.d0
      do j = 1, num_dd
         do i = 1, num_parm
            sol_obs(j) = sol_obs(j) - apart(i,j)*sol_vec(i)
         end do 
         rms_dd = rms_dd + sol_obs(j)**2
*        Site dependent calculation
         is = ow_des(1,dd_des(1,j))

         ssqr(is) = ssqr(is) + sol_obs(j)**2/obvar(j)
         swgh(is) = swgh(is) + 1/obvar(j)
         num_dd_site(is) = num_dd_site(is) + 1
  
      end do
      var_sum = rms_dd

***   Finish up calculation
      do is = 1, num_site
         if( num_dd_site(is).gt. 0 ) then
            wrms_dd_site(is) = sqrt(ssqr(is)/swgh(is))*1000.d0
         else
            wrms_dd_site(is) = 0
         end if
      end do


*     Accumulate the average RMS DD residuals
      if( num_dd-ndf.gt.0 ) then
          rms_dd_avg = rms_dd_avg + rms_dd
          num_dd_avg = num_dd_avg + (num_dd-ndf)
      end if

      if( num_dd-ndf.gt.0 ) then
          rms_dd = sqrt(rms_dd/(num_dd-ndf))*1000.d0
      else 
          rms_dd = sqrt(rms_dd)*1000.d0
      end if

****  See if RMS is bad enough that we should edit the data 
      if( ep.ge.debug_start .and. ep.le.debug_end ) then
         write(*,1503) ep, rms_dd, rms_edtol, sqrt(ow_var(1))*1.d3,
     .                 num_dd
1503     format('DD_RMS ',I7,2F12.2,1x,' OW Sig ',F10.3,' Numdd ',i4)
      endif

*     Save the RMS in case we re-set because too large
      var_sum = rms_dd

      if( rms_dd.gt.100.d0 .and. 
     .    rms_dd.gt.rms_edtol*sqrt(ow_var(1))*1.d3 ) then
         rms_dd = 100.d0
      end if

      if( rms_dd.gt.rms_edtol*sqrt(ow_var(1))*1.d3 ) then 
          max_err = 0.d0
          do i = 1, num_dd
             if( abs(sol_obs(i)).ge.max_err ) then
                 max_err = abs(sol_obs(i))
                 max_ind = i
             end if
          end do

*         Now decompose which term has the largest error
C         print *,'Checking ',ep,' RMS_DD/Max_error ',rms_dd, max_err,
C    .             max_ind, num_dd, ndf, num_parm, 
C    .             rms_edtol*sqrt(ow_var(1))
          if( max_err*1000.d0.gt.rms_edtol*rms_dd ) then
             ns =  ow_des(1,dd_des(1,max_ind))
             lv =  ow_des(2,dd_des(1,max_ind))
             call mjd_to_ymdhms((ep-1)*usr_interval/86400.d0
     .                      +ref_start,date, sectag)

             write(*,720) ep, date, sectag, site_names(ns), 
     .                    lv, sol_obs(max_ind)*1.d3, 
     .                    site_names(ow_des(1,dd_des(2,max_ind))),
     .                    ow_des(2,dd_des(3,max_ind)),
     .                    var_sum
 720         format('+ Post-fit Edit at Epoch ',i6,1x,
     .               I4,4(1x,I2.2),1x,f6.2, ' Site ',a4,' PRN ',
     .              I3.2,' Residual ',F12.1,' mm; DD ',a4,
     .             ' PRN ',I3.2,' RMS ',F8.1,' mm')
             call mark_slip(ep, ns, lv, 2, sol_obs(max_ind))
          end if
      end if



****  Save the kinematic position only if the site appears to non-static
*     (For static sites, the position changes are retained in the 
*     sol_vec solution vector, for kinematic sites the position is
*     updated.  In the first run using P1 the static array is set 
*     false).      
      converged = .true.                 
      do i = 1, num_kine
          ks = kine_to_site(i)
          do j = 1, 3
             if( .not.static(i) ) then
                 kine_xyz(j,i,ep) = curr_site_xyz(j,ks)
             end if
             kine_out(j,i) = curr_site_xyz(j,ks)
          end do
          if( ep.ge.debug_start .and. ep.le.debug_end ) then
             write(*,995) ep, i, ks, (kine_out(j,i),j=1,3)
 995         format('KINEOUT: ',i7, 2i3, 3F14.4)
          end if

*         Make position as OK if sigma < 1
          np = pos_parn(1,i)
          if( cov_parm(np,np)+cov_parm(np+1,np+1)+
     .        cov_parm(np+2,np+2).lt.1 ) then
              call sbit(kine_OK(1,i),ep,1)
          else
              call sbit(kine_OK(1,i),ep,0)
          end if
      end do

****  See if we should write the covariance matrix and solution
      if( index(options,'SMOOTH').gt.0  ) then
          if( an.lt.0 ) then
              call wr_cov_parm(ep,'W',data_type)
          end if
      end if

      return
      end

CTITLE PARTIALS

      subroutine partials(ep, ns, na, nf, site_pos, svs_ear, range,  
     .                    elev, wet_map, parts)

      implicit none
             
*     Routine to form the partial derivatives for the current observation.

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* ns -- Site number
* na -- Ambiquity number.  If non-biased observable (P1/P2/PC) then
*       na is set to 0.
* nf -- Frequency number (nf=1 is L1 or LC; nf=2 is L2)
* ep -- Epoch number

      integer*4 ns, na, nf, ep

* site_pos(3) -- Position of the site
* svs_ear(3)  -- Position of satellite
* range       -- Range to satellite (m)
* elev        -- Elevation angle to satellite (deg)
* wet_map     -- Dry mapping function 
* parts(max_parm) -- Partials derivatives for oneways

      real*8 site_pos(3), svs_ear(3), range, elev, wet_map,
     .       parts(max_parm)

      real*8 rot_mat(3,3), geod(3)  ! Rotation matrix and Geod coordinates

* LOCAL VARIABLES
* j   -- Loop counter
* ks  -- Kinematic site number
* np  -- Parameter number
* ks  -- Kinematic site numbber

      integer*4 j,ks, np

****  Loop over the parameters to be estimated and compute
*     partial if needed
*     Clear all the partials first
      do j = 1, num_parm
         parts(j) = 0.d0
      end do

*     See if this a kinematic site and therefore will have
*     a site position partial      
      ks = 0
      do j = 1, num_kine
         if( kine_to_site(j).eq.ns ) ks = j
      end do 

*     If ks is set, then this is a kinematic site with postion partials
      if( ks.gt.0 ) then 
          do j = 1, 3
              np = pos_parn(j,ks)
              if( np.gt.0 ) then
                  parts(np) = (site_pos(j)-svs_ear(j))/range
              end if
          end do
      end if

*     Now see if atmospheric delay estimated
      np = atm_parn(ns)
      if( np.gt.0 ) then
*         Use just an approximate formula for the moment
c          parts(np) = 1.d0/(sin(elev*pi/180.d0)+1.232d-3)
          if( .not. atm_scale(ns) ) then
             parts(np) = wet_map
          else  ! Use Delta-Height scale factor
!            Compute geodetic coordinates of site
             call XYZ_to_GEOD( rot_mat, site_pos, geod) 
             parts(np) = (geod(3)-site_geod(3,1))*1.d-3*wet_map
          end if
!         if( nint(ep/1000.)*1000-ep.eq.0 ) then
!             print *,'ATM_PART ',ep, ns, np, atm_scale(ns),
!     .            parts(np), (geod(3)-site_geod(3,1))*1.d-3, elev
!         endif
      end if

****  See if we have any ambiquity partials to add
      if( na.gt.0 ) then
         np = amb_parn(nf,na)
         if( np.gt.0 ) then
             parts(np) = 1.d0
         end if
      end if

****  Thats all
      return
      end

CTITLE FORM_DD

      subroutine form_dd(ep)

      implicit none

*     Routine to form the double differences and the double difference
*     partials based on the dd_des pointers and the ow o-minus-c and
*     partials.

      include 'track_com.h'

* PASSED VARIABLES
      integer*4 ep    ! Epoch number (used for debug)

* LOCAL VARIABLES
* i,j,k,l  -- Loop counters
* dd_sgn(4) -- Signs with which the one-ways are applied (1 -1 -1 1)
* ref_owdt  -- First site description of data type.

      integer*4 i,j,k,l, dd_sgn(4), ref_owdt

      data dd_sgn / 1, -1, -1, 1 /            


****  Loop over the number of double differences we have
      do i = 1, num_dd
         sol_obs(i) = (ow_vec(dd_des(1,i)) - ow_vec(dd_des(2,i)))-
     .                (ow_vec(dd_des(3,i)) - ow_vec(dd_des(4,i)))
         if( ep.ge.debug_start .and. ep.le.debug_end ) then
             write(*,100) ep, i, num_dd, dd_des(:,i),
     .           ow_vec(dd_des(:,i)), sol_obs(i)  
 100         format('DDALL Ep ',i6,' DD ',2I4,' DES ',4i4,
     .              ' OW_VEC ',4F16.4,' DDiff ',F10.4) 
         endif
*        Set and check the data type for this double difference
         ref_owdt = ow_dt(dd_des(1,i))
         do j = 2,4
            if( ref_owdt.ne.ow_dt(dd_des(j,i))) then
               write(*,120) i, j, dd_des(j,i), ref_owdt,
     .         ow_dt(dd_des(j,i))
 120           format('** DATA TYPE ERROR ** for DD ',i4,
     .                ' Comp index ',i2,' DD_DES ',i4,
     .                ' Data Types ',2o4)
            end if
         end do
         dd_dt(i) = ref_owdt 
c        write(*,*) i, sol_obs(i),ow_vec(dd_des(1,i)),
c     .      ow_vec(dd_des(2,i)), ow_vec(dd_des(3,i)),
c     .      ow_vec(dd_des(4,i)) 

*        Now form the partials
         do j = 1, num_parm
            apart(j,i)=(ow_part(j,dd_des(1,i))-ow_part(j,dd_des(2,i)))-
     .                 (ow_part(j,dd_des(3,i))-ow_part(j,dd_des(4,i)))
c            write(*,*) 'Parts ',i,j,ow_part(j,dd_des(1,i)),
c     .                  ow_part(j,dd_des(2,i)),
c     .                  ow_part(j,dd_des(3,i)),ow_part(j,dd_des(4,i)) 
         end do

*        Now form the double difference covariance matrix.  Here we add 
*        a contribution if the one-ways used in the forming the double
*        difference are the same.  The dd_sng accounts for the different
*        signs that the oneways may be used in.
         do j = 1, num_dd
            cov_obs(i,j) = 0.d0
            do k = 1,4
               do l = 1,4
                  if( dd_des(k,i).eq.dd_des(l,j) ) then
                      cov_obs(i,j) = cov_obs(i,j) + 
     .                   dd_sgn(k)*ow_var(dd_des(k,i))*
     .                   dd_sgn(l)
                  end if
               end do
            end do
         end do
      end do

****  Thats all (unless we form covariance matrix at the same time?)

      return
      end

CTITLE GET_OW_OMC

      subroutine get_ow_omc(ep, i, j, pn, data_type, no, nc, range, 
     .             phase, elev, dry_map, wet_map, svs_ear)

      implicit none

*     Routine to get the one-way o-minus-c values at epoch ep, for
*     site i, channel j.

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* ep   -- Epoch number
* i, j -- Site and channel numbers
* pn   -- Satellite PRN number
* no   -- Running sum of number of oneways
* nc   -- Running sum of number of oneways at the current site
* ne   -- Number of estimated data types 

      integer*4 ep, i,j, pn, no, nc, ne

* range(2) -- Theoretical range to the satellite at L1 and L2
* phase(2) -- Theoretical phase (L1/L2 cycles)
* elev  -- Elevation angle (deg)
* dry_map -- Dry mapping function.  Used for atm partial
* wet_map -- Wet mapping function 
* svs_ear(3) -- Position of satellite being used (m)

* elev_wght -- Weight to given to elevation dependence

      real*8 range(2), phase(2), elev, dry_map, wet_map, svs_ear(3), 
     .       elev_wght

* data_type -- Character string with type of data to use

      character*(*) data_type

* LOCAL VARIABLES
* na   -- Number of current ambiquity
* k    -- Loop counter

      integer*4 na, k

      integer*4 PtoL   ! Generate the satellite list number from the PRN
     .,         lv     ! Satellite list number 

      elev_wght = 0.0d0

***   OK, loop over each data type
*     OK, based on data types, compute the o-minus-c values
      k = 0
      ne = 0 
      na = amb_point_cse(j,i,ep)
      if( na.le.0 .and. index(data_type,'L').gt.0 ) then
          write(*,120) ctop_cse(j,i,ep),j,i,ep
 120      format('** ERROR ** No ambiquity number for PRN ',i3.2,
     .           ' Chan ',i2,' Site ',i2,' Epoch ',i6)
          RETURN
      end if  
      if( index(data_type,'L1').gt.0 ) then
          na = amb_point_cse(j,i,ep)
          if( na.le.0 ) then
             write(*,*) '**WARNING** Not ambiguity number at epoch ',
     .                  ep,' Site ',i,' PRN ', ctop_cse(j,i,ep)
             RETURN
          end if
          no = no + 1
          nc = nc + 1
          ne = ne + 1
          ow_vec(no) = (L1o_all_cse(j,i,ep) +
     .                  ambiq_all(1,na))*vel_light/fR1 - 
     .                  phase(1)*vel_light/fR1
C     .                  ambiq_all(1,na))*vel_light/fR1 - range(1)
          ow_des(1,no) = i
          ow_des(2,no) = pn
*         Set data type bit to indicate L1 phase (similar steps are done
*         other observables).
          ow_dt(no) = 1

*         Need to change later
          lv = PtoL(pn)

          ow_var(no) = data_var(1,lv)*
     .                 (1.d0+(data_var(5,lv)/sin(elev*pi/180.d0)**2)) 
*         If the bias is not fixed, give the data to range sigma
C         if( .not.kbit(bf_ents(5,na),2) ) ow_var(no) = data_var(3)
          call partials(ep, i, na, ne, curr_site_xyz(1,i), svs_ear, 
     .                  range(1), elev, wet_map, ow_part(1,no))
      end if

      if( index(data_type,'L2').gt.0 ) then
          na = amb_point_cse(j,i,ep)
          if( na.le.0 ) then
             write(*,*) '**WARNING** Not ambiguity number at epoch ',
     .                  ep,' Site ',i,' PRN ', ctop_cse(j,i,ep)
             RETURN
          end if
          no = no + 1
          nc = nc + 1
          ne = ne + 1
          ow_vec(no) = (L2o_all_cse(j,i,ep) +
     .                  ambiq_all(2,na))*vel_light/fR2 - 
     .                  phase(2)*vel_light/fR2
C     .                  ambiq_all(2,na))*vel_light/fR2 - range(2)
          ow_des(1,no) = i
          ow_des(2,no) = pn
          ow_dt(no) = 2

*         Need to change later
          lv = PtoL(pn)
          ow_var(no) = data_var(2,lv)*
     .                 (1.d0+(data_var(5,lv)/sin(elev*pi/180.d0)**2)) 
*         If the bias is not fixed, give the data to range sigma
C          if( .not.kbit(bf_ents(5,na),2) ) ow_var(no) = data_var(4)
          call partials(ep, i, na, ne, curr_site_xyz(1,i), svs_ear, 
     .                  range(1), elev, wet_map, ow_part(1,no))
      end if

      if( index(data_type,'LC').gt.0 ) then
          na = amb_point_cse(j,i,ep)
          if( na.le.0 ) then
             write(*,*) '**WARNING** Not ambiguity number at epoch ',
     .                  ep,' Site ',i,' PRN ', ctop_cse(j,i,ep)
             RETURN
          end if
          no = no + 1
          nc = nc + 1
          ne = ne + 1
          ow_vec(no) = (lcf1*(L1o_all_cse(j,i,ep)+ambiq_all(1,na)) +
     .                  lcf2*(L2o_all_cse(j,i,ep)+ambiq_all(2,na)) )*
     .                  vel_light/fR1 - 
     .                 (pcf1*phase(1)*vel_light/fR1 + 
     .                  pcf2*phase(2)*vel_light/fR2)
          if( ep.ge.debug_start .and. ep.le.debug_end ) then
             write(*,220) ep, j,i, na, L1o_all_cse(j,i,ep), 
     .               ambiq_all(1,na),  L2o_all_cse(j,i,ep),
     .               ambiq_all(2,na), phase(1), phase(2)
 220         format('OWOMC ',I6,' CSA ',3i4,' PhaseObs ',4F14.3,
     .              ' Theory ',2F14.4)
          end if
C     .                  vel_light/fR1 - (pcf1*range(1)+pcf2*range(2))
*         Compute ion free delay accounting for the two ranges at 
*         L1 and L2. Note: due to definition of lcf2 we divide the
*         L2 phase by the L1 frequency.
C         ow_vec(no) =  lcf1*(L1o_all_cse(j,i,ep)+ambiq_all(1,na))*
C    .                  vel_light/fR1 - range(1) +
C    .                  lcf2*(L2o_all_cse(j,i,ep)+ambiq_all(2,na))*
C    .                  vel_light/fR1 - range(2)

          ow_des(1,no) = i
          ow_des(2,no) = pn
          ow_dt(no) = 3
          lv = PtoL(pn)

*         Need to change later
          ow_var(no) = (data_var(1,lv)*pcf1**2 + 
     .                  data_var(2,lv)*pcf2**2)*
     .                 (1.d0+(data_var(5,lv)/sin(elev*pi/180.d0)**2)) 
*         If the bias is not fixed, give the data to range sigma
C          if( .not.kbit(bf_ents(5,na),2) ) ow_var(no) = data_var(1)
          call partials(ep, i, na, ne, curr_site_xyz(1,i), svs_ear, 
     .                  range(1), elev, wet_map, ow_part(1,no))
      end if

      if( index(data_type,'P1').gt.0 ) then
          no = no + 1
          nc = nc + 1
          ow_vec(no) = P1o_all_cse(j,i,ep) - range(1)
          ow_des(1,no) = i
          ow_des(2,no) = pn
          ow_dt(no) = 4
          lv = PtoL(pn)

*         Need to change later
          ow_var(no) = data_var(3,lv)*
     .                 (1.d0+elev_wght/sin(elev*pi/180.d0)) 
          call partials(ep, i, 0, 0, curr_site_xyz(1,i), svs_ear,  
     .                  range(1), elev, wet_map, ow_part(1,no))
      end if 
      if( index(data_type,'P2').gt.0 ) then
          no = no + 1
          nc = nc + 1
          ow_vec(no) = P2o_all_cse(j,i,ep) - range(2)
          ow_des(1,no) = i
          ow_des(2,no) = pn
          ow_dt(no) = 8
*         Need to change later
          lv = PtoL(pn)
         ow_var(no) = data_var(3,lv)*
     .                 (1.d0+elev_wght/sin(elev*pi/180.d0)) 
          call partials(ep, i, 0, 0, curr_site_xyz(1,i), svs_ear,  
     .                  range(1), elev, wet_map, ow_part(1,no))
      end if

      if( index(data_type,'PC').gt.0 ) then
          no = no + 1
          nc = nc + 1
          ow_vec(no) = pcf1*(P1o_all_cse(j,i,ep)-range(1)) +
     .                 pcf2*(P2o_all_cse(j,i,ep)-range(2))
          ow_des(1,no) = i
          ow_des(2,no) = pn
          ow_dt(no) = 12

*         Need to change later
          lv = PtoL(pn)
          ow_var(no) = (data_var(3,lv)*pcf1**2 + 
     .                  data_var(4,lv)*pcf2**2)*
     .                 (1.d0+elev_wght/sin(elev*pi/180.d0)) 
          call partials(ep, i, 0, 0, curr_site_xyz(1,i), svs_ear,  
     .                  range(1), elev, wet_map, ow_part(1,no))
      end if

****  Thats all for the moment, add more later

      return
      end

CTITLE FORM_IONDD

      subroutine form_iondd

      implicit none

*     Routine to form the ionospheric delay double differences and 
*     partials based on the dd_des pointers and the ow o-minus-c and
*     partials.

      include 'track_com.h'

* PASSED VARIABLES
* NONE

* LOCAL VARIABLES
* i,j,k,l  -- Loop counters
* dd_sgn(4) -- Signs with which the one-ways are applied (1 -1 -1 1)

      integer*4 i,j,k,l,m,n, dd_sgn(4)

      data dd_sgn / 1, -1, -1, 1 /            


****  Loop over the double difference pointers from the main
*     data analyis.  Since we use only one observable here skip
*     any entries due to additional data
      i = 0
      do k = 1, num_dd, num_dtype
         i = i + 1
         ion_obs(i) = (ow_ion(dd_dsi(1,i)) - ow_ion(dd_dsi(2,i)))-
     .                (ow_ion(dd_dsi(3,i)) - ow_ion(dd_dsi(4,i))) 
C        write(*,*) 'ION ', k,i, ion_obs(i), (dd_dsi(j,i),j=1,4)

*        Now form the partials
         do j = 1, num_site
            apion(j,i)=(ow_ionp(j,dd_dsi(1,i))-ow_ionp(j,dd_dsi(2,i)))-
     .                 (ow_ionp(j,dd_dsi(3,i))-ow_ionp(j,dd_dsi(4,i)))
         end do

*        Now form the double difference covariance matrix.  Here we add 
*        a contribution if the one-ways used in the forming the double
*        difference are the same.  The dd_sng accounts for the different
*        signs that the oneways may be used in.
         j = 0
         do l = 1, num_dd, num_dtype
            j = j + 1
            ion_cov(i,j) = 0.d0
            do m = 1,4
               do n = 1,4
                  if( dd_dsi(m,i).eq.dd_dsi(n,j) ) then
                      ion_cov(i,j) = ion_cov(i,j) + 
     .                               dd_sgn(m)*ow_iov(dd_dsi(m,i))*
     .                               dd_sgn(n)*ow_iov(dd_dsi(n,j))
                  end if
               end do
            end do
         end do
      end do

****  Thats all (unless we form covariance matrix at the same time?)

      return
      end

CTITLE POS_ANAL

      subroutine pos_anal(type)

      implicit none

*     Routine to look at the mean and RMS of the kinematic postions to
*     see if the site is truly moving.  If it appear stationary then
*     use a fixed value

      include 'track_com.h'


* PASSED VARIABLES 
* type -- Type of analysis being run

      character*(*) type

* LOCAL VARIABLES
* ks, ns  -- Kinematic and main site number
* ep      -- Epoch counter
* i       -- Loop counter
* num_sum -- Number of values in sum
* gep     -- Epoch number of first good epoch

      integer*4 ks, ns, ep, i, num_sum, gep

* sum_mean(3), sum_sqr(3) -- Sum of mean offset (XYZ) and res**2 (XYZ)
* mean(3), rms(3)  -- Mean and RMS in XYZ of trajectory
* dpos  -- Position diffference

      real*8 sum_mean(3), sum_sqr(3), mean(3), rms(3), dpos

* good_apr -- Set true if the apriori values look good.

      logical good_apr, kbit


****  OK, Loop over all the kinematic sites
      write(*,205) type 
 205  format('POSITION ANALYSIS BASED ON ',a)
      do ks = 1, num_kine

*         Clear the mean and RMS arrays
          ns = kine_to_site(ks)
          do i = 1, 3
             sum_mean(i) = 0.d0
             sum_sqr(i)  = 0.d0
             num_sum     = 0
          end do

****      Now loop over all the epochs
          do ep = 1, num_epochs
             if( num_chan_se(ns,ep).gt.3 ) then
*                Range position should be OK
                 num_sum = num_sum + 1
                 do i = 1, 3
                    dpos = kine_xyz(i,ks,ep)-site_apr(i,ns)
                    sum_mean(i) = sum_mean(i) + dpos
                    sum_sqr(i)  = sum_sqr(i)  + dpos**2
                 end do
             end if
          end do

****      Now see what we have
          if( num_sum.gt.0 ) then
             static(ks) = .true.
             good_apr = .true.
             do i = 1,3 
                 mean(i) = sum_mean(i)/num_sum
                 rms(i)  = sqrt(abs((sum_sqr(i) - 
     .                           mean(i)**2*num_sum)/num_sum)) 
                 if( rms(i).gt. dynamic_tol .and.
     .               mar_site(i,ns).gt. 0  ) then
                     static(ks) = .false.
                 end if
                 if( static(ks).and. abs(mean(i)).gt.3.d0 .and.
     .               apr_site(i,ns).ge.3.0 ) then
* MOD TAH 100518: Only mark bad if apriori sigmas were loose.
                      good_apr = .false.
                 end if
             end do 

****         Tell user the status:
             if( static(ks) ) then
                 if( good_apr ) then 
                     write(*,210) site_names(ns), rms, mean
 210                 format('Kinematic site ',a4,' appears static ',
     .                      ' Coordinate RMS XYZ ',3F6.2,' m, ',
     .                      ' Apriori coordinates good: Diff XYZ ',
     .                      3F6.2,' m')
                 else
                     write(*,220) site_names(ns), rms, mean
 220                 format('Kinematic site ',a4,' appears static ',
     .                      ' Coordinate RMS XYZ ',3F6.2,' m, ',
     .                      ' Apriori coordinates bad: Diff XYZ ',
     .                      3F6.2,' m')
                     write(*,225) site_names(ns), 
     .                      (site_apr(i,ns)+mean(i),i=1,3)
 225                 format('Initial coordinates of ',a4,' are XYZ ',
     .                       3F15.3,' m')

                 end if
*                If we appear static: Save one site of coordinates
*                either mean or apriori
                 do ep = 1, num_epochs
                    if( good_apr ) then
                        do i = 1,3
                           kine_xyz(i,ks,ep) = site_apr(i,ns)
                        end do
                    else
                        do i = 1,3
                           kine_xyz(i,ks,ep) = site_apr(i,ns)+mean(i)
                        end do
                    end if
                 end do  
             else
                 write(*,230) site_names(ns), rms
 230             format('Kinematic site ',a4,' appears dynamic ',
     .                      ' Coordinate RMS XYZ ',3F12.2,' m. ')
                 gep = 0
                 do ep = 1, num_epochs
                    if ( gep .eq.0 .and. kbit(kine_OK(1,ks),ep) ) then
                         gep = ep
                    endif
                 end do
		 
		 if( gep.eq.0 ) gep = 1

                 write(*,235) site_names(ns), (kine_xyz(i,ks,gep),
     .                        i = 1,3), gep
 235             format('Initial coordinates of ',a4,' are XYZ ',
     .                   3F15.3,' m at Epoch ',i8)

            end if

         else
*           We seem to have no data:
            write(*,240) site_names(ns)
 240        format('Kinematic site ',a4,' appears to have no data')
         end if

      end do

***** Thats all
      return
      end

CTITLE RESOLVE_FLOAT

      subroutine resolve_float(iter, iter_res, inc_res)

      implicit none

*     Routine to resolve the floating point ambiguities.
* MOD TAH 071028: Changed calculation order plus updated solutin
*     as values are resolved.
      
      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* iter  -- Iteration number
      integer*4 iter
      integer*4 iter_res   ! Number of iterations in the resolve code
     ,,         inc_res    ! Incremental number of resolved ambiguities

* LOCAL VARIABLES
* ne -- Number of observables (1 for LC, 2 for L1+L2)
* na -- Ambiguity number
* j  -- Loop counter
* NL1, NL2 -- Final estimates of number of ambiquities to be removed
* TL1, TL2 -- Trial values for the L1 and L2 ambqiuities
* RL2  -- Ranges of L2 values that need trying (based on MW-wl sig)
* L1_min, L1_max -- Range of L1 to try for the LC solution

      integer*4 ne, na, j, NL1, NL2, TL1, TL2, RL2, L1_min, L1_max

* dL1, dL2  -- Differences between estimate and integer for:
*       L1 and L2 when L1+L2 used for float; and
*       LC and WL when LC used for float
* dLC  -- Difference between estimate and integer based value of LC
* TLC  -- Trial value of LC ambiquity
* best_fit, next_fit -- Best and next fit to LC for the LC analysis
* relrank -- Relative rank of best choice
* dWL  -- Error in widelane
* fit  -- Compostite fit of LC and WL for judging best fit.

      real*8 dL1, dL2, dLC, TLC, best_fit, next_fit,
     .       relrank, relrank_L1, relrank_L2, dWL, fit, 
     .       dLg, TLg

      real*8 relrank_luse  ! Used relrank limit (1/(1-exp(-it/ref))

* fit_conts(3) -- Contributions to the factors for the best fit case
* fit_fin(3)   -- Contributions for final values
      real*8 fit_conts(3), fit_fin(3)

* Fcode -- Character string code to indicate why a bias is not fixed
* Mapping of the code is:
* Char #  Code Reason
*   1      S   Floating point uncertainity
*   2      W   Wide lane uncertainity
*   3      R   Relative rank
*   4      C   Chi**2 increment too Large
*   5      O   Other biases not fixed
* For an L1+L2 type solution the Fcode is:
* 

      character*5 Fcode

* resolved(2) -- Logical for L1 and L2 resolved
* kbit  -- Check bit status
* done  -- Set true when iteration finished

      logical resolved(2), kbit, done, err
 
* amb_tested(int((max_ambs-1)/32)+1)  ! set when tested
      integer*4 amb_tested(((max_ambs-1)/32)+1)

* best_sig -- Best float sigma for estimate
      real*8 best_sig

***** Loop over the ambiquities and see which ones we can fix.
*     First see how many types there are (ie LC or L1+L2)
      ne = 0
      if( index(float_type,'LC').gt.0 )  ne = 1
      if( index(float_type,'L1').gt.0 .and.
     .    index(float_type,'L2').gt.0 )  ne = 2
      lambda(1) = (vel_light/fR1)
      lambda(2) = (vel_light/fR2) 

* MOD TAH 100711: Compute the relrank_limit for this iteration
      relrank_luse = relrank_limit/(1.d0 - exp(-iter/relrank_iter))
      write(*,110) iter, iter_res, relrank_luse , relrank_iter
 110  format('* Iteration ',I4,' Sub-iter ',i4,
     .       ' RelRank ',F10.2,' RelRank iter ',f6.3)

      if( ne.eq.1 ) then
         write(*,120) 
 120     format('* BF Site  S  PRN       Epoch Range     F     ',
     .          'Estimate dLC       SigLim Relative Rank  ',
     .          'Fix Fcode Change   L1   L2 Residual  L1       ',
     .          'L2   Fits   Best    LC    WL    LG')
      else
         write(*,140)
 140     format('* BF  Site S  PRN      Epoch Range    F     ',
     .          'Estimate L1         Estimate  L2    Sig Limit ',
     .          'Relative Rank  Fix Fcode Change  L1   L2  ',
     .          'Residual  L1       L2')
      end if

*     OK, now loop over ambiquities and see which ones we have
*     estimated and can resolve
* MPD TAH 071028: Loop in order of best determined.
      done = .false.
      inc_res = 0 
      do j = 1, int((max_ambs-1)/32)+1
         amb_tested(j) = 0
      end do


      do while ( .not. done )

*        Find best sigma and any that are left.
         best_sig = 1.d10
         na = 0
         do j = 1, num_ambs
            if( .not.kbit(bf_ents(5,j),2) .and. 
     .          .not. kbit(amb_tested,j)   .and.
     .          wls_ref(1,j).gt.0 ) then
*               Test sigma on L1 or LC depending on ne. 
                if( ambiq_flsig(1,j).lt.best_sig ) then
                   na = j 
                   best_sig = ambiq_flsig(1,j)
                end if
            end if
         end do
*        If nothing found then we are done for this iteration
         if( na.eq.0 ) done = .true.

         IF( .NOT. DONE ) THEN

*           Show we have tested. Could be improved but second iteration
*           should catch
            call sbit(amb_tested,na,1)
* MOD TAH 100214: Set the 'O' bit in bf_ents(5,na) to zero.  If can't be
*           resolved due to other biases, this bit will be set to 1 below
            call sbit(bf_ents(5,na),4,0)
*
*           Set the code depending on LC or L1+L2 float type
            if( ne.eq.1 ) then
               Fcode = 'SWRCO'
            else
               Fcode ='SSRRO'
            endif 

            if( .not.kbit(bf_ents(5,na),2) .and. 
     .           wls_ref(1,na).gt.0 ) then

*               OK, free bias.  See what the adjustment is and if
*               we can resolve.  Do this separately for L1+L2 versus
*               LC
                resolved(1) = .false.
                resolved(2) = .false.
                relrank = 0
                nL1 = 0
                nL2 = 0
                dL1 = 0.d0
                dL2 = 0.d0

*               See if we will check ambiquiuties.  Here if we are within
*               a factor of 4 of the tolerance we try (this way the user
*               can see how far off they are).  We don't want to make the
*               factor of 4 too large or else some very poor values could
*               used to set the initial estimates.
                if( ne.eq.1 .and. 
     .             ambiq_flsig(1,na)/lambda(1)
     .                                 .le.4*float_limit(1) ) then
*                  LC resolution.  Check the L1/L2 combination which
*                  if consistent with the MW wide lane and the LC bias.
*                  OK, see how large our limits may be based on the noise
*                  in the MW-widelane (NOTE: lcf2 is negative).
*                  Wide lane sigma OK

                   resolved(1) = .true.
                   resolved(2) = .true.
                   RL2 = nint(abs(wls_sig(1,na))+
     .                        abs(wls_all(1,na))) + 1
c                   RL2 = nint(abs(wls_sig(1,na)))
                   L1_max = nint(ambiq_float(1,na)/lambda(1)-lcf2*RL2)/
     .                      (lcf1+lcf2) + 1
                   L1_min = nint(ambiq_float(1,na)/lambda(1)+lcf2*RL2)/
     .                      (lcf1+lcf2) - 1
                   best_fit = 9999.d0
                   next_fit = 9999.d0
                   fit_fin(1) = 9999.d0
                   fit_fin(2) = 9999.d0
                   fit_fin(3) = 9999.d0

                   if( na.eq.0 )
     .             print *, 'RESFL ',na, L1_min, L1_max, RL2,
     .                        ambiq_float(1,na) 
                   do TL1 = L1_min, L1_max
                      do TL2 = TL1-RL2, TL1+RL2
*                        Compute the expected value of LC
                         TLC = (lcf1*TL1 + lcf2*TL2)*lambda(1)
                         dLC = ambiq_float(1,na) - TLC
                         dWL = wls_all(1,na) + (TL1-TL2)
                         TLG = exf1*TL1+exf2*TL2
                         dLG = wls_all(2,na) + TLG

* MOD TAH 010306: Place a lower bound on the sigma of the
*                        dLC estimate.
*                        Save the contributions
                         fit_conts(1) = 
     .                         (dLC/(ambiq_flsig(1,na) + min_lcsig))**2
                         fit_conts(2) = 
     .                         (dWL/wls_sig(1,na))**2
                         fit_conts(3) = 
     .                         (dLG/wls_sig(2,na))**2
                         fit = fit_conts(1)+
     ,                         wl_fact*fit_conts(2)+
     .                         lg_fact*fit_conts(3)
*                        fit = (dLC/(ambiq_flsig(1,na)+0.01))**2 +
*    .                         wl_fact*(dWL/wls_sig(1,na))**2 +
*    .                         lg_fact*(dLG/wls_sig(2,na))**2
C                        if( abs(dLC).lt.best_fit ) then
                         if( abs(fit).lt.best_fit ) then
                             next_fit = best_fit
C                             best_fit = abs(dLC)
                             fit_fin = fit_conts
                             best_fit = abs(fit)
                             NL1 = TL1
                             NL2 = TL2

C                        else if( abs(dLC).lt.next_fit ) then
C                            next_fit = abs(dLC)
                         else if( abs(fit).lt.next_fit ) then
                             next_fit = abs(fit)
                         end if
                         if( na.eq.0 )
     .                   write(*,998) 'RESFL ', na, tl1, tl2, tlc, tlg, 
     .                              dlc,ambiq_flsig(1,na),
     .                              dlg,best_fit, next_fit,
     .                              nl1, nl2, dWL, fit_conts
 998                    format(a,3i4,' FLC/G ',2F10.4,' dLC +- ',3f12.3,' Best/Next ',
     .                        ' Best/Next ', 2F12.2, ' NL1/2 ',2i4,
     .                        ' dWL ',F10.3,' FitC ',3f12.3)
                      end do
                   end do

****               OK, now check the constrast
                   relrank = (next_fit/best_fit)**2

*                  Now compute how close we came to getting it correct
                   dL1 = (ambiq_float(1,na) - 
     .                   (lcf1*NL1+lcf2*NL2)*lambda(1))/
     .                   ((lcf1+lcf2)*lambda(1))
                   dL2 = wls_all(1,na) + (nL1-nL2)
* MOD TAH 100515: Added check on second site being fixed.
                   if( relrank.lt.relrank_luse .or.
     .                 ambiq_flsig(1,na)/lambda(1).gt.
     .                 float_limit(1)  .or. best_fit.gt.max_fit .or.
     .                 wls_sig(1,na).gt. float_limit(2) .or.
     .                 .not.kbit(bf_ents(5,WLS_ref(1,na)),2,1).or.
     .                 .not.kbit(bf_ents(5,WLS_ref(2,na)),2,1) ) then
                       resolved(1) = .false.
                       resolved(2) = .false.
                   
*                      See why not resolved
                       if( relrank.ge.relrank_luse ) Fcode(3:3) = '-'
                       if( ambiq_flsig(1,na)/lambda(1).le.
     .                     float_limit(1)) Fcode(1:1) = '-'
                       if( best_fit.le.max_fit ) Fcode(4:4) = '-'
                       if( kbit(bf_ents(5,WLS_ref(1,na)),2,1).and.
     .                     kbit(bf_ents(5,WLS_ref(2,na)),2,1)) then
                           Fcode(5:5) = '-'
                       else   ! Set the 'O' bit on to show this is cause
                           call sbit(bf_ents(5,na),4,1)
                       end if
                        
                       if( wls_sig(1,na).le. float_limit(2)) 
     .                                           Fcode(2:2)='-'
                   else
                       Fcode = '-----' 
                   end if

                else if ( ne.eq.2 ) then

*                  Resolve L1 and L2 separately (although both have to
*                  be resolvable)
                   resolved(1) = .true.
                   resolved(2) = .true.
                   dL1 = ambiq_float(1,na)/lambda(1) - 
     .                    nint(ambiq_float(1,na)/lambda(1))
                   best_fit = dL1
                   NL1 =  nint(ambiq_float(1,na)/lambda(1))
                   relrank_L1 = ((1.d0-abs(dL1))/dL1)**2

                   if( relrank_L1.lt.relrank_luse .or.
     .                 ambiq_flsig(1,na)/lambda(1).gt.float_limit(1).or.
     .                 .not.kbit(bf_ents(5,WLS_ref(1,na)),2,1) .or.
     .                 .not.kbit(bf_ents(5,WLS_ref(2,na)),2,1) ) then
                       resolved(1) = .false.
*                      See which peices were resolved OK
                       if( relrank_L1.ge.relrank_luse )
     .                           Fcode(3:3) = '-'
                       if( ambiq_flsig(1,na)/lambda(1).le.
     .                      float_limit(1) ) Fcode(1:1) = '-'
                       if( kbit(bf_ents(5,WLS_ref(1,na)),2,1) .and.
     .                     kbit(bf_ents(5,WLS_ref(2,na)),2,1)) then
                          Fcode(5:5) = '-'
                       else  ! Set O bit on
                          call sbit(bf_ents(5,na),4,1)
                       end if

                   else
                       Fcode(1:1) = '-'
                       Fcode(3:3) = '-'
                       Fcode(5:5) = '-'
                   end if

*                  Now do L2
                   dL2 = ambiq_float(2,na)/lambda(2) - 
     .                    nint(ambiq_float(2,na)/lambda(2))
                   NL2 =  nint(ambiq_float(2,na)/lambda(2))
                   relrank_L2 = ((1.d0-abs(dL2))/dL2)**2
                   if( relrank_L2.lt.relrank_luse .or.
     .                 ambiq_flsig(2,na)/lambda(2).gt.float_limit(1).or.
     .                 .not.kbit(bf_ents(5,WLS_ref(1,na)),2,1).or.
     .                 .not.kbit(bf_ents(5,WLS_ref(2,na)),2,1) ) then
                       resolved(2) = .false.
*                      See which peices were resolved OK
                       if( relrank_L2.ge.relrank_luse )
     .                           Fcode(4:4) = '-'
                       if( ambiq_flsig(2,na)/lambda(2).le.
     .                      float_limit(1) ) Fcode(2:2) = '-'
                       if( kbit(bf_ents(5,WLS_ref(1,na)),2,1).and.
     .                     kbit(bf_ents(5,WLS_ref(2,na)),2,1) ) then
                          Fcode(5:5) = '-'
                       else
                          call sbit(bf_ents(5,na),4,1)
                       end if
                   else
                       Fcode(2:2) = '-'
                       Fcode(4:4) = '-'
                       Fcode(5:5) = '-'
                   end if

                   relrank = min(relrank_L1,relrank_L2)
                end if

****            Now tell user what is happening. 
                if( resolved(1) .and. resolved(2) ) then
                    call sbit(bf_ents(5,na),2,1)
                    num_tot_resolved = num_tot_resolved + 1
                    inc_res = inc_res + 1
                end if

                if( ne.eq.1 ) then
                    write(*,410) na, site_names(bf_ents(1,na)), 
     .                  bf_ents(1,na), (bf_ents(j,na),j=2,5),
     .                  ambiq_float(1,na)/lambda(1),
     .                  ambiq_flsig(1,na)/lambda(1),
     .                  float_limit(1),
     .                  relrank, resolved, Fcode, nL1, nL2, 
     .                  dL1, dL2,  best_fit, fit_fin
 410                format(i4,1x,a4,1x,i2,1x,' PRN ',i3.2, 2I8,3x,o2.2,
     .                     (F9.2,' +-',  F8.2),' SL ',F6.2,
     .                     ' RR ',f10.2,1x,2L2,1x,a5,' dL1,2 ',2I7,
     .                     ' dL12 ',2F9.2,' Fits ',4F6.1)
                else
                    write(*,420) na, site_names(bf_ents(1,na)), 
     .                   bf_ents(1,na), (bf_ents(j,na),j=2,5),
     .                   (ambiq_float(j,na)/lambda(j),
     .                    ambiq_flsig(j,na)/lambda(j), j=1,ne),
     .                   float_limit(1),
     .                   relrank, resolved, Fcode,
     .                   nL1, nL2, dL1, dL2
 420                format(i4,1x,a4,1x,i2,1x,' PRN ',i3.2, 2I8,3x,o2.2,
     .                     2(F9.2,' +-',  F8.2),' SL ',F6.2,
     .                     ' RR ',f10.2,1x,2L2,1x,a5, ' dL1,2 ',2I7,
     .                     ' dL12 ',2F9.2)
                end if

                if ( iter.gt.1 ) then
*                   After iteration 1 updated
!                    ambiq_all(1,na) = ambiq_all(1,na) - NL1
!                    ambiq_all(2,na) = ambiq_all(2,na) - NL2
                    call update_wls_all(na, -nL1, -nL2, .true. )
                end if
                if( resolved(1) .and. resolved(2) ) then
* MOD TAH 071028:   Moved update to only when values resolved.
*                   Also added progation into other terms.  If 
*                   first iteration would not have been done above
                    if( iter.le.1 ) then 
 !                      ambiq_all(1,na) = ambiq_all(1,na) - NL1
 !                      ambiq_all(2,na) = ambiq_all(2,na) - NL2
                       call update_wls_all(na, -nL1, -nL2, .true. )
                    endif
*
!                   Update_wls_all does update the float estimates
!                   here we use 0 0 as the change in L1/L2 cycles.
!                   call update_estsol(na, ne, NL1, NL2, err)
                    call update_estsol(na, ne, 0, 0, err)
                    if( err ) then
C                       do j = 1, num_ambs
C                          print *,'AMP_PARN ',j,  amb_parn(1,j), 
C     .                            amb_save(1,j)
C                       end do
c                       resolved(1) = .false.
c                       resolved(2) = .false.
c                       call sbit(bf_ents(5,na),2,0)
c                       num_tot_resolved = num_tot_resolved - 1
                    end if
               endif
 
            end if
         END IF
      end do 

****  Thats all
      return
      end
 
CTITLE FIX_MS

      subroutine fix_ms(ep, nd, ns, lv, msi)

      implicit none

*     Routine to apply millisecond corrections to phase and/or range
*
      include '../includes/const_param.h'
      include 'track_com.h'


* PASSED VARIABLES
      integer*4 ep   ! Epoch
     .,         nd   ! Double difference number
     .,         ns   ! Site number
     .,         lv   ! Satellite PRN Number
     .,         msi  ! Correction (ms) to be applied

* LOCAL VARIABLES
      integer*4 ch   ! Channel number
     .,         j    ! Loop counter

      logical kbit   ! Test bit status (true if set)

*     First map the PRN number back to a channel for this site at this epoch
      ch = -1
      do j = 1,num_chan_se(ns,ep)
         if( ctop_cse(j,ns,ep).eq.lv ) ch = j
      end do
      if( ch.le.0 ) then
         write(*,120) ep, nd, lv 
 120     format('** ERROR ** Could not map site/prn to channel',
     .          ' at epoch ',I7,' Site # ',i3,' PRN ',i3.2)
         return
      end if

****  Now see which data types are affected 
      if( kbit(dd_dt(nd),1) ) then      ! L1 phase
          L1o_all_cse(ch,ns,ep) = L1o_all_cse(ch,ns,ep) - 
     .                            msi*1.d-3*fR1
      end if

      if( kbit(dd_dt(nd),2) ) then      ! L2 phase
          L2o_all_cse(ch,ns,ep) = L2o_all_cse(ch,ns,ep) - 
     .                            msi*1.d-3*fR2
      end if
      if( kbit(dd_dt(nd),3) ) then      ! P1 range
          P1o_all_cse(ch,ns,ep) = P1o_all_cse(ch,ns,ep) - 
     .                            msi*1.d-3*vel_light
      end if

      if( kbit(dd_dt(nd),4) ) then      ! P2 range
          P2o_all_cse(ch,ns,ep) = P2o_all_cse(ch,ns,ep) - 
     .                            msi*1.d-3*vel_light
      end if

****  Thats all 
      return
      end

CTITLE CHECK_CLK

      subroutine check_clk(ep, ns, omc, avclk)

      implicit none

*     Routine to check that the clock epoch calculation looks OK 

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES 

      integer*4 ep   ! Epoch number
     .,         ns   ! Site number

      real*8 omc(max_chan)   ! OMC values
     .,      avclk           ! Input value for the clock error.

* LOCAL VARIABLES
      integer*4 j
     .,         nc    ! count of number of values in clock estimate
     .,         max_ind   ! Index of worst residual

      logical omc_OK(max_chan)  ! Set true to indicate OMC OK,
     .,       data_OK     ! Function to check data

      real*8 rms, sum  ! Sums for rms and mean value
     .,      max_err   ! Worst residual

      real*8 noise_tol ! Clock noise tol based on apriori position
                       ! simga
      real*8 absomc(max_chan)  ! Absolute value of error, used to find 
                               ! median
     .,      median_absomc     ! Median abs omc
     .,      median_omc        ! Median OMC as estimate of clock error

***** Loop over omc and make RMS looks OK 
      if( num_chan_se(ns,ep).le.1 ) RETURN

* MOD TAH 101030: Modified error calculation to include process
*     noise on first past 
C     noise_tol = max(sqrt(apr_site(1,ns))*10.d0,100.d0)
      noise_tol = max(sqrt(apr_site(1,ns)**2+mar_site(1,ns)**2)*10.d0,
     .                40.d0)
!      if( .not. kine_known ) then  ! When kine not known; make the values 
!                    ! large to allow motion of vechile; allow 300 m/s 
!                    ! motions
!        noise_tol = max(noise_tol,300.d0*usr_interval)
!      end if

! MOD TAH 110520: Find the median abs error to get a measure of the noise
! MOD TAH 120611: Refined this use of the median error; remove the median
!     OMC first and then find median max difference.
      call mdian2(omc, num_chan_se(ns,ep), median_omc)
      do j = 1,num_chan_se(ns,ep)
         absomc(j) = abs(omc(j)-median_omc)
      end do
!     Find the median and remove; then see max difference 
      call mdian2( absomc, num_chan_se(ns,ep), median_absomc )

      noise_tol = min(max(noise_tol, 6*median_absomc),999.0)


      rms = 1000.d0
      do j = 1, num_chan_se(ns,ep)
         omc_ok(j) = data_OK(data_flag_cse(j,ns,ep),data_mask)
      end do

      do while ( rms.gt.noise_tol )
         rms = 0.d0
         sum = 0.d0
         max_err = 0.d0
         nc = 0
         do j = 1, num_chan_se(ns,ep)
            if ( omc_ok(j) ) then 
                nc = nc + 1
                if( abs(omc(j)-avclk).ge.max_err ) then
                    max_err = abs(omc(j)-avclk)
                    max_ind = j
                endif
                rms = rms + (omc(j)-avclk)**2
                sum = sum + omc(j)
            end if
         end do 
         if( nc.gt.1 ) then
             rms = sqrt(rms/(nc-1))
             if( rms.gt.noise_tol ) then
*               Adjust clock for BAD value
                avclk = (sum-omc(max_ind))/(nc-1)
                write(*,120) ep, site_names(ns), 
     .                       ctop_cse(max_ind,ns,ep),
     .                       max_err, rms, avclk, noise_tol
 120            format('+ Bad Clock estimate: Epoch ',i7,1x,a4,' PRN ',
     .                 I3.2,' Error ',f10.2,' m; RMS ',F10.2,' m',
     .                 ' Clock ',F10.2,' m, Tol ',F10.2,' m')
                omc_ok(max_ind) = .false.
*               Mark the data point as bad.
                call sbit(data_flag_cse(max_ind,ns,ep),2,1)

             end if
         else
             write(*,160) ep, site_names(ns),
     .                    (ctop_cse(j,ns,ep), omc(j), 
     .                        j=1,num_chan_se(ns,ep))
 160         format('+ Clock Warning: All data deleted in clock ', 
     .              'estimate at Epoch ',i7,' Site ',a4,' OMC: ',
     .              32('PRN',i3.2,1x,F10.2,' m '))

             rms = 0.d0
         end if
      end do

****  Thats all
      return
      end


CTITLE UPDATE_ESTSOL

      subroutine update_estsol(na, ne, dNL1, dNL2, err)

      implicit none

*     Routine to update the sol_vec, cov_parm and amgiq values
*     when a term has been resolved.

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
      integer*4 na   ! Ambiguity number that was just resolved
     .,         ne   ! Number of equations (1 == LC, 2 = L1+L2)
     .,         dNL1, dNL2 ! Number of cycles change at L1 and L2

      logical err

* LOCAL
      integer*4 i,j   ! Loop counters
     .,         ipivot(2)  ! Used in inversion
     .,         np(2) ! Parameters to be forced.

*   cov_col(max_parm, 2)    - The columns of the
*               -  covariance matrix for the parameters
*               - being forced.
*   sol_col(2)      - The change in the parameter
*               - estimates needed to get the forced values.
*   var_force(2) - Variances of focced parameters
*   dsol(2)     - Variance with which value should be forced.
*   avat(2,2)   - The multiplication of
*               - the partials matrix of (00001000,000001000)
*               - and the covaraiance matrix.
*   equ_gn(max_parm)  - Kalman gain matrix for
*               - getting solution.
*   scale(2)        - Scaling vector for solution. (passed
*               - to invert_vis)
*   dchi        - Change in Chi**2
*   sol_save(2) - Saved values of solution vector
 
      real*8 cov_col(max_parm, 2), sol_col(2),
     .    dsol(2), var_force(2),
     .    avat(2,2),
     .    equ_gn(max_parm,2), scale(2), dchi, sol_save(2)

* Saved values of the solution vector and sigme before
*     update.  These are compared to post-values and large
*     changes are reported.
* dsolp, dsigp -- Change in values
      real*8 sol_prior(max_parm), sig_prior(max_parm),
     .       dsolp, dsigp
 

****  Compute value to which sol_vec will need to be forced
C      np(1) = amb_parn(1,na)
      np(1) = amb_save(1,na)
      err = .false.
      if( np(1).le.0 ) then
          print *,'ERROR in update_estsol Amb ',na,np
          err = .true.
          RETURN
      endif

****  Save values before update
      do i = 1, num_parm
         sol_prior(i) = sol_vec(i)
         sig_prior(i) = sqrt(cov_parm(i,i))
      end do

****  OK get change value
      if( ne.eq. 1 ) then
          dsol(1) = (lcf1*dNL1+lcf2*dNL2)*lambda(1)
          var_force(1) = 1.d-9
          sol_save(1) = sol_vec(np(1))

          call force_track( cov_parm, sol_vec, max_parm, num_parm,
     .         cov_col, sol_col, np, 1,  dsol,
     .         var_force, dchi, avat, equ_gn, ipivot, scale)
c          print *,'Force1 ',na, np(1),sol_save(1), sol_vec(np(1)),
c     .                     dsol(1),sol_save(1)-sol_vec(np(1)),dchi

*         Now resave ambiq values 
          do i = 1,num_ambs
             if( amb_parn(1,i).gt.0 ) then
                ambiq_float(1,i) = sol_vec(amb_parn(1,i))
                ambiq_flsig(1,i) = sqrt(cov_parm(amb_parn(1,i), 
     .                                           amb_parn(1,i))) 
c                print *,'NEW EST Amb ',i,amb_parn(1,i),
c     .                  ambiq_float(1,i),ambiq_flsig(1,i) 
             endif
          end do
*     else do L1 an dL2
      else if( amb_parn(2,na).gt.0 ) then
          np(2) = amb_parn(2,na)
          dsol(1) = dNl1*lambda(1)
          dsol(2) = dNL2*lambda(2)
          var_force(1) = 1.d-9
          var_force(2) = 1.d-9
          sol_save(1) = sol_vec(np(1))
          sol_save(2) = sol_vec(np(2))

          call force_track( cov_parm, sol_vec, max_parm, num_parm,
     .         cov_col, sol_col, np, 2,  dsol,
     .         var_force, dchi, avat, equ_gn, ipivot, scale)

*         Now resave ambiq values 
          do i = 1,num_ambs
             do j = 1,ne
                if( amb_parn(j,i).gt.0 ) then
                   ambiq_float(j,i) = sol_vec(amb_parn(j,i))
                   ambiq_flsig(j,i) = sqrt(cov_parm(amb_parn(j,i), 
     .                                              amb_parn(j,i))) 
                endif
             end do
          end do
      else
          print *,'ERROR in update_estsol amb ',na,' Expected L2'
      endif

****  Now compare values and output big values
C     write(*,510) np(1), dsol(1), sol_prior(np(1)), sol_vec(np(1)),
C    .             sig_prior(np(1)),cov_parm(np(1),np(1))
C510  format('FORCE # ',i4,' Dsol/prior/fin ',3F12.5,' Sig/Var ',
C    .       F12.5,1x,E20.4)
C     do i = 1, num_parm
C        if ( (i.ne. np(1) .and. ne.eq.1) .or.
C    .        (ne.eq.2 .and. i.ne.np(1) .and. i.ne.np(2)) ) then
C           dsolp = sol_vec(i) - sol_prior(i)
C           dsigp = sqrt(sig_prior(i)**2-cov_parm(i,i))
C           if( abs(dsolp/dsigp).gt. 3.d-2 ) then
C              write(*,520) i, np, sol_prior(i), sol_vec(i),
C    .             sig_prior(i), sqrt(cov_parm(i,i)), dsigp, 
C    .             abs(dsolp/dsigp)
C520           format('Large Param Change ',i3,' from dP ',2i4,
C    .                ' Prior/Fin Ests ',2F12.5,' Sigmas ',3F12.6,
C    .                ' N-sig ',F6.2)
C           end if
C        end if
C     end do

****  Thats all
      return
      end

CTILE FORCE_TRACK

      subroutine force_track( cov_parm, sol_parm, dim_parm, num_parm,
     .    cov_col, sol_col, nc, num_force, force_values,
     .    force_var, dchi, avat, equ_gn, ipivot, scale )

      implicit none
 
*     This routine will force parameters in a covariance matrix to have
*     specific values.  The parameters to be forced are given in the
*     NC array.  There are num_force of them and the values to be forced
*     to are given in the force_values array.  If local array mapping
*     is to be then map_force should be called before this routine
*     setup the mapping of the arrays needed.
 
*   dum_parm    - Dimensions of cov_parm
*   num_parm        - Number of parameters in cov_parm.  Cov_parm
*               - is assumed to be dimensioned to this size.
*   num_force   - Number of parameters to be focesd to specific
*               - values
*   nc(num_force)   - The list of parameters to be made equal.
*   ipivot(num_force)   - The pivot elements for the matrix
*               - inversion
 
      integer*4 dim_parm, num_parm, num_force, 
     .          nc(num_force), ipivot(num_force)
 
*   cov_parm(num_parm,num_parm) - Full covariance matrix to
*               - be modified.
*   sol_parm(num_parm)      - The solution vector.  The nc
*               - elements of this vector will be set equal.
*   cov_col(num_parm, num_force)    - The columns of the
*               -  covariance matrix for the parameters
*               - being forced.
*   sol_col(num_force)      - The change in the parameter
*               - estimates needed to get the forced values.
*   force_values(num_force) - The values the forced parameters
*               - should take on.
*   force_var(num_force) - Variance with which value should be forced.
*   avat(num_force,num_force)   - The multiplication of
*               - the partials matrix of (00001000,000001000)
*               - and the covaraiance matrix.
*   equ_gn(num_parm,num_force)  - Kalman gain matrix for
*               - getting solution.
*   scale(num_force)        - Scaling vector for solution. (passed
*               - to invert_vis)
*   dchi        - Change in Chi**2
 
      real*8 cov_parm(dim_parm,dim_parm), sol_parm(dim_parm),
     .    cov_col(dim_parm, num_force), sol_col(num_force),
     .    force_values(num_force), force_var(num_force),
     .    avat(num_force,num_force),
     .    equ_gn(dim_parm,num_force), scale(num_force), dchi
 
* LOCAL PARAMETERS
 
*   i,j,k       - Loop counters
 
      integer*4 i,j,k
 
*   dsol, dcov  - Summation variables for computing corrections
*               - to the solution vector and covariance matrix.
 
 
      real*8 dsol, dcov
 
*                        T
***** Start, form the AVA  matrix, where A is of the form:
*     y = Ax (y is vector of zero observations)
*     and
*     A =  0 0 0 1 0 0 0 0 0 0....... to num_parm
*          0 0 0 0 0 0 0 1 0 0
*          0 0 0 0 0 0 0 0 1 0
*     where the above form would set parameters 4,8, and 9.
*     For num_force parameters being equated, there are num_force
*     rows in A.
*
*     The above form is much simpler to compute is we save the
*     columns of the covariance matrix for thhose parameters to
*     be forced first.
 
      do i = 1, num_parm
          do j = 1, num_force
              cov_col(i,j) = cov_parm(i,nc(j))
          end do
      end do
 
*     Save the forced parts of the solution vector as well
      do j = 1, num_force
          sol_col(j) = sol_parm(nc(j))
      end do
 
*               T
****  Now do AVA
*
      do i = 1, num_force
          do j = 1, num_force
              avat(i,j) = cov_col(nc(i),j)
* MOD TAH 000302: Added variance to diagonal for constrained 
*             forcing.              
              if( i.eq.j ) avat(i,i) = avat(i,i) + force_var(i) 
          end do
      end do
 
****  Now invert this matrix.  (If "sort of equal" was desired we could
*     add value to diagonal now representing variance of y above). Pass
*     zero as number in solution vector, we dont want to multiply.
*     kgain below is dummy argument.
 
      call invert_vis(avat, equ_gn, scale, ipivot, num_force,
     .               num_force, 0 )

**** Before continuing compute the change in Chi**2 due to condition
      dchi = 0.d0
      do i = 1, num_force  
         do j = 1, num_force 
            dchi = dchi + (force_values(i)-sol_col(i))*avat(i,j)*
     .                    (force_values(j)-sol_col(j))
         end do
      end do
      dchi = dchi/num_force  

 
*     Now form the Kalman gain, equ_gn given by
*                T     T -1
*     equ_gn = VA  (AVA )
*
      do i = 1, num_parm
          do j = 1, num_force
 
*             Do the multiply (could use VIS but stick to straight)
*             call dwmul(equ_gn(i,j), col_col(i,1), num_force,
*    .                    avat(1,j),1, num_force)
 
              equ_gn(i,j) = 0.d0
              do k = 1, num_force
                  equ_gn(i,j) = equ_gn(i,j) + cov_col(i,k)*
     .                                     avat(k,j)
              end do
          end do
      end do

****  Now get the change to the solution vector
*
*     x  = x  -  equ_gn*(force_values-sol_col)
*      n    o
*
      do i = 1,num_parm
          dsol = 0.d0
          do j = 1, num_force
              dsol = dsol + equ_gn(i,j)*(force_values(j)-sol_col(j))
          end do
          sol_parm(i) = sol_parm(i) + dsol
      end do

*     Now update the covariance matrix
*
*     V  = V - equ_gn*cov_col
*      n    o
 
      do i = 1, num_parm
          do j = 1, num_parm
 
*             Do summation loop
              dcov = 0.d0
              do k = 1, num_force
                  dcov = dcov + equ_gn(i,k)*cov_col(j,k)
              end do
              cov_parm(i,j) = cov_parm(i,j) - dcov
          end do
      end do
 
****  Thats all
      return
      end

CTITLE DUMP_OWS
     
      subroutine dump_ows(i)

      implicit none

*     Routine to output the one way phase and range residuls
*     for ionospheric delay analysis.

      include '../includes/const_param.h'
      include '../includes/xfile_def.h'
      include 'track_com.h'

* PASSED VARIABLES
* i   -- Site number
      integer*4 i

* LOCAL VARIABLES 
* ep -- Epoch number at which to compute position

      integer*4 ep

* LOCAL VARIABLES
* no, ns, nd, nc -- Short name counters for number of one-ways, single differences,
*     double differences, and number of sd at each site.
* pn        -- Current PRN being processed
* i,j        -- Loop counters

      integer*4 no, nc,pn, j, k 

* rcv_time  -- MJD of time given by the receiver clock
* trn_time  -- MJD of tranmission time as given by satellite clock (needs to be
*              corrected for satellite clock error).
* rcv_sec   -- Seconds from start of SP3 file for time
* trn_sec   -- Seconds from start of SP3 file for tme
* svs_ear(3) -- Satellite coordinates for current PRN and site
* range(2)   -- Range computed to satellite at L1 and L2 (m). (Difference due to
*               antenna offsets at the two frequencies)
* elev       -- Real*8 version of elevation angle (deg)
* az         -- Azimuth (real*8) (deg)

      real*8 rcv_time, trn_time, svs_ear(3), range(2), phase(2),
     .       elev, az, trn_sec


* omc(max_chan) -- Range obs-minus-computed
* prev_range(max_chan)  -- Range to satellite from previous iteration.
* dry_map   -- Dry mapping function, passed from theory into partials
* wet_map   -- Wet mapping function (use for partials)
* max_err   -- Max dd error when checking for outliers (mm)

      real*8 site_clk, dry_map, wet_map

* na -- Ambiguity number
      integer*4 date(5), na, ierr
      real*8 sectag
      real*8 tec_los    ! LOS of tec


      character*512 file  ! Name of output file
      character*16 all_data
      integer*4 und   ! Unit number for dump files

 
****  Open the output file name
      all_data = 'L1L2P1P2'
 
      write(file,110) trim(posit_root), site_names(i)
 110  format(a,'_',a,'_dump.txt')
      und = 119 
      open(und,file=file,iostat=ierr,status='unknown')
      write(und,120)
 120  format('*  Epoch   PRN          dL1 (m)        ',
     .       ' dL2 (m)      dP1 (m)',
     .       '      dP2 (m)     Amb L1 (cyc)    Amb L2 (cyc)     ',
     .       'El (dg)  Az (deg)    DataF  LoS TECU',
     .       '   SI Long    SI Lat ')

      do ep = 1, num_epochs

*        Loop over the satellites at this epoch
         rcv_time = ((ep-1)*usr_interval+sec_offset(i,ep))/86400.d0+
     .                ref_start

*        Compute all the satellite clock corrections
         call comp_svs_clk( rcv_time )

*        Get site coordinates
         call get_kine_apr(ep)
         
****     Now we are ready to start final calculation                 
         do j = 1, num_chan_se(i,ep)

*           Compute the theoretical range (use transmitt time as the
*           receive time minus the psuedorange)
            pn = ctop_cse(j,i,ep)
            trn_time = rcv_time - 20000.d3/vel_light/86400.d0
            trn_sec  = ((ep-1)*usr_interval+sec_offset(i,ep)) +
     .                 ref_sec - 20000.d3/vel_light
       
            call theory(trn_time, trn_sec, i, j, curr_site_xyz(1,i), 
     .                  svs_ear,  range, phase, elev, az, dry_map, 
     .                  wet_map, 'SN', ep) 

            trn_time = rcv_time - range(1)/vel_light/86400.d0
            trn_sec  = ((ep-1)*usr_interval+sec_offset(i,ep)) +
     .                 ref_sec - range(1)/vel_light
       
            call theory(trn_time, trn_sec, i, j,
     .                  curr_site_xyz(1,i), 
     .                  svs_ear,  range, phase, elev, az, dry_map, 
     .                  wet_map, 'FN', ep) 

*           TEST CODE: Compute ion tec_los
            tec_los = 0

            na = amb_point_cse(j,i,ep)
            if( na.gt.0 ) then 
               no = 0 
               nc = 0
               call get_ow_omc(ep, i, j, pn, all_data, no, nc, range, 
     .                  phase, elev, dry_map, wet_map, svs_ear)

               if( use_ionex ) then
                   call interp_ionex(ep,i,j, curr_site_xyz(1,i),
     .                         az, elev,  rcv_time, tec_los )
*                  Since we removed the ion delay from O-C; add back now
*                  (For range we have o-(c+tec) = o-c-tec to to get original
*                  o-c need to add tec value (substract for phase)
                   ow_vec(1) = ow_vec(1) - l1tecu*tec_los  ! Phase (subtract)
                   ow_vec(2) = ow_vec(2) - l2tecu*tec_los
                   ow_vec(3) = ow_vec(3) + l1tecu*tec_los  ! Range (add)
                   ow_vec(4) = ow_vec(4) + l2tecu*tec_los
               end if
 

*              Write out values
               write(und,220) ep, pn, ow_vec(1:4)+svs_clk(pn)*vel_light, 
     .                ambiq_all(:,na),  elev, az, data_flag_cse(j,i,ep),
     .                tec_los, tec_subI_cse(:,j,i,ep)
 220           format(I8,4x,I3,1x,2F16.4,2F13.4,1x,2F16.1,1x,
     .                F10.3,1x,F10.3,1x,O8,1x,F9.3,1x, F9.3,1x,F9.3)
            end if

         end do

      end do
      close(und)

      return
      end

