      subroutine update_outfs(ep)

      implicit none

*     Routine to open the output files and update the names if we have
*     reached the correct output interval

      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep            ! Counter for number of epochs of data processed.

* LOCAL
      integer*4 i   ! Loop over sites
      integer*4 date(5)  ! YMD HM
      integer*4 ierr     ! IOSTAT error
      integer*4 un       ! Generic unit number
      integer*4 trimlen  ! Length of string
      integer*4 run_time(5)  ! Time when files are created.

      real*8 mjd_file ! MJD for naming files
      real*8 sec      ! seconds tag

      character*256 timetag  ! Test time tag to see if matches
      character*256 file     ! Generic file name
    


****  Based on the current time; generate the time tag name
      mjd_file = int(RT_MJD_obs(sblk)/file_updint)*file_updint
      call mjd_to_ymdhms(mjd_file,date,sec)
      if( file_updint.ge.1 ) then ! ge one-day, generate name to 1-day
         write(timetag,110) date(1:3)
 110     format(I4.4,I2.2,I2.2)
      else
         write(timetag,130) date
 130     format(I4.4,I2.2,I2.2,'_',i2.2,I2.2)
      endif

****  See if tag has changed
      if( timetag.eq.curr_timetag ) RETURN  
      curr_timetag = timetag

***** Before resetting the epoch, make the last epoch of all the bias 
*     flags consistent with new counter
! MOD TAH 120613: Loop through as set through the last good epoch
!     values so that ep-bf_ents(6,na) has correct value.  Do this 
!     before resetting epoch
      do i = 1, num_ambs 
         if( bf_ents(4,i).eq.ep ) bf_ents(4,i) = 0 ! consistent with one
                                                   ! before current
         if( bf_ents(6,i).gt.0 )  bf_ents(6,i) = bf_ents(6,i)- ep
      end do

***** Reset the epoch counter to 1 when new file opened (do not set zero
*     since zero epoch is used for initilization).
      ep = 1

****  OK: New set of file names are needed
      do i = 1, num_site
         if( index(posit_type,'GEOD').gt.0 .and.
     .       write_pos(i) ) then 
             file = posit_root(1:trimlen(posit_root)) // 
     .             '.GEOD.' // site_names(i) //
     .             '.' // timetag(1:trimlen(timetag)) // 
     .             '.' // anal_types(1)
             un = 200 + i
             close(un)
             open(un,file=file,status='unknown',iostat=ierr)
             write(*,*) 'Creating ', file(1:trimlen(file))
             write(un,120)
 120         format('*YY  DOY   Seconds     Latitude      ',
     .              'Longitude    Height   SigN  SigE  SigH',
     .              '   RMS  #     Atm     +-       Fract',
     .              ' DOY     Epoch  #BF NotF',/,
     .              '*                       (deg)        ',
     .              '  (deg)       (m)     (cm)  (cm)  (cm)',
     .              '  (mm)  DD    (mm)   (mm)')
          end if
          if( index(posit_type,'NEU').gt.0 .and.
     .        write_pos(i) ) then 
              file = posit_root(1:trimlen(posit_root)) // 
     .             '.NEU.' // site_names(i) //
     .             '.' // timetag(1:trimlen(timetag)) // 
     .             '.' // anal_types(1)
              un = 200 + num_site + i
              open(un, file=file, status='unknown',iostat=ierr)
              write(*,*) 'Creating ',file(1:trimlen(file))
              write(un, 140)
 140          format('* YY  MM DD HR MIN  Sec          dNorth',
     .               '      +-            dEast        +-    ',
     .               '      dHeight      +-        RMS    #  ',
     .               '    Atm     +-       Fract DOY     Epoch',
     .               '  #BF NotF  Rho_UA',/,
     .               '*                                 (m)  ',
     .               '      (m)            (m)        (m)    ',
     .               '     (m)           (m)      (mm)   DD  ',
     .               '    (mm)    (mm)')
          end if
          if( index(posit_type,'XYZ').gt.0 .and.
     .        write_pos(i) ) then 
              file = posit_root(1:trimlen(posit_root)) // 
     .             '.XYZ.' // site_names(i) //
     .             '.' // timetag(1:trimlen(timetag)) // 
     .             '.' // anal_types(1)
              un = 200 + 2*num_site + i
              open(un, file=file, status='unknown',iostat=ierr)
              write(*,*) 'Creating ',file(1:trimlen(file))
              write(un, 160)
 160          format('* YY  MM DD HR MIN  Sec             X  ',
     .               '      +-               Y         +-    ',
     .               '          Z        +-        RMS    #  ',
     .               '    Atm     +-       Fract DOY     Epoch',
     .               '  #BF NotF  Rho_UA',/,
     .               '*                                 (m)  ',
     .               '      (m)            (m)        (m)    ',
     .               '     (m)           (m)      (mm)   DD  ',
     .               '    (mm)    (mm)')
          end if
          if( index(posit_type,'DHU').gt.0 .and.
     .        write_pos(i) ) then 
              file = posit_root(1:trimlen(posit_root)) // 
     .             '.DHU.' // site_names(i) //
     .             '.' // timetag(1:trimlen(timetag)) // 
     .             '.' // anal_types(1)
              un = 200 + 3*num_site + i
              open(un, file=file, status='unknown',iostat=ierr)
              write(*,*) 'Creating ',file(1:trimlen(file))
              write(un, 180)
 180          format('* YY  MM DD HR MIN  Sec          dNorth',
     .               '      +-            dEast        +-    ',
     .               '      dHeight      +-        RMS    #  ',
     .               '    Atm     +-       Fract DOY     Epoch',
     .               '  #BF NotF  Rho_UA',/,
     .               '*                                 (mm) ',
     .               '      (mm)           (mm)       (mm)   ',
     .               '     (mm)          (mm)     (mm)   DD  ',
     .               '    (mm)    (mm)')
          end if
          if( trimlen(csv_root).gt.0 .and.
     .        write_pos(i) ) then
              file = csv_root(1:trimlen(csv_root)) // 
     .             '_' // site_names(i) //
     .             '_' // timetag(1:trimlen(timetag)) // 
     .             '.csv'
              un = 200 + 4*num_site + i
              open(un, file=file, status='unknown',iostat=ierr)
              write(*,*) 'Creating ',file(1:trimlen(file))
              write(un,200) site_names(i), init_geod(2,i),
     .            init_geod(1,i), init_geod(3,i), ant_name(i), 
     .            rcv_type(i)
 200          format('Site ',a4,1x,F9.5,1x,F8.5,1x,F10.3,1x,a20,1x,
     .               ' DCB Type ',a1)
          end if

      end do

****  Now see if new summary file is needed
      if( trimlen(sum_file).gt.0 ) then
         file = sum_file(1:trimlen(sum_file)) // 
     .             '.' // timetag(1:trimlen(timetag)) // 
     .             '.sum'
         lus = 98
         close(lus)
         write(*,*) 'Creating ',file(1:trimlen(file))
         open(lus,file=file,status='unknown',iostat=ierr)
         call report_error('IOSTAT',ierr,'open',file,0,
     .                     'TRACK/SUM_FILE')
         if( ierr.ne.0 ) lus = 6
         call systime(run_time,sec)
         write(lus,320) trackRT_version, run_time
 320     format('SUMMARY FILE: Track Vers ',a,' Run ',
     .          i4,'/',i2.2,'/',i2.2,' ',i2.2,':',i2.2)

         call report_setupRT(lus)
 
      endif

****  Now see if new summary file is needed
      if( trimlen(prt_root).gt.0 ) then
         file = prt_root(1:trimlen(prt_root)) // 
     .             '.' // timetag(1:trimlen(timetag)) // 
     .             '.out'
         write(*,*) 'Creating ',file(1:trimlen(file))
         open(6,file=file,status='unknown',iostat=ierr)
         call report_error('IOSTAT',ierr,'open',file,0,
     .                     'TRACK/prt_root')
         call systime(run_time,sec)
         write(*,420) trackRT_version, run_time
 420     format('OUTPUT FILE: Track Vers ',a,' Run ',
     .          i4,'/',i2.2,'/',i2.2,' ',i2.2,':',i2.2)

         call report_setupRT(6)
 
      endif

****  Thats all
      return
      end


      subroutine output_outfs(ep)

      implicit none

*     Routine to output the position and atmospheric delay estimates at this
*     epoch.
      include '../includes/const_param.h'

      include 'trackRTObs.h'      ! Real-time data structures 
      include 'trackRT.h'         ! Common block

* PASSED 
      integer*4 ep            ! Counter for number of epochs of data processed.

* LOCAL

* geod_coord(3) -- Geodetic coordinates
* xyz_coord(3)  -- XYZ coordinates of current site 
* geod_ref(3)   -- Geodetic coordinates of reference site
* dNEU(3)       -- NEU difference in position (m)
* cov_xyz(3,3)  -- Position XYZ covarinance matrix 
* cov_neu(3,3)  -- Position NEU covariance
* cov_tmp(3,3)  -- Temporary matrix
* sig_neu(3)    -- Sigmas for NEU 
* sectag        -- Seconds tag
* fday          -- Fractional day
* sec_in_day    -- Number of seconds since start of day
* av_static(3,num_site) -- Average GEOD lat/long/and height during static times
* av_square(3,num_site) -- Squares of values to compute RMS
* av_sig(3,num_site) -- Sigam of each component based on RMS
* av_init(3,num_site) -- Initial coordinates of points (improves accuracy in RMS
*            calculation)
* adj_xyz(3), adj_neu(3) -- Adjustment to XYZ and NEU to apriori final
*            output.
      real*8 geod_coord(3), geod_ref(3), dNEU(3), cov_xyz(4,4),
     .       cov_neu(4,4), sig_neu(3), cov_tmp(4,4), sectag, fday, 
     .       sec_in_day, av_static(3,max_site),av_square(3,max_site),
     .       av_sig(3,max_site), av_sectag, av_mjd, 
     .       av_init(3,max_site), adj_xyz(3), adj_neu(3), xyz_coord(3), 
     .       loc_coord(3) 

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
      integer*4 i, j, k, na, date(5), pj, pk, doy, pnt, num_total,
     .          num_fixd, ks, av_num(max_site), av_date(5)

      integer*4 trimlen

      logical kbit  ! Tests if bit is set


****  Get the time tag of this output 
      kf_curr_mjd = RT_MJD_obs(sblk)
      call mjd_to_ymdhms(kf_curr_mjd,date,sectag)
      call ymd_to_doy( date, doy)
      sec_in_day = date(4)*3600.d0 + date(5)*60.d0 + sectag


***   Output results
      do i = 1, num_site

*        Count the number of bias flags and how many are fixed
*        at the this sites
         num_total = 0
         num_fixd = 0
         do j = 1, num_ambs
            if( bf_ents(1,j).eq.i .and. bf_ents(4,j).eq.ep ) then
*              Ambiquity points to this station
               num_total = num_total + 1
               if( kbit(bf_ents(5,j),2) .or.
     .             wls_ref(2,j).le.0    ) then
                  num_fixd = num_fixd + 1
               end if
            endif
         end do

*        Compute the position of site
         do j = 1,3
            if( site_parn(j,i).gt.0 ) then
               pj = site_parn(j,i)
               xyz_coord(j) = curr_site_xyz(j,i) + sol_vecp(pj)
               adj_xyz(j) = xyz_coord(j) - (site_int(j,i) + 
     .             site_vel(j,i)*(RT_MJD_obs(sblk)-site_ep(i))/365.25d0)
            else
               xyz_coord(j) = curr_site_xyz(j,i)
               adj_xyz(j) = 0.0d0
            end if
         end do
                    
*        Compute geodetic position so that atmospherc delay can
*        be determined
         call XYZ_to_GEOD( rot_mat, xyz_coord(1), geod_coord )

         lat = pi/2.d0-geod_coord(1)
         long = geod_coord(2)
         ht = geod_coord(3)
         if( atm_parn(i).gt.0 ) then
             na = atm_parn(i)
             atm_delay = sol_vecp(na)  ! Set from sol_vecm once OK
             atm_sig = sqrt(cov_parmp(na,na))
         else
             atm_delay = 0.0d0  ! default
             atm_sig   = 0.0d0  ! Not estimated
         end if

*        Now compute apriori delay
         rel_hum = ref_rel_hum

*        Get the time: For output we want GPS time (not time on
*        the receiver).
         call mjd_to_ymdhms(kf_curr_mjd, date, sectag)
          
         if ( atm_mtt ) then
            call met_seasonal( temp_c, press, rel_hum, tbias, 
     .            kf_curr_mjd,lat, ht ) 
                           
*           Now get the zenith delay terms and map to the elevation angle
            if( .not.use_atm_ref ) then 
                call dry_saas_zen( temp_C, lat, ht, press,  
     .                             dry_zen)
                call wet_saas_zen( temp_C, lat, ht, rel_hum,
     .                             wet_zen)
            else
                call get_atmRT( i, kf_curr_mjd, dry_zen, wet_zen, 
     .                     lat, long, ht )
            end if
         else 
*           Use GPT and GMF
            call gpt(kf_curr_mjd,lat,long,ht, press, temp_C, undu) 
            if( .not.use_atm_ref ) then 
                call dry_saas_zen( temp_C, lat, ht, press,   
     .                             dry_zen)
                call wet_saas_zen( temp_C, lat, ht, rel_hum, 
     .                             wet_zen)
            else
                call get_atmRT( i, kf_curr_mjd, dry_zen, wet_zen, 
     .                        lat, long, ht )
            end if
         endif

         atm_total = dry_zen + wet_zen + atm_delay

         fday = doy*1.d0 + date(4)/24.d0 + date(5)/1440.d0 +
     .          sectag/86400.d0 

****     Now get the position estimates
*        Pull of the covariance matrix elements
         cov_xyz(:,:) = 0.d0   ! Clear the cov_XYZ matrix
         do j = 1,3
            pj = site_parn(j,i)
            if( pj.gt.0 ) then
               do k = 1, 3
                  pk = site_parn(k,i)
                  if( pk.gt.0 ) then
                     cov_xyz(j,k) = cov_parmp(pj,pk)
                  end if
               end do
               pk = atm_parn(i)
               if( pk.gt.0 ) then
                   cov_xyz(j,4) = cov_parmp(pj,pk)
                   cov_xyz(4,j) = cov_parmp(pj,pk)
               endif
            end if
         end do
         pk = atm_parn(i)
         if( pk.gt.0 ) then
            cov_xyz(4,4) = cov_parmp(pk,pk)
         end if
*        Make the trans matrix
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

*        Get the NEU variance-covariance
         call var_comp(tran_mat, cov_xyz, cov_neu, cov_tmp, 4,4,1)
         if ( pk.gt.0 .and. pj.gt.0 ) then
            rho_ha = cov_neu(3,4)/sqrt(cov_neu(3,3)*cov_neu(4,4))
         else   ! Either position or atm not estimated so zero
            rho_ha = 0.d0
         end if

*        Save sigmas for NEU 
         do j = 1,3
            sig_neu(j) = sqrt(cov_neu(j,j))
         end do

*        Update atm delay estimate
         if( pk.gt.0 ) then
            atm_delay = sol_vecp(pk)
            atm_sig = sqrt(cov_parmp(pk,pk))
         else
            atm_delay = 0
            atm_sig = 0
         end if

****     Now see which types will be output
         if( index(posit_type,'GEOD').gt.0 .and.
     .       write_pos(i) ) then


            atm_total = dry_zen + wet_zen + atm_delay

            write(200+i, 220) date(1), doy, sec_in_day, 
     .        90.d0-geod_coord(1)*180/pi, 
     .        geod_coord(2)*180/pi, geod_coord(3),
     .        (sig_neu(j)*100.0,j=1,3), 
     .          wrms_dd_site(i), num_dd_site(i),
     .        atm_total*1000., atm_sig*1000.0, fday, ep,
     .        num_total, num_total-num_fixd
 220        format(i5,i5,f11.2,2f14.8,f10.3,3f6.1,f6.1,i3,
     .           2F9.2,1x,f15.9,1x,i6,2I4)
         endif

         if( index(posit_type,'NEU').gt.0  .and.
     .       write_pos(i)) then
*           Get the NEU difference in position 
* MOD TAH 0901112: Added use of ref_xyz instead of first site
            if( ref_xyz(1).eq.0 ) then
                do j = 1,3
                   ref_xyz(j) = site_apr(j,1)
                end do
            end if

            call XYZ_to_geod(rot_mat, ref_xyz, geod_ref)
            dNEU(1) = -earth_rad*(geod_coord(1)-geod_ref(1))
            dNEU(2) =  earth_rad*sin(geod_ref(1))*
     .                           (geod_coord(2)-geod_ref(2))
            dNEU(3) = (geod_coord(3)-geod_ref(3))
            write(200+num_site+i, 420) date, sectag, 
     .          (dNEU(j), sig_neu(j),j=1,3), 
     .          wrms_dd_site(i), num_dd_site(i),
     .          atm_delay*1000., atm_sig*1000.0, fday, ep,
     .          num_total, num_total-num_fixd, rho_ha  
 420        format(i5,4i3,1x,F7.3, 3(F15.4,1x,F9.4,1x),
     .           F8.2,1x,i3,1x, 2F9.2,1x,f15.9,1x,i6,2i4,
     .           1x,F6.3)
         end if

         if( index(posit_type,'XYZ').gt.0  .and.
     .       write_pos(i)) then
            write(200+2*num_site+i, 620) date, sectag, 
     .          (xyz_coord(j), sqrt(cov_xyz(j,j)),j=1,3), 
     .          wrms_dd_site(i), num_dd_site(i),
     .          atm_delay*1000., atm_sig*1000.0,
     .          fday, ep, num_total, num_total-num_fixd, rho_ha  
 620        format(i5,4i3,1x,F7.3, 3(F15.4,1x,F9.4,1x),
     .           F8.2,1x,i3,1x, 2F9.2,1x,f15.9,1x,i6,2i4,
     .           1x,F6.3)
         end if

         if( index(posit_type,'DHU').gt.0  .and.
     .       write_pos(i)) then
*           Rotate the dXYZ from the initial coordimates into
*           dNEU
            call rotate_geod(adj_xyz, adj_neu, 'XYZ','NEU',
     .             xyz_coord, loc_coord, rot_mat)
            write(200+3*num_site+i, 720) date, sectag, 
     .          (adj_NEU(j)*1000, sig_neu(j)*1000,j=1,3), 
     .          wrms_dd_site(i), num_dd_site(i),
     .          atm_delay*1000., atm_sig*1000.0, fday, ep,
     .          num_total, num_total-num_fixd, rho_ha 
 720        format(i5,4i3,1x,F7.3, 3(F15.1,1x,F9.1,1x),
     .           F8.2,1x,i3,1x, 2F9.2,1x,f15.9,1x,i6,2i4,
     .           1x,F6.3)
         endif

****     See of CSV format needed
         if( trimlen(csv_root).gt.0 .and.
     .       write_pos(i)) then
             call rotate_geod(adj_xyz, adj_neu, 'XYZ','NEU',
     .             xyz_coord, loc_coord, rot_mat)

             write(200+4*num_site+i,820) date,nint(sectag),
     .         ((adj_NEU(j)-sig_neu(j))*1000,
     .          (adj_NEU(j))*1000,
     .          (adj_NEU(j)+sig_neu(j))*1000,j=1,3),
     .          (atm_delay-atm_sig)*1000,atm_delay*1000,
     .          (atm_delay+atm_sig)*1000,wrms_dd_site(i),
     .          num_dd_site(i)
 820        format(I4.4,'-',I2.2,'-',I2.2,1x,I2.2,':',I2.2,':',I2.2,
     .          ',',9(F8.2,','),3(F8.2,','),F6.2,',',I2.2,',')
          end if

      end do

      return
      end

