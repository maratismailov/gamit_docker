CTITLE THEORY

      subroutine theory( trn_time, trn_sec, ns, ch, site_pos,
     .           svs_ear,
     .           range, phase, elev, az, dry_map, wet_map, opt, ep )
* MOD AZ 190305: use built Module to store VMF grids
      use vmf_type_build

      implicit none

*     Routine to compute the theoretical delay from station to satellite
*     for PRN prn at transmitt time trn_time (as given by satellite
*     clock).

      include '../includes/const_param.h'   
      include 'track_com.h'
* MOD AZ 190305: include header for VMF grid path
      include 'vmf_com.h'

* PASSED VARIABLES
* ns       -- Site number 
* ch       -- Channel (use ctop_cse to get pn)
* ep       -- Epoch number.  Used when times gets too far off due to bad data.

      integer*4 ns, ch, ep

* pn       -- PRN of satellite
* MOD AZ 190305: additioncal variables for tidal correction
* date_et,date_otl -- date arrary for tidal correction
      integer*4 pn, date_et(5), date_otl(5)

* trn_time  -- Transmitt time as given by satellite clock (MJD)
* trn_sec   -- Transmit time as seconds from start of SP3 file
* site_pos(3) -- Approximate position of the site (XYZ, m) Earth fixed
* svs_ear(3)  -- Position of the satellite in earth fixed frame,
* range(2)    -- Range to the satellite at the L1 and L2 frequencies (m)
* phase(2)    -- Phase for the satellite at L1 and L2 freqs (cycles at each freq)
* elev        -- Elevation angle to satellite (deg)
* dry_map     -- Dry mapping function
      real*8 trn_time, trn_sec, site_pos(3), svs_ear(3), range(2), 
     .       phase(2), elev, az, dry_map
* MOD AZ 190305: additional variables for IERS 2010 Tidal Changes
* fhr -- fractional hour in the day in UTC
* second_e -- seconds of signal sent time for earth tide
* jd_cent -- julian centuries of signal sent time
* fhr_sec -- fractional hour in seconds
* seconds_otl -- seconds of signal sent time for OTL in UTC
* pjd_etide -- PEP Julian date of signal sent time for earth tide
* dXYZ_etide -- earth tide coorection in XYZ
* taiutc -- difference between TAI and UTC (leap seconds)
* pjd_eph -- PEP Julian date of TT time for Sun&Moon position
* jd_etide_itoe -- timestamp used for conversion from inertial to earth
* sun -- sun position array interpolated
* moon -- moon position array interpolated
* sunf -- geocentric position of sun
* moonf -- geocentric position of moon
* suni -- inertial position of sun
* mooni -- inertial position of moon
* pjd_etide_int -- floored value of pjd_etide
* seconds_otl_int -- floored value of seconds_otl
* blq_unit -- open unit tag for BLQ files
* dNEU_otl -- OTL correction in NEU
* dXYZ_otl -- OTL correction in XYZ
* jd_utc -- julian date of signal sent time for OTL in UTC
      real*8 fhr, seconds_e, jd_cent, fhr_sec, seconds_otl,
     .       pjd_etide, dXYZ_etide(3), taiutc, pjd_eph, jd_etide_itoe
      real*8 sun(6), moon(6),sunf(3),moonf(3),suni(3),mooni(3)
      integer*4 pjd_etide_int, seconds_otl_int, blq_unit
      real*8 dNEU_otl(3), dXYZ_otl(3),jd_utc
* MOD AZ 190305: additional variable for Relativistic Effect - Propagation
* sitrad -- geocentric range of the station
* satrad -- geocentric range of the satellite
* relcon -- relativistic constant for conversion
* reldel -- relativistic effect correction
* range_sit_sat -- geocentric range from station to satellite
* dot -- dot product of vectors
* max_rel -- the thereotical maximum value for relativistic effect
      real*8 sitrad, satrad, relcon, reldel, range_sit_sat,
     .       dot, max_rel
* opt   -- Option passed: S for short calculation (no atm etc) or
*                         F for full calculation, N means no debug

      character*(*) opt

* LOCAL VARIABLES
* i, j  -- Loop counters
       integer*4 i, j
* send_mjd  -- Actual transmission time corrected for satellite clock
*         errors
* atm_range  -- Atmospheric delay correction to the range

      real*8 send_mjd, atm_range
                                    
* dNEU_L1(3), dNEU_L2(3) -- Changes site coordinates due to antenna
*     height and phase center offset (local)
* dXYZ_L1(3), dXYZ_L2(3) -- Changes site coordinates due to antenna
*     height and phase center offset (Cartesean)
* site_L1e(3), site_L2e(3) -- Carttesan coordinates of phase centers of
*     L1 and L2.
* loc_coord(3) -- Co-lat, long and height of site (rad,rad,m)
* rot_matrix(3,3) -- Rotation matrix from XYZ to NEU

* dry_zen,  -- Hydrostatic zenith delay (m) and mapping function
* wet_zen, wet_map -- Wet zenith delay (m) and mapping function

* site_L1i(3), site_L2i(3)    -- Site coordinates in non-rotating frame
* svs_inert(3)     -- Satellite coordinates in non-rotating frame
* site_t, svs_t    -- Times at site and satellite for rotating into
*                     non-rotating time
* send_sec         -- Send time in seconds
* dnadir           -- Nadir angle from satellite to ground (deg)
* zend        -- Zenith distance (rads) 
* undu  -- Geoid height from GPT model.

      real*8 dNEU_L1(3), dNEU_L2(3), dXYZ_L1(3), dXYZ_L2(3),
     .       site_L1e(3), site_L2e(3), loc_coord(3), rot_matrix(3,3),
     .       dry_zen,  wet_zen, wet_map, send_sec, dnadir, zend,
     .       undu

      real*8 site_L1i(3), site_L2i(3), svs_inert(3), site_t, svs_t

      real*8 svs_phs(3) ! Position of satellite phase center in inertial frame

      real*8 dXYZ_tide(3), jd_tide

* dsitL1, dsitL2 -- Site antanna phase center model (m)
* dsvsL1, dsvsL2 -- Satellite antenna phase center model (m)
* IntAntMod      -- Function to interpolate phase center table

      real*8 dsitL1, dsitL2, dsvsL1, dsvsL2, IntAntMod 

      real*8 dcbap(2)  ! Applied dcb (m) for debug

* lat, long, ht  -- Latitude and height of site (used for atmospheric delay
*     calculations)
* temp_c, press, rel_hum -- Temperature, press, and rel_hum at site
*     based on seasonal model (currently)
* tbias    -- Bias in surface temperature due to inversions (not used)

      real*8 lat, long, ht, temp_C, press, rel_hum, tbias 

      real*8 tec_los   ! LOS TEC value from IONEX file

* ob(2)  -- Set true if ZD/Nadir angle out of bounds for table (Values
*     for site and satellite)
      logical ob(2)

      integer*4 PtoL   ! Generate the satellite list number from the PRN
     .,         lv     ! Satellite list number 


***** Make sure it is OK

      pn = ctop_cse(ch,ns,ep)
      lv = PtoL(pn)
* MOD AZ 190305: In some Linux system, when lv = -1, the program will stop
*     from an error showing that "index cannot be -1", so introduce a "if"
*     structure to avoid this.
      if ( lv.le.0 ) then
        return
      end if

      if( pn.eq.0 .or. lv.eq. 0 ) then
         write(*,120) ep, site_names(ns), ns, pn, lv 
 120     format('*** NO PRN at Epoch ',i6,' Site ',a,1x,i3,' PL ',2I4)
         return
      end if

****  Get the satellite clock correction
      send_mjd = trn_time - svs_clk(lv)/86400.d0
      send_sec = trn_sec - svs_clk(lv)
* MOD AZ 190305: theoretical max value of 19 mm
      max_rel = 19.d0/1000.d0

      call eph_to_xyz( send_mjd, send_sec, lv, 'E', ep )

      if( opt(1:1).eq.'S' ) then
          range(1) = sqrt( (site_pos(1)-svs_xyz(1,lv))**2 +
     .                 (site_pos(2)-svs_xyz(2,lv))**2 +
     .                 (site_pos(3)-svs_xyz(3,lv))**2 )

c          if( ep.gt.debug_start .and. ep.le.debug_end .and.
c     .        opt(2:2).ne.'N'  ) then
c              write(*,140) ep, send_mjd, send_sec, pn, 
c     .              site_pos, svs_xyz(:,lv)
c 140          format('SRange EP MJD Sec PRN ',I6,F14.6,F14.6,1x,i3,
c     .               ' Site SVS Pos ',3F14.3,1x,3F16.3)
c          end if
          RETURN
      end if
 
****  Save the Earth-fixed coordinates of the satellite
      do j = 1, 3
         svs_ear(j) = svs_xyz(j,lv)
      end do

****  Get the elevation angle to the satellite
      call get_elev(site_pos, svs_ear, elev, az)
* MOD TAH 100517: Check elevation angle to make sure not less than zero
      if( elev.lt. 0.1  ) then
          if( debug_start.gt.0 ) 
     .    write(*,150) ep, site_names(ns), pn, send_mjd, elev
 150      format('WARNING: Small elvation angle Ep ',i6,' Site ',a,
     .           ' PRN ',I3,' MJD ',F15.8,' Elev ',F6.2,' dg')
          elev = 0.1d0
      end if
      zend = pi/2 - elev*pi/180.d0
      call get_nadir(site_pos, svs_ear, dnadir)

****  Compute the antenna offset corrections for the site
      do i =1, 3
         dNEU_L1(i)= site_offarp(i,ns) + sit_L12(i,1,ns)
         dNEU_L2(i)= site_offarp(i,ns) + sit_L12(i,2,ns)
      enddo
      
      call rotate_geod(dNEU_L1,dXYZ_L1, 'NEU', 'XYZ',
     .              site_pos, loc_coord, rot_matrix) 


*     Save the latitude and height for computing the atmospheric
*     delay
      lat = (pi/2.d0-loc_coord(1))
      long = loc_coord(2)
      ht  = loc_coord(3)

      call rotate_geod(dNEU_L2,dXYZ_L2, 'NEU', 'XYZ',
     .              site_pos, loc_coord, rot_matrix)


****  Now for L1 and L2 compute the phase center coordinates of the
*     antenna
      do i = 1,3
          site_L1e(i) = site_pos(i)+dXYZ_L1(i)
          site_L2e(i) = site_pos(i)+dXYZ_L2(i)
      end do

* MOD AZ 190305: Start of IERS2010 Solid Earth Tide calculation
* MOD TAH 200519: Allow option for IERS2010 model or older analytic
*     model.
      if ( etide_2010 ) then
      
        jd_tide = send_mjd + 2400000.5d0
        pjd_etide = jd_tide + 0.5d0
        jd_cent = ( jd_tide - DJ2000 )/36525.d0
        jd_cent = jd_cent + (19.d0 + 32.184d0)/(3600.d0*
     .            24.d0*36525.d0) 

        call JD_to_YMDHMS(jd_tide, date_et, seconds_e)
        
        fhr_sec = nint(seconds_e) + date_et(4)*3600.d0 
     .            + date_et(5)*60.d0
        pjd_etide_int = floor(pjd_etide)  

*        write(*,*) pjd_etide_int      
        fhr = (fhr_sec - taiutc(pjd_etide_int) + 19.d0)/3600.d0
        pjd_eph = jd_tide + 0.5d0 + ( 32.184d0 + 19.d0 )/86400.d0

        call ephtrp(pjd_eph,3,0,sun)     
        call ephtrp(pjd_eph,10,0,moon)
        
        do i = 1,3
          suni(i) = -sun(i)*1000.d0
          mooni(i) = moon(i)*1000.d0
        enddo
     
        jd_etide_itoe = send_mjd
        call earth_to_inert( jd_etide_itoe, suni, sunf, 'I','E')      
        call earth_to_inert( jd_etide_itoe, mooni, moonf, 'I','E')

        call DEHANTTIDEINEL(site_pos,jd_cent,fhr,sunf,moonf,dXYZ_etide)

      else
* MOD TAH 200519: Allow old model to be used as well.
    
*       MOD MAF 210624: Added definition of jd_tide, as above when etide_2010 is true
        jd_tide = send_mjd + 2400000.5d0
        call earth_tide( jd_tide, site_pos, dXYZ_etide )

      end if

* MOD AZ 190305: End of IERS2010 Solid Earth Tide Calculation
* MOD 090116: Added option to turn off tides
      if( noetide ) then
         do i = 1,3
            dXYZ_etide(i) = 0.0d0
         end do
      end if
      do i = 1,3
          site_L1e(i) = site_L1e(i)+dXYZ_etide(i)
          site_L2e(i) = site_L2e(i)+dXYZ_etide(i)
      end do         
* MOD AZ 190305: Start of IERS2010 OTL Calculation
c      if ( pin_holder.ge.1 ) then
* MOD TAH 200225: Only call ocean tide if requested
      if( use_blq ) then 

         blq_unit = 170 + ns

*        MOD MAF 210624: Added definition of jd_tide, pjd_etide and pjd_etide_int, as above when etide_2010 is true
         jd_tide = send_mjd + 2400000.5d0
         pjd_etide = jd_tide + 0.5d0
         pjd_etide_int = floor(pjd_etide)  
         jd_utc = jd_tide - ( taiutc(pjd_etide_int) - 19.d0 )/86400.d0

         call JD_to_YMDHMS(jd_utc,date_otl,seconds_otl)

         seconds_otl_int=nint(seconds_otl)
     
!         write(*,*) 'Okay till now before call hardisp'

!         write(*,*) date_otl(1),' ', date_otl(2),' ',date_otl(3),
!        .         ' ',date_otl(4),' ',date_otl(5),' ',seconds_otl_int,
!        .         ' ',blq_unit,' ',dNEU_otl
         call hardisp(date_otl(1),date_otl(2),date_otl(3),
     .            date_otl(4),date_otl(5),seconds_otl_int,
     .            blq_unit,dNEU_otl)
         rewind blq_unit
!         write(*,*) 'Okay till now before call rotate_geod'
         call rotate_geod(dNEU_otl, dXYZ_otl, 'NEU', 'XYZ',
     .                site_pos,loc_coord,rot_matrix)
         
         do i = 1,3
             site_L1e(i) = site_L1e(i)+dXYZ_otl(i)
             site_L2e(i) = site_L2e(i)+dXYZ_otl(i)
         end do     
      endif    

c      write(*,*) 'IERS2010 OTL Time Tag: ',jd_utc,
c     .           ' Site: ',site_names(ns),
c     .           ' dX(mm): ',dXYZ_otl(1)*1000.d0,
c     .           ' dY(mm): ',dXYZ_otl(2)*1000.d0,
c     .           ' dZ(mm): ',dXYZ_otl(3)*1000.d0,
c     .           ' dN(mm): ',dNEU_otl(1)*1000.d0,
c     .           ' dE(mm): ',dNEU_otl(2)*1000.d0,
c     .           ' dU(mm): ',dNEU_otl(3)*1000.d0
c      end if

* MOD AZ 190305: End of IERS2010 OTL Calculation
****  Now compute the atmospheric delay correction.  For the moment
*     use the seasonal model.
      if ( atm_mtt ) then
         call met_seasonal( temp_c, press, rel_hum, tbias, trn_time,
     .                   lat, ht ) 
                                  
*        Now get the zenith delay terms and map to the elevation angle
         if( .not.use_atm_ref ) then 
             call dry_saas_zen( temp_C, lat, ht, press,   dry_zen)
             call wet_saas_zen( temp_C, lat, ht, rel_hum, wet_zen)
         else
             call get_atm( ns, trn_time, dry_zen, wet_zen, 
     .                     lat, long, ht )
         end if
 
         call dry_mtt_map( temp_C, lat, ht, elev, dry_map )
         call wet_mtt_map( temp_C, lat, ht, elev, wet_map )
      elseif ( atm_gmf ) then
*        Use GPT and GMF
         call gpt(trn_time,lat,long,ht, press, temp_C, undu) 
         if( .not.use_atm_ref ) then 
             call dry_saas_zen( temp_C, lat, ht, press,   dry_zen)
             rel_hum = ref_rel_hum
             call wet_saas_zen( temp_C, lat, ht, rel_hum, wet_zen)
         else
             call get_atm( ns, trn_time, dry_zen, wet_zen, 
     .                     lat, long, ht )
         end if
         call gmf(trn_time,lat,long,ht, zend, dry_map, wet_map)

      elseif ( atm_vmf ) then
* MOD AZ 190305: Use VMF
         if( .not.use_atm_ref ) then
             call vmf1_grid( VMF1_grid_file, trn_time, lat,
     .                      long, ht, zend, dry_map, wet_map,
     .                      dry_zen, wet_zen )
        else
             call get_atm( ns, trn_time, dry_zen, wet_zen,
     .                     lat, long, ht )
        end if
      else
             write(*,155)
 155         format('ERROR: ATM Model wrongly read',
     .               'ATM_MODELC command')
             stop 'TRACK: Commands out of order'
      endif

*     Now compute the atmospheric delay

      atm_range = (dry_zen+atm_offset(ns))*dry_map +
     .            wet_zen*wet_map
       

****  Now compute the approximate range  at L1 to get times difference
*     for inertial rotation
      range(1) = sqrt( (site_L1e(1)-svs_ear(1))**2 +
     .                 (site_L1e(2)-svs_ear(2))**2 +
     .                 (site_L1e(3)-svs_ear(3))**2 ) 
* MOD AZ 190305: Preparation of Relativistic Effect Calculation
      range_sit_sat = range(1)

      svs_t  = trn_time
      site_t = trn_time + (range(1)/vel_light)/86400.d0  

      call earth_to_inert( site_t, site_L1e, site_L1i, 'E','I')
      call earth_to_inert( site_t, site_L2e, site_L2i, 'E','I')      
      call earth_to_inert( svs_t, svs_ear, svs_inert,  'E','I')

****  Correct the satellite for the phase center position
      call svs_cm_to_phs( svs_t, lv, svs_inert, svs_phs )


      range(1) = sqrt( (site_L1i(1)-svs_phs(1))**2 +
     .                 (site_L1i(2)-svs_phs(2))**2 +
     .                 (site_L1i(3)-svs_phs(3))**2 ) + atm_range
      range(2) = sqrt( (site_L2i(1)-svs_phs(1))**2 +
     .                 (site_L2i(2)-svs_phs(2))**2 +
     .                 (site_L2i(3)-svs_phs(3))**2 ) + atm_range
      if( ep.ge.debug_start .and. ep.le.debug_end ) then
          write(*,160) ep, send_mjd, send_sec, lv, 
     .          site_L1i, svs_phs(:)
 160      format('SRFull EP MJD Sec PRN ',I6,F14.6,F14.6,1x,i3,
     .           ' Site L1 SVS Pos ',3F14.3,1x,3F16.3)
      end if
* MOD AZ 190305: Start of Relativity Correction - Propagation

      sitrad = dsqrt(dot(site_L1e,site_L1e))

      satrad = dsqrt(dot(svs_ear,svs_ear))

      relcon = (2.d0*GM_earth/vel_light**2)

      reldel = relcon
     .        *log((sitrad+satrad+range_sit_sat)
     .        /(sitrad+satrad-range_sit_sat))
* MOD AZ 190305: If larger than the theoretical max,
*               mark as bad
      if ( reldel.ge.max_rel ) then

c        write(*,*) '**WARNING** Bad Satellite: PRN ',pn

        reldel = 0.d0

        call sbit(data_flag_cse(ch,ns,ep),2,1)

c        write(*,*) '**WARNING** Not Using PRN ',pn
      
      end if
    
c      if ( pin_holder.ge.1 ) then
      
c      write(*,*) 'Time Tag: ', jd_tide,
c     .           'PRN: ', pn,
c     .           'Rel(mm): ', reldel*1000.d0
c     
c      end if

* MOD AZ 190305: Endof Relativity Correction - Propagation

* MOD TAH 061229: Version 1.14: Compute the phase center model
*     contributions
      dsitL1 = IntAntMod(num_sit_dph(1,ns), sit_dzn(1,1,ns),
     .                   sit_dphs(1,1,1,ns), 90-elev, az, 
     .                   max_zen, max_az, ob(1))
      dsitL2 = IntAntMod(num_sit_dph(1,ns), sit_dzn(1,1,ns),
     .                   sit_dphs(1,1,2,ns), 90-elev, az, 
     .                   max_zen, max_az, ob(1))

*     Compute the satellite nadir angle
      
      dsvsL1 = IntAntMod(num_svs_dph(1,lv), svs_dna(1,1,lv),
     .                   svs_dphs(1,1,1,lv), dnadir, az,
     .                   max_dna, 1, ob(2))
      dsvsL2 = IntAntMod(num_svs_dph(1,lv), svs_dna(1,1,lv),
     .                   svs_dphs(1,1,2,lv), dnadir, az,
     .                   max_dna, 1, ob(2))

      if( ob(1) ) then
cd        write(*,160) ns, pn, ep, elev, dsitL1, dsitL2
cd 160     format('PVC SITE Out-of-range: Site ',i2,' PRN ',I2,
cd     .          ' EP ',i6,' Elev  ',F6.2,' PCV ',2F9.4,' m')
      endif

      if( ob(2) ) then
cd        write(*,165) ns, pn, ep, dnadir, dsvsL1, dsvsL2
cd 165     format('PVC SVS  Out-of-range: Site ',i2,' PRN ',I2,
cd     .          ' EP ',i6,' Nadir ',F6.2,' PCV ',2F9.4,' m')
      endif
* MOD AZ 190305: Apply Relativistic Effect Correction to ranges
      range(1) = range(1) + dsitL1 + dsvsL1 + reldel
      range(2) = range(2) + dsitL2 + dsvsL2 + reldel

****  Convert range to phase values and then apply DCB coorections
*     to ranges: Corrections add to data and therefore removed
*     from theory module here 
      phase(1) = range(1)*fR1/vel_light
      phase(2) = range(2)*fR2/vel_light

*     Get DCB corrections
      call get_dcb(i, pn, dcbap) 
      range(1) = range(1) - dcbap(1)
      range(2) = range(2) - dcbap(2)

****  If IONEX file read; Add correction to model
      if( use_ionex ) then
         call interp_ionex( ep, ns, ch, site_pos, az, elev,  
     .                      send_mjd, tec_los )
         range(1) = range(1) + l1tecu*tec_los
         range(2) = range(2) + l2tecu*tec_los
         phase(1) = phase(1) - l1tecu*tec_los*fR1/vel_light
         phase(2) = phase(2) - l2tecu*tec_los*fR2/vel_light
      end if


      if( ep.ge.debug_start .and. ep.le.debug_end ) then
         write(*,210) site_names(ns), site_pos,
     .        (dNEU_L1(i), i=1,3), (dNEU_L2(i),i=1,3)
 210     format('Site ',a4,' Position ',3F14.4,/,
     .          'L1 NEU phase center ',3F8.4,
     .         ' L2 NEU phase center ',3F8.4)
         write(*,220) 
     .        (dXYZ_L1(i), i=1,3), (dXYZ_L2(i),i=1,3)
 220     format('L1 XYZ phase center ',3F8.4,
     .         ' L2 XYZ phase center ',3F8.4)
         write(*,230) site_names(ns), range(1), site_L1i,
     .           site_names(ns), range(2), site_L2i
 230     format('Inertial coords ',a4,' L1 ', F18.4, 3F14.4,/,
     .          'Inertial coords ',a4,' L1 ', F18.4, 3F14.4)

         write(*,240) (dXYZ_etide(i),i=1,3)
 240     format('Earthtide dXYZ ',3F8.4)

         write(*,250) pn, rcv_type(ns),dcbap
 250     format('APPLIED DCB for PRN ',i3.2,' Type ',a1,
     .          ' P1/P2 ',2F8.3,' m')

         write(*,260) pn, (svs_phs(i)-svs_inert(i),i=1,3)
 260     format('SVS phase center PRN ',i3,' dXYZ ',3F8.4)

         write(*,280) ns, trn_time, temp_c, press, rel_hum, lat, ht,
     .             dry_zen, wet_zen, dry_map, wet_map, atm_range,
     .             elev
 280     format('ATM ',i2,F14.6,' TPH ',3F9.2,' LAT,HT ',2F10.3,
     .       ' ZEN ',2F9.4, ' MAP ',2F9.4, ' TOT ',2f9.4)

         write(*,300) ns, pn, elev, dnadir, az, 
     .           dsitL1, dsitL2, dsvsL1, dsvsL2
 300     format('PVC: Site ',i2,' PRN ',i3,' Elev, Nadir, AZ ',
     .          3F8.2,' Site L1/L2 ', 2F9.4,' SVS L1/2 ',2F9.4)
         if( use_ionex ) then 
            write(*,340) ns, pn, elev, tec_los, 
     .                   l1tecu*tec_los, l2tecu*tec_los
 340        format('IONEX: Site ',i3,' PRN ',i3,' Elev ',F8.2,
     .             ' TEC_LOS ',F6.2,' L1/L2 dR ',2F8.4,' m') 
         endif 
         write(*,400) ep, site_names(ns),pn, range, phase
 400     format('Final Range Ep ',i6,1x,a,1x,' PRN ',I3,1x,
     .          'Range ',2F16.4,' m Phase ',2F16.4, ' cyc')

      end if
 
****  Thats all
      return
      end

      subroutine earth_tide( epoch, site_pos, dXYZ_tide )

      implicit none
 
*     Routine to compute the solid Earth tide.  Based on
*     Formualtion in DSR thesis with extension to the number of
*     coefficients.
 
 
      include '../includes/const_param.h'
 
* PASSED Variables
 
* epoch         - Epoch (JD or MJD)
* site_pos(3) - Site position (XYZ)
* dXYZ_tide(3) - Tidal contribution to site position (m)
 
 
      real*8 epoch, site_pos(3), dXYZ_tide(3)
 
* LOCAL VARIABLES
 
*    dn, de, dh   - Tide displacements in North, East and Height
*    dNEU_tide(3)  - Tides displacements as a vector
*    latr, longr - Latitude and longitude
*    U0, Udlat, Udlong  - Potential (in mm) at point and displaced
*                 - in lat and long
*    dUdlat, dUdLong    - Deriviatives in lat and long
*    loveh, lovel - Love numbers
*    jd           - Julian date
*    dh_K1   - K1 FCN contribution to tide
*    gast    - Greenwich sidereal time
 
      real*8 loc_coord(3), rot_matrix(3,3),
     .    dn, de, dh, dNEU_tide(3), latr, longr, U0, Udlat, Udlong,
     .    dUdlat, dUdLong, loveh, lovel, jd, dh_K1, gast
 
      DATA LOVEH/.6090d0/, LOVEL/.0852d0/
 
****  See type of date passed
      if( epoch.gt.2000000.d0 ) then
          jd = epoch
      else
          jd = epoch + 2400000.5d0 - 94554.d0 + 142350.d0
      end if 
****  Get the lat and long of the site
      call XYZ_to_GEOD( rot_matrix, site_pos, loc_coord )
 
      latr = pi/2 - loc_coord(1)
      longr = loc_coord(2)
      call tide_u(jd,  latr, longr, U0 )
 
*     Now do the derivatives
      call tide_u(jd,  latr+1.d-6, longr, Udlat)
      dudlat = (Udlat-U0)/1.d-6
 
      call tide_u(jd,  latr, longr+1.d-6, Udlong)
      dudlong = (Udlong-U0)/1.d-6
 
*     Now compute the displacements
      dh = loveh*U0
      dn = lovel*dUdlat
* MOD TAH 950626; Added cos(latr) scaling to de
      de = lovel*dUdlong/cos(latr)

*     Compute the K1 correction
      call gst_jd( jd, gast )
      dh_K1 = -25.3d0*dsin(latr)*dcos(latr)*
     .       dsin(gast+longr)

      dh = dh + dh_K1

*     Save and convert values from mm to meters
      dNEU_tide(1) = dn/1000.d0
      dNEU_tide(2) = de/1000.d0
      dNEU_tide(3) = dh/1000.d0
 
*     Now rotate to XYZ didplacements
      call rotate_geod( dNEU_tide, dXYZ_tide, 'NEU', 'XYZ',
     .        site_pos, loc_coord, rot_matrix)

C      write(*,300) (epoch-int(epoch))*86400.d0, dNEU_tide,
C     .              dXYZ_tide, jd
C 300  format('TIDE:',f16.6,6F12.5,1x,F15.6 )
 
****  Thats all
      end
 
CTITLE TIDE_U
 
      subroutine tide_u(jd, latr, longr, U)

      implicit none
 
*     Routine to compute the tidal potential in mm.
 
      include '../includes/const_param.h'
 
*         num_lp, num_di, num_se   - Number of terms in
*                    - Long period, diurnal and semidiurnal
*                    - series
 
      integer*4 num_lp, num_di, num_se
 
      parameter ( num_lp = 13 )
      parameter ( num_di = 18 )
      parameter ( num_se = 12 )
 
 
*      latr, longr   - Lat and Long in rads
*      jd            - JD for determination
*      U             - Potenital (mm)
*      lm, ls, w, gst, tc   - Long of moon, of sun,
*                    - Argument of lunar perigee,
*                    - Greenwich sidreal time, and
*                    - time in centuries since 1900.
*      lha, lst      - Local hour angle of moon and sun
 
*      fund_arg(6)   - Browns fundamental arguments
*      dood_arg(6)   - Doodson's arguments in following
*                    - order:
*                    - tau - Time angle in lunar days
*                    - s   - Mean longitude of Moon
*                    - h   - Mean longitude of Sun
*                    - p   - Long of Moon's perigee
*                    - N'  - Negative of long of Moon's Node
*                    - p1  - Longitude of Sun's Perigee.
*      arg           - Argument of tide (rads)
*      A_lp, A_di, A_se   - Long period, diurnal and semi-diurnal
*                    - ampltituds (not quite potential since stills
*                    - need to be multiplied 268.8 mm.
 
      real*8 latr, longr, jd, U, lm, ls, w, gst, tc, lha, lst,
     .    fund_arg(6), dood_arg(6), arg, A_lp, A_di, A_se
 
*         lp_tides(7, num_lp)   - Long period Doodson arguments
*                    - and amplitude by 1d-5
*         di_tides(7, num_di)   - Diurnal args and amp
*         se_tides(7, num_se)   - Semidiurnal args and amp
*         i,j        - Loop counters
 
 
      integer*4 lp_tides(7, num_lp), di_tides(7, num_di),
     .    se_tides(7, num_se), i,j
 
      data lp_tides /  0,  0,  0,  0,  0,  0, 73807,
     .                 0,  0,  0,  0,  1,  0, -6556,
     .                 0,  0,  1,  0,  0, -1,  1156,
     .                 0,  0,  2,  0,  0,  0,  7266,
     .                 0,  1, -2,  1,  0,  0,  1579,
     .                 0,  1,  0, -1,  0,  0,  8255,
     .                 0,  2, -2,  0,  0,  0,  1366,
     .                 0,  2,  0, -2,  0,  0,   676,
     .                 0,  2,  0,  0,  0,  0, 15645,
     .                 0,  2,  0,  0,  1,  0,  6482,
     .                 0,  2,  0,  0,  2,  0,   605,
     .                 0,  3,  0, -1,  0,  0,  2995,
     .                 0,  3,  0, -1,  1,  0,  1241  /
 
      data di_tides /  1, -3,  0,  2,  0,  0,   954,
     .                 1, -3,  2,  0,  0,  0,  1151,
     .                 1, -2,  0,  1, -1,  0,  1359,
     .                 1, -2,  0,  1,  0,  0,  7214,
     .                 1, -2,  2, -1,  0,  0,  1370,
     .                 1, -1,  0,  0, -1,  0,  7105,
     .                 1, -1,  0,  0,  0,  0, 37690,
     .                 1,  0,  0, -1,  0,  0, -1066,
     .                 1,  0,  0,  1,  0,  0, -2963,
     .                 1,  1, -3,  0,  0,  1,  1028,
     .                 1,  1, -2,  0,  0,  0, 17546,
     .                 1,  1,  0,  0, -1,  0,  1050,
     .                 1,  1,  0,  0,  0,  0,-53009,
     .                 1,  1,  0,  0,  1,  0, -7186,
     .                 1,  1,  2,  0,  0,  0,  -755,
     .                 1,  2,  0, -1,  0,  0, -2963,
     .                 1,  3,  0,  0,  0,  0, -1623,
     .                 1,  3,  0,  0,  1,  0, -1039 /
 
      data se_tides /  2, -3,  2,  1,  0,  0,   669,
     .                 2, -2,  0,  2,  0,  0,  2298,
     .                 2, -2,  2,  0,  0,  0,  2774,
     .                 2, -1,  0,  1, -1,  0,  -649,
     .                 2, -1,  0,  1,  0,  0, 17380,
     .                 2, -1,  2, -1,  0,  0,  3301,
     .                 2,  0,  0,  0, -1,  0, -3390,
     .                 2,  0,  0,  0,  0,  0, 90805,
     .                 2,  1,  0, -1,  0,  0, -2567,
     .                 2,  2, -2,  0,  0,  0, 42248,
     .                 2,  2,  0,  0,  0,  0, 11495,
     .                 2,  2,  0,  0,  1,  0,  3424 /
 
      call gst_jd( jd, gst )
 
      tc = (jd - 2415020.5d0)/ 36525.d0
 
      lm =  4.719967d0 + 8399.709d0*tc
      ls =  4.881628d0 + 628.3319d0*tc
      w  =  5.835152d0 + 71.01803d0*tc
 
      lha = gst - lm + longr
      lst = gst - ls + longr
 
***** Get the fundamental arguments at this time and then
*     convert to Doodson argument
      call fund_angles( jd, fund_arg )
 
*     Now computed Doodson's angles: NOTE:
*     fund_arg(6) is gst+pi (so we don't need to add pi below)
      dood_arg(2) = fund_arg(3) + fund_arg(5)
      dood_arg(1) = fund_arg(6) - dood_arg(2) + longr
      dood_arg(3) = dood_arg(2) - fund_arg(4)
      dood_arg(4) = dood_arg(2) - fund_arg(1)
      dood_arg(5) = -fund_arg(5)
      dood_arg(6) = dood_arg(2) - fund_arg(4) - fund_arg(2)
 
****  Now compute the potential
*     Start with Long Period
      A_lp = 0
      do i = 1, num_lp
         arg = 0
         do j = 1,6
            arg = arg + lp_tides(j,i)*dood_arg(j)
         end do
         A_lp = A_lp + lp_tides(7,i)*1.d-5*cos(arg)
      end do
 
*     Do the diurnal terms
      A_di = 0.0d0
      do i = 1, num_di
         arg = -pi/2
         do j = 1,6
            arg = arg + di_tides(j,i)*dood_arg(j)
         end do
         A_di = A_di + di_tides(7,i)*1.d-5*cos(arg)
      end do
 
*     Do the Semidiurnal tides.
      A_se = 0.0d0
      do i = 1, num_se
         arg =  0
         do j = 1,6
            arg = arg + se_tides(j,i)*dood_arg(j)
         end do
         A_se = A_se + se_tides(7,i)*1.d-5*cos(arg)
      end do
 
****  Now we can add up all the peices
      U = 268.8d0*( cos(latr)**2*A_se +
     .              sin(2*latr) *A_di -
     .             (1.5d0*sin(latr)**2-0.5d0)*A_lp )
 
****  Thats all
      return
      end

CTITLE SVS_CM_TO_PHS
 
      subroutine svs_cm_to_phs( mjd, pn, svs_non, svs_phs )

      implicit none
 
*     This routine computes the phase center position of satellite
*     given the position of the center of mass
 
      include '../includes/const_param.h'
      include 'track_com.h' 
 
*   mjd  - Modified Julian date (not used yet but may be needed later)
*   svs_non(3)  - Position of CM of satellite
*   svs_phs(3)  - Position of phase center
 
      real*8 mjd, svs_non(3), svs_phs(3)
 
*   pn      - PRN number of satellite
 
      integer*4 pn
 
* LOCAL VARIABLES
 
*   ehat(3) - Unit vector in direction of Earth
*   shat(3) - Unit vector in direction of sun
*   yhat(3) - Satellite y axis vector
*   mag     - Magnitude of svs_non vector
*   tdt     - Dynamical time for Sun's position
*   sun_dist - Distance to Sun (AU)
*   sun_long, sun_lat - Longitude and Latitude off sun
*   sv_cm(3)  - Position of Satellite phase center in local
*              system.
 
      real*8 ehat(3), shat(3), yhat(3), mag, tdt, sun_dist, sun_long,
     .       sun_lat, sv_cm(3), xhat(3)

c     real*8 old_tdt
 
*   i       - Loop counter
 
 
      integer*4 i

      logical warn_out  ! Set true once warning on sun angle output

      data warn_out / .false. /

c     save old_tdt
 
****  Compute unit vector in direction on Earth
      mag = sqrt(svs_non(1)**2+svs_non(2)**2+svs_non(3)**2)
      do i = 1,3
          ehat(i) = -svs_non(i)/mag
      end do

****  Get TDT and compute vector to sun
      tdt =  mjd - 94554.d0 + 142350.d0 +(19.d0+32.184d0)/86400.d0

      call SUN20(tdt ,shat ,sun_dist, sun_long, sun_lat)


****  Now get Y and X-axis
      call xprod(ehat, shat, yhat )
      mag = sqrt(yhat(1)**2+yhat(2)**2+yhat(3)**2)
      if( mag.lt.0.01d0 .and. .not.warn_out ) write(*,110) mag
 110  format('WARNING: S/C z-axis close to Sun direction. Mag ', F9.6)
      if( mag.lt.0.01d0 ) warn_out = .true.
      if( yhat(1)**2+yhat(2)**2+yhat(3)**2.lt.0 ) then
         print *,'BAD YHAT Mag ',yhat(1)**2+yhat(2)**2+yhat(3)**2, 
     .      ' yhat ', yhat, ' MJD ', mjd, ' PRN ',pn
         stop 'BAD YHAT Mag '
      end if

      do i = 1,3
          yhat(i) = yhat(i)/mag
      end do
     
      call xprod(yhat, ehat, xhat )


*     Use antenna offsets from antex file. Assume L1 and L2 are the
*     same and so just apply L1 values
      sv_cm(1) = svs_L12(1,1,pn)
      sv_cm(2) = svs_L12(2,1,pn)
      sv_cm(3) = svs_L12(3,1,pn)

 
*     Be really slack and only do radial on Block II
      do i = 1,3
C         svs_phs(i) = svs_non(i) + 1.0259d0*ehat(i)
          svs_phs(i) = sv_cm(1)*xhat(i) + sv_cm(2)*yhat(i) + 
     .                 sv_cm(3)*ehat(i) + svs_non(i) 
      end do
 
****  That all
      return
      end
 
 
C*
      SUBROUTINE SUN20(XMJD,X,R,L,B)

      implicit none

CC
CC NAME       :  SUN20
CC
CC    CALL SUN20(XMJD,X,R,L,B)
CC
CC PURPOSE    :  COMPUTATION OF POSITION OF THE SUN AT TIME XMJD
CC               (MODIFIED JULIAN DATE). THIS SR WAS WRITTEN USING
CC               SIMON NEWCOMB'S "TABLES OF THE SUN".
CC
CC               PRECISION  :        MAXIMUM        MEAN
CC                                       DIFFERENCES
CC                              L      .08"         .03"
CC                              B      .03"         .005"
CC                              R     3.E-7        .5E-7
CC
CC PARAMETERS :
CC         IN :  XMJD   : EPOCH IN MODIFIED JULIAN DATE IN    R*8
CC                        BARYCENTRIC DYNAMICAL TIME
CC                        CORRESPONDING TO EPHEMERIS TIME
CC        OUT :  X(K),K=1,2,3 : RECTANGULAR COORDINATES OF    R*8
CC                        THE SUN IN EQUATORIAL SYSTEM
CC                        J2000.0 (IN AU)
CC               R      : DISTANCE EARTH-SUN (IN AU)          R*8
CC               L , B  : ECLIPTICAL LONGITUDE, LATITUDE IN   R*8
CC                        MEAN SYSTEM OF EPOCH XMJD (FK4!)
CC
CC SR CALLED  :  PRAE, FK4FK5, SPROD, PREN20
CC
CC REMARKS    :  ---
CC
CC AUTHOR     :  G.BEUTLER, U.HUGENTOBLER
CC
CC VERSION    :  3.4  (JAN 93)
CC
CC CREATED    :  31-MAY-92                  LAST MODIFIED :  23-SEP-93
CC
CC CHANGES    :  23-SEP-93 : DECLARATION OF "PR" AND "VR" AS REAL*8
CC
CC COPYRIGHT  :  ASTRONOMICAL INSTITUTE
CC      1992      UNIVERSITY OF BERNE
CC                    SWITZERLAND
CC
C*
        REAL*8 XMJD,T,L,B,R,X(3),Y(3),PID,GD,EKL,TH
c* not used real*8 v(3),praez(3,3),vr,pr

        real*8 mag
        real*8 pi, tt, ge, gme, gv, gj, gs, gm, xl1, w1, w2 
        real*8 c, xlme, rme, w, xlv, rv, wh, xlm, rm, xlj, rj 
        real*8 xls, rs 
        real*8 xbv, xbm, xbj, xbs, xm, d, f, us, xlmo, rmo, bmo 

        integer*4 i,k

        INTEGER*4 JME(4),IME(4),S1ME(4),K1ME(4),S2ME(4),K2ME(4),
     1            JV(39),IV(39),S1V(39),K1V(39),S2V(39),K2V(39),
     2            JM(45),IM(45),S1M(45),K1M(45),S2M(45),K2M(45),
     3            JJ(21),IJ(21),S1J(21),K1J(21),S2J(21),K2J(21),
     4            JS(11),IS(11),S1S(11),K1S(11),S2S(11),K2S(11)
        INTEGER*4 JBV(22),IBV(22),SBV(22),KBV(22),
     1            JBM(3),IBM(3),SBM(3),KBM(3),
     2            JBJ(7),IBJ(7),SBJ(7),KBJ(7),
     3            JBS(2),IBS(2),SBS(2),KBS(2)
        DATA JME/4*-1/
        DATA IME/1,2,3,4/
        DATA S1ME/13,5,15,23/
        DATA K1ME/243,225,357,326/
        DATA S2ME/28,6,18,5/
        DATA K2ME/335,130,267,239/
        DATA JV/4*-1,5*-2,5*-3,5*-4,4*-5,4*-6,4*-7,5*-8,2*-9,-10/
        DATA IV/0,1,2,3,0,1,2,3,4,2,3,4,5,6,3,4,5,6,7,5,6,7,8,
     1          6,7,8,9,7,8,9,10,8,9,12,13,14,9,10,10/
        DATA S1V/75,4838,74,9,3,116,5526,2497,44,13,666,1559,
     1           1024,17,3,210,144,152,6,84,37,123,154,38,14,
     2           10,14,20,6,3,0,11,0,42,0,32,6,0,3/
        DATA K1V/2962,2991,2076,2490,1620,1488,1483,3159,3123,
     1           1760,1777,3453,3182,3150,1980,2062,1954,3438,
     2           3220,2356,2218,1953,3596,2641,2530,2300,
     3           120,2940,2790,2880,0,3220,0,2592,0,488,
     4           3510,0,180/
        DATA S2V/94,2359,69,16,4,160,6842,869,52,21,1045,1497,194,
     1           19,6,376,196,94,6,163,59,141,26,80,25,14,12,42,12,
     2           4,4,24,6,44,12,33,13,4,8/
        DATA K2V/2050,2091,3485,3300,900,584,583,2267,388,900,
     1           876,2552,495,430,900,1163,1052,2548,590,1454,
     2           1322,1054,2700,1743,1640,1350,2840,2035,1940,
     3           1660,1350,2340,2180,1697,2220,1387,2610,2560,
     4           2930/
        DATA JM/3*1,4*2,4*3,4*4,4*5,4*6,3*7,4*8,3*9,3*10,
     1          2*11,12,2*13,2*15,2*17/
        DATA IM/2,1,0,3,2,1,0,4,3,2,1,4,3,2,1,5,4,3,2,6,5,4,3,
     1          6,5,4,7,6,5,4,7,6,5,7,6,5,7,6,7,8,7,
     2          9,8,10,9/
        DATA S1M/6,273,48,41,2043,1770,28,4,129,425,8,34,500,585,
     1           9,7,85,204,3,0,20,154,101,6,49,106,3,10,52,21,4,
     2           28,62,5,19,5,17,44,6,13,45,21,0,4,26/
        DATA K1M/2180,2177,2603,3460,3439,2004,1480,2840,2942,
     1           3389,70,710,1052,3341,3250,1720,546,1008,180,
     2           0,1860,2274,963,3010,1765,2227,720,3070,3489,
     3           2152,570,2980,3460,680,1110,3380,590,1059,2320,
     4           1840,2278,3090,0,2430,1130/
        DATA S2M/8,150,28,52,2057,151,31,6,168,215,6,49,478,105,
     1           10,12,107,89,3,5,30,139,27,10,60,38,5,15,45,8,
     2           6,34,17,8,15,0,20,9,5,15,5,22,6,4,0/
        DATA K2M/1300,1277,3470,2554,2538,2950,2343,1800,2035,
     1           2490,900,3397,152,659,530,900,3246,110,1080,
     2           2170,957,1373,1880,2090,862,1329,3490,2170,2597,
     3           3100,3290,2081,2570,3370,230,0,3300,210,1430,
     4           940,1430,2200,2610,1530,0/
        DATA JJ/5*1,4*2,4*3,4*4,4*5/
        DATA IJ/-3,-2,-1,0,1,-3,-2,-1,0,-4,-3,-2,-1,-4,-3,-2,-1,
     1           -5,-4,-3,-2/
        DATA S1J/3,163,7208,2600,73,69,2731,1610,73,5,164,556,210,
     1          16,44,80,23,0,5,7,9/
        DATA K1J/1980,1986,1795,2632,2763,808,871,1095,2526,1580,
     2          1705,827,985,2590,1682,777,930,0,2590,1640,710/
        DATA S2J/5,208,7067,244,80,103,4026,1459,8,9,281,803,174,29,
     1           74,113,17,3,10,12,14/
        DATA K2J/1120,1120,895,3386,65,3505,3571,195,2630,
     2           690,812,3526,86,1700,799,3477,30,2520,1690,760,
     3     3430/
        DATA JS/4*1,4*2,2*3,4/
        DATA IS/-2,-1,0,1,-3,-2,-1,0,-2,-1,-2/
        DATA S1S/11,419,320,8,0,108,112,17,21,17,3/
        DATA K1S/1050,1006,2695,2700,0,2906,2936,2770,2890,
     1          2910,2880/
        DATA S2S/15,429,8,8,3,162,112,0,32,17,4/
        DATA K2S/110,106,3530,0,1980,2006,2031,0,2001,2010,
     1           1940/
        DATA JBV/4*-1,4*-2,5*-3,3*-4,2*-5,3*-6,-8/
        DATA IBV/0,1,2,3,1,2,3,4,2,3,4,5,6,3,5,6,6,7,
     1           5,7,8,12/
        DATA SBV/29,5,92,7,23,12,67,14,14,8,210,7,4,6,31,
     1           12,9,19,6,4,4,10/
        DATA KBV/1450,3230,937,2620,1730,1490,1230,1110,2010,
     1           1870,1518,1530,2960,2320,18,1800,270,180,2880,
     2           570,570,610/
        DATA JBM/2*2,4/
        DATA IBM/-2,0,-3/
        DATA SBM/8,8,7/
        DATA KBM/900,3460,1880/
        DATA JBJ/4*1,2,2*3/
        DATA IBJ/-2,-1,0,1,-1,-2,-1/
        DATA SBJ/7,17,16,23,166,6,18/
        DATA KBJ/1800,2730,1800,2680,2655,1710,2670/
        DATA JBS/2*1/
        DATA IBS/-1,1/
        DATA SBS/6,6/
        DATA KBS/2600,2800/
        PID=4*DATAN(1.D0)
        PI=PID
        T=(XMJD-15019.5D0)/36525.D0
        TT=T
        L=279.D0+41.D0/60+48.04D0/3600
        L=L+(129602768.13D0*T+1.089D0*T**2)/3600
        L=DMOD(L/180*PID,2*PID)
        GD=358.D0+28.D0/60+33.D0/3600
        GD=GD+(129596579.1D0*T-.54D0*T**2-.012D0*T**3)/3600
        GE=DMOD(GD/180*PID,2*PID)
        TH=(XMJD+3242.297D0)/365.25D0
        GME=DMOD((248.07D0+1494.7235D0*TH)/180*PID,2*PID)
        GV =DMOD((63.07037D0+22518.442986D0*T)/180*PID,2*PID)
        GJ =DMOD((221.64742D0+32964.466939D0*T)/180*PID,2*PID)
        GS =DMOD((193.13230D0+34777.259042D0*T)/180*PID,2*PID)
        GM =DMOD((165.94905D0+16859.069667D0*T)/180*PID,2*PID)
        XL1=6.4*SIN((231.19+20.2*TT)/180*PI)
     1      +(1.882-.016*TT)*SIN((57.24+150.27*TT)/180*PI)
     2      +.266*SIN((31.8+119.0*TT)/180*PI)
     3      +.202*SIN((315.6+893.3*TT)/180*PI)
        L=L+XL1/3600*PI/180
        GE=GE+XL1/3600*PI/180
        GV=GV-XL1/3600/180*PI
        GJ=GJ+XL1/3600/180*PI
        GS=GS+XL1/3600/180*PI
        W1=299.1+(GV-GE)/PI*180
        GD=63.07037D0+22518.442986D0*T
        W2=90.+DMOD(GD,360.D0)
        C=(6910.057-17.24*TT-.052*TT**2)*SIN(GE)
     1  +(72.338-.361*TT)*SIN(2*GE)
     2    +(1.054-.001*TT)*SIN(3*GE)+.018*SIN(4*GE)
        R=3057-15*T+COS(GE)*(-727412.D0+1814*T+5*T**2)
     1    +COS(2*GE)*(-9138+46*T)+COS(3*GE)*(-145+T)
     2    +COS(4*GE)*(-2)
        XLME=0
        RME=0
        DO 10 K=1,4
        W=-(JME(K)*GME+IME(K)*GE)
        W1=K1ME(K)/180.*PI
        W2=K2ME(K)/180.*PI
        XLME=XLME+S1ME(K)*COS(W+W1)
10      RME=RME+S2ME(K)*COS(W+W2)
        XLV=0
        RV=0
        DO 20 K=1,39
        W=-(JV(K)*GV+IV(K)*GE)
        W1=K1V(K)/1800.*PI
        W2=K2V(K)/1800.*PI
        IF(K.EQ.2)W1=299.1017/180*PI
        IF(K.EQ.7)W1=148.3133/180*PI
        IF(K.EQ.8)W1=315.9433/180*PI
        IF(K.EQ.12)W1=345.2533/180*PI
        IF(K.EQ.11)W1=177.71/180*PI
        IF(K.EQ.13)W1=318.12/180*PI
        IF(K.EQ.2)W2=209.08/180*PI
        IF(K.EQ.7)W2=58.3183/180*PI
        IF(K.EQ.11)W2=87.57/180*PI
        IF(K.EQ.12)W2=255.25/180*PI
        IF(K.EQ.16)W2=116.28/180*PI
        WH=-JV(K)*330.9017/180*PI
     1    -(IV(K)+JV(K))*GE-JV(K)*(GV+PI)
        XLV=XLV+S1V(K)*COS(WH+W1)
20      RV=RV+S2V(K)*COS(WH+W2)
        XLM=0
        RM=0
        DO 30 K=1,45
        W=(-JM(K)*GM+IM(K)*GE)
        W1=K1M(K)/1800.*PI
        W2=K2M(K)/1800.*PI
        IF(K.EQ.5)W1=343.8883/180*PI
        IF(K.EQ.6)W1=200.4017/180*PI
        IF(K.EQ.10)W1=338.88/180*PI
        IF(K.EQ.13)W1=105.18/180*PI
        IF(K.EQ.14)W1=334.05/180*PI
        IF(K.EQ.5)W2=253.8283/180*PI
        IF(K.EQ.13)W2=15.17/180*PI
        WH=-JM(K)*127.0633/180*PI+(IM(K)-JM(K))*GE+JM(K)*GM
        XLM=XLM+S1M(K)*COS(WH+W1)
30      RM=RM+S2M(K)*COS(WH+W2)
        XLJ=0
        RJ=0
        DO 40 K=1,21
        W=-(JJ(K)*GJ+IJ(K)*GE)
        W1=K1J(K)/1800.*PI
        W2=K2J(K)/1800.*PI
        IF(K.EQ.3)W1=179.5317/180*PI
        IF(K.EQ.4)W1=263.2167/180*PI
        IF(K.EQ.7)W1=87.145/180*PI
        IF(K.EQ.8)W1=109.4933/180*PI
        IF(K.EQ.12)W1=82.65/180*PI
        IF(K.EQ.3)W2=89.545/180*PI
        IF(K.EQ.7)W2=357.1083/180*PI
        IF(K.EQ.8)W2=19.4667/180*PI
        IF(K.EQ.12)W2=352.56/180*PI
        WH=-JJ(K)*88.4450/180*PI-(JJ(K)+IJ(K))*GE+JJ(K)*GJ
        XLJ=XLJ+S1J(K)*COS(WH+W1)
40      RJ=RJ+S2J(K)*COS(WH+W2)
        XLS=0
        RS=0
        DO 50 K=1,11
        W=-(JS(K)*GS+IS(K)*GE)
        W1=K1S(K)/1800.*PI
        W2=K2S(K)/1800.*PI
        IF(K.EQ.2)W1=100.58/180*PI
        IF(K.EQ.3)W1=269.46/180*PI
        WH=-JS(K)*10.2417/180*PI-(IS(K)+JS(K))*GE+JS(K)*GS
        XLS=XLS+S1S(K)*COS(WH+W1)
50      RS=RS+S2S(K)*COS(WH+W2)
        XBV=0
        DO 60 K=1,22
        W=(KBV(K)/10.-JBV(K)*330.9017)/180*PI
     1    -(IBV(K)+JBV(K))*GE-JBV(K)*(GV+PI)
60      XBV=XBV+SBV(K)*COS(W)
        XBM=0
        DO 70 K=1,3
        W=(KBM(K)/10.-JBM(K)*127.0633)/180*PI
     1     -(IBM(K)+JBM(K))*GE+JBM(K)*GM
70      XBM=XBM+SBM(K)*COS(W)
        XBJ=0
        DO 80 K=1,7
        W=(KBJ(K)/10.-JBJ(K)*88.445)/180*PI
     1    -(IBJ(K)+JBJ(K))*GE+JBJ(K)*GJ
80      XBJ=XBJ+SBJ(K)*COS(W)
        XBS=0
        DO 90 K=1,2
        W=(KBS(K)/10.-JBS(K)*10.2417)/180*PI
     1    -(IBS(K)+JBS(K))*GE+JBS(K)*GS
90      XBS=XBS+SBS(K)*COS(W)
C MONDSTOERUNGEN
        GD=296.104608D0+477198.849108D0*T+.9192D-2*T**2+14.D-6*T**3
        XM=DMOD(GD/180*PID,2*PID)
        GD=350.737486D0+445267.114217D0*T-.1436D-2*T**2+2.D-6*T**3
        D=DMOD(GD/180*PID,2*PID)
        GD=11.D0+15.D0/60+3.2D0/3600
     1     +(1739527290.54D0*T-11.56D0*T**2-.12D-2*T**3)/3600
        F=DMOD(GD/180*PID,2*PID)
        GD=259+10.D0/60+59.79D0/3600-6962911.23D0*T/3600
     2     +(7.48*T**2+.0086*T**3)/3600
        US=L-DMOD(GD*PID/180,2*PID)
        XLMO=6.454*SIN(D)+.013*SIN(3*D)+.177*SIN(D+XM)
     1       -.424*SIN(D-XM)+.039*SIN(3*D-XM)-.064*SIN(D+GE)
     2       +.172*SIN(D-GE)-.013*SIN(D-XM-GE)-.013*SIN(2*US)
        RMO=1336*COS(D)+3*COS(3*D)+37*COS(D+XM)-133*COS(D-XM)
     1      +8*COS(3*D-XM)-14*COS(D+GE)+36*COS(D-GE)
     2     -3*COS(D-GE-XM)+3*COS(2*US)
        BMO=.576*SIN(F)+.016*SIN(F+XM)-.047*SIN(F-XM)
     1      +.021*SIN(F-2*US)+.005*SIN(F-2*US-XM)
     2      +.005*SIN(F+GE)+.005*SIN(F-GE)
        L=L+(XLME+XLV+XLM+XLJ+XLS)/3600000./180*PI
        L=L+(C+XLMO)/3600*PI/180
        B=(XBV+XBM+XBJ+XBS)/3600000./180*PI
        B=-B+BMO/3600*PI/180
        R=(R+(RME+RV+RM+RJ+RS)/10+RMO)*1.D-8
        R=10.D0**R
        EKL=(23.D0+27.D0/60+8.26D0/3600)
        EKL=EKL-(46.845*T+.59D-2*T**2-.181D-2*T**3)/3600
        EKL=EKL/180*PID
        X(1)=DCOS(L)*DCOS(B)
        X(2)=DSIN(L)*DCOS(B)
        X(3)=DSIN(B)
        Y(1)=X(1)
        Y(2)=X(2)*DCOS(EKL)-X(3)*DSIN(EKL)
        Y(3)=X(2)*DSIN(EKL)+X(3)*DCOS(EKL)

*       Do not transform.  Leave values as system of date.

C PRECESSION TO 1950.0
C       CALL PRAE(XMJD,PRAEZ)
C       DO 77 I=1,3
C         X(I)=0.D0
C         V(I)=0.D0
C         DO 77 K=1,3
C           X(I)=X(I)+PRAEZ(K,I)*Y(K)
C           
C77       CONTINUE
C TRANSFORMATION FROM B1950.0 TO J2000.0
C       CALL FK4FK5(X,V,0D0,0D0,X,V,PR,VR)
C       DEPOCH=(XMJD-51544.5D0)/36525.D0
C       DO 79 I=1,3
C         X(I)=R*(X(I)+V(I)*DEPOCH)
C79      CONTINUE
C
        mag = sqrt(y(1)**2+y(2)**2+y(3)**2)
        do i = 1,3
           x(i) = y(i)/mag
        end do
        RETURN
        END

CTITLE XPROD

        subroutine xprod(a,b,c)

      implicit none

*       rotuine to take c = a x b

        real*8 a(3), b(3), c(3)

****    Directly compute cross product

        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)

        return
        end




 

