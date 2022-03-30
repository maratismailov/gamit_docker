CTITLE COMPUTE_OC
 
      subroutine compute_OC( sol_obs, type, indx )

      implicit none  
 
*     This routine actually does the o minus c calcualtion
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/sd_common.h'
 
      include '../includes/glb_hdr_def.h'
 
*   indx        - Index for the parameter (eg site/source number)
*   i           - Loop counter
*   type        - Parameter type from code.
 
      integer*4 indx, type
 
*   sol_obs     - Parameter estimates from solution.  These values
*               - will be replaced with O_minus_C.
 
      real*8 sol_obs

*   sv_num      - satellite number (local)
*   orb_el      - Orbital element
*   m1, m2      - Split of multi-day PMU index

      integer*4 sv_num, orb_el, m1, m2, mul_pmu_ent
 
 
***** Use a compute goto to get to correct calcualtion for this
*     parameter

* NOTE: Should include types 52-55 which are translations, translation
*     rate, scale and scale rate but at the moment zero is the apriori
*     and so we don' need to update.
* MOD TAH 030906: Increase type to 64 to account for log term estimation
 
      if( type.lt.7 .or. type.gt.64 ) return
 
      goto(  100,  200,  300,  400,  500,  600,  700,  800,  900, 1000,
     .      1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
     .      2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
     .      3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
     .      4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000,
     .      5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800   ) type-6
 
***** X-site coordinate  (In the future if we change PM aprioris we
*     may need to include that here as well)
 100  continue
          call site_OC(1, ltog_sites(indx), sol_obs)
          return
 
***** Y-site coordinate
 200  continue
          call site_OC(2, ltog_sites(indx), sol_obs)
          return
 
***** Z-site coordinate
 300  continue
          call site_OC(3, ltog_sites(indx), sol_obs)
          return
 
***** Axis offset
  400 continue
          sol_obs = sol_obs - apr_val_axo(1,ltog_sites(indx)) -
     .              apr_val_axo(2,ltog_sites(indx))*
     .              (gepoch_expt-axo_epoch(ltog_sites(indx)))/365.25d0
          return
 
***** Source RA
  500 continue
          call source_OC(1, ltog_sources(indx), sol_obs)
          return
 
***** Source Dec
  600 continue
          call source_OC(2, ltog_sources(indx), sol_obs)
          return
 
***** Wobble components (the component is given by index)
  700 continue
          sol_obs = sol_obs - apr_val_wob(indx)
 
*         For the moment also remove the wobble which the data were
*         processed with
* MOD TAH 9405017: Update first four values

          if( indx.le.4 ) then
              sol_obs = sol_obs - gwob_apr(indx)
          end if
*         print *,'DEBUG: SOL_OBS GWOB ', sol_obs, indx, 
*    .                 gwob_apr(indx)
         return
 
***** UT1-AT
  800 continue
          sol_obs = sol_obs - apr_val_UT1(indx)
* MOD TAH 9405017: Update first two  values
          if( indx.le.2 ) then
              sol_obs = sol_obs - gut1_apr(indx)
* MOD TAH 990222: Remove any leap seconds
              sol_obs = sol_obs - nint(sol_obs/15.d3)*15.d3
          end if
          return
 
***** Nutation angle values
  900 continue
          sol_obs = sol_obs - apr_val_nut_ang(indx)
          if( indx.le.2 ) then
              sol_obs = sol_obs - vnut_ang_apr(indx)
          end if
          return
 
***** Tide Love l
 1000 continue
          sol_obs = sol_obs - apr_val_tid(1,ltog_sites(indx))
          return
 
***** Tide Love h
 1100 continue
          sol_obs = sol_obs - apr_val_tid(2,ltog_sites(indx))
          return
 
***** Tide Lag angle
 1200 continue
          sol_obs = sol_obs - apr_val_tid(3,ltog_sites(indx))
          return
 
***** Extended Earth tide model
 1300 continue
 1400 continue
 1500 continue
 1600 continue
 1700 continue
 1800 continue
 1900 continue
 2000 continue
 2100 continue
 2200 continue
 2300 continue
 2400 continue

* MOD TAH 910130: Quick fix on not enough tides
         if( ltog_sites(indx).gt.max_sites) RETURN
         sol_obs = sol_obs - sd_tides_mid(type-18, ltog_sites(indx))
         return
***** UT1 diurnal and semidiurnal coefficents
 2500 continue
 2600 continue
 2700 continue
 2800 continue
         sol_obs = sol_obs - sd_ut1_mid(type-30)
         return

***** Six values for the polar diurnal and semidiurnal
*     signals.
 2900 continue
 3000 continue
 3100 continue
 3200 continue
 3300 continue
 3400 continue
          sol_obs = sol_obs - sd_xy_mid(type-34)
          return
 
***** Gamma
 3500 continue
          sol_obs = sol_obs - apr_val_gamma
          return
 
***** X position rate
 3600 continue
          call site_OC(4, ltog_sites(indx), sol_obs)
          return
 
***** Y position rate
 3700 continue
          call site_OC(5, ltog_sites(indx), sol_obs)
          return
 
***** Z position rate
 3800 continue
          call site_OC(6, ltog_sites(indx), sol_obs)
          return
 
***** RA rate
 3900 continue
          call source_OC(3, ltog_sources(indx), sol_obs)
          return
 
***** Dec rate
 4000 continue
          call source_OC(4, ltog_sources(indx), sol_obs)
          return
 
***** Nutation series coefficients
 4100 continue
          call coeff_OC(type, indx, apr_val_nut_coeff, sol_obs)
          return
 
 
***** EXTENDED Earth tide coefficiants (Implemement later)
 4200 continue
          return

***** UT1 diurnal and semidiurnal coefficients
 4300 continue
          return
 
***** Polar motion diurnal and semidiurnal coefficients
 4400 continue
          return

***** SV Orbits 
 4500 continue

*         decode the index into sv number and orbital element
          call decode_code( indx, orb_el, sv_num )

          sol_obs = sol_obs - apr_val_svs(orb_el, ltog_svs(sv_num) )
          return

***** Translation parameters
 4600 continue
          return

***** Translation rate parameters
 4700 continue
          return

***** Scale parameters       
 4800 continue
          return

***** Scale parameters       
 4900 continue
          return

***** Multiday X-pole position
 5000 continue

*         Split the indx to get the offset/rate index and the epoch
*         number that this one referrs to
          call decode_code( indx, m1, m2)
          call get_mul_pmu_ent ( cmul_pmu_ep(m2), mul_pmu_ent, 'Y' )
          sol_obs = sol_obs - apr_val_mul_pmu(m1,type-55,mul_pmu_ent )
*         print *,'DEBUG: SOL_OBS MUL ', sol_obs, m1, type-55,
*    .           mul_pmu_ent,apr_val_mul_pmu(m1,type-55,mul_pmu_ent )
 
          return

***** Multiday Y-pole position
 5100 continue

*         Split the indx to get the offset/rate index and the epoch
*         number that this one referrs to
          call decode_code( indx, m1, m2)
          call get_mul_pmu_ent ( cmul_pmu_ep(m2), mul_pmu_ent, 'Y' )
          sol_obs = sol_obs - apr_val_mul_pmu(m1,type-55,mul_pmu_ent )
          return

***** Multiday  UT1 values
 5200 continue

*         Split the indx to get the offset/rate index and the epoch
*         number that this one referrs to
          call decode_code( indx, m1, m2)
          call get_mul_pmu_ent ( cmul_pmu_ep(m2), mul_pmu_ent, 'Y' )
          sol_obs = sol_obs - apr_val_mul_pmu(m1,type-55,mul_pmu_ent )

* MOD TAH 990222: Remove any leap seconds
          sol_obs = sol_obs - nint(sol_obs/15.d3)*15.d3

          return


***** Generic rotations about the XYZ axes
 5300 continue

          return

***** Generic rotation rates about the XYZ axes
 5400 continue

          return

***** Atmospheric Zenith delay parameters (also multi-valued)
 5500 continue
          sol_obs = sol_obs - apr_val_atm(ltog_sites(indx))
          return

***** Log N component
 5600 continue
          sol_obs = sol_obs - apr_val_log(1,ltog_sites(indx))
          return

***** Log E component
 5700 continue
          sol_obs = sol_obs - apr_val_log(2,ltog_sites(indx))
          return

***** Log U component
 5800 continue
          sol_obs = sol_obs - apr_val_log(3,ltog_sites(indx))
          return

***** For the moment that is all
      end
 
 
 
