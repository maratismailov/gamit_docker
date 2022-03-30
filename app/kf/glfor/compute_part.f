CTITLE COMPUTE_PARTIALS
 
      subroutine compute_partials( pn, type, indx, part_pnt, a_part,
     .                             amp_deriv )

      implicit none  
 
*     This rotuine actually cpmputes the partial derivatives.  A
*     compute GOTO is used to separate the different types.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
 
*   coeff       - Coefficient number for nutation series coefficients
*   indx        - index of parameter eg site/source number
*   i,j,k       - Loop counters
*   part_pnt(2,1)   - pointer to partials (gives start parameter,
*               - and number for each parameter to be used)
*   pn          - Parameter number from current solution being
*               - processed
*   phase       - Phase for nutation coefficients (1-in phase,
*               - 2-out of phase)
*   type        - Type of parameter

*   sv_num      - Satellite number
*   orb_el      - Orbital element

*   ts          - Tide site.  Actual number if local tides, and
*                 1 if glb tides.
*   loc_type    - Local type relative to type (eg extended tides 
*                 go from 1 to 12
*   amp_type    - Type of amplitude paramtere (for tides values from
*                 1 to 6 (diurn, rad/N/E and semi Rad/N/E )
*   cs_comp     - 1 if cos component, 2 if sin component
*   num_amp     - Generated while we add partials for different amplitude
*                 coefficients wrt estimated dialy value.

 
      integer*4 coeff, indx, i,j,k, part_pnt(2,*), pn, phase, type,
     .          sv_num, orb_el

*   mul_pmu_ent -- Epoch entry number for multipmu estimates
*   m1, m2      -- Rate and XY indices for mul_pmu parn entries

      integer*4 mul_pmu_ent, m1, m2

      integer*4 ts, loc_type, amp_type, cs_comp, num_amp

*   kbit        - Sees if a bit is set

      logical kbit
 
*   a_part(1)   - The partials for this parameter
*   coeff_deriv(2,max_nut_coeff)    - Nutation series coefficient
*               - derivatives
*   namp_deriv(2,4,max_nut_coeff)    - Nutation amplitude derivatives
*   deriv(2)    - Two dummy derivatives for nut_coeff_part
*   epoch       - Main memory epoch
*   fcn_period  - Main memory copy of FCN period
*   amp_deriv(2,4,max_coeff_periods) - Amplitude derivatives for a+ and a-
*                 in-phase and out-of-phase for all periods (computed
*                 before hand)
 
      real*8 a_part(*), coeff_deriv(2,max_nut_coeff), deriv(2), epoch,
     .    fcn_period, amp_deriv(2,4,max_coeff_periods)
 
 
***** Use the GOTO
 
      if ( type.lt.7 .or. type.gt.65 ) return
 
      goto(  100,  200,  300,  400,  500,  600,  700,  800,  900, 1000,
     .      1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
     .      2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000,
     .      3100, 3200, 3300, 3400, 3500, 3600, 3700, 3800, 3900, 4000,
     .      4100, 4200, 4300, 4400, 4500, 4600, 4700, 4800, 4900, 5000,
     .      5100, 5200, 5300, 5400, 5500, 5600, 5700, 5800  ) type-6
 
***** X_site coordinate
 100  continue
          if( kbit(guse_site,ltog_sites(indx))) then
              call site_partial( 1, pn, part_pnt, a_part, 
     .                           ltog_sites(indx))
          end if
          return
 
***** Y_site coordinate
 200  continue
          if( kbit(guse_site,ltog_sites(indx))) then
              call site_partial( 2, pn, part_pnt, a_part, 
     .                           ltog_sites(indx))
          end if
          return
 
 
***** Z_site coordinate
 300  continue
          if( kbit(guse_site,ltog_sites(indx))) then
              call site_partial( 3, pn, part_pnt, a_part, 
     .                           ltog_sites(indx))
          end if
          return
 
***** Axis offset partial
 400  continue
 
*         See if turned on
          call add_partial( indx_pnt(pn), part_pnt,
     .         parn_axo(1,ltog_sites(indx)), a_part, 1.d0)
 
          return
 
***** Source RA
  500 continue
          call source_partial( 1, pn, part_pnt, a_part,
     .                         ltog_sources(indx))
          return
 
***** Source Dec
  600 continue
          call source_partial( 2, pn, part_pnt, a_part,
     .                         ltog_sources(indx))
          return
 
***** Wobble components
  700 continue
 
*         See if turned on and if we are using multiday or single estimates
          if( num_mul_pmu.eq.0 ) then
             call add_partial( indx_pnt(pn), part_pnt, parn_wob(indx),
     .                         a_part, 1.d0)
          else
*            Find which day this should be referred to.  Decode the epoch
*            pointer
             call get_mul_pmu_ent( cepoch_expt, mul_pmu_ent, 'Y' )
             call map_wob_indx('STOM', indx, m1, m2)
             if( mul_pmu_ent.gt.0 ) then
                call add_partial( indx_pnt(pn), part_pnt, 
     .                            parn_mul_pmu(m1, m2, mul_pmu_ent),
     .                            a_part, 1.d0)
             else
*               Warn user
                if( kbit(mul_pmu_opt,5) )
     .          write(*,720) cepoch_expt
 720            format('** WARNING ** No match to MUL_PMU_EP for',
     .                 ' Experiment epoch ',f12.3,' WOBBLE')
             end if
          end if

          return
 
***** UT1components
  800 continue
 
*         See if turned on
          if( num_mul_pmu.eq.0 ) then
              call add_partial( indx_pnt(pn), part_pnt,
     .            parn_ut1(indx), a_part, 1.d0)
          else
*            Find which day this should be referred to.  Decode the epoch
*            pointer
             call get_mul_pmu_ent( cepoch_expt, mul_pmu_ent, 'Y' )
             if( mul_pmu_ent.gt.0 ) then
                call add_partial( indx_pnt(pn), part_pnt, 
     .                            parn_mul_pmu(indx, 3, mul_pmu_ent),
     .                            a_part, 1.d0)
             else
*               Warn user
                if( kbit(mul_pmu_opt,5) )
     .          write(*,820) cepoch_expt
 820            format('** WARNING ** No match to MUL_PMU_EP for',
     .                 ' Experiment epoch ',f12.3,' UT1/LOD')
             end if
          end if
****      See if we are ignoring UT1 values
          if( indx.eq.1 .and. kbit(mul_pmu_opt,4) ) then
               indx_pnt(pn) = 0
          end if
          return
 
***** Nutation angles
  900 continue
 
*         See if turned on
          call add_partial( indx_pnt(pn), part_pnt,
     .         parn_nut_ang(indx), a_part, 1.d0)
 
*         Now check for series coefficients
          if( indx.eq.1 ) then
              deriv(1) = 1.d0
              deriv(2) = 0.d0
          end if
          if( indx.eq.2 ) then
              deriv(1) = 0.d0
              deriv(2) = 1.d0
          end if
 
          fcn_period = nut_period
          epoch = gepoch_expt
 
*                                              ! Contribution to series
          if( indx.eq.1 .or. indx.eq.2 ) then
              call nut_coeff_part(deriv,epoch,fcn_period,coeff_deriv)
c             call gamp_parts(epoch, fcn_period, namp_deriv)
              do i = 1, max_nut_coeff
*                                 ! In and Out of phase
                  do j = 1,2
                      call add_partial( indx_pnt(pn), part_pnt,
     .                  parn_nut_coeff(j,i), a_part, coeff_deriv(j,i))
                  end do
c                 do j = 1, 4
c                     call add_partial( indx_pnt(pn), part_pnt,
c    .                  parn_nut_coeff(j,i), a_part, 
c    .                  namp_deriv(indx,j,i))
c                 end do
              end do
          end if
 
          return
 
****  Earth tide l
 1000 continue
          call add_partial( indx_pnt(pn), part_pnt,
     .        parn_tid(1,ltog_sites(indx)), a_part, 1.d0)
          return
 
 
****  Earth tide h
 1100 continue
          call add_partial( indx_pnt(pn), part_pnt,
     .        parn_tid(2,ltog_sites(indx)), a_part, 1.d0)
          return
 
****  Earth tide Lag angle
 1200 continue
          call add_partial( indx_pnt(pn), part_pnt,
     .        parn_tid(3,ltog_sites(indx)), a_part, 1.d0)
          return
 
***** Now we need the code for extended earth tide parmeters to coefficients
*     All of these can be done at once.  The correct partial is computed from
*     from the type.      
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

*****	  All the tides are done at once.  First get the tide site
          if( glb_glb_tides ) then
              ts = 1
          else
              ts = ltog_sites(indx)
          end if
              
*         Get the tide type and cos/sin index
          loc_type = type - 18 
          amp_type = (loc_type-1)/2 + 1
          cs_comp  = mod(loc_type+1,2) + 1

*         Add partial for itself
          call add_partial( indx_pnt(pn), part_pnt,
     .         parn_eor_etd(loc_type,ts), a_part, 1.d0)

*         Now do the coefficient partial.  Only do if ts is
*         id less than max number of tidal sites.
          if( ts.le.max_etd_sites ) then
              do i = 1, max_coeff_periods
                  num_amp = (i-1)*max_etd_names + (amp_type-1)*2
*                 Loop for +- period with j, and in/out phase with k
                  do j = 1,2
                     do k = 1,2 
                        call add_partial(indx_pnt(pn), part_pnt,
     .                     parn_etd_coeff(k,num_amp+j,ts), a_part,
     .                     amp_deriv(cs_comp,(j-1)*2+k,i))
                     end do
                  end do
              end do
          end if

          return

****  UT1 amplitude partials.  Code is similiar to Extended Earth tides
*     but not as complicated since we do not have site dependence
 2500 continue
 2600 continue
 2700 continue
 2800 continue

*         Get the ampltude type and cos/sin index
          loc_type = type - 30
          amp_type = (loc_type-1)/2 + 1
          cs_comp  = mod(loc_type+1,2) + 1

*         Add partial for itself
          call add_partial( indx_pnt(pn), part_pnt,
     .         parn_eor_ut1(loc_type), a_part, 1.d0)

*         Now do the coefficient partial
          do i = 1, max_coeff_periods
              num_amp = (i-1)*max_ut1_names + (amp_type-1)*2  
*             Loop for +- period with j, and in/out phase with k
              do j = 1,2
                 do k = 1,2
                    call add_partial(indx_pnt(pn), part_pnt,
     .                 parn_ut1_coeff(k,num_amp+j), a_part,
     .                 amp_deriv(cs_comp,(j-1)*2+k,i))
                 end do
              end do
          end do

          return

***** XY Pole diurnal and semidiurnal partialsd
 2900 continue
 3000 continue
 3100 continue
 3200 continue
 3300 continue
 3400 continue

*         Get the ampltude type and cos/sin index
          loc_type = type - 34
          amp_type = (loc_type-1)/2 + 1
          cs_comp  = mod(loc_type+1,2) + 1

*         Add partial for itself
          call add_partial( indx_pnt(pn), part_pnt,
     .         parn_eor_xy(loc_type), a_part, 1.d0)

*         Now do the coefficient partial
          do i = 1, max_coeff_periods
              num_amp = (i-1)*max_xy_names + (amp_type-1)*2
*             Loop for +- period with j, and in/out phase with k
              do j = 1,2
                 do k = 1,2
                    call add_partial(indx_pnt(pn), part_pnt,
     .                 parn_xy_coeff(k,num_amp+j), a_part,
     .                 amp_deriv(cs_comp,(j-1)*2+k,i))
                 end do
              end do
          end do

          return

***** Gamma partial
 3500 continue 
          call add_partial( indx_pnt(pn), part_pnt,
     .                      parn_gamma, a_part, 1.d0)
          return
 
***** X coordinate rate
 3600 continue
          if( kbit(guse_site,ltog_sites(indx))) then
              call site_rate_part( 1, pn, part_pnt, a_part, 
     .                           ltog_sites(indx))
C             call add_partial( indx_pnt(pn), part_pnt,
C    .                parn_site(1,2,ltog_sites(indx)), a_part, 1.d0)
          end if
          return
 
***** Y coordinate rate
 3700 continue
          if( kbit(guse_site,ltog_sites(indx))) then
              call site_rate_part( 2, pn, part_pnt, a_part, 
     .                           ltog_sites(indx))
C             call add_partial( indx_pnt(pn), part_pnt,
C    .                parn_site(2,2,ltog_sites(indx)), a_part, 1.d0)
          end if
          return
 
***** Z coordinate rate
 3800 continue
          if( kbit(guse_site,ltog_sites(indx))) then
              call site_rate_part( 3, pn, part_pnt, a_part, 
     .                           ltog_sites(indx))
C             call add_partial( indx_pnt(pn), part_pnt,
C    .                parn_site(3,2,ltog_sites(indx)), a_part, 1.d0)
          end if
          return
 
***** RA rate
 3900 continue
          call add_partial( indx_pnt(pn), part_pnt,
     .            parn_source(1,2,indx), a_part, 1.d0)
 
***** Dec rate
 4000 continue
          call add_partial( indx_pnt(pn), part_pnt,
     .            parn_source(2,2,indx), a_part, 1.d0)
 
***** Nutation series coefficients
 4100 continue
          coeff = mod(indx,128)
          phase = (indx-coeff)/128 + 1
          call add_partial( indx_pnt(pn), part_pnt,
     .            parn_nut_coeff(phase,coeff), a_part, 1.d0)
          return

***** Extended Earth tide coefficients.  There is a problem here.
*     if periods above 10 are used then then the decoding will 
*     not work.  Maybe we should exchange the position of site and
*     coefficient number.  ONLY A PROBLEM WHEN COEFFICIENTS FROM
*     PREVIOUS GLOBAL ARE USED.
 4200 continue
          call decode_code( indx, amp_type, ts )
          if( kbit(amp_type,16) ) then
              cs_comp = 2
          else
              cs_comp = 1
          end if
          amp_type = mod(amp_type,128)
          ts = ltog_sites(ts)
          if( glb_glb_tides ) ts = 1
          if( ts.le.max_etd_sites ) then
              call add_partial( indx_pnt(pn), part_pnt,
     .               parn_etd_coeff(cs_comp, amp_type, ts), 
     .               a_part, 1.d0)
          end if

          return

****  UT1 coefficients
 4300 continue
          amp_type = indx
          if( kbit(amp_type,16) ) then
              cs_comp = 2
          else
              cs_comp = 1
          end if
          amp_type = mod(amp_type,128)
          call add_partial(indx_pnt(pn), part_pnt,
     .           parn_ut1_coeff(cs_comp, amp_type), a_part, 1.d0)

          return

****  XY coefficients
 4400 continue
          amp_type = indx
          if( kbit(amp_type,16) ) then
              cs_comp = 2
          else
              cs_comp = 1
          end if
          amp_type = mod(amp_type,128)
          call add_partial(indx_pnt(pn), part_pnt,
     .           parn_xy_coeff(cs_comp, amp_type), a_part, 1.d0)

          return

***** Satellite orbits 
 4500 continue

*        split the sv number from the orbital element
         call decode_code( indx, orb_el, sv_num )
         call add_partial( indx_pnt(pn), part_pnt, 
     .        parn_svs(orb_el, ltog_svs(sv_num)), a_part, 1.d0 )

         return

***** Translation parameters
 4600 continue
          call add_partial(indx_pnt(pn), part_pnt, parn_tran(indx,1),
     .                     a_part, 1.d0)
          return

***** Translation rate parameters
 4700 continue
          call add_partial(indx_pnt(pn), part_pnt, parn_tran(indx,2),
     .                     a_part, 1.d0)
          return

***** Scale parameters       
 4800 continue
          call add_partial(indx_pnt(pn), part_pnt, parn_scale(1),
     .                     a_part, 1.d0)
          return

***** Scale parameters       
 4900 continue
          call add_partial(indx_pnt(pn), part_pnt, parn_scale(2),
     .                     a_part, 1.d0)
          return

***** Multiday X-pole position
 5000 continue

*         Split the indx to get the offset/rate index and the epoch
*         number that this one referrs to.
          call decode_code( indx, m1, m2)
*         Now get the epoch that this corresponds to in current solution
          call get_mul_pmu_ent ( cmul_pmu_ep(m2), mul_pmu_ent, 'Y' )
          if ( mul_pmu_ent.gt.0 ) then
              call add_partial( indx_pnt(pn), part_pnt, 
     .                          parn_mul_pmu(m1,1, mul_pmu_ent),
     .                          a_part, 1.d0)
          else
*             Warn user
              if( kbit(mul_pmu_opt,5) ) 
     .        write(*,5020) cmul_pmu_ep(m2),'X wob',m1,m2
 5020         format('** WARNING ** No match to MUL_PMU_EP for',
     .                 ' entry epoch ',f12.3,1x,a,' Type/Ent ',2i4)
          end if

          return

***** Multiday Y-pole position
 5100 continue

*         Split the indx to get the offset/rate index and the epoch
*         number that this one referrs to.
          call decode_code( indx, m1, m2)
*         Now get the epoch that this corresponds to in current solution
          call get_mul_pmu_ent ( cmul_pmu_ep(m2), mul_pmu_ent, 'Y' )
          if ( mul_pmu_ent.gt.0 ) then
              call add_partial( indx_pnt(pn), part_pnt, 
     .                          parn_mul_pmu(m1,2, mul_pmu_ent),
     .                          a_part, 1.d0)
          else
*             Warn user
              if( kbit(mul_pmu_opt,5) )
     .        write(*,5020) cmul_pmu_ep(m2),'Y wob',m1,m2
          end if

          return

***** Multiday  UT1 values
 5200 continue

*         Split the indx to get the offset/rate index and the epoch
*         number that this one referrs to.
          call decode_code( indx, m1, m2)
*         Now get the epoch that this corresponds to in current solution
          call get_mul_pmu_ent ( cmul_pmu_ep(m2), mul_pmu_ent, 'Y' )
          if ( mul_pmu_ent.gt.0 ) then
              call add_partial( indx_pnt(pn), part_pnt, 
     .                          parn_mul_pmu(m1, 3, mul_pmu_ent),
     .                          a_part, 1.d0)
          else
*             Warn user
              if( kbit(mul_pmu_opt,5) )
     .        write(*,5020) cmul_pmu_ep(m2),'UT1',m1,m2
          end if 

****      See if we are ignoring UT1 values
          if( m1.eq.1 .and. kbit(mul_pmu_opt,4) ) then
               indx_pnt(pn) = 0
          end if
          return


***** Generic rotations about the XYZ axes
 5300 continue
          call add_partial(indx_pnt(pn), part_pnt, parn_rot(indx,1),
     .                     a_part, 1.d0)

          return

***** Generic rotation rates about the XYZ axes
 5400 continue
          call add_partial(indx_pnt(pn), part_pnt, parn_rot(indx,2),
     .                     a_part, 1.d0)
          return

***** Atmospheric Zenith delay parameters (also multi-valued)
 5500 continue
          call add_partial(indx_pnt(pn), part_pnt, 
     .                     parn_atm(ltog_sites(indx)),
     .                     a_part, 1.d0)
          return

****  LOG Esimates in North
 5600 continue
          call add_partial(indx_pnt(pn), part_pnt, 
     .                     parn_log(1,ltog_sites(indx)),
     .                     a_part, 1.d0)
          return

****  LOG Esimates in East
 5700 continue
          call add_partial(indx_pnt(pn), part_pnt, 
     .                     parn_log(2,ltog_sites(indx)),
     .                     a_part, 1.d0)
          return

****  LOG Esimates in Height
 5800 continue
          call add_partial(indx_pnt(pn), part_pnt, 
     .                     parn_log(3,ltog_sites(indx)),
     .                     a_part, 1.d0)
          return
        
 
***** Thats all for the moment we will extend later on
      end
 
 
