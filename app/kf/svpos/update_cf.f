CTITLE UPDATE_CF
 
      subroutine update_cf

      implicit none 
 
*     This routine will update the cfile using theorectical range
*     and phase computed in Earth fixed frame from an SP3 file.

      include '../../libraries/includes/const_freq.h'
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
      include 'modear.h'
 
* LOCAL variables
 
*   ierr    - IOSTAT error
*   i,j,k   - Loop counters
*   date(5) - Date
*   pn  - PRN number of current satellite
 
      integer*4 ierr, i,j,k, date(5), pn
 
*   sectag  - Seconds tag
*   mjd  - Julian date
 
 
      real*8 sectag, mjd
 
****  Open the cfile
      call open_cf( 100, in_cf, ierr )
      call report_error('IOSTAT',ierr,'open', in_cf,1,'update_cf')
 
****  Now create the output cfile
      call create_cf(200, out_cf, ierr )
      call report_error('IOSTAT',ierr,'creat', out_cf,1,'update_cf')
 
*     Read the first cfile block
      call read_cf1( 100, 'ALL', ierr )
 
*     Update the history entries on the file
      if( cf_ntext.lt.cf_maxtxt) then
          k = cf_ntext + 1
          cf_ntext = k
          call systime( date, sectag)
          write(cf_text(k),100) date, sectag
 100      format('MODEAR: Run ',i4,'/',i2,'/',i2,
     .            1x,i2,':',i2,1x,F5.2)
      end if
 
      call write_cf1(200, 'ALL', ierr )
 
***   Now read and write the type 2 and type 3 records
      call read_cf2(100, 'ALL',ierr )
      call write_cf2(200, 'ALL', ierr)
      call report_error('IOSTAT',ierr,'cf2 writ', out_cf,1,'update_cf')
      call read_cf3(100, 'ALL',ierr )

****  Now do the pre-ample computations.  These include:
*     (1) computing the position of the L1 and L2 phase centers of reciever
*     (2) computing the average atmospheroc delay for this site.
 
      call init_comp
 
*     Clear the clock parameters
      do i = 4,6
          cf_preval(i) = 0.d0
      end do
      if( cf_nextra.eq.1 ) cf_extra(1) = 0.d0

*     Save the value in meters
      cf_preval(7) = dry_zen
 
      call write_cf3(200, 'ALL', ierr)
      call report_error('IOSTAT',ierr,'cf3 writ', out_cf,1,'update_cf')
 
 
****  Now loop over the cfile
      do i = 1, cf_nepoch
          call read_cf4(100,'ALL', ierr )
 
*         Set the value of cf_rclock to zero (since all the clocks
*         are accounted for in this analysis.  For matching with
*         autcln)
          if( abs(cf_rclock-site_clock(i))*1.d6.gt.1.0 ) then
              write(*,410) i, cf_rclock*1.d9, site_clock(i)*1.d9
 410          format('CLKEP ',i5,' rclock ',f12.2,' site clock ',
     .               f12.2,' ns')
          end if
          
          cf_rclock = 0.d0
 
          call write_cf4(200, 'ALL',ierr )
          call report_error('IOSTAT',ierr,'cf4 read', in_cf,
     .                1,'update_cf')
 
****          Get the nominal time for this epoch
          call ydsd_to_mjd( cf_iyr, cf_idoy, cf_sod, mjd )
 
          data_epoch = mjd
 
****      Now do the epoch dependent computions:
*         (1) Compute the tidally displaced position of the site
 
          call epdep_comp
 
          do j = 1, cf_msat
              call read_cf5(100,'ALL', ierr )
 
*             Now compute the theoretical range and phase
*             and the omc values (only the latter will be used
*             by solve for the estimation).
                        pn = cf_iprn
              call theory_comp( i, pn)
 
              call write_cf5(200,'ALL',ierr)
              call report_error('IOSTAT',ierr,'cf5 writ',1,
     .                'update_cf')
 
          end do
      end do
 
****  Thats all
      close(200)
      close(100)
 
      return
      end
 
CTITLE INIT_COMP
 
      subroutine init_comp

      implicit none 
 
*     This routine will initialize the quantities needed in processing
*     Earth fixed c-files.
*     The main role of this routine is
*     (1) compute the XYZ coordinates of the phase centers of
*         the L1 and L2 antennas and
*     (2) the atmospheric delay
 
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
 
      include 'modear.h'
 
* LOCAL Variables
 
*   loc_coord(3)    - Geodetic Lat, long and height of site
*   rot_matrix(3,3) - Rotation maxtrix from NEU to XYZ
*   dNEU            - Vector formed into standard NEU order
 
 
      real*8 loc_coord(3), rot_matrix(3,3), dNEU(3)
 
****  Start, rotate the L1 and L2 phase centers into XYZ cartessian coordinates
*     Reorder the Phase center offsets
      dNEU(1) = cf_offsL1(2)
      dNEU(2) = cf_offsL1(3)
      dNEU(3) = cf_offsL1(1)
      call rotate_geod(dNEU,dXYZ_L1, 'NEU', 'XYZ',
     .              site_xyz, loc_coord, rot_matrix)

      write(*,120) cf_sitnam(1:16), site_xyz
 120  format(' Site ',a16,' Coordinates ',3F14.4,' m')
      write(*,140) dNEU
 140  format(' L1 Phase center offset ',3F14.4,' m')
     
      dNEU(1) = cf_offsL2(2)
      dNEU(2) = cf_offsL2(3)
      dNEU(3) = cf_offsL2(1)
      call rotate_geod(dNEU,dXYZ_L2, 'NEU', 'XYZ',
     .              site_xyz, loc_coord, rot_matrix)
 
      write(*,160) dNEU
 160  format(' L2 Phase center offset ',3F14.4,' m')

***** Compute the height dependent zenith delay
      call dry_hgtcon( loc_coord(3), loc_coord(1), dry_zen )
      dry_zen = cf_preval(7) 
 
      write(*,180) dry_zen
 180  format(' Zenith delay           ',F14.4,' m')

****  Thats all
      return
      end
 
CTITLE DRY_HGTCON
 
      subroutine dry_hgtcon (ellip_hgt, colatitude, dry_zen )

      implicit none 
 
*     Rotuine to compute zenith delay based on pressure computed from
*     height alone.  Hydrostatic equilrium is used to get the pressure
*     This rouitine takes the site number.
*
 
      include '../includes/const_param.h'
 
*  ellip_hgt  - Ellipsodial height (m)
*  colatitude - Geodetic colatitude radians
*  dry_zen    - Dry zenith delay (m)
 
 
      real*8 ellip_hgt, colatitude, dry_zen
 
* LOCAL VARIABLES
*   fphih       - Gravity function
*   latr        - real*4 Latitude
*   hgtr        - Real*4 value of height
 
 
      real*4 fphih, latr, hgtr
 
*   axis_press  - Pressure computed at axis of telescope
*   axis_temp   - Temperature at axis (computed assumming
*               - 293.15K at MSL)
 
 
      real*8 axis_press, axis_temp
 
***** Compute temperature and pressure at axis
 
      axis_temp = 293.15d0 - 6.5d-3*ellip_hgt
 
      axis_press = 1013.25d0*(axis_temp/293.15d0)**5.26d0
 
      latr = pi - colatitude
      hgtr = ellip_hgt
 
      call phi_function(latr, hgtr, fphih)
 
      dry_zen = (0.0022768d0*axis_press)/fphih

***** Thats all
      return
      end
 
CTITLE EPDEP_COMP
 
      subroutine epdep_comp

      implicit none 
 
*     This routine computes the epoch dependent model components
*     These are:
*     (1) Solid Earth tide
 
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
 
      include 'modear.h'
 
****  We only compute the tidal displacement
 
      call earth_tide( data_epoch, site_xyz, dXYZ_tide )
 
****  Thats all
      return
      end
 
CTITLE EARTH_TIDE
 
 
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
 
 
      real*8 loc_coord(3), rot_matrix(3,3),
     .    dn, de, dh, dNEU_tide(3), latr, longr, U0, Udlat, Udlong,
     .    dUdlat, dUdLong, loveh, lovel, jd
 
      DATA LOVEH/.6090d0/, LOVEL/.0852d0/
 
****  See type of date passed
      if( epoch.gt.2000000.d0 ) then
          jd = epoch
      else
          jd = epoch + 2400000.d0
      end if 
****  Get the lat and long of the site
      call XYZ_to_GEOD( rot_matrix, site_pos, loc_coord )
 
      latr = pi - loc_coord(1)
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
      de = lovel*dUdlong

*     Save and convert values from mm to meters
      dNEU_tide(1) = dn/1000.d0
      dNEU_tide(2) = de/1000.d0
      dNEU_tide(3) = dh/1000.d0
 
*     Now rotate to XYZ didplacements
      call rotate_geod( dNEU_tide, dXYZ_tide, 'NEU', 'XYZ',
     .        site_pos, loc_coord, rot_matrix)
 
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
 
CTITLE THEORY_COMP
 
      subroutine theory_comp( ep, pn )

      implicit none 
 
*     This routine does the computation of the theoretical range
*     to the satellite and computes the obs-minus-comp values for
*     storage in the cfiles.
 
 
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/cfile_def.h'
 
      include 'modear.h'
 
* PASSED VARIABLES
 
*   ep  - Epoch number of this measurement
*   pn  - PRN number to be computed for
 
      integer*4 ep, pn
 
* LOCAL VARIABLES
*   range   - Computed range from satellite to site
*   send_epoch  - Transmission epoch from SV
*   svs_phs(3)  - Coordinates of the satellite
*           - phase center
*   elev, map_fn    - Elevation angle (deg) and mapping
*           - function value
*   atm_range   - Atmospheric delay contribution to
*           - range
*   dpos_L1(3), dpos_L2(3)  - Difference in position from
*           - Site to SV for phase centers at L1 and L2
*   dclk        - Difference between ground and SV clock to
*           - be added to theoretical (sec)
*   theo        - Theoretical value for obsevable being processes
*   L1r, L2r    - L1 and L2 range with atmospheric delay (m)
*   prev_range  - Previous rnage estimate
*   drange      - Change in range between iterations
 
      real*8 range, send_epoch, svs_phs(3), elev, map_fn, atm_range,
     .    dpos_L1(3), dpos_L2(3), dclk, theo, L1r, L2r, prev_range,
     .    drange 
     
      integer*4 iter, i
 
****  Start by getting the transmitt time from the satelite
      range = 20000.d3
      drange = 10.0
      iter = 0
      do while ( abs(drange).gt.1.d-4 .and. iter.lt.5 )
          iter = iter + 1
*         Compute the range accounting for the propagation
*         delay.
*         Compute the transmission time
 
          send_epoch = data_epoch - site_clock(ep)/86400.d0 -
     .                range/vel_light/84000.d0
          call eph_to_xyz( send_epoch, pn, 'E')

          prev_range = range
 
          range = sqrt( (site_xyz(1)-svs_xyz(1,pn))**2+
     .                  (site_xyz(2)-svs_xyz(2,pn))**2+
     .                  (site_xyz(3)-svs_xyz(3,pn))**2)
          drange = range - prev_range
      end do

***** Check we did not exceed iteration bound
      if( iter.ge.5 ) then
          write(*,220) ep, pn, range, prev_range
 220      format(' WARNING: 5 or more range iterations at Epoch ',i5,
     .           ' PRN,'i2.2,' Last two ',2F16.4,' m')
      end if
 
****  Now compute the final position of the satellite center
*     of mass
      send_epoch = data_epoch - site_clock(ep)/86400.d0 -
     .            range/vel_light/84000.d0
 
      call eph_to_xyz( send_epoch, pn, 'E')
 
***** Now get the phase center of satellite
      call svs_cm_to_phs( data_epoch, pn, svs_xyz(1,pn), svs_phs )
 
***** Get the elevation angle to satellite
      call get_elev(site_xyz, svs_xyz(1,pn), elev )
 
***** Get the atmospheric delay contribution
      call mit_dry( elev, map_fn )
      atm_range = map_fn*dry_zen
 
****  Now for L1 and L2 compute the theorectical range
 
      do i = 1,3
          dpos_L1(i) = site_xyz(i)+dXYZ_L1(i)+dXYZ_tide(i) -
     .            svs_phs(i)
          dpos_L2(i) = site_xyz(i)+dXYZ_L2(i)+dXYZ_tide(i) -
     .            svs_phs(i)
      end do
 
****  Now get the L1 and L2 ranges
      L1r = sqrt(dpos_L1(1)**2+dpos_L1(2)**2+dpos_L1(3)**2) + atm_range
      L2r = sqrt(dpos_L2(1)**2+dpos_L2(2)**2+dpos_L2(3)**2) + atm_range

      write(*,200) ep, pn, L1r-atm_range, cf_tau*vel_light,
     .            (L1r-atm_range-cf_tau*vel_light)
 200  format('DELAY: ',2i4,3f15.4)
 
      cf_tau = L1r/vel_light
      cf_atmdel = atm_range/vel_light
      dclk = site_clock(ep) - svs_clk(pn)
 
*     Now compute OMC for each observable type
      do i = 1, cf_ndat
 
*         Check each of the observable types
*                                         ! L1 Carrier phase
          if( cf_dattyp(i).eq.1 ) then
              theo = (L1r/vel_light + dclk)*fL1
              cf_omcs(i) = cf_obsv(i) - theo
*                                             ! L2 Carrier phase
          else if( cf_dattyp(i).eq.2 ) then
              theo = (L2r/vel_light + dclk)*fL2
              cf_omcs(i) = cf_obsv(i) - theo
*                                             ! L1 P-code Range
          else if( cf_dattyp(i).eq.3 ) then
              theo = (L1r/vel_light + dclk)*fL1
              cf_omcs(i) = cf_obsv(i)/vel_light*fL1 - theo
*                                             ! L2 P-code Range
          else if( cf_dattyp(i).eq.4 ) then
              theo = (L2r/vel_light + dclk)*fL2
              cf_omcs(i) = cf_obsv(i)/vel_light*fL2 - theo
*                                             ! L1 CA-code Range
          else if( cf_dattyp(i).eq.5 ) then
              theo = (L1r/vel_light + dclk)*fL1
              cf_omcs(i) = cf_obsv(i)/vel_light*fL1 - theo
          end if
      end do
 
****  Thats all
      return
      end
 
CTITLE SVS_CM_TO_PHS
 
      subroutine svs_cm_to_phs( mjd, pn, svs_xyz, svs_phs )

      implicit none 
 
*     This routine computes the phase center position of satellite
*     given the position of the center of mass
 
 
*   mjd  - Modified Julian date (not used yet but may be needed later)
*   svs_xyz(3)  - Position of CM of satellite
*   svs_phs(3)  - Position of phase center
 
      real*8 mjd, svs_xyz(3), svs_phs(3)
 
*   pn      - PRN number of satellite
 
      integer*4 pn
 
* LOCAL VARIABLES
 
*   ehat(3) - Unit vector in direction of Earth
*   shat(3) - Unit vector in direction of sun
*   yhat(3) - Satellite y axis vector
*   mag     - Magnitude of svs_xyz vector
 
      real*8 ehat(3), shat(3), yhat(3), mag
 
*   i       - Loop counter
 
 
      integer*4 i
 
****  Compute unit vector in direction on Earth
      mag = sqrt(svs_xyz(1)**2+svs_xyz(2)**2+svs_xyz(3)**2)
      do i = 1,3
          ehat(i) = -svs_xyz(i)/mag
      end do
 
*     Be really slack and only do radial on Block II
      do i = 1,3
          svs_phs(i) = svs_xyz(i) + 0.9519d0*ehat(i)
      end do
 
****  That all
      return
      end
 
 
