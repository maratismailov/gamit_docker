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
C     cf_preval(7) = dry_zen
 
      call write_cf3(200, 'ALL', ierr)
      call report_error('IOSTAT',ierr,'cf3 writ', out_cf,1,'update_cf')
 
 
****  Now loop over the cfile
      do i = 1, cf_nepoch
          call read_cf4(100,'ALL', ierr )
 
*         Set the value of cf_rclock to zero (since all the clocks
*         are accounted for in this analysis.  For matching with
*         autcln)
          if( abs(cf_rclock-site_clock(i))*1.d6.gt.1.0 .and.
     .        cf_rclock.ne.0.d0      ) then
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
              call report_error('IOSTAT',ierr,'cf5 writ',in_cf,1,
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

C     do i = 1,3
C        dXYZ_tide(i) = 0
C     end do
 
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

      write(*,300) (epoch-int(epoch))*86400.d0, dNEU_tide,
     .              dXYZ_tide, jd
 300  format('TIDE:',f16.6,6F12.5,1x,F15.6 )
 
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
      include '../../libraries/includes/const_freq.h'    
 
      include 'modear.h'
 
* PASSED VARIABLES
 
*   ep  - Epoch number of this measurement
*   pn  - PRN number to be computed for
 
      integer*4 ep, pn
 
* LOCAL VARIABLES
*   range   - Computed range from satellite to site
*   send_epoch  - Transmission epoch from SV
*   recv_epoch  - Receive time for the measurement
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
*   site_int(3) - Site position in inertial space
*   site_L1i(3), site_L2i(3) - L1 and L2 phase center in inertial space
*   site_L1e(3), site_L2e(3) - L1 and L2 phase center in Earth fixed
*   svs_ear(3)  - Satelitte in Earth fixed frame (for debug)

*   dphs_l1_mm, dphs_l2_mm  - Phase center controbutions (mm)
*   dphs_l1, dphs_l2  - Phase center controbutions (cycles)
 
      real*8 range, send_epoch, svs_phs(3), elev, map_fn, atm_range,
     .    dpos_L1(3), dpos_L2(3), dclk, theo, L1r, L2r, prev_range,
     .    drange, site_int(3), site_L1i(3), site_L2i(3), recv_epoch,
     .    site_L1e(3), site_L2e(3), recv_sec, send_sec, svs_ear(3),
     .    dphs_l1_mm, dphs_l2_mm, dphs_l1, dphs_l2
     
      integer*4 iter, i
 
****  Start by getting the transmitt time from the satelite
      range = 20000.d3
      drange = 10.0
      iter = 0
      recv_epoch = data_epoch - site_clock(ep)/86400.d0
      call earth_to_inert(recv_epoch, site_xyz, site_int,'E','I')

*     Get the satellite clock errors
      call comp_svs_clk( data_epoch )
            
      do while ( abs(drange).gt.1.d-4 .and. iter.lt.5 )
          iter = iter + 1
*         Compute the range accounting for the propagation
*         delay.
*         Compute the transmission time
         
          send_epoch = recv_epoch - range/vel_light/86400.d0
          call eph_to_xyz( send_epoch, pn, 'I')

          prev_range = range
 
          range = sqrt( (site_int(1)-svs_xyz(1,pn))**2+
     .                  (site_int(2)-svs_xyz(2,pn))**2+
     .                  (site_int(3)-svs_xyz(3,pn))**2)
          drange = range - prev_range
      end do

***** Check we did not exceed iteration bound
      if( iter.ge.5 ) then
          write(*,220) ep, pn, range, prev_range
 220      format(' WARNING: 5 or more range iterations at Epoch ',i5,
     .           ' PRN ',i2.2,' Last two ',2F16.4,' m')
      end if
 
****  Now compute the final position of the satellite center
*     of mass
      send_epoch = recv_epoch - range/vel_light/86400.d0
      
*     Get the satellite clock errors
      call comp_svs_clk( send_epoch )
 
      call eph_to_xyz( send_epoch, pn, 'I')
 
***** Now get the phase center of satellite
      call svs_cm_to_phs( data_epoch, pn, svs_xyz(1,pn), svs_phs )
 
***** Get the elevation angle to satellite
      call get_elev(site_int, svs_xyz(1,pn), elev )
 
***** Get the atmospheric delay contribution
      call mit_dry( elev, map_fn )

      write(*,400) elev, map_fn, cf_tmpart(4)*(vel_light/fL1),
     .          ( map_fn -  cf_tmpart(4)*(vel_light/fL1))*dry_zen 
 400  format('ATM: ',f10.4, 2F10.4, F10.4 )
      
C     map_fn = cf_tmpart(4)*(vel_light/fL1)
      
      atm_range = map_fn*dry_zen

****  Get phase center controbution
      call phase_mod(elev, dphs_l1_mm, dphs_l2_mm)
      dphs_l1 = dphs_l1_mm/1000.d0/vel_light*fL1
      dphs_l2 = dphs_l2_mm/1000.d0/vel_light*fL2
 
****  Now for L1 and L2 compute the theorectical range
 
      do i = 1,3
          site_L1e(i) = site_xyz(i)+dXYZ_L1(i)+dXYZ_tide(i)
          site_L2e(i) = site_xyz(i)+dXYZ_L2(i)+dXYZ_tide(i)
      end do
                
      call earth_to_inert(recv_epoch, site_L1e, site_L1i,'E','I')
      call earth_to_inert(recv_epoch, site_L2e, site_L2i,'E','I')
      
      do i = 1,3
          dpos_L1(i) = site_L1i(i) - svs_phs(i)
          dpos_L2(i) = site_L2i(i) - svs_phs(i)
      end do
 
****  Now get the L1 and L2 ranges
      L1r = sqrt(dpos_L1(1)**2+dpos_L1(2)**2+dpos_L1(3)**2) + atm_range
      L2r = sqrt(dpos_L2(1)**2+dpos_L2(2)**2+dpos_L2(3)**2) + atm_range

* DEBUG: TAH Get earth fixed for output
C     call eph_to_xyz( send_epoch, pn, 'E')
 
***** Now get the phase center of satellite
C     call svs_cm_to_phs( data_epoch, pn, svs_xyz(1,pn), svs_phs )

*     Rotate satellite back to Earth fixed frame
      call earth_to_inert( send_epoch, svs_phs, svs_ear, 'I','E')

      recv_sec= (recv_epoch-int(recv_epoch))*86400.d0
      send_sec= (send_epoch-int(send_epoch))*86400.d0

      write(*,200) ep, pn, recv_sec, send_sec, site_clock(ep),
     .             (l1r-atm_range), cf_tau*vel_light,
     .             (l1r-atm_range)- cf_tau*vel_light,
     .             pn, svs_ear, site_L1e
 200  format('DELAY: ',i4,i3,3f17.10, 2f15.4,F8.4,1x,/,
     .       ' PRN ',i2,1x,6f15.4)

 
      cf_tau = (L1r-atm_range)/vel_light
      cf_atmdel = atm_range/vel_light
      dclk = site_clock(ep) - svs_clk(pn)
 
*     Now compute OMC for each observable type
      do i = 1, cf_ndat
 
*         Check each of the observable types
*                                         ! L1 Carrier phase
          if( cf_dattyp(i).eq.1 ) then
              theo = (L1r/vel_light + dclk)*fL1 + dphs_l1
              cf_omcs(i) = cf_obsv(i) - theo
*                                             ! L2 Carrier phase
          else if( cf_dattyp(i).eq.2 ) then
              theo = (L2r/vel_light + dclk)*fL2 + dphs_l2
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
 
      include '../includes/const_param.h' 
 
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
*   tdt     - Dynamical time for Sun's position
*   sun_dist - Distance to Sun (AU)
*   sun_long, sun_lat - Longitude and Latitude off sun
*   sv_cm(3)  - Position of Satellite phase center in local
*              system.
 
      real*8 ehat(3), shat(3), yhat(3), mag, tdt, sun_dist, sun_long,
     .       sun_lat, sv_cm(3), xhat(3)

      real*8 old_tdt
 
*   i       - Loop counter
 
 
      integer*4 i

      save old_tdt
 
****  Compute unit vector in direction on Earth
      mag = sqrt(svs_xyz(1)**2+svs_xyz(2)**2+svs_xyz(3)**2)
      do i = 1,3
          ehat(i) = -svs_xyz(i)/mag
      end do

****  Get TDT and compute vector to sun
      tdt =  mjd - 94554.d0 + 142350.d0 +(19.d0+32.184d0)/86400.d0

      call SUN20(tdt ,shat ,sun_dist, sun_long, sun_lat)


****  Now get Y and X-axis
      call xprod(ehat, shat, yhat )
      mag = sqrt(yhat(1)**2+yhat(2)**2+yhat(3)**2)
      if( mag.lt.0.01d0 ) write(*,110) mag
 110  format('WARNING: S/C z-axis close to Sun direction. Mag ', F9.6)
      do i = 1,3
          yhat(i) = yhat(i)/mag
      end do
     
      call xprod(yhat, ehat, xhat )

      if( abs(tdt-old_tdt).gt.0.02 ) then
          write(*,100) (mjd-int(mjd))*86400.d0, shat,
     .                 xhat, yhat, ehat, 
     .                (mjd-int(mjd))*86400.d0, svs_xyz

 100      format('SUN: ',f10.4,3F10.6,' Xhat ',3F7.4,' Yhat ',3F7.4,
     .           ' Zhat ',3F7.4,/,'SUN SVS: ',f10.4,3F16.4)
          old_tdt = tdt
      end if

*     See which type of sateliitet (valid Jan 1994).
      if( pn.eq.3 .or. pn.eq.12 .or. pn.eq.13 ) then
          sv_cm(1) = 0.2100d0
          sv_cm(2) = 0.d0
          sv_cm(3) = 0.8540d0
      else
          sv_cm(1) = 0.2794d0
          sv_cm(2) = 0.d0
          sv_cm(3) = 1.0259d0
      end if
 
*     Be really slack and only do radial on Block II
      do i = 1,3
C         svs_phs(i) = svs_xyz(i) + 1.0259d0*ehat(i)
          svs_phs(i) = sv_cm(1)*xhat(i) + sv_cm(2)*yhat(i) + 
     .                 sv_cm(3)*ehat(i) + svs_xyz(i) 
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

