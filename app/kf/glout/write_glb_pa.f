CTITLE WRITE_GLB_PARAMS
 
      subroutine write_glb_params( iout, options, cov_parm, sol_parm )

      implicit none  
 
*     Routine to write out the estimates of the parameters from the
*     global solution.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
*                                           ! To get nutation names
      include '../includes/globk_cmds.h'
      include '../includes/sd_common.h'
      include '../includes/const_param.h'
 
*   icnt,jcnt       - Counts number of components of site position
*                   - estimated
*   icnr            - Counter for number of rate components
*                   - estimated
*   iel, jel, kel   - Lookup positions in the solution vector
*                   - and covariance matrix
*   i,j,k           - loop counters
*   iout            - Output lu
*   options         - Options for the output
*   nt              - Number of tide stations
*   np              - Generic parameter number
*   num_amp         - Ampitude coefficient number orders by period
*                     and type 
 
      integer*4 icnt, jcnt, icnr, iel, i,j,k,n, iout, options, nt, np,
     .          num_amp, jel
 
*   covar(36)       - Six by six matrix needed for covariance
*                   - calculations
*   cov_parm(num_glb_parn,num_glb_parn) - Covariance matrix from
*                   - the solution
*   dt              - Time difference between site/source epoch
*                   - and current experiment
*   loc_coord(3)    - Local coordinates of the site (latitude,
*                   - longitude and radius)
*   NEU_covar(3,3)  - Covariance matrix for the NEU coordinates
*   pole_pos        - Temporary storage for pole positions
*   pos_axo         - Final estimate of the Axis offset or rate
*   pos_atm         - Final esimtate of atm delay
*   pos_xyz_adj(3)  - adjustment to XYZ position
*   pos_NEU_adj(3)  - adjustment to NEU position
*   pos_xyz_fin(3)  - Final XYZ position
*   pos_NEU_fin(3)  - Final NEU position
*   pos_radec(2)    - Final positions for RA and Dec
*   rat_radec(2)    - Final rates for RA and Dec
*   rat_xyz_adj(3)  - adjustment to XYZ rate
*   rat_NEU_adj(3)  - adjustment to NEU rate
*   rat_xyz_fin(3)  - Final XYZ rate
*   rat_NEU_fin(3)  - Final NEU rate
*   rot_matrix(3,3) - Rotation matrix between XYZ and NEU
*   scr_real(10)    - Scratch area
*   sol_parm(num_glb_parn)  - Solution vector from
*                   - the solution
*   temp_covar(36)  - Temporay storage for COMP_VAR
*   rat_xyz_sig(3)  - XYZ rate sigmas
 
      real*8 covar(36), cov_parm(num_glb_parn,num_glb_parn), dt, 
     .    loc_coord(3), NEU_covar(3,3), pole_pos, pos_axo, pos_atm,
     .    pos_xyz_adj(3), pos_NEU_adj(3), pos_xyz_fin(3),
     .    pos_NEU_fin(3), pos_radec(2), rat_radec(2), rat_xyz_adj(3),
     .    rat_NEU_adj(3), rat_xyz_fin(3), rat_NEU_fin(3),
     .    rot_matrix(3,3), scr_real(10), sol_parm(num_glb_parn),
     .    temp_covar(36), rat_xyz_sig(3)

* dsol_av(3) -- Change in position at non-correlated epoch due to
*     non-secular terms
      real*8 dsol_av(3)

* xyz_std(6) -- XYZ sigma (m) and corelations XY, XZ, YZ
* neu_std(6) -- NEU sigma (m) and corelations NE, NU, EU
* rxyz_std(6) -- Rate XYZ sigma (m) and corelations XY, XZ, YZ
* rneu_std(6) -- Rate NEU sigma (m) and corelations NE, NU, EU
* llu_std(6) -- Sigma Latitude, longitude (10^9 deg) and Height (m)

      real*8 xyz_std(6), neu_std(6), llu_std(6), 
     .       rxyz_std(6), rneu_std(6)

* Declarations for mapping to uncorrelated site position epcoh
* dtx(3), dtav - Uncorrelated epochs for XYZ and the average (years
*         from experiment epoch.
* jdav, decyrsav - Jd of avergage time and decimimal years.
* unc_fin(3) - Uncorrelated position value
* unc_adj(3) - Adjustment to position at uncorrelaed time
* unc_sig(3) - Sigma at uncorrelated time.
* conv       - Conversion factor

* log_fin(3) - Log final estimates
* log_adj(3) - Adjustment to apriori
* log_sig(3) - Sigmas of log estimates

      real*8 dtx(3), dtav, jdav, decyrsav, unc_fin(3), unc_adj(3),
     .       unc_sig(3), conv, log_fin(3), log_adj(3), log_sig(3)

* Orbit conversion:
* dXYZ_orb(6)  - Adjustments to XYZ XYZdot of orbit (m and mm/s)
* FXYZ_orb(6)  - Final orbit XYZ XYZdot with velocity in m/s
* aeinpa_rad(6) - Orbit a, e, i, RA of Node, Arg Perigee and Mean
*                 anomaly (angles in rads)
* dKdXYZ(6,6)   - Transformation from XYZ XYZdot to Keplerian elements
*                 (Scaled so that Keplerian elements are all in meters)
* aei_covar(6,6) - Covariance matrix from Keplerian elments in meters
* daei_orb(6)    - Adjustements to orbital elments in meters

      real*8 dXYZ_orb(6),FXYZ_orb(6), aeinpa_rad(6), dKdXYZ(6,6),
     .       aei_covar(6,6), daei_orb(6)

* aei_lab(6)  - Names of the orbital elements

      character*23 aei_lab(6)

* SVS file output variables

* yr, doy, secod, sec, hrs, mns - Variables for IC epoch format
* date(5)  - Yr, mon, day, hr, min for experiment epoch
      integer*4 yr, doy, secod, sec, hrs, mns, date(5), dats(5),
     .          datf(5)

* svs_fin(max_svs_elem) - Satellite elements
* sectag - Seconds tag of day

       real*8 svs_fin(max_svs_elem), sectag, secf, secs


* Nutation ampltiude values
* nut_amp_est(4) - Estimated nutation amplitudes
* nut_amp_adj(4) - Adjustments to nutation amplitudes
* nut_amp_sig(4) - Sigmas for the amplitudes
* nut_coeffs(4)   - Nutation series coefficients (when series used)

* SD coefficient values and arrays
* bs_ut_sd(4,max_ut1_names) - Band and +- designator for SD ut1 values
*                 - plus the cosine and sine signs          
* bs_xy_sd(4,max_xy_names ) Same type of arguuments for the polar motion
*                   terms

* sd_coeffs(2)    - Apriori coefficients for the UT1/XY diurnal
*                  semidiurnal terms

      integer*4 bs_ut_sd(4, max_ut1_names), bs_xy_sd(4, max_xy_names )

      real*8 nut_amp_est(4), nut_amp_adj(4), nut_amp_sig(4), 
     .       nut_coeffs(4), sd_coeffs(2)


* output      - Logical which indicates thats we should outpuy
*               the nutation amplitudes

      logical output, kbit, logout

* amp_types(4) - Types of nutation amplitudes

      character*8 amp_types(4)

* ntav   - Number of values in uncorrelated time estimate
* ne     - Earthquake number for log output

      integer*4 ntav, ne

* unc_label(3) - Labels for uncorrelated position estimates
* log_label(3) - Labels for log estimates

      character*14 unc_label(3), log_label(3)


*   etd_site_name       - Holds site name or global

      character*8 etd_site_name
 
*   st_parm_label(24)   - Labels for different types of paramters
*                       - which have a site or source label
*   orb_lab(max_svs_elem)          - Labels for orbital parameters
 
      character*14 st_parm_label(25), orb_lab(max_svs_elem)
 
*   nl_parm_label(27)   - Labels for parameters which do not have
*                       - a site or source name
*   temp_label          - Label used for nutation coefficients
 
 
      character*22 nl_parm_label(27), temp_label
      character*17 temp_short
      character*24 lab

*   eor_outl            - Label used for eor and etd parameters

      character*26 eor_outl
 
*   unit_label(10)      - Labels for the units
*   temp_unit           - Nutation coefficient label
*   orb_unit(max_svs_elem)         - UNits for orbital parameters
*   mp_lab(2,3)         - Labels for multi polar motion
 
      character*6 unit_label(10), temp_unit, orb_unit(max_svs_elem),
     .            mp_lab(2,3)

*   coeff_phase(2)      - Gives in or out of phase for coefficients
*   coeff_trig(2)       - Trig function for the coefficients.

      character*4 coeff_phase(2), coeff_trig(2)
 
*   nl                  - Format for skip a line
 
      character*4 nl

      character*12 full_ver, hsver

* MOD TAH 030730: Variables needed for GEOD and UTM output
      real*8 unc_geod(3)  ! Uncorrelated epoch Geod co-lat, lng and
                          ! and height (rads, rads, m)
      real*8 dgdt(3)      ! Geodetic rates of changes (deg/yr and m/yr)
      real*8 unc_utm(3)   ! UTM coordimates at uncorrelated time
      integer*4 zone      ! Zone for UTM coordinates
      character*1 hemi    ! Hemisphere for UTM coordinates
 
 
      common covar, temp_covar
 
      data st_parm_label / 'X coordinate  '
     .,                    'Y coordinate  '
     .,                    'Z coordinate  '
     .,                    'X rate        '
     .,                    'Y rate        '
     .,                    'Z rate        '
     .,                    'N coordinate  '
     .,                    'E coordinate  '
     .,                    'U coordinate  '
     .,                    'N rate        '
     .,                    'E rate        '
     .,                    'U rate        '
     .,                    'Axis offset   '
     .,                    'Axis rate     '
     .,                    'Right ascens  '
     .,                    'Declination   '
     .,                    'RA rate       '
     .,                    'Dec rate      '
     .,                    'Horiz. Love # '
     .,                    'Radial Love # '
     .,                    'Tidal Lag     '
     .,                    'X-axis        '
     .,                    'Y-axis        '
     .,                    'Z-axis        '
     .,                    'Zenith Delay  ' /

      data unc_label      /'X uncorr pos. ' 
     .,                    'Y uncorr pos. '
     .,                    'Z uncorr pos. ' /

      data  log_label     /'N Log         '
     .,                    'E Log         '
     .,                    'U Log         ' /
 
      data nl_parm_label / 'X-pole position       '
     .,                    'Y-pole position       '
     .,                    'X-pole rate           '
     .,                    'Y-pole rate           '
     .,                    'X-pole seas. offset   '
     .,                    'X-pole seas. rate     '
     .,                    'Y-pole seas. offset   '
     .,                    'Y-pole seas. rate     '
     .,                    'UT1-AT                '
     .,                    'UT1-AT rate           '
     .,                    'UT1-AT annual offset  '
     .,                    'UT1-AT annual rate    '
     .,                    'UT1-AT semian. offset '
     .,                    'UT1-AT semian. rate   '
     .,                    'Nutation in longitude '
     .,                    'Nutation in obliquity '
     .,                    'Long. random walk     '
     .,                    'Obliq. random walk    '
     .,                    'Long. seas. offset    '
     .,                    'Long. seas. rate      '
     .,                    'Obliq. seas. offset   '
     .,                    'Obliq. seas. rate     '
     .,                    'Gamma                 '
     .,                    'ETD In                '
     .,                    'ETD Out               '
     .,                    'NUT In                '
     .,                    'NUT Out               ' /

* MOD TAH 190610: Changed labels to accomadate ECOMC model.
C     data orb_lab       / 'Inert.  X     '
C    .,                    'Inert.  Y     '
C    .,                    'Inert.  Z     '
C    .,                    'Inert.  dX/dT '
C    .,                    'Inert.  dY/dT '
C    .,                    'Inert.  dZ/dT '
C    .,                    'Direct Rad    '
C    .,                    'Y Axis Bias   '
C    .,                    'Z Axis Bias   '
C    .,                    'B Axis Bias   '
C    .,                    'X Axis Bias   '
C    .,                    'Cos Direct    '
C    .,                    'Sin Direct    '
C    .,                    'Cos Y Bias    '
C    .,                    'Sin Y Bias    '
C    .,                    'Cos B Bias    '
C    .,                    'Sin B Bias    '
C    .,                    'Sin X1 Bias   '
C    .,                    'Sin X3 Bias   '
C    .,                    'Sin Z1 Bias   '
C    .,                    'AntOffest X   ' 
C    .,                    'AntOffest Y   '
C    .,                    'AntOffest Z   ' /

      data orb_lab       / 'Inert.  X     '
     .,                    'Inert.  Y     '
     .,                    'Inert.  Z     '
     .,                    'Inert.  dX/dT '
     .,                    'Inert.  dY/dT '
     .,                    'Inert.  dZ/dT '
     .,                    'Direct Rad    '
     .,                    'Y Axis Bias   '
     .,                    'B Axis Bias   '
     .,                    'Cos Direct    '
     .,                    'Sin Direct    '
     .,                    'Cos Y Bias    '
     .,                    'Sin Y Bias    '
     .,                    'Cos B Bias    '
     .,                    'Sin B Bias    '
     .,                    'Cos 2U Direct '
     .,                    'Sin 2U Direct '
     .,                    'Cos 4U Direct '
     .,                    'Sin 4U Direct '
     .,                    'UNKNOWN       '
     .,                    'AntOffest X   ' 
     .,                    'AntOffest Y   '
     .,                    'AntOffest Z   ' /


      data orb_unit      / '(m)   '
     .,                    '(m)   '
     .,                    '(m)   '
     .,                    '(mm/s)'
     .,                    '(mm/s)'
     .,                    '(mm/s)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(none)'
     .,                    '(m)   '
     .,                    '(m)   '
     .,                    '(m)   ' / 

      data aei_lab  /  'Semimajor axis      (m)',
     .                 'Eccentricity (none & m)',
     .                 'Inclination  (degs & m)',
     .                 'RA Node      (degs & m)',
     .                 'Arg. Perigee (degs & m)',
     .                 'M + w        (degs & m)'  /

      data coeff_phase   / 'IN  ', 'OUT ' /
      data coeff_trig    / 'Sin ', 'Cos'  /
 
      data unit_label    / '(m)   '
     .,                    '(m/yr)'
     .,                    '(mas) '
     .,                    '(ms)  '
     .,                    '(ms/y)'
     .,                    '(deg) '
     .,                    '(mm)  '
     .,                    '(uas) '
     .,                    '(uts) '
     .,                    '(ms/d)'    /

      data amp_types     / 'A+ in',
     .                     'A+ out',
     .                     'A- in',
     .                     'A- out'   /

      data bs_ut_sd  /   -1, -1, -1, -1,
     .                   -1,  1, -1, -1,
     .                   -2, -1,  1,  1,
     .                   -2,  1,  1,  1  /

      data bs_xy_sd  /    1, -1,  1, -1,
     .                    1,  1,  1, -1,
     .                   -2, -1, -1,  1,
     .                   -2,  1, -1,  1,
     .                    2, -1, -1,  1,
     .                    2,  1, -1,  1   / 

      data mp_lab    /  'X-pole','X-rate','Y-pole','Y-rate',
     .                  'UT1-AT','LOD   '  /
 
      data nl            / '(1x)'   /
 
* MOD TAH 030730: See if we are output GEOD or UTM coordinates
      full_ver = hsver(globk_version)
      if( kbit(options,25) ) then  ! GEOD coordinates to be output
         write(iout,50) full_ver
 50      format('GEOD* Geodetic (WGS84) coordinates from GLOBK ',
     .          'Vers ',a,/,
     .          'GEOD*  Site        Latitude        Longitude    ',
     .          'Height        dlat/dt       dlng/dt      dH/dt    ',
     .          'Epoch',/,
     .          'GEOD*               (deg)            (deg)       ',
     .          '(m)          (deg/yr)     (deg/yr)      (m/yr)    ',
     .          'Year')
      endif
      if( kbit(options,26) ) then  ! UTM coordinates to be output
         write(iout,60) full_ver
 60      format('UTM* UTM (WGS84) coordinates from GLOBK ',
     .          'Vers ',a,/,
     .          'UTM*  Site        Northing       Easting           ',
     .          'Height  Zn H      dN/dt     dE/dt     dH/dt     ',
     .          'Epoch',/,
     .          'UTM*                 (m)            (m)             ',
     .          '(m)              (m/yr)    (m/yr)    (m/yr)    Year')
      endif

***** Output the header line
      write(iout,100) full_ver 
  100 format( /,' PARAMETER ESTIMATES FROM GLOBK Vers ',a,/,
     .          '  #',t10,'PARAMETER',t47,'Estimate',t62,'Adjustment',
     .                t77,'Sigma')

      do j = 1,3
         rat_xyz_sig(j) = -1.0d0
      enddo
     
*     Site positions
      do i = 1, gnum_sites
*                             ! Indicates if any componets estimated
          neu_std = 0 
          rneu_std = 0
          icnt = 0
          dt   = (gepoch_out-site_epoch(i))/365.25d0
          if( kbit(guse_site,i) .and. site_epoch(i).ne. 0 )
     .    write(iout,120) gsite_names(i),(apr_val_site(j,1,i),j=1,3),
     .                    (apr_val_site(j,2,i),j=1,3), 
     .                    (site_epoch(i)- 2451544.5d0)/365.25d0+2000
          if( kbit(guse_site,i) .and. site_epoch(i).eq. 0 )
     .    write(iout,120) gsite_names(i),(apr_val_site(j,1,i),j=1,3),
     .                    (apr_val_site(j,2,i),j=1,3), 
     .                    0.00d0

 120      format('Int. ',a8,1x,3(f14.5,1x),3(f10.5,1x),f8.3,
     .                1x,3(f7.4,1x) )
 
*         Do values first
*                             ! Loop over XYZ
          do j = 1,3
*                                              ! estimated
* MOD TAH 150829: Initilize pos_xyz_fin so that if this no output
*             pbo line still has values.
              pos_xyz_fin(j) = apr_val_site(j,1,i)    +
     .                         apr_val_site(j,2,i)*dt +
     .                         cont_nonsec(j,i)
              if( parn_site(j,1,i).ne.0 ) then
                  iel  = parn_site(j,1,i)

*                 Save the sigma for output on the Unc. line
*                 incase rates are not estimated.
                  if( iel.gt.0 ) then
                      unc_sig(j) = sqrt(cov_parm(iel,iel))
                  else
                      unc_sig(j) = 0.d0
                  end if

                  icnt = icnt + 1
 
*                 Get final position of site
                  pos_xyz_adj(j) = sol_parm(iel)
                  pos_xyz_fin(j) = apr_val_site(j,1,i)    +
     .                             apr_val_site(j,2,i)*dt +
     .                             cont_nonsec(j,i)       +
     .                             pos_xyz_adj(j)
                  unc_fin(j) = pos_xyz_fin(j)

                  call out_glbp(iout,parn_site(j,1,i), pos_xyz_fin(j),
     .                 sol_parm, cov_parm, 1.d0, gsite_names(i),
     .                 st_parm_label(j), unit_label(1), 0)
*                             ! estimated
* MOD TAH 050927: Save sigmas for PBOP output
                  xyz_std(j) = sqrt(cov_parm(iel,iel))
                  if( j.eq.2 ) then
                      xyz_std(4) = cov_parm(iel-1,iel)/ 
     .                            (xyz_std(1)*xyz_std(2))
                  else if ( j.eq.3 ) then
                      xyz_std(5) = cov_parm(iel-2,iel)/ 
     .                            (xyz_std(1)*xyz_std(3))
                      xyz_std(6) = cov_parm(iel-1,iel)/ 
     .                            (xyz_std(2)*xyz_std(3))
                  endif
*                             ! estimated
              else
                  unc_fin(j) = apr_val_site(j,1,i)    +
     .                         apr_val_site(j,2,i)*dt +
     .                         cont_nonsec(j,i) 
              end if
*                             ! Looping over XYZ
          end do
 
*         Now do rates
*                             ! Rate counter
          icnr = 0
*                             ! Loop over XYZ dot
          do j = 1,3
*                                              ! estimated
              if( parn_site(j,2,i).ne.0 ) then
                  iel  = parn_site(j,2,i)
                  icnr = icnr + 1
 
*                 Get final position of site
                  rat_xyz_adj(j) = sol_parm(iel)
                  rat_xyz_fin(j) = apr_val_site(j,2,i) + rat_xyz_adj(j)
                  rat_xyz_sig(j) = sqrt(cov_parm(iel,iel))
 
                  call out_glbp(iout,parn_site(j,2,i), rat_xyz_fin(j),
     .                 sol_parm, cov_parm, 1.d0, gsite_names(i),
     .                 st_parm_label(3+j), unit_label(2), 0)
* MOD TAH 050927: Save sigmas for PBOP output
*                 Save the Rate XYZ sigmas
                  rxyz_std(j) = sqrt(cov_parm(iel,iel))
                  if( j.eq.2 ) then
                      rxyz_std(4) = cov_parm(iel-1,iel)/ 
     .                            (rxyz_std(1)*rxyz_std(2))
                  else if ( j.eq.3 ) then
                      rxyz_std(5) = cov_parm(iel-2,iel)/ 
     .                            (rxyz_std(1)*rxyz_std(3))
                      rxyz_std(6) = cov_parm(iel-1,iel)/ 
     .                            (rxyz_std(2)*rxyz_std(3))
                  endif
*                             ! estimated
              else
                  rat_xyz_fin(j) = apr_val_site(j,2,i)
                  call jd_to_decyrs(gepoch_out, decyrsav)
              end if
*                             ! Looping over XYZ dot
          end do

*         Now see if we should adjust the position epoch to
*         make the rate uncorrelated with the position.

*         first compute the average epoch shoft needed.
          dtav = 0.d0
          ntav = 0.d0
          do j = 1,3
             iel = parn_site(j,2,i)
             dtx(j) = 0.d0
             if( iel.gt.0 ) then
                 jel = parn_site(j,1,i)

*                Make sure that rate variance is not effectively zero.  This
*                way even those sites being forced to some rate will appear in
*                the output for when site positions are extracted.
                 if( cov_parm(iel,iel).gt.1.d-8 ) then
                     dtx(j) = -cov_parm(jel,iel)/cov_parm(iel,iel)
                 else
                     dtx(j) = 0.d0
                 end if
                 dtav   = dtav + dtx(j)
                 ntav   = ntav + 1
             end if
          end do

*         Now compute average epoch offset
          if( ntav.gt.0 ) then
              dtav = dtav/ntav
              jdav = gepoch_out + dtav*365.25d0
              call jd_to_decyrs( jdav, decyrsav)
              write(iout,400) gsite_names(i), decyrsav, dtx
 400          format(' Postion of ',a8,' referred to ',F9.4,
     .               '   XYZ offsets ', 3f9.4,' years')
              do j = 1,3
*                                              ! estimated
                  if( parn_site(j,2,i).ne.0 ) then
                      iel  = parn_site(j,1,i)
                      jel  = parn_site(j,2,i)
     
*                     Get final position of site
                      unc_adj(j) = sol_parm(iel) + 
     .                             sol_parm(jel)*dtav
* MOD TAH 050926: Add back in non-secular terms
                      unc_fin(j) = apr_val_site(j,1,i)    +
     .                             apr_val_site(j,2,i)*(dt+dtav) +
     .                             cont_nonsec(j,i) +
     .                             unc_adj(j)
*                     Compute the nonsecular terms at the uncoorelated
*                     epoch
* MOD TAH 050926: Do not need to do this calculation.  Done below for
*                 Apr lines only
                      call eval_nonsec(i, jdav, num_nonsec, 
     .                     param_nonsec, apr_val_nonsec, dsol_av, 0)
C                      unc_fin(j) = unc_fin(j) + dsol_av(j)

*                     Make sure that rate variance is not effectivelt zero
                      if( cov_parm(jel,jel).gt.1.d-8 ) then
                          unc_sig(j) = sqrt(cov_parm(iel,iel)+
     .                                2*cov_parm(iel,jel)*dtav + 
     .                                  cov_parm(jel,jel)*dtav**2)
                      else
                          unc_sig(j) = sqrt(cov_parm(iel,iel))
                      end if
 
                      call write_line( iout,-1, gsite_names(i),
     .                    unc_label(j), unit_label(1),
     .                    unc_fin(j), unc_adj(j), unc_sig(j),0)
*                             ! estimated
                  end if
               end do
          end if

*         Now write out the one line postion if the site was used
          if( kbit(guse_site,i) ) then
               write(iout,420) gsite_names(i), unc_fin,
     .               rat_xyz_fin, decyrsav, unc_sig
 420           format('Unc. ',a8,1x,3(f14.5,1x),3(f10.5,1x),f8.3,
     .                1x,3(f7.4,1x) )
* MOD TAH 050925: When non-secular terms have been used in the apriori
*              we should remove these from the output apiori values.
               write(iout,425) gsite_names(i), 
     .               (unc_fin(j)-cont_nonsec(j,i),j=1,3),
     .               rat_xyz_fin, decyrsav, unc_sig, rat_xyz_sig
 425           format('Apr. ',a8,1x,3(f14.5,1x),3(f10.5,1x),f8.3,
     .                1x,3(f7.4,1x), 1x, 3(f7.4,1x) )

          end if

****      See if we estimated log functions
          if( kbit(guse_site,i) ) then
               logout = .false.
*              See if apriori non-secular log term for this site
               call get_nonlog(i,apr_val_log(1,i))
               do j = 1,3
                  if( parn_log(j,i).ne.0 ) then
                      iel = parn_log(j,i)
                      log_adj(j) = sol_parm(iel)
                      log_fin(j) = apr_val_log(j,i) + log_adj(j)
                      log_sig(j) = sqrt(cov_parm(iel,iel))
                      call write_line( iout,iel, gsite_names(i),
     .                    log_label(j), unit_label(1),
     .                    log_fin(j), log_adj(j), log_sig(j),0)
                      logout = .true.
                  else
                      log_adj(j) = 0
                      log_fin(j) = apr_val_log(j,i)
                      log_sig(j) = 0
                  end if
               enddo
               if( logout ) then
*                  find the earthquake this log estimate is associated 
*                  with
                   do j = 1, num_eq
                      if( gsite_names(i)(7:8).eq.
     .                    eq_codes(j)(1:2) ) ne = j
                   end do
                   call jd_to_ymdhms(eq_epoch(ne),date,sectag)
                   write(iout,450) gsite_names(i),date, eq_log_tau(ne),
     .                  log_fin
 450               format('Apr. EXTENDED ',a8,' LOG ',i4,4i3,1x, F10.2,
     .                    3(1x,F10.5),' ! From GLOBK')
               end if

          end if
 
*****     Now see if we should do NEU components
*                                 ! Convert to NEU coordinates
          if( icnt.gt.0 ) then
              call rotate_geod(pos_xyz_adj, pos_neu_adj, 'XYZ','NEU',
     .                         pos_xyz_fin, loc_coord, rot_matrix)
 
*             Convert the local coordinates
              call loc_to_geod( loc_coord, pos_neu_fin )
 
*             Now compute the sigmas. Firstly get the elements of
*             covariance matrix we need
              call mov_cov(parn_site(1,1,i), parn_site(1,1,i), 3,
     .                     covar,3, cov_parm, num_glb_parn)
              call var_comp(rot_matrix, covar, NEU_covar, temp_covar,
     .                      3,3,1)
 
*             Now output the results
              neu_std = 0 
              rneu_std = 0
              do j = 1,3
                  scr_real(j) = sqrt(abs(NEU_covar(j,j)))
 
                  call write_line( iout,-1, gsite_names(i),
     .                st_parm_label(6+j), unit_label(1),
     .                pos_neu_fin(j), pos_neu_adj(j), scr_real(j),0)

* MOD TAH 050927: Save sigmas for PBOP output
                  neu_std(j) = scr_real(j)
                  if( j.eq.2 .and. neu_std(1)*neu_std(2).gt.0 ) then
                      neu_std(4) = NEU_covar(1,2)/ 
     .                            (neu_std(1)*neu_std(2))
                  elseif ( j.eq.3 .and. 
     .                neu_std(1)*neu_std(3).gt.0 ) then
                      neu_std(5) = NEU_covar(1,3)/ 
     .                            (neu_std(1)*neu_std(3))
                      neu_std(6) = NEU_covar(2,3)/ 
     .                            (neu_std(2)*neu_std(3))
                  endif

              end do
 
              call out_loc_correl(iout, NEU_covar, 3,3, NEU_covar,
     .             'NE,NU,EU position correlations')
 
*                         ! Any local coordinates to be output
          end if
 
*         Now do the rates in local coordinates
*                                 ! Convert to NEU rates
          if( icnr.gt.0 ) then
 
*             Convert adjustment to rates to NEU
              call rotate_geod(rat_xyz_adj, rat_neu_adj, 'XYZ','NEU',
     .                         pos_xyz_fin, loc_coord, rot_matrix)
 
*             Convert final values of rate to NEU
              do j = 1,3
                  call dvdot(rat_neu_fin(j), rot_matrix(j,1),3,
     .                       rat_xyz_fin,1, 3)
              end do
 
*             Now compute the sigmas. Firstly get the elements of
*             covariance matrix we need
              call mov_cov(parn_site(1,2,i), parn_site(1,2,i), 3,
     .                     covar,3, cov_parm, num_glb_parn)
              call var_comp(rot_matrix, covar, NEU_covar, temp_covar,
     .                      3,3,1)
 
*             Now output the results
              do j = 1,3
              scr_real(j) = sqrt(abs(NEU_covar(j,j)))
 
                  call write_line( iout,-1, gsite_names(i),
     .                st_parm_label(9+j), unit_label(2),
     .                rat_neu_fin(j), rat_neu_adj(j), scr_real(j),0)
* MOD TAH 050927: Save sigmas for PBOP output
                  rneu_std(j) = scr_real(j)
                  if( j.eq.2 .and. rneu_std(1)*rneu_std(2).gt.0 ) then
                      rneu_std(4) = NEU_covar(1,2)/ 
     .                            (rneu_std(1)*rneu_std(2))
                  elseif ( j.eq.3 .and. 
     .                rneu_std(1)*rneu_std(3).gt.0 ) then
                      rneu_std(5) = NEU_covar(1,3)/ 
     .                            (rneu_std(1)*rneu_std(3))
                      rneu_std(6) = NEU_covar(2,3)/ 
     .                            (rneu_std(2)*rneu_std(3))
                  endif

               end do
 
              call out_loc_correl(iout, NEU_covar, 3,3, NEU_covar,
     .             'NE,NU,EU rate correlations')

 
*                         ! Any local rates to be output
          end if

* MOD TAH 030730: See if we are outputing geod coordinates
          if( kbit(guse_site,i) ) then
               if( kbit(options,25) ) then
*                  Convert the unc_fin values to geodetic lat/long
                   call geod_to_geod(unc_fin, unc_geod, 'XYZ', 'GEOD',
     .                   'WGS84','WGS84',zone,hemi)
*                  Convert the rates to angular rates
                   dgdt(1) = rat_neu_fin(1)/Earth_rad*180/pi
                   dgdt(2) = rat_neu_fin(2)/Earth_rad*180/pi/
     .                       sin(unc_geod(1))
                   dgdt(3) = rat_neu_fin(3)
                   write(iout,460) gsite_names(i), 
     .                  (pi/2-unc_geod(1))*180/pi,unc_geod(2)*180/pi,
     .                  unc_geod(3), dgdt, decyrsav
 460               format('GEOD',2x,a8,2F15.9,F10.4,1x,2F15.9,F10.5,
     .                    1x,F10.4)
               endif
               if( kbit(options,26) ) then
*                  Convert the unc_fin values to UTM coordinates
                   call geod_to_geod(unc_fin, unc_utm, 'XYZ', 'UTM',
     .                   'WGS84','WGS84',zone,hemi)
                   write(iout,470) gsite_names(i), unc_utm, zone, hemi,
     .                  rat_neu_fin, decyrsav
 470               format('UTM',2x,a8,3F15.4,1x,i2.2,1x,a1,2x,
     .                     3F10.5,1x,F10.4)
               end if 

****           See if PBOP (position timeseries output)
               if( kbit(options,28) ) then
                   call jd_to_ymdhms(gepoch_out,date,sectag)
                   call geod_to_geod(pos_xyz_fin, unc_geod, 
     .                   'XYZ', 'GEOD','WGS84','WGS84',zone,hemi)
                   llu_std(1) = neu_std(1)/Earth_rad*180/pi*1.d9
                   llu_std(2) = neu_std(2)/Earth_rad*180/pi/
     .                       sin(unc_geod(1))*1.d9
                   llu_std(3) = neu_std(3)

*                  Write out the very long line
                   write(iout,480) gsite_names(i), gsite_full(i)(1:16),
     .                  date, gepoch_out-2400000.5d0, 
     .                 (pos_xyz_fin(j), j=1,3), (xyz_std(j),j=1,6),
     .                 (pi/2-unc_geod(1))*180/pi,unc_geod(2)*180/pi,
     .                  unc_geod(3), (llu_std(j),j=1,3),
     .                 (pos_neu_fin(j),j=1,3), (neu_std(j),j=1,6)
 480               format('pbo. ',a8,1x,a16,1x,i5,4(1x,i2.2),1x,F10.4,
     .                  1x,3F15.5,3F8.5,3F7.3,' | ',
     .                  2F16.10,1x,F10.5,1x,2F8.1,1x,F10.5,' | ',
     .                  2F15.5,1x,F10.5,3F8.5,3F7.3)
               endif

****           See if we should write the velocity line
               if( kbit(options,28) .and. icnr.gt.0 ) then
                   call mjd_to_ymdhms(gdata_st(i)+1.d-6,dats,secs)
                   call mjd_to_ymdhms(gdata_en(i)+1.d-6,datf,secf)
                   write(iout,490) gsite_names(i), gsite_full(i)(1:16),
     .                  date, gepoch_out-2400000.5d0, 
     .                 (pos_xyz_fin(j), j=1,3),
     .                 (pi/2-unc_geod(1))*180/pi,unc_geod(2)*180/pi,
     .                  unc_geod(3), 
     .                 (rat_xyz_fin(j),j=1,3), (rxyz_std(j),j=1,6),
     .                 (rat_neu_fin(j),j=1,3), (rneu_std(j),j=1,6),
     .                  dats,nint(secs), datf,nint(secf)
 490               format('pbr. ',a8,1x,a16,1x,i4,4i2.2,'00',1x,F10.4,
     .                  1x,3F15.5,2F16.10,1x,F10.5,2x,
     .                  3F9.5,3F8.5,3F7.3,2x, 3F9.5,3F8.5,3F7.3,2x,
     .                  i4,5i2.2,1x,i4,5i2.2)
               endif 

          end if
 
*         Now if we output any thing, write a line
          if( icnt.ne.0 .or. icnr.ne.0 ) write(iout,nl)
 
*                         ! Looping over the sites
      end do

* MOD TAH 040703: Now output atmopsheric delay
      icnt = 0
      do i = 1,gnum_sites
         iel = parn_atm(i)
*        			 ! Value estimated
         if( iel.gt.0 ) then
             icnt = icnt + 1
             pos_atm = apr_val_atm(i) + sol_parm(iel)
 
 
             call out_glbp( iout, parn_atm(i), pos_atm,
     .   	  sol_parm, cov_parm, 1.d0, gsite_names(i),
     .   	  st_parm_label(25), unit_label(1), 0)
         end if
      end do
      if( icnt.gt.0 ) write(iout,nl)

 
***** Axis offset
 
      icnt = 0
      do i = 1,gnum_sites
 
          dt = (gepoch_out-axo_epoch(i))/365.25d0
 
*                             ! Value and rate
          do j = 1,2
              iel = parn_axo(j,i)
*                                     ! Value estimated
              if( iel.gt.0 ) then
                  icnt = icnt + 1
                  pos_axo = apr_val_axo(j,i) + sol_parm(iel)
 
*                                     ! Add velocity effect as well
                  if( j.eq.1 ) then
                      pos_axo = pos_axo + apr_val_axo(2,i)*dt
                  end if
 
                  call out_glbp( iout, parn_axo(j,i), pos_axo,
     .                 sol_parm, cov_parm, 1.d0, gsite_names(i),
     .                 st_parm_label(12+j), unit_label(j), 0)
              end if
          end do
      end do
      if( icnt.gt.0 ) write(iout,nl)


***** Now do source positions
 
      do i = 1, gnum_sources
 
          dt = (gepoch_out-source_epoch(i))/365.25d0
 
*                             ! loop over RA and Dec
          do j = 1,2
*                                                ! estimated
              if( parn_source(j,1,i).ne.0 ) then
                  iel = parn_source(j,1,i)
                  pos_radec(j) = apr_val_source(j,1,i)    +
     .                           apr_val_source(j,2,i)*dt +
     .                           sol_parm(iel)
 
*                 If RA convert to mts
                  if( j.eq.1 ) then
                      pos_radec(j) = pos_radec(j)/15.d0
                  end if
 
                  call out_glbp( iout, parn_source(j,1,i), pos_radec(j),
     .                 sol_parm, cov_parm, 1.d0, gsource_names(i),
     .                 st_parm_label(14+j),unit_label(3),1)
 
                  if( j.eq.2 ) write(iout,nl)
              end if
          end do
      end do
 
*     Now do rates
 
      do i = 1, gnum_sources
*                             ! loop over RA and Dec rates
          do j = 1,2
*                                                ! estimated
              if( parn_source(j,2,i).ne.0 ) then
                  iel = parn_source(j,2,i)
                  rat_radec(j) = apr_val_source(j,2,i) + sol_parm(iel)
*                 If RA convert to mts
                  if( j.eq.1 ) then
                      rat_radec(j) = rat_radec(j)/15.d0
                  end if
 
                  call out_glbp( iout, parn_source(j,2,i), rat_radec(j),
     .                 sol_parm, cov_parm, 1.d0, gsource_names(i),
     .                 st_parm_label(16+j),unit_label(5),1)
 
                  if( j.eq.2 ) write(iout,nl)
              end if
          end do
      end do
 
***** Now do the satillite orbit parameters.  
* MOD TAH 971012: Write out the header line with IC Epoch and Type
      if( gnum_svs.gt.0 ) then
          call jd_to_yds( svs_epoch(1)+1.d-9, yr, doy, secod)
          if( yr.gt.500) yr = yr-1900
          hrs = secod/3600
          mns = (secod-hrs*3600)/60
          sec = (secod-hrs*3600-mns*60)
          write(iout,500) yr, doy, hrs, mns, sec, ggtime, ggframe,
     .                    ggprec, ggsrpmod
 500      format('Eph. #IC ',i2,1x,i3,1x,i2,1x,i2,1x,i2,20x,a4,1x,a5,1x,
     .           a5,1x,a5)

*         Convert to Solution time to yr,day,mon, and hour 
*         (Use the last experiment date)
          call jd_to_ymdhms( gepoch_expt, date, sectag)
          
      end if 


      do i = 1, gnum_svs
         icnt = 0
         do j = 1,max_svs_elem
             np = parn_svs(j,i)

             if( np.eq.0 ) svs_fin(j) = apr_val_svs(j,i)

             if( np.gt.0 ) then
                 scr_real(1) = apr_val_svs(j,i) + sol_parm(np)
                 svs_fin(j) = scr_real(1)

*                Save main part of orbital elements (ie. XYZ and
*                XYZdot.  For the velocity terms convert to m/s)
                 if( j.le.6 ) then 
                     dXYZ_orb(j) = sol_parm(np)
                     if( j.le.3 ) then
                         FXYZ_orb(j) = scr_real(1)
                     else
                         FXYZ_orb(j) = scr_real(1)/1000.d0
                     end if
                 end if
*                
                 call out_glbp(iout, np, scr_real(1), sol_parm,
     .               cov_parm, 1.d0, gsvs_names(i), orb_lab(j),
     .               orb_unit(j), 0 )
                 icnt = icnt + 1
             else
                 if( j.le.6 ) then
                     dXYZ_orb(j) = 0.d0
                     FXYZ_orb(j) = apr_val_svs(j,i)
                     if( j.gt.3 ) FXYZ_orb(j) = FXYZ_orb(j)/1000.d0
                 end if                 
             end if
         end do

*        Now output one-line svs_file format line
         if( icnt.gt.0 ) then
             write(iout,505) (date(j),j=1,4), gsvs_names(i), 
     .                       (svs_fin(j),j=1,max_svs_elem)
 505         format('Eph. ',i5,3i3,1x,a8,3(1x,F14.4),3(1x,F14.5),
     .              17(1x,F9.5))
         end if

****     Now convert the adjustments and cartesian elements to 
*        Keplerian elememts.
         if( icnt.gt.0 .and. abs(FXYZ_orb(1)).gt.1.d0 ) then
             call elem(FXYZ_orb(1), FXYZ_orb(2), FXYZ_orb(3),
     .                 FXYZ_orb(4), FXYZ_orb(5), FXYZ_orb(6),
     .                 aeinpa_rad, 3, dKdXYZ, GM_earth, 0) 

*            Add the argument of perigee to the Mean Anomaly
             aeinpa_rad(6) = aeinpa_rad(5) + aeinpa_rad(6)
             aeinpa_rad(6) = mod(aeinpa_rad(6), 2*pi)

*            Convert the angular arguments to meters (multiplu by
*            a and convert the velocities to be mm/s)
             do k = 4,6
                 dKdXYZ(1,k) = dKdXYZ(1,k) / 1000.d0
             end do
             do j = 2,6
                do k = 1,3
                   dKdXYZ(j,k) = dKdXYZ(j,k)*aeinpa_rad(1)
                end do
                do k = 4,6
                   dKdXYZ(j,k) = dKdXYZ(j,k)*aeinpa_rad(1)/1000.d0
                end do
             end do

*            Finally add the Mean Anomaly to the argument of 
*            perigee
             do j = 1,6
                dKdXYZ(6,j) = dKdXYZ(5,j) + dKdXYZ(6,j)
             end do

*            Now get the changes in the elements and compute the 
*            covariance matrix
             call mov_cov(parn_svs(1,i), parn_svs(1,i), 6, 
     .                    covar,  6, cov_parm, num_glb_parn)
             call var_comp(dKdXYZ, covar, aei_covar, temp_covar,
     .                     6,6,1)

*            Now compute the adjustments to the values
             do j = 1,6
                call dvdot(daei_orb(j), dKdXYZ(j,1),6, 
     .                     dXYZ_orb,1,6) 
             end do

*            Now write results
             do j = 1,1
                write(iout,511) gsvs_names(i), aei_lab(j), 
     .                aeinpa_rad(j), daei_orb(j), 
     .                sqrt(aei_covar(j,j))
 511            format(' Loc. ',a,1x,a,t39,F18.4,2(1x,f12.4))
             end do
            
             do j = 2,2
                write(iout,512) gsvs_names(i), aei_lab(j), 
     .                aeinpa_rad(j), daei_orb(j), 
     .                sqrt(aei_covar(j,j))
 512            format(' Loc. ',a,1x,a,t39,F18.10,2(1x,f12.4))
             end do
             do j = 3,6 
                write(iout,513) gsvs_names(i), aei_lab(j), 
     .                aeinpa_rad(j)*180.d0/pi, daei_orb(j), 
     .                sqrt(aei_covar(j,j))
 513            format(' Loc. ',a,1x,a,t39,F18.10,2(1x,f12.4))
             end do

         end if

         if( icnt.ne.0 ) write(iout, nl)
      end do

* MOD TAH 131223: Moved translations and rotations to be near PMU estimates
****  Translation parameter
      icnt = 0
      do i = 1,3
         if(  parn_tran(i,1).gt.0 ) then 
             icnt = icnt + 1
             call out_glbp(iout, parn_tran(i,1), 0.0d0, 
     .                     sol_parm, cov_parm,
     .                     1.d0, 'TRANSLTN', st_parm_label(i), 
     .                     unit_label(1),0 )
         end if
      end do
      do i = 1,3
         if(  parn_tran(i,2).gt.0 ) then 
             icnt = icnt + 1
             call out_glbp(iout, parn_tran(i,2), 0.0d0, 
     .                     sol_parm, cov_parm,
     .                     1.d0, 'TRANRATE', st_parm_label(i), 
     .                     unit_label(2),0 )
         end if
      end do

****  Rotation parameter
      icnt = 0
      do i = 1,3
         if(  parn_rot(i,1).gt.0 ) then 
             icnt = icnt + 1
             call out_glbp(iout, parn_rot(i,1), 0.0d0, 
     .                     sol_parm, cov_parm,
     .                     1.d0, 'ROTATION', st_parm_label(i+21), 
     .                     unit_label(3),0 )
         end if
      end do
      do i = 1,3
         if(  parn_rot(i,2).gt.0 ) then 
             icnt = icnt + 1
             call out_glbp(iout, parn_rot(i,2), 0.0d0, 
     .                     sol_parm, cov_parm,
     .                     1.d0, 'ROT_RATE', st_parm_label(i+21), 
     .                     unit_label(5),0 )
         end if
      end do

****  Scale parameters
      if( parn_scale(1).gt.0 ) then
          icnt = icnt + 1
          call out_glbp(iout, parn_scale(1), 0.0d0, 
     .                  sol_parm, cov_parm,
     .                  1.d0, 'SCALE   ', '(ppb)', ' ', 0 )
      end if
      if( parn_scale(2).gt.0 ) then
          icnt = icnt + 1
          call out_glbp(iout, parn_scale(2), 0.0d0, 
     .                  sol_parm, cov_parm,
     .                  1.d0, 'SCALRATE', '(ppb/yr)', ' ', 0 )
      end if
      if( icnt.gt.0 ) write(iout,nl)

***** Now do polar motion and UT1
* MOD TAH 070823: If we are outputing at mid-point convert the
*     dailty polar motion/UT1 to multiday with the correct 
*     epoch (this needs to be the last value)
      if( kbit(options,29) ) then
          call pmu_midpt
      end if

*     Continue with polar motion/ut1.
      icnt = 0
* MOD TAH 200806: Updated comment
*                     ! Loop over the 4 compoent model
      do j = 1,4
          if( parn_wob(j).gt.0 ) then
              icnt = icnt + 1
              pole_pos = apr_val_wob(j) + gwob_apr(j) +
     .                                    sol_parm(parn_wob(j))
              n = 3
*             Change label for rate terms
              if( j.eq.3 .or. j.eq.4 ) n = 10 
              call out_glbp(iout, parn_wob(j), pole_pos,
     .             sol_parm, cov_parm, 1.d0, nl_parm_label(j),
     .             unit_label(n), ' ',0)
          end if
      end do
 
      if( icnt.gt.0 ) write(iout,nl)
 
      icnt = 0
*                     ! Loop over the 2 compoent model
      do j = 1,2
          if( parn_ut1(j).gt.0 ) then
              icnt = icnt + 1
              pole_pos = (apr_val_ut1(j) + gut1_apr(j)
     .                                   + sol_parm(parn_ut1(j)))/15.d0
              n = 4
              if( j.eq.2 ) n = 10 
              call out_glbp(iout, parn_ut1(j), pole_pos,
     .             sol_parm, cov_parm, 1/15.d0, nl_parm_label(8+j),
     .             unit_label(n), ' ',0)
          end if
      end do

****  Now write out the correlations
      call write_pmu_corr( iout, cov_parm, num_glb_parn, parn_wob(1),
     .                     parn_wob(2), parn_ut1(1) )
      if( icnt.gt.0 ) write(iout,nl)

****  Now write out the multi-epoch polar motion/UT1 values
*     Do X/Y Pole first
      do j = 1, 2
         do i = 1, 2
            icnt = 0
            do k = 1, num_mul_pmu
               if(  parn_mul_pmu(i,j,k).ne.0 ) then
                    icnt = icnt + 1
                    pole_pos = apr_val_mul_pmu(i,j,k)+
     .                         sol_parm(parn_mul_pmu(i,j,k))
                    n = 3
                    if( i.eq.2 ) n = 10

                    call jd_to_ymdhms(gmul_pmu_ep(k)+1.d-4, 
     .                                date, sectag)                 
                    write(lab,600) mp_lab(i,j), date
 600                format(a6,1x,I4,'/',I2.2,'/',I2.2,1x,I2.2,':',I2.2)          
                    call out_glbp(iout, parn_mul_pmu(i,j,k), pole_pos,
     .                  sol_parm, cov_parm, 1.d0, lab,
     .                  unit_label(n), ' ',0)
               endif
            end do
            if( icnt.gt.0 .and. num_mul_pmu.gt.1 ) write(iout,nl)
         end do
      end do 

*     Now do the UT1 and LOD values 
      do j = 3, 3
         do i = 1, 2
            icnt = 0
            do k = 1, num_mul_pmu
               if(  parn_mul_pmu(i,j,k).ne.0 ) then
                    icnt = icnt + 1
                    pole_pos = apr_val_mul_pmu(i,j,k)+
     .                         sol_parm(parn_mul_pmu(i,j,k))
                    conv = 1.d0/15.d0
                    if( i.eq.2 ) conv = -conv
                    pole_pos = pole_pos*conv
                    n = 4
                    call jd_to_ymdhms(gmul_pmu_ep(k)+1.d-4, 
     .                                date, sectag)                 
                    write(lab,600) mp_lab(i,j), date
                    call out_glbp(iout, parn_mul_pmu(i,j,k), pole_pos,
     .                  sol_parm, cov_parm, conv, lab,
     .                  unit_label(n), ' ',0)
               endif
            end do
            if( icnt.gt.0 .and. num_mul_pmu.gt.1 ) write(iout,nl)
         end do
      end do 

      if( icnt.gt.0 ) write(iout,nl)
***** Nutation angles
      icnt = 0
*                     ! Loop over the 8 compoent model
      do j = 1,8
          if( parn_nut_ang(j).gt.0 ) then
              icnt = icnt + 1

              pole_pos = apr_val_nut_ang(j) + gnut_ang_apr(j)
     .                                      + sol_parm(parn_nut_ang(j))
 
              call out_glbp(iout, parn_nut_ang(j), pole_pos,
     .             sol_parm, cov_parm, 1.d0, nl_parm_label(14+j),
     .             unit_label(3), ' ',0)
          end if
      end do
      if( icnt.gt.0 ) write(iout,nl)
      call write_iers( iout, cov_parm, sol_parm )

***** Extended Earth tide coefficicents
      if( glb_glb_tides ) then
          nt = 1
      else 
          nt = gnum_sites
      end if

      do i = 1, nt
         do j = 1, max_etd_names
             np = parn_eor_etd(j,i)
             if( np.ne.0 ) then
                 eor_outl = etd_names(j)
                 call sub_char( eor_outl, 'A+', 'COS' )
                 call sub_char( eor_outl, 'A-', 'SIN' )
                 if( glb_glb_tides ) then
                     call out_glbp(iout, np, 0.d0, sol_parm,
     .                    cov_parm, 1.d3, 'GLOBAL', eor_outl(1:14),
     .                    unit_label(7),0)
                 else
                     call out_glbp(iout, np, 0.d0, sol_parm,
     .                    cov_parm, 1.d3, gsite_names(i), 
     .                    eor_outl(1:14), unit_label(7),0)
                 end if
             end if
         end do
      end do

*     UT1 angle estimates 
      do i = 1, max_ut1_names
          np = parn_eor_ut1(i)
          if( np.ne.0 ) then
              eor_outl = ut1_names(i)
              call sub_char( eor_outl, 'A+', 'COS' )
              call sub_char( eor_outl, 'A-', 'SIN' )
              pole_pos = (sd_ut1_mid(i) + sol_parm(np))/15.d-3
              call out_glbp(iout, np, pole_pos, sol_parm,
     .             cov_parm, 1/15.d-3, eor_outl(1:14), '        ',
     .             unit_label(9),0)
          end if
      end do

*     XY angle estimates
      do i = 1, max_xy_names
          np = parn_eor_xy(i)
          if( np.ne.0 ) then
              eor_outl = xy_names(i)
              call sub_char( eor_outl, 'A+', 'COS' )
              call sub_char( eor_outl, 'A-', 'SIN' )
              pole_pos = (sd_xy_mid(i) + sol_parm(np))*1000.d0
              call out_glbp(iout, np, pole_pos, sol_parm,
     .             cov_parm, 1.d3, eor_outl(1:14), '        ',
     .             unit_label(8),0)
          end if
      end do

***** Gamma

      if( parn_gamma.ne.0 ) then
          write(iout,nl)
          scr_real(1) = apr_val_gamma + sol_parm(parn_gamma)

          call out_glbp(iout, parn_gamma, scr_real, sol_parm, cov_parm,
     .         1.d0, nl_parm_label(23),' ',' ',0)
      end if

***** Nutation series coefficients
 
      icnt = 0
      do i = 1,4
         nut_coeffs(i) = 0.0
      end do
      do i = 1, max_nut_coeff
 
*         Get any apririo series coefficents
          if ( i.gt.4) call get_nut_coeff((i+1)/2, nut_coeffs)
 
*                       ! In and Out of phase
          do j = 1,2
 
              temp_label = nl_parm_label(25+j)(1:10) // nut_names(i)
              if ( j.eq.1 .and. i.le.2 ) then
                  temp_unit = '(ms/y)'
              else
                  temp_unit = '(mas)'
              end if

              iel = parn_nut_coeff(j,i)
*                                                            ! Estimated
              if( iel.gt.0 ) then
*                                                            ! estimate
                  pole_pos = apr_val_nut_coeff(j,i) +
     .                       sol_parm(iel)
*                                                ! Set estimate to zero
*                 Add the series coefficients
                  if( i-(i/2)*2.eq.0 ) then
                     pole_pos = pole_pos + nut_coeffs(2+j)
                  else
                     pole_pos = pole_pos + nut_coeffs(j)
                  end if
 
                  call out_glbp(iout, parn_nut_coeff(j,i), pole_pos,
     .                 sol_parm, cov_parm, 1.d0, temp_label,temp_unit,
     .                 ' ',0)
              end if
          end do
      end do

      write(iout,nl)

****  Now do the ampltitudes 
      icnt = 0
*     Start with the FCN mode (skip offsets and rate), Loop two
*     at a time to get all long and obliquityA
      do  i = 1,4
         nut_coeffs(i) = 0.0
      end do

      do i = 3, max_nut_coeff, 2
          if( i.gt.4 ) call get_nut_coeff((i+1)/2, nut_coeffs)
          call comp_nut_amp(parn_nut_coeff(1,i), nut_coeffs, 
     .         apr_val_nut_coeff(1,i),
     .         cov_parm, sol_parm, num_glb_parn, nut_amp_est,
     .         nut_amp_adj, nut_amp_sig, output )

*         See if values to output
          if( output ) then
              icnt = icnt + 1
              do j = 1,4
                 temp_short = amp_types(j) // coeff_periods((i+1)/2)

                 call write_line( iout,-1, 'NUT AMP ',
     .                 temp_short, temp_unit, nut_amp_est(j),
     .                 nut_amp_adj(j), nut_amp_sig(j), 0 )
              end do
          end if
      end do

      if( icnt.ne.0 ) write(iout,nl)

*     Now do the extended Earth parameter amplitudes
*     Loop over sites (or just 1 if global)
      do n = 1, min(nt,max_etd_sites)
         if( glb_glb_tides ) then
             etd_site_name = 'GLOBAL'
         else
             etd_site_name = gsite_names(n)
         end if

         jcnt = 0

*        Loop over the different types of coefficients
         do i = 1, max_etd_names
            icnt = 0

*           Now loop over the periods
            do j = 1, max_coeff_periods
               num_amp = (j-1)*max_etd_names + i

*              Loop over in and outof phase
               do k = 1,2
                  np = parn_etd_coeff(k,num_amp,n)
                  if( np.gt.0 ) then
                     if( icnt.eq.0 ) then
                         write(iout,710) etd_names(i), etd_site_name
 710                     format(5x,'Extended Earth Tides ',a,1x,a)
                         icnt = icnt + 1
                     end if
                     call out_glbp(iout, np,
     .                   apr_val_etd_coeff(k,num_amp,n), sol_parm,
     .                   cov_parm, 1000.d0, coeff_periods(j),
     .                      coeff_phase(k), unit_label(7),0)
                  end if
               end do
               jcnt = jcnt + icnt
            end do
          end do 
      end do

***** Now do the UT1 coeffiecients
      do i = 1, max_ut1_names
         icnt = 0
         do j = 1, max_coeff_periods
             num_amp = (j-1)*max_ut1_names + i

*            get the apriori coefficient used.
             call get_sd_coeff(j, bs_ut_sd(1,i), sd_ut1_val, sd_ut1_arg,
     .                         sd_ut1_num, sd_coeffs )

             do k = 1,2
                np = parn_ut1_coeff(k,num_amp)
                if( np.ne.0 ) then
                    eor_outl = ut1_names(i)(1:12) // coeff_periods(j)

                    pole_pos = (sd_coeffs(k)+
     .                          sol_parm(np)*bs_ut_sd(k+2,i))/15.d-3
                    call out_glbp(iout, np, pole_pos, sol_parm,
     .                  cov_parm, bs_ut_sd(k+2,i)*1.d3/15.d0,
     .                  eor_outl(1:20), 
     .                  coeff_trig(k), unit_label(9), 0 )
                    icnt = icnt + 1
                end if
             end do
          end do 
          if( icnt.gt.0 ) write(iout,nl)
      end do

***** Now the XY polar position coefficients
      do i = 1, max_xy_names
         icnt = 0
         do j = 1, max_coeff_periods
             num_amp = (j-1)*max_xy_names + i

*            get the apriori coefficient used.
             call get_sd_coeff(j, bs_xy_sd(1,i), sd_xy_val, sd_xy_arg,
     .                         sd_xy_num, sd_coeffs )
             do k = 1,2
                np = parn_xy_coeff(k,num_amp)
                if( np.ne.0 ) then
                    eor_outl = xy_names(i)(1:15) // coeff_periods(j)
                    
                    pole_pos = (sd_coeffs(k)+
     .                          sol_parm(np)*bs_xy_sd(k+2,i))*1.d3  
                    call out_glbp(iout, np, pole_pos, sol_parm,
     .                  cov_parm, bs_xy_sd(k+2,i)*1.d3, eor_outl(1:21), 
     .                  coeff_trig(k), unit_label(8), 0 )
                    icnt = icnt + 1
                end if
             end do
          end do
          if( icnt.gt.0 ) write(iout,nl)
      end do

****  Write the out the translation/rate correlations
C     do i = 1, 3
C        np = parn_tran(i,1)
C        if( np.gt.0 ) then
C            if( i.eq.1 ) then
C                write(iout,810)
C810             format(/,'MAXIMUM Correlations',/,
C    .                    ' NP    Type  Comp   N    Rho  ....')
C            end if 

C            call max_corr(iout,np, 'TRANSLTN',st_parm_label(i),
C    .                     cov_parm)
C        end if
C     end do

C     do k = 1, num_mul_pmu
C         np = parn_mul_pmu(2,3,k)
C         if( np.gt.0 ) then
C            call max_corr(iout,np, 'LOD     ',' ',
C    .                     cov_parm)
C         end if
C     enddo
              
 
***** Thats all
      return
      end

CTITLE MAX_CORR 

      subroutine max_corr(iout,np, label1,label2 , cov_parm)

      implicit none 

*     Routine to look for max correlations and write them out

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
*                                           ! To get nutation names
      include '../includes/globk_cmds.h'
      include '../includes/sd_common.h'
      include '../includes/const_param.h' 

* PASSED PARAMETERS
* np  -- Parameter number to be tested
* iout -- Output unit number
      integer*4 np, iout

* cova_parm(num_glb_parn,num_glb_parn) -- Covariance matrix
      real*8 cov_parm(num_glb_parn,num_glb_parn)

* label1, label2 -- Parameter labels

      character*(*) label1, label2

* LOCAL VARIABLES
      integer*4 max_cs   ! Maximum number of correlations to list
      parameter ( max_cs = 10 )

      real*8 max_correl(max_cs)   ! Maximum correlations
      integer*4 n_max(max_cs)     ! parameter numbers with max

      integer*4 i, j, k  ! Loop counters
      logical done
      real*8 rho      ! Current correlations

*     Clear correlations
      do i = 1, max_cs
         max_correl(i) = 0.d0
         n_max(i) = 0
      end do

****  Now loop over all parameters
      do i = 1, num_glb_parn
         if( i.ne.np ) then
            rho = cov_parm(i,np)/sqrt(cov_parm(i,i)*cov_parm(np,np))
*           See where this fits in list
            j = max_cs+1
            done = .false.
            do while ( j.gt.2 .and. .not. done )
               j = j - 1
               if( abs(rho).gt.abs(max_correl(j)) ) then
*                  Push this correlation onto list
                   do k = 2,j
                      max_correl(k-1) = max_correl(k)
                      n_max(k-1) = n_max(k)
                   end do
                   max_correl(j) = rho
                   n_max(j) = i
                   done = .true.
               end if
            end do
          end if
      end do

*     Now write the results
      write(iout,140 ) np, label1, label2(1:1), 
     .      (n_max(j), max_correl(j),j=1,max_cs)
 140  format(I5,1x,a,1x,a,10(1x,i4,1x,F6.3))

      return
      end

CTITLE GET_LONG_NAMES

      subroutine get_long_names

      implicit none 
      
*     Routine to read head.snx and extra the long names for
*     sites

      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'

 
* unitc  -- Comments file unit numbers

      integer*4 unitc
      
* LOCAL VARIABLES      
*   ierr        - IOSTAT error
*   trimlen     - Length of string
 
      integer*4 ierr, trimlen, ns, i

 
*   start_found - Set true when start line found
*   stop_found  - Set true when stop line found
 
      logical start_found, stop_found
 
*   line        - Line read from comments
 
      character*128  line
            
*   domes       - DOMES number for site (9 characters needed)
*   type        - P=GPS, R=VLBI, S=SLR
*   pt          - Point number.  This is added to the end of the four character
*                 code.  Later if it is not A we will make _GP[Pt]
*   globk_extent - Last four characters of globk name

      character*12 domes
      character*4 type, pt, globk_extent, code
      character*8 full_name
* MOD TAH 050927: Added long name to be read from head.snx
      character*22 long_name
      character*128 snx_comfile 

      data unitc / 110 /
      
****  Check to see if comment unit
      snx_comfile = 'head.snx'
      open(unitc, file=snx_comfile, iostat=ierr, status='old')

      if( ierr.ne.0 ) then

*         Try to open version in HELP_DIR (use line for name)
          call getenv('HELP_DIR',line)
          line(trimlen(line)+1:) = '/' // snx_comfile
          open(unitc, file=line, iostat=ierr, status='old')
      endif
      if( ierr.ne.0 ) then
         write(*,110) snx_comfile(1:trimlen(snx_comfile))
 110     format('**WARNING** Unable to open ',a,
     .          ' NO UPDATED LONG SITE NAMES')
         do i = 1,gnum_sites
             call sub_null(gsite_full(i))
         end do
         RETURN
      endif

***** If file opened OK, read through names and assign to sites

      line = ' '
      ierr = 0
      start_found = .false.
      stop_found = .false.
      do while ( ierr.eq.0 .and. .not.stop_found )
 
          read(unitc, '(a)', iostat=ierr) line
          if( index(line,'-DOMES').eq.1 ) stop_found = .true.
          if( ierr.eq.0 .and. .not.stop_found .and. start_found .and.
     .        trimlen(line).gt.0 .and. line(1:1).eq.' ' ) then
              
              read(line, 120) code, pt, domes, type, globk_extent,
     .                        long_name
 120          format(1x,a4,1x,a2,1x,a9,1x,a1,1x,a4,1x,a16)
              full_name = code(1:4) // globk_extent

*             Match the 4-char code
              do ns = 1, gnum_sites
                  if( code(1:4).eq.gsite_names(ns)(1:4) ) then
                     if( trimlen(long_name).gt.0 ) then
                         gsite_full(ns)(1:22) = long_name
                     endif
                     gsite_full(ns)(23:31) = domes(1:9)
                     gsite_full(ns)(32:32) = type
                  end if
              end do
          end if
*         Check for start afterwards so that the start line
*         is not echoed.
          if( index(line,'+DOMES').eq.1 ) start_found = .true.
      end do

****  Now make sure none of the long names are nulls or blanks
      do ns = 1, gnum_sites
        if( trimlen(gsite_full(ns)).eq.0 ) then
            gsite_full(ns) = gsite_names(ns)
        endif
      end do

      close(unitc)

****  Thats all
      return 
      end

            
    
 
