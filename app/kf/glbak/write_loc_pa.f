CTITLE WRITE_LOC_PAR
 
      subroutine write_loc_par( iout, options, type, parn, ndimp, ltog,
     .                          gnum, cov_parm, sol_parm,  ndimc,
     .                          obs_corr, ndimo )

      implicit none  
 
*     Routine to write out the estimates of the parameters from the
*     global back solution solution.
*     The user passes the type of parameter that the current parn
*     array is for, the parn array (which looks like the standard
*     parn arrays for the global solutions), the dimension of parn
*     array, and the number of values to check in parn.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/const_param.h'
 
*   gnum            - Number of values to check. eg. number of sites
*                   - or number of sources
*   icnt            - Counts number of components of site position
*                   - estimated
*   icnr            - Counter for number of rate components
*                   - estimated
*   ident           - Indentifies use of scratch common (not used)
*   iel, jel, kel   - Lookup positions in the solution vector
*                   - and covariance matrix
*   i,j,k           - loop counters
*   iout            - Output lu
*   is              - Local site or source number (found from the
*                   - ltog_ arrays.
*   ltog(1)         - Allows mappning of parn index to site or
*                   - source number.  If first element is -1 then
*                   - return use index number
*   ltog_map        - integer*4 function to convert ltog to site
*                   - or source number.  If first element of ltog
*                   - is -1 then, index is used
*   ndimc           - Dimension of the covariance matrix passed
*   ndimp           - First dimension of the parn array (also
*                   - number values to check over)
*   ndimo           - Dimension of the obs_corr vector.  If the
*                   - dimension is zero then these corrections
*                   - are not aplied.
*   options         - Options for the output
*   parn(ndimp,1)   - Parameter numbers for current output
*   type            - start type of parameters to br output
*                   - eg. site positions would be 7, nutation
*                   - angles 16 etc.
 
      integer*4 gnum, icnt, icnr, ident, iel, i,j, k, iout,
     .    is, ltog(*), ltog_map, ndimc, ndimp, ndimo, options,
     .    parn(ndimp,1), type, dumglbak, date(5), n
 
*   wgh             - Used for weight in NEU statistics accumultions
 
      real*4 wgh
 
*   covar(36)       - Six by six matrix needed for covariance
*                   - calculations
*   cov_parm(ndimc,ndimc) - Covariance matrix from
*                   - the solution
*   dt              - Difference in time between current experiment
*                   - and reference time for parameters which
*                   - have rates (years)
*   loc_coord(3)    - Local coordinates of the site (latitude,
*                   - longitude and radius)
*   obs_corr(1)     - Correction to be applied to the observed
*                   - values to remove translations and rotations
*                   - of the coordinate system
*   pole_pos        - Temporary storage for pole positions
*   pos_axo         - Final estimate of the Axis offset or rate
*   pos_atm         - Final estimate of atmospheric delay
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
*   sol_parm(ndimc) - Solution vector from
*                   - the solution
*   temp_covar(36)  - Temporay storage for COMP_VAR
 
 
      real*8 covar(36), cov_parm(ndimc,ndimc), dt, loc_coord(3),
     .    obs_corr(1), pole_pos, pos_axo, pos_xyz_adj(3),
     .    pos_NEU_adj(3), pos_xyz_fin(3), pos_NEU_fin(3), pos_radec(2),
     .    rat_radec(2), rat_xyz_adj(3), rat_NEU_adj(3), rat_xyz_fin(3),
     .    rat_NEU_fin(3), rot_matrix(3,3), scr_real(10),
     .    sol_parm(ndimc), temp_covar(36), conv, pos_atm
 
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
     
* svs_fin(max_svs_elem) - Satellite elements
* sectag - Seconds tag of day

       real*8 svs_fin(max_svs_elem), sectag

* VARIABLES needed for PBOP format
* xyz_std(6) -- XYZ sigma (m) and corelations XY, XZ, YZ
* neu_std(6) -- NEU sigma (m) and corelations NE, NU, EU
* llu_std(6) -- Sigma Latitude, longitude (10^9 deg) and Height (m)
* NEU_covar(3,3)  - Covariance matrix for the NEU coordinates

      real*8 xyz_std(6), neu_std(6), llu_std(6), NEU_covar(3,3)
      real*8 unc_geod(3)  ! Uncorrelated epoch Geod co-lat, lng and
                          ! and height (rads, rads, m)

      integer*4 yr, doy, secod, sec, hrs, mns  ! Time tag for orbit IC

      logical kbit   ! tes Bit
* aei_lab(6)  - Names of the orbital elements

       character*23 aei_lab(6)

*   st_parm_label(21)   - Labels for different types of paramters
*                       - which have a site or source label
*   orb_lab(max_svs_elem)          - Labels for orbital parameters
 
      character*14 st_parm_label(22), orb_lab(max_svs_elem)
 
*   nl_parm_label(27)   - Labels for parameters which do not have
*                       - a site or source name 
*   lab                 - Label for multipmu epochs
 
      character*22 nl_parm_label(27)
      character*24 lab
 
*   unit_label( 6)      - Labels for the units
*   orb_unit(max_svs_elem)         - Units for orbit parameters
*   mp_lab(2,3)         - Labels for multi polar motion
 
      character*6 unit_label( 7), orb_unit(max_svs_elem),
     .            mp_lab(2,3)
 
*   nl                  - Format for skip a line
 
      character*4 nl

      character*1 hemi    ! Hemisphere for UTM coordinates 
      integer*4 zone      ! Zone for UTM coordinates

 
      common ident, dumglbak, covar, temp_covar
 
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
     .,                    'Zenith Delay  ' /

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

* MOD TAH 200806: Change random walk to rate (same as write_glb_pa.f) 
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
 
      data unit_label    / '(m)   '
     .,                    '(m/yr)'
     .,                    '(mas) '
     .,                    '(mts) '
     .,                    '(ms/y)'
     .,                    '(deg) '
     .,                    '(ms/d)' /
 
      data nl            / '(1x)'   /
 
      data aei_lab  /  'Semimajor axis      (m)',
     .                 'Eccentricity (none & m)',
     .                 'Inclination  (degs & m)',
     .                 'RA Node      (degs & m)',
     .                 'Arg. Perigee (degs & m)',
     .                 'M + w        (degs & m)'  / 

      data mp_lab    /  'X-pole','X-rate','Y-pole','Y-rate',
     .                  'UT1-AT','LOD   '  /
 
 
***** START checking the type
 
*     See if site positions
      IF( TYPE.EQ.7 ) THEN
 
*                             ! Loop over sites
      do i = 1, gnum
*                             ! Indicates if any componets estimated
          icnt = 0
          is   = ltog_map( ltog,i)
 
          dt   = (gepoch_expt-site_epoch(is))/365.25d0
* MOD TAH 190528: Output apriori coordinate
          if ( parn(1,i).ne.0 )
     .    write(iout,120) gsite_names(is),(apr_val_site(j,1,is),j=1,3),
     .                    (apr_val_site(j,2,is),j=1,3), 
     .                    (site_epoch(is)- 2451544.5d0)/365.25d0+2000
 120      format('Int. ',a8,1x,3(f14.5,1x),3(f10.5,1x),f8.3,
     .                1x,3(f7.4,1x) )
*         Do values first
*                             ! Loop over XYZ
          do j = 1,3
*                                       ! estimated
              if( parn(j,i).ne.0 ) then
                  iel  = parn(j,i)
                  icnt = icnt + 1
 
*                 Get final position of site
                  pos_xyz_adj(j) = sol_parm(iel)
* MOD TAH 190528: Addin non-secular terms.
                  pos_xyz_fin(j) = apr_val_site(j,1,is)    +
     .                             apr_val_site(j,2,is)*dt +
     .                             cont_nonsec(j,is) +
     .                             pos_xyz_adj(j)
 
                  call correct_obs( pos_xyz_fin(j), iel, obs_corr,ndimo)
 
                  call out_glbl(iout,parn(j,i), pos_xyz_fin(j),
     .                 sol_parm, cov_parm, ndimc, 1.d0, gsite_names(is),
     .                 st_parm_label(j), unit_label(1), 0)
*                             ! estimated
* MOD TAH 190526 Save sigmas for PBOP output
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
              end if
*                             ! Looping over XYZ
          end do
 
****      Now see if we should do NEU coordinates
 
*                                 ! Convert to NEU coordinates
          if( icnt.gt.0 ) then
              call rotate_geod(pos_xyz_adj, pos_neu_adj, 'XYZ','NEU',
     .                         pos_xyz_fin, loc_coord, rot_matrix)
 
*             Convert the local coordinates
              call loc_to_geod( loc_coord, pos_neu_fin )
 
*             Now compute the sigmas. Firstly get the elements of
*             covariance matrix we need
              call mov_cov(parn(1,i), parn(1,i), 3,
     .                     covar,3, cov_parm, ndimc)
              call var_comp(rot_matrix, covar, NEU_covar, temp_covar,
     .                      3,3,1)
 
*             Now output the results
              neu_std = 0 
              do j = 1,3
C                 scr_real(j) = sqrt(scr_real(j))
                  scr_real(j) = sqrt(abs(NEU_covar(j,j)))
                  call write_line( iout,-1, gsite_names(is),
     .                st_parm_label(6+j), unit_label(1),
     .                pos_neu_fin(j), pos_neu_adj(j), scr_real(j),0)


* MOD TAH 090926: Save sigmas for PBOP output
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

*                 Accumualte the statitics here as well, if these are
*                 are residuals we are doing.
* MOD TAH 190526: Residuals don't really work on loose solutions.
*                                             ! These are residuals
C                 if( ltog(1).ne.-1) then
C                     wgh = 1.d0/scr_real(j)**2
C                     iel = max_chi_types - 3 + j
C                     sum_type_wgh(iel) = sum_type_wgh(iel) + wgh
C                     sum_type_res(iel) = sum_type_res(iel) +
C    .                                    pos_neu_adj(j)**2*wgh
C                     sum_type_num(iel) = sum_type_num(iel) + 1
C                 end if
              end do

****          See if PBOP (position timeseries output)
              if( kbit(bak_opts,28) ) then
                  call jd_to_ymdhms(gepoch_expt,date,sectag)
                  call geod_to_geod(pos_xyz_fin, unc_geod, 
     .                  'XYZ', 'GEOD','WGS84','WGS84',zone,hemi)
                  llu_std(1) = neu_std(1)/Earth_rad*180/pi*1.d9
                  llu_std(2) = neu_std(2)/Earth_rad*180/pi/
     .                      sin(unc_geod(1))*1.d9
                  llu_std(3) = neu_std(3)

*                 Write out the very long line
                  write(iout,480) gsite_names(is),gsite_full(is)(1:16), 
     .                 date, gepoch_expt-2400000.5d0, 
     .                (pos_xyz_fin(j), j=1,3), (xyz_std(j),j=1,6),
     .                (pi/2-unc_geod(1))*180/pi,unc_geod(2)*180/pi,
     .                 unc_geod(3), (llu_std(j),j=1,3),
     .                (pos_neu_fin(j),j=1,3), (neu_std(j),j=1,6)
 480              format('pbo. ',a8,1x,a16,1x,i5,4(1x,i2.2),1x,F10.4,
     .                 1x,3F15.5,3F8.5,3F7.3,' | ',
     .                 2F16.10,1x,F10.5,1x,2F8.1,1x,F10.5,' | ',
     .                 2F15.5,1x,F10.5,3F8.5,3F7.3)
              endif
 
              write(iout,nl)
 
*                         ! Any local coordinates to be output
          end if
*                         ! Looping over the number of entries in parn
      end do
 
*                         ! TYPE matchs site positions
      ENDIF
 
***** See if we should do rates
 
*                             ! Do the rates
      if( type.eq.28 ) then
 
      do i = 1, gnum
 
*         Now do rates
*                             ! Rate counter
          icnr = 0
          is   = ltog_map( ltog, i )
 
*                             ! Loop over XYZ dot
          do j = 1,3
*                                       ! estimated
              if( parn(j,i).ne.0 ) then
                  iel  = parn(j,i)
                  icnr = icnr + 1
 
*                 Get final position of site
                  rat_xyz_adj(j) = sol_parm(iel)
                  rat_xyz_fin(j) = apr_val_site(j,2,is) + rat_xyz_adj(j)
 
                  call correct_obs( rat_xyz_fin(j), iel, obs_corr,ndimo)
 
                  call out_glbl(iout,parn(j,i), rat_xyz_fin(j),
     .                 sol_parm, cov_parm, ndimc, 1.d0, gsite_names(is),
     .                 st_parm_label(3+j), unit_label(2), 0)
*                             ! estimated
              end if
*                             ! Looping over XYZ dot
          end do
 
*****     Now do the rates in local coordinates
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
              call mov_cov(parn(1,i), parn(1,i), 3,
     .                     covar,3, cov_parm, ndimc)
              call var_comp(rot_matrix, covar, scr_real, temp_covar,
     .                      3,3,0)
 
*             Now output the results
              do j = 1,3
                  scr_real(j) = sqrt(scr_real(j))
                  call write_line( iout,-1, gsite_names(is),
     .                st_parm_label(9+j), unit_label(2),
     .                rat_neu_fin(j), rat_neu_adj(j), scr_real(j),0)
              end do
 
              write(iout,nl)
 
*                         ! Any local rates to be output
          end if
*                         ! Looping over the sites
      end do
 
*                         ! TYPE was for rates
      END IF

* MOD TAH 040703: Check for atmospheric delay estimates
      IF( TYPE.EQ.61 ) then
         icnt = 0
         do i = 1,gnum
            is = ltog_map(ltog,i)
            iel = parn(1,i)
            if( iel.gt.0 ) then
               pos_atm = apr_val_atm(is) + sol_parm(iel)
               call correct_obs( pos_atm, iel, obs_corr, ndimo)
               call out_glbl( iout, parn(1,i), pos_atm,
     .                 sol_parm, cov_parm, ndimc, 1.d0, gsite_names(is),
     .                 st_parm_label(22), unit_label(1), 0)
            end if
         end do
      end if

 
***** Axis offset ( values only )
 
*                             ! Axis offsets
      IF ( TYPE.EQ.10 ) then
 
          icnt = 0
          do i = 1,gnum
              is  = ltog_map( ltog,i )
 
              dt  = (gepoch_expt-axo_epoch(is))/365.25d0
 
              iel = parn(1,i)
*                                     ! Value estimated
              if( iel.gt.0 ) then
                  icnt = icnt + 1
                  pos_axo = apr_val_axo(1,is)    +
     .                      apr_val_axo(2,is)*dt +
     .                      sol_parm(iel)
 
                  call correct_obs( pos_axo, iel, obs_corr, ndimo)
 
                  call out_glbl( iout, parn(1,i), pos_axo,
     .                 sol_parm, cov_parm, ndimc, 1.d0, gsite_names(is),
     .                 st_parm_label(13), unit_label(1), 0)
              end if
          end do
 
          if ( icnt.gt.0 ) write(iout,nl)
 
      END IF
 
***** Now do source positions
 
      IF ( TYPE.EQ.11 ) then
 
      do i = 1, gnum
 
          is = ltog_map( ltog, i)
          dt = (gepoch_expt - source_epoch(is))/365.25d0
 
*                             ! loop over RA and Dec
          do j = 1,2
*                                       ! estimated
              if( parn(j,i).ne.0 ) then
                  iel = parn(j,i)
                  pos_radec(j) = apr_val_source(j,1,is)    +
     .                           apr_val_source(j,2,is)*dt +
     .                           sol_parm(iel)
 
                  call correct_obs( pos_radec(j), iel, obs_corr,ndimo)
 
*                 If RA convert to mts
                  if( j.eq.1 ) then
                      pos_radec(j) = pos_radec(j)/15.d0
                  end if
 
                  call out_glbl( iout, parn(j,i), pos_radec(j),
     .                 sol_parm, cov_parm, ndimc, 1.d0,
     .                 gsource_names(is),
     .                 st_parm_label(14+j),unit_label(3),1)
 
                  if( j.eq.2 ) write(iout,nl)
              end if
          end do
      end do
 
      END IF
 
*     Now do rates
 
*                             ! do rates
      IF ( TYPE.EQ. 31 ) THEN
 
      do i = 1, gnum
          is = ltog_map( ltog,i )
*                             ! loop over RA and Dec rates
          do j = 1,2
*                                       ! estimated
              if( parn(j,i).ne.0 ) then
                  iel = parn(j,i)
                  rat_radec(j) = apr_val_source(j,2,is) + sol_parm(iel)
 
                  call correct_obs( rat_radec(j), iel, obs_corr,ndimo)
 
*                 If RA convert to mts
                  if( j.eq.1 ) then
                      rat_radec(j) = rat_radec(j)/15.d0
                  end if
 
                  call out_glbl( iout, parn(j,i), rat_radec(j),
     .                 sol_parm, cov_parm, ndimc, 1.d0,
     .                 gsource_names(is),
     .                 st_parm_label(16+j),unit_label(5),1)
 
                  if( j.eq.2 ) write(iout,nl)
              end if
          end do
      end do
 
*                 ! Source rates
      END IF

***** Now do the satellite orbit parameters
      if( type.eq.51 ) then
         do i = 1, gnum
            icnt = 0
            call jd_to_yds( svs_epoch(1)+1.d-9, yr, doy, secod)
            if( yr.gt.500) yr = yr-1900
            hrs = secod/3600
            mns = (secod-hrs*3600)/60
            sec = (secod-hrs*3600-mns*60)
            write(iout,500) yr, doy, hrs, mns, sec, ggtime, ggframe,
     .                      ggprec, ggsrpmod
 500        format('Eph. #IC ',i2,1x,i3,1x,i2,1x,i2,1x,i2,20x,a4,1x,
     .             a5,1x, a5,1x,a5)

            do j = 1,max_svs_elem
                iel = parn(j,i)
                if( iel.gt.0 ) then
                    scr_real(1) = apr_val_svs(j,i) + sol_parm(iel)
                    
* MOD TAH 980226: Save svs_fin for output on the Eph. line.                
                    svs_fin(j) = scr_real(1)
                    
                    call out_glbl(iout, iel, scr_real(1), sol_parm,
     .                  cov_parm, ndimc, 1.d0, gsvs_names(i), 
     .                  orb_lab(j), orb_unit(j), 0 )
                    icnt = icnt + 1

*                   Save main part of orbital elements (ie. XYZ and
*                   XYZdot.  For the velocity terms convert to m/s)
                    if( j.le.6 ) then 
                        dXYZ_orb(j) = sol_parm(iel)
                        if( j.le.3 ) then
                            FXYZ_orb(j) = scr_real(1)
                        else
                            FXYZ_orb(j) = scr_real(1)/1000.d0
                        end if
                    end if
*                
                else
* MOD TAH 980226: Save svs_fin for output on the Eph. line.                
                    svs_fin(j) = apr_val_svs(j,i)
                    if( j.le.6 ) then
                        dXYZ_orb(j) = 0.d0
                        FXYZ_orb(j) = apr_val_svs(j,i)
                        if( j.ge.3 ) FXYZ_orb(j) = FXYZ_orb(j)/1000.d0
                    end if                 
                end if
            end do
            
* MOD TAH 980226: Now output one-line svs_file format line
            if( icnt.gt.0 ) then
                call jd_to_ymdhms(gepoch_expt, date, sectag)
                
                write(iout,505) (date(j),j=1,4), gsvs_names(i), 
     .                          (svs_fin(j),j=1,max_svs_elem)
 505            format('Eph. ',i5,3i3,1x,a8,3(1x,F14.4),3(1x,F14.5),
     .                 17(1x,F9.5))
            end if
            
****        Now convert the adjustments and cartesian elements to 
*           Keplerian elememts.
            if( icnt.gt.0 ) then
                call elem(FXYZ_orb(1), FXYZ_orb(2), FXYZ_orb(3),
     .                    FXYZ_orb(4), FXYZ_orb(5), FXYZ_orb(6),
     .                    aeinpa_rad, 3, dKdXYZ, GM_earth, 1) 

*               Add the argument of perigee to the Mean Anomaly
                aeinpa_rad(6) = aeinpa_rad(5) + aeinpa_rad(6)
                aeinpa_rad(6) = mod(aeinpa_rad(6), 2*pi)

*               Convert the angular arguments to meters (multiplu by
*               a and convert the velocities to be mm/s)
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

*               Finally add the Mean Anomaly to the argument of 
*               perigee
                do j = 1,6
                   dKdXYZ(6,j) = dKdXYZ(5,j) + dKdXYZ(6,j)
                end do

*               Now get the changes in the elements and compute the 
*               covariance matrix
                call mov_cov(parn(1,i), parn(1,i), 6, 
     .                       covar,  6, cov_parm, ndimc)
                call var_comp(dKdXYZ, covar, aei_covar, temp_covar,
     .                        6,6,1)

*               Now compute the adjustments to the values
                do j = 1,6
                   call dvdot(daei_orb(j), dKdXYZ(j,1),6, 
     .                        dXYZ_orb,1,6) 
                end do

*               Now write results
                do j = 1,1
                   write(iout,501) gsvs_names(i), aei_lab(j), 
     .                   aeinpa_rad(j), daei_orb(j), 
     .                   sqrt(aei_covar(j,j))
 501               format(' Loc. ',a,1x,a,t39,F18.4,2(1x,f12.4))
                end do
            
                do j = 2,2
                   write(iout,502) gsvs_names(i), aei_lab(j), 
     .                   aeinpa_rad(j), daei_orb(j), 
     .                   sqrt(aei_covar(j,j))
 502               format(' Loc. ',a,1x,a,t39,F18.10,2(1x,f12.4))
                end do
                do j = 3,6 
                   write(iout,503) gsvs_names(i), aei_lab(j), 
     .                   aeinpa_rad(j)*180.d0/pi, daei_orb(j), 
     .                   sqrt(aei_covar(j,j))
 503               format(' Loc. ',a,1x,a,t39,F18.10,2(1x,f12.4))
                end do

            end if

            if( icnt.ne.0 ) write(iout, nl)
         end do
      end if

 
***** Now do polar motion and UT1
 
*                             !   do wobble
      IF ( TYPE.EQ.13 ) THEN
 
      icnt = 0
*                     ! Loop over the 8 compoent model
* MOD TAH 200806: Change loop boundaries to match just offet+rate
      do j = 1,4
          if( parn(j,1).gt.0 ) then
              icnt = icnt + 1
              iel = parn(j,1)
              pole_pos = apr_val_wob(j) + gwob_apr(j) +
     .                                    sol_parm(iel)
 
              call correct_obs( pole_pos, iel, obs_corr,ndimo)
              n = 3
*             Change label for rate terms
* MOD TAH 210119: Fixed indexing to unit label (should be 7 not 10)
              if( j.eq.3 .or. j.eq.4 ) n = 7 

              call out_glbl(iout, iel , pole_pos,
     .             sol_parm, cov_parm, ndimc, 1.d0, nl_parm_label(j),
     .             unit_label(n), ' ',0)
          end if
      end do
      if( icnt.gt.0 ) write(iout,nl)
 
*                 ! wobble output
      END IF
 
****  Check UT1
*                             ! UT1 to be output
      IF ( TYPE.EQ.14 ) THEN
 
      icnt = 0
*                     ! Loop over the 2 compoent model
      do j = 1,2
          if( parn(j,1).gt.0 ) then
              icnt = icnt + 1
              iel = parn(j,1)
 
              pole_pos = apr_val_ut1(j) + gut1_apr(j) +
     .                                     sol_parm(iel)
              call correct_obs( pole_pos, iel, obs_corr,ndimo)
*                                          ! Convert to mts
              pole_pos = pole_pos / 15.d0
              n = 4
* MOD TAH 210119: Fixed indexing to unit label (should be 7 not 10)
              if( j.eq.2 ) n = 7 
 
              call out_glbl(iout, iel, pole_pos,
     .             sol_parm, cov_parm, ndimc, 1/15.d0,
     .             nl_parm_label(8+j), unit_label(n), ' ',0)
          end if
      end do
      if( icnt.gt.0 ) write(iout,nl)
 
*                     ! UT1 output
      END IF


****  Now write out the multi-epoch polar motion/UT1 values
*     Do X/Y Pole first 
      IF ( TYPE.EQ.56 .OR. TYPE.EQ.57 ) THEN
         j = type - 55
         do i = 1, 2
            icnt = 0
            do k = 1, num_mul_pmu
               if(  parn_mul_pmu(i,j,k).ne.0 ) then
                    icnt = icnt + 1
                    pole_pos = apr_val_mul_pmu(i,j,k)+
     .                         sol_parm(parn_mul_pmu(i,j,k))
                    if( j.eq.3 ) pole_pos = pole_pos/15.d0
                    n = 3
                    if( i.eq.2 ) n = 7

                    call jd_to_ymdhms(gmul_pmu_ep(k)+1.d-4, 
     .                                date, sectag)                 
                    write(lab,600) mp_lab(i,j), date
 600                format(a6,1x,I4,'/',I2.2,'/',I2.2,1x,I2.2,':',I2.2)          
                    call out_glbp(iout, parn_mul_pmu(i,j,k), pole_pos,
     .                  sol_parm, cov_parm, 1.d0, lab,
     .                  unit_label(n), ' ',0)
               endif
            end do
            if( icnt.gt.0 ) write(iout,nl)
         end do
      END IF 

*     Now do the UT1 and LOD values 
      IF ( TYPE.EQ.58 ) then
         j = 3
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
            if( icnt.gt.0 ) write(iout,nl)
         end do
      END IF
 
***** Nutation angles
 
*                             ! Output nutation angle
      IF ( TYPE.EQ. 15 ) THEN
 
      icnt = 0
*                     ! Loop over the 8 compoent model
      do j = 1,8
          if( parn(j,1).gt.0 ) then
              icnt = icnt + 1
              iel  = parn(j,1)
              pole_pos = apr_val_nut_ang(j) + gnut_ang_apr(j) +
     .                                        sol_parm(iel)
 
              call correct_obs( pole_pos, iel, obs_corr,ndimo)
 
              call out_glbl(iout, iel, pole_pos,
     .             sol_parm, cov_parm, ndimc, 1.d0,
     .             nl_parm_label(14+j), unit_label(3), ' ',0)
 
          end if
      end do
      if( icnt.gt.0 ) write(iout,nl)
 
*                 ! Nutations angles to be output
      END IF
 
*     Next add tides
 
***** Gamma
 
*                              ! Output gamma
      IF ( TYPE.EQ.27 ) THEN
 
      if( parn(1,1).ne.0 ) then
          iel = parn(1,1)
 
          scr_real(1) = apr_val_gamma + sol_parm(parn(1,1))
 
          call correct_obs( scr_real, iel, obs_corr,ndimo)
 
          call out_glbl(iout, parn, scr_real, sol_parm, cov_parm, ndimc,
     .         1.d0, nl_parm_label(23),' ',' ',0)
      end if
 
*                              ! Output gamma
      END IF
 
***** Next add ETD coefficients and NUT coefficients
 
 
***** Thats all
      return
      end
 
