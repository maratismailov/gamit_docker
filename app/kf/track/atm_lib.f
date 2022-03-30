
CTITLE DRY_MTT_MAP
 
      subroutine dry_mtt_map (temp_C, lat, ht, elev, dry_map)

      implicit none
 
*     Routine to compute the new dry_mit mapping function which
*     depends only on temperature (and position)

      include '../includes/const_param.h'

* PASSED VARIABLES
* dry_map  -- Dry mapping function
* elev     -- Elevation angle (deg)

      real*8 dry_map, elev

* temp_C   -- Temperature (C)
* lat, ht  -- Latitude (rad) and height of site (m)

      real*8 temp_C, lat, ht

* LOCAL VARIABLES 
*   A,B,C,D     - the A,B,C, and D coeffiecents in mit-2.2
*   beta        - Term in mit-2.2
*   cose        - Cos of elevation angle
*   gamma       - Term in mit-2.2
*   sine        - Sine of elevation angle
*   topcon      - Constant of top of mapping fuinction to ensure
*                 that value is 1.0000 at zenith 

      real*8 A,B,C, beta, cose,  gamma, sine, topcon

*   cosphi     - Cosine of latitude
*   hs_km      - Height of site in km

      real*8 cosphi, hs_km


***** Start compute the coefficients in the mapping function.
      cosphi = cos(lat)
      hs_km  = ht/1000.d0
 
*     Compute the coefficienst of the mapping function
      A  = 1.23200d-3 + 0.01391d-3*cosphi - 0.02089d-3*hs_km
     .                + 0.002154d-3*(temp_C - 10.d0)
      B  = 3.16116d-3 - 0.16004d-3*cosphi - 0.03306d-3*hs_km
     .                + 0.002064d-3*(temp_C - 10.d0)
      C  =71.24372d-3 - 4.29342d-3*cosphi - 0.14908d-3*hs_km
     .                - 0.002098d-3*(temp_C - 10.d0)
 
      sine  = sin( elev*pi/180 )
      cose  = cos( elev*pi/180 )
      beta  = B/( sine + C )
      gamma = A/( sine + beta)

      topcon = (1.d0 + A/(1.d0 + B/(1.d0 + C)))

      dry_map = topcon / ( sine + gamma )

 
***** Thats all
      return
      end
 
CTITLE WET_MTT_MAP
 
      subroutine wet_mtt_map ( temp_C, lat, ht, elev, wet_map ) 

      implicit none
 
*     Routine to compute the new wet_mit mapping function which
*     depends only on temperature (and position) 

      include '../includes/const_param.h'
 
* PASSED VARIABLES
* wet_map  -- Dry mapping function
* elev     -- Elevation angle (deg)

      real*8 wet_map,  elev 

* temp_C   -- Temperature (C)
* lat, ht  -- Latitude (rad) and height of site (m)

      real*8 temp_C, lat, ht

* LOCAL VARIABLES 
*   A,B,C,D     - the A,B,C, and D coeffiecents in mit-2.2
*   beta        - Term in mit-2.2
*   cose        - Cos of elevation angle
*   gamma       - Term in mit-2.2
*   sine        - Sine of elevation angle
*   topcon      - Constant of top of mapping fuinction to ensure
*                 that value is 1.0000 at zenith 

      real*8 A,B,C, beta, cose,  gamma, sine, topcon

*   cosphi     - Cosine of latitude
*   hs_km      - Height of site in km

      real*8 cosphi, hs_km

***** Start compute the coefficients in the mapping function.
      cosphi = cos(lat)
      hs_km  = ht/1000.d0
 
*     Compute the coefficienst of the mapping function
      A  = 0.58266d-3 - 0.01105d-3*cosphi - 0.05181d-3*hs_km
     .                + 0.001442d-3*(temp_C - 10.d0)
      B  = 1.40218d-3 + 0.10249d-3*cosphi - 0.10128d-3*hs_km
     .                + 0.002046d-3*(temp_C - 10.d0)
      C  =45.85450d-3 - 1.91277d-3*cosphi - 1.28787d-3*hs_km
     .                + 0.015136d-3*(temp_C - 10.d0)
 
      sine  = sin( elev*pi/180.d0 )
      cose  = cos( elev*pi/180.d0 )
      beta  = B/( sine + C )
      gamma = A/( sine + beta)

      topcon = (1.d0 + A/(1.d0 + B/(1.d0 + C)))

      wet_map = topcon / ( sine + gamma )

***** Thats all
      return
      end

CTITLE WET_SAAS_ZEN
 
      subroutine wet_saas_zen(temp_C, lat, ht, rel_hum, wet_zen )
 

      implicit none
 
*     Routine to compute wet zenith delay using Saastamoinen formula
*

* PASSED VARIABLES
* wet_zen  -- Wet Zenith delay (m)

      real*8 wet_zen
 
* temp_C   -- Temperature (C)
* lat, ht  -- Latitude (rad) and height of site (m)
* rel_hum  -- Relative humidity

      real*8 temp_C, lat, ht, rel_hum

* LOCAL VARIABLES
* fphih  -- Function which accounts for mean changes in gravity
* wet_pp -- Partial pressure of water vapor

      real*8 fphih, wet_pp
 
***** Compute the gravity function
      fphih =  1 - 0.00266d0 *cos(2*lat) - 0.00028D-03*ht 

***** Calculate the saturation pressure of water vapor (mbars)
      call wet_press(temp_C, rel_hum, wet_pp)
 
      wet_zen = 0.0022768d0* (1255.d0/(temp_C+273.15d0)+0.05d0)*
     .          wet_pp/fphih
 
****  Thats all
      return
      end

CTITLE DRY_SAAS_ZEN
 
      subroutine dry_saas_zen( temp_C, lat, ht, press_mb, dry_zen )

      implicit none
 
*     Routine to compute dry zenith delay using Saastamoinen formula
*

* PASSED VARIABLES
* dry_zen  -- Dry Zenith delay (m) 

      real*8 dry_zen

* temp_C   -- Temperature (C)
* lat, ht  -- Latitude (rad) and height of site (m)
* press_mb -- Pressure in mbars

      real*8 temp_C, lat, ht, press_mb

* LOCAL VARIABLES
* fphih  -- Function which accounts for mean changes in gravity

      real*8 fphih
 
***** Compute the gravity function
      fphih =  1 - 0.00266d0 *cos(2*lat) - 0.00028D-03*ht 

****  Compute the dry zenith delay 
      dry_zen = 0.0022768d0*press_mb/fphih
 
****  Thats all
      return
      end

CTITLE MET_SEASONAL
 
      subroutine met_seasonal( T_C, P, rh, tbias, epoch, lat, hgt)

      implicit none
 
*     Routine to return an estimate of temperature, pressure, and
*     relative humidity based on time, latitude, and height
 
      include '../includes/const_param.h'
 
* PASSED VARIABLES
 
*   t_c     - Temperature in C
*   P       - Pressure in mbar
*   rh      - Relative humidity in 0-1
*   lat     - Latitude of site in radians
*   hgt     - Height of site in m.
*   tbias   - Bias in surface temperature (Should not be added to
*             t_c returned from this routine).
*   t_k     - Temperature Kelvin
 
      real*8 t_c, P, rh, lat, hgt, tbias, t_k
 
*   epoch   - JD of current measurement.
*   dt      - Argument for seasonal term
 
      real*8 epoch, dt
 
*****
*     Get and estimate of temperature first.  Get seasonal argument
      dt = mod ((epoch-dj2000)/365.25,1.d0)*2*pi

      t_c   = ( -20.5 + 48.4*cos(lat)  - 3.1*hgt/1000.d0 ) +
     .         (-14.3 + 3.3*hgt/1000.d0)*sin(lat)*cos(dt)  +
     .         ( -4.7 + 1.1*hgt/1000.d0)*sin(lat)*sin(dt)  
 
      tbias = (  -3.6 - 0.3*cos(lat)  + 2.3*hgt/1000.d0 ) +
     .         ( -1.7 - 1.0*hgt/1000.d0)*sin(lat)*cos(dt)  +
     .         ( -0.2 + 0.8*hgt/1000.d0)*sin(lat)*sin(dt)  
 
*     Now based on standard lapse rate compute the pressure
      t_k = t_c + 273.15
      P = 1013.25d0 * (t_k/(t_k + 6.5d-3*hgt))**5.26d0
 
*     Set the relative humidity
      rh = 0.5
 
****  Thata all
      return
      end
 


CTITLE WET_PRESS
 
      subroutine wet_press(temp_C, rel_hum, wet_pp)

      implicit none
 
*     J.L. Davis                   3:02 PM  MON., 13  JULY, 1987
*
*     Calculates the wet partial pressure in mbars iven the
*     temperature and relative humidity
 

* PASSED VARIABLES
* wet_pp  -- Computed paritial pressure of water vapor

      real*8 wet_pp

* temp_C  -- Temperature C
* rel_hum -- Relative humidity (0-1)

      real*8 temp_C, rel_hum
 
* LOCAL VARIABLES
* temp_0  -- 0C in K
* temp_K  -- Temperature in K

      real*8 temp_0, temp_K
 
      data temp_0 / 273.15d0 /
 
***** Calculate the absolute temperature
      temp_K = temp_C + temp_0
 
***** Calculate the saturation pressure in mbars
* MOD TAH 880120: multiplied the saturated wet pressure by the rel_hum to
*     to get actual wet pressure.
 
      wet_pp = 6.11d0 * (temp_K / temp_0) ** (-5.3d0)
     .                *  exp(25.2d0 * temp_C / temp_K) * rel_hum
      
      return
      END
 

CTITLE HYDRO_PRESS
 
      subroutine hyrdo_press(press_surf, temp_C_surf, ht, press_ht)

      implicit none
 
*     Routine to return pressure at altitude given a surface pressure
*     temperature and the height

* PASSED VARIABLES
* press_surf  -- Pressure at the surface (mbar)
* temp_c_surf -- Temperature at the surface (C)
* ht          -- Height above the surface (m)
* press_ht    -- Hydrostatic pressure at height ht

      real*8 press_surf, temp_C_surf, ht, press_ht

* LOCAL VARIABLES
* temp_0    -- 0C in degrees Kelvin
* temp_K_ht -- Temperature at height assumming a 6.5 K/km lapse rate.
* temp_K_surf -- Temperatue in K at the surface  

      real*8 temp_0, temp_K_ht, temp_K_surf
 
 
      data temp_0 / 273.15d0 /
 
***** Calculate the absolute temperature at surface and altitude
      temp_K_surf = temp_C_surf + temp_0
      temp_K_ht   = temp_K_surf - 6.5d-3*ht

*     Compute the pressure at altitude
      press_ht = press_surf*(temp_K_ht/temp_K_surf)**5.26d0  

      return
      end

CTITLE GET_ATM

      subroutine get_atm( ns, mjd, dry_zen, wet_zen, lat, long, ht )

      implicit none

*     Routine to get the zenith atmospheric delays from the files
*     output by GAMIT.  These are peice-wise linear functions.

      include 'track_com.h'

* PASSED VARIABLES
* ns  -- Station number

      integer*4 ns

* mjd -- Time at which values are needed (MJD+fractional day)
* dry_zen -- Dry zenith delay from tabular interpolation
* wet_zen -- Wet zenith delay. (Returned as zero here)

      real*8 mjd, dry_zen, wet_zen

* LOCAL VARIABLES
* ent -- Entry number in the table of zenith delays 

      integer*4 ent

* dtd -- Time difference in days between the current epoch and
*        the previous tabular point
* dzen -- Difference in zenith delays between the tabule points

      real*8 dtd, dzen

      real*8 lat, long, ht, temp_C, press, rel_hum, tbias, undu 

***** Ok, set the wet delay to 0.0
      wet_zen = 0.0

*     Now find where we are in the tabular points
      ent = (mjd-atm_ref_mjd)/atm_ref_sep + 1
      if( ent.lt.0 .or. ent.gt. atm_ref_num ) then
         write(*,120) mjd, atm_ref_mjd, atm_ref_sep, atm_ref_num, ent
 120     format('**ERROR** Need to extrapolate atmosphere too far',
     .          ' Current MJD ',F10.4,' Atm Ref MJD ',F10.4,
     .          ' Separation ',f6.4,' days, Num ',i4,' Ent ',i4)
 
*        Set the use back to false so that standard model will be used
         use_atm_ref = .false.
         dry_zen = 2.3d0
         RETURN
      end if

*     If the entry is zero, start at first
      if( ent.lt.1 ) ent = 1
      if( ent.ge.atm_ref_num ) ent = atm_ref_num-1

*     OK, get the time separation
      dtd =  mjd - (atm_ref_mjd+(ent-1)*atm_ref_sep)
      dzen = atm_ref_zen(ent+1,ns)-atm_ref_zen(ent,ns)

*     Compute the zenith delay
      dry_zen = atm_ref_zen(ent,ns) + (dtd/atm_ref_sep)*dzen

      if( dry_zen .le.1.0 ) then
*         Values not valid; use regular model.
          if( atm_mtt ) then 
             call met_seasonal( temp_c, press, rel_hum, tbias, 
     .             mjd,lat, ht ) 
                            
*            Now get the zenith delay terms and map to the elevation angle
             call dry_saas_zen( temp_C, lat, ht, press,  dry_zen)
             call wet_saas_zen( temp_C, lat, ht, rel_hum,wet_zen)
          else
             call gpt(mjd,lat,long,ht, press, temp_C, undu) 
             call dry_saas_zen( temp_C, lat, ht, press,  dry_zen)
             rel_hum = ref_rel_hum
             call wet_saas_zen( temp_C, lat, ht, rel_hum,wet_zen)
          endif

      endif

****  Thats all
      return 
      end

CTITLE READ_ATM_FILE

      subroutine read_atm_file

      implicit none

*     Routine to read the atmosphere zenith delay file.

      include 'track_com.h'

* LOCAL VARIABLES
* ierr, jerr -- IOSTAT error
* indx -- Pointer to position in string
* ns, prev_ns  -- Site number and previous site number
* na   -- Number of atm entries
* i    -- Loop counter
* date  -- Date read from file 

      integer*4 ierr, jerr, indx, ns, prev_ns, na, i, date(5)

* sectag  -- Seconds tag of date
* read_mjd -- MJD read from file
* read_zen -- Zenith delat read from file

      real*8 sectag, read_mjd, read_zen

* line   -- Line read from file
* word   -- Word extracted from line
* name   -- Site name read from file
* test   -- Casefolded version of our site names
      character*256 line
      character*16  word
      character*4   name, test

****  OK, open the file
      open(99, file = atm_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',atm_file,0,
     .                  'READ_ATM_FILE')
      if( ierr.ne.0 ) then
          write(*,'(a)') 'Not able to apply values, default model used'
          use_atm_ref = .false.
          RETURN
      end if

****  OK, loop over the file
      prev_ns = 0
      do while ( ierr.eq.0 )
         read(99,'(a)',iostat=ierr) line

         if( ierr.eq.0 .and. line(1:7).eq.'ATM_ZEN' ) then 
*            OK, see what we have.  Get the first two words and
*            then the site name
             indx = 0
             call getword(line, word, indx)
             call getword(line, word, indx)
             call getword(line, name, indx)

*            See if name matches one of our sites
             ns = 0
             do i = 1, num_site
                test = site_names(i)
                call casefold(test)
                call casefold(name)
                if( test.eq.name ) then
                    ns = i
                endif
             end do 

*            If the atm site number matches the process the rest
*            of the line.
             if( ns.gt.0 ) then
*                Decode the rest of the line.  Skip next entry
                 call getword( line, word, indx)
                 call multiread(line, indx, 'I4', jerr, date, word, 5)
                 sectag = 0
                 call ymdhms_to_mjd( date, sectag, read_mjd)

*                Skip next three entries
                 call getword( line, word, indx)
                 call getword( line, word, indx)
                 call getword( line, word, indx)

*                Now read zenith delay
                 call read_line(line,indx,'R8', jerr, read_zen,word)

*                Is this the first or second entry?
                 if( ns.ne.prev_ns ) then
*                    This is a new entry.  Set the intial values
                     prev_ns = ns
                     na = 1
                     atm_ref_zen(1,ns) = read_zen
                     atm_ref_mjd = read_mjd 
                 else
                     na = na + 1
                     atm_ref_zen(na,ns) = read_zen
                     if( na.eq.2 ) then
                         atm_ref_sep = read_mjd - atm_ref_mjd
                     end if
                 end if
             end if
         end if
      end do

*     Save the number of atm values
      atm_ref_num = na
      use_atm_ref = .true.

      close(99)

***** Thats all
      return
      end

CTITLE REPORT_ATM_REF

      subroutine report_atm_ref(unit)

      implicit none

*     Routine to report the values being used for the zenith
*     delays

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
* unit  -- Unit number for output

      integer*4 unit

* LOCAL VARIABLES
* date(5)  -- Date
* i, j     -- Loop counterss

      integer*4 date(5), i, j, ns, trimlen

* sectag   -- Seconds tag
* mjd      -- Compute time of zenith delay

      real*8 sectag, mjd

* int_atm(max_site)
* lat, ht  -- Latitude and height of site (used for atmospheric delay
* long -- Lonitude (rads) 

      real*8 lat(max_site), long(max_site), ht(max_site),
     .       int_atm(max_site) 

      real*8 loc_coord(3), rot_matrix(3,3),dry_zen,  wet_zen, dt

****  OK, Report the values
      if( .not. use_atm_ref  ) RETURN

      write(unit, 120) atm_file(1:trimlen(atm_file)),
     .                 (site_names(i),i=1,num_site)
 120  format('ZENITH DELAYS FROM ',a,/,
     .       'YYYY MM DD HR MIN ',20(3x,a4,1x))
      do i = 1, atm_ref_num
         mjd = atm_ref_mjd + (i-1)*atm_ref_sep + 1.d-6
         call mjd_to_ymdhms(mjd, date,sectag )
         write(unit,160) date, (atm_ref_zen(i,j),j=1,num_site)
 160     format(i4,4i3,2x,20(1x,F7.4))
      end do

***   Test out put test values
      do i = 1, num_site
         call XYZ_to_GEOD( rot_matrix, site_apr(1,i), loc_coord)

*        Save the latitude and height for computing the atmospheric
*        delay
         lat(i) = (pi/2.d0-loc_coord(1))
         long(i) = loc_coord(2)
         ht(i)  = loc_coord(3)
      end do

      write(unit, 220) atm_file(1:trimlen(atm_file)),
     .                 (site_names(i),i=1,num_site)
 220  format('INTERPOLATED ZENITH DELAYS FROM ',a,/,
     .       'YYYY MM DD HR MIN ',20(3x,a4,1x))
      dt = (data_end(1)-ref_start)/12
      do i = 1,12
         mjd = ref_start + (i-1)*dt
         call mjd_to_ymdhms(mjd, date,sectag )
         do ns = 1, num_site
             call get_atm( ns, mjd, dry_zen, wet_zen, lat(ns), 
     .                    long(ns), ht(ns) )
             int_atm(ns) = dry_zen + wet_zen
         enddo

         write(unit,160) date, (int_atm(j),j=1,num_site)
      enddo

****  Thats all
      return 
      end

      subroutine gpt (dmjd,dlat,dlon,dhgt,pres,temp,undu)

      implicit none

C     This subroutine determines Global Pressure and Temperature
C     based on Spherical Harmonics up to degree and order 9
C
C     input data
C     ----------
C     dmjd: modified julian date
C     dlat: latitude in radians
C     dlon: longitude in radians
C     dhgt: ellipsoidal height in m
C
C     output data
C     -----------
C     pres: pressure in hPa
C     temp: temperature in Celsius
C     undu: Geoid undulation in m (from a 9x9 EGM based model)
C
C     Johannes Boehm, 2006 June 12
C     rev 2006 June 16: geoid undulation is accounted for
C
C     Reference:
C     J. Boehm, and H. Schuh,
C     Global Pressure and Temperature (GPT): A spherical harmonic expansion 
C     of annual pressure and temperature variations for geodetic applications,
C     to be submitted to Journal of Geodesy, 2006
c
c     PT060615: variables declared explicitly 
c
c      implicit double precision (a-h,o-z)
c
c      dimension dfac(20),P(10,10),aP(55),bP(55),
c     .          ap_mean(55),bp_mean(55),ap_amp(55),bp_amp(55),
c     .          at_mean(55),bt_mean(55),at_amp(55),bt_amp(55)

       real*8  dfac(20),P(10,10),aP(55),bP(55),
     .          ap_mean(55),bp_mean(55),ap_amp(55),bp_amp(55)
     .          ,at_mean(55),bt_mean(55),at_amp(55),bt_amp(55)
     .          ,a_geoid(55),b_geoid(55)
     .          ,pi,doy,dmjd,pres,temp,dlat,dlon,dhgt,t,apm,apa
     .          ,pres0,atm,ata,temp0,sum,hort,undu

       integer i,j,k,m,n,ir



      pi = 3.14159265359d0

C     reference day is 28 January
C     this is taken from Niell (1996) to be consistent
      doy = dmjd  - 44239.d0 + 1 - 28

      data (a_geoid(i),i=1,55)/
     .-5.6195d-001,-6.0794d-002,-2.0125d-001,-6.4180d-002,-3.6997d-002,
     .+1.0098d+001,+1.6436d+001,+1.4065d+001,+1.9881d+000,+6.4414d-001,
     .-4.7482d+000,-3.2290d+000,+5.0652d-001,+3.8279d-001,-2.6646d-002,
     .+1.7224d+000,-2.7970d-001,+6.8177d-001,-9.6658d-002,-1.5113d-002,
     .+2.9206d-003,-3.4621d+000,-3.8198d-001,+3.2306d-002,+6.9915d-003,
     .-2.3068d-003,-1.3548d-003,+4.7324d-006,+2.3527d+000,+1.2985d+000,
     .+2.1232d-001,+2.2571d-002,-3.7855d-003,+2.9449d-005,-1.6265d-004,
     .+1.1711d-007,+1.6732d+000,+1.9858d-001,+2.3975d-002,-9.0013d-004,
     .-2.2475d-003,-3.3095d-005,-1.2040d-005,+2.2010d-006,-1.0083d-006,
     .+8.6297d-001,+5.8231d-001,+2.0545d-002,-7.8110d-003,-1.4085d-004,
     .-8.8459d-006,+5.7256d-006,-1.5068d-006,+4.0095d-007,-2.4185d-008/

      data (b_geoid(i),i=1,55)/
     .+0.0000d+000,+0.0000d+000,-6.5993d-002,+0.0000d+000,+6.5364d-002,
     .-5.8320d+000,+0.0000d+000,+1.6961d+000,-1.3557d+000,+1.2694d+000,
     .+0.0000d+000,-2.9310d+000,+9.4805d-001,-7.6243d-002,+4.1076d-002,
     .+0.0000d+000,-5.1808d-001,-3.4583d-001,-4.3632d-002,+2.2101d-003,
     .-1.0663d-002,+0.0000d+000,+1.0927d-001,-2.9463d-001,+1.4371d-003,
     .-1.1452d-002,-2.8156d-003,-3.5330d-004,+0.0000d+000,+4.4049d-001,
     .+5.5653d-002,-2.0396d-002,-1.7312d-003,+3.5805d-005,+7.2682d-005,
     .+2.2535d-006,+0.0000d+000,+1.9502d-002,+2.7919d-002,-8.1812d-003,
     .+4.4540d-004,+8.8663d-005,+5.5596d-005,+2.4826d-006,+1.0279d-006,
     .+0.0000d+000,+6.0529d-002,-3.5824d-002,-5.1367d-003,+3.0119d-005,
     .-2.9911d-005,+1.9844d-005,-1.2349d-006,-7.6756d-009,+5.0100d-008/

      data (ap_mean(i),i=1,55)/ 
     .+1.0108d+003,+8.4886d+000,+1.4799d+000,-1.3897d+001,+3.7516d-003,
     .-1.4936d-001,+1.2232d+001,-7.6615d-001,-6.7699d-002,+8.1002d-003,
     .-1.5874d+001,+3.6614d-001,-6.7807d-002,-3.6309d-003,+5.9966d-004,
     .+4.8163d+000,-3.7363d-001,-7.2071d-002,+1.9998d-003,-6.2385d-004,
     .-3.7916d-004,+4.7609d+000,-3.9534d-001,+8.6667d-003,+1.1569d-002,
     .+1.1441d-003,-1.4193d-004,-8.5723d-005,+6.5008d-001,-5.0889d-001,
     .-1.5754d-002,-2.8305d-003,+5.7458d-004,+3.2577d-005,-9.6052d-006,
     .-2.7974d-006,+1.3530d+000,-2.7271d-001,-3.0276d-004,+3.6286d-003,
     .-2.0398d-004,+1.5846d-005,-7.7787d-006,+1.1210d-006,+9.9020d-008,
     .+5.5046d-001,-2.7312d-001,+3.2532d-003,-2.4277d-003,+1.1596d-004,
     .+2.6421d-007,-1.3263d-006,+2.7322d-007,+1.4058d-007,+4.9414d-009/ 

      data (bp_mean(i),i=1,55)/ 
     .+0.0000d+000,+0.0000d+000,-1.2878d+000,+0.0000d+000,+7.0444d-001,
     .+3.3222d-001,+0.0000d+000,-2.9636d-001,+7.2248d-003,+7.9655d-003,
     .+0.0000d+000,+1.0854d+000,+1.1145d-002,-3.6513d-002,+3.1527d-003,
     .+0.0000d+000,-4.8434d-001,+5.2023d-002,-1.3091d-002,+1.8515d-003,
     .+1.5422d-004,+0.0000d+000,+6.8298d-001,+2.5261d-003,-9.9703d-004,
     .-1.0829d-003,+1.7688d-004,-3.1418d-005,+0.0000d+000,-3.7018d-001,
     .+4.3234d-002,+7.2559d-003,+3.1516d-004,+2.0024d-005,-8.0581d-006,
     .-2.3653d-006,+0.0000d+000,+1.0298d-001,-1.5086d-002,+5.6186d-003,
     .+3.2613d-005,+4.0567d-005,-1.3925d-006,-3.6219d-007,-2.0176d-008,
     .+0.0000d+000,-1.8364d-001,+1.8508d-002,+7.5016d-004,-9.6139d-005,
     .-3.1995d-006,+1.3868d-007,-1.9486d-007,+3.0165d-010,-6.4376d-010/ 

      data (ap_amp(i),i=1,55)/ 
     .-1.0444d-001,+1.6618d-001,-6.3974d-002,+1.0922d+000,+5.7472d-001,
     .-3.0277d-001,-3.5087d+000,+7.1264d-003,-1.4030d-001,+3.7050d-002,
     .+4.0208d-001,-3.0431d-001,-1.3292d-001,+4.6746d-003,-1.5902d-004,
     .+2.8624d+000,-3.9315d-001,-6.4371d-002,+1.6444d-002,-2.3403d-003,
     .+4.2127d-005,+1.9945d+000,-6.0907d-001,-3.5386d-002,-1.0910d-003,
     .-1.2799d-004,+4.0970d-005,+2.2131d-005,-5.3292d-001,-2.9765d-001,
     .-3.2877d-002,+1.7691d-003,+5.9692d-005,+3.1725d-005,+2.0741d-005,
     .-3.7622d-007,+2.6372d+000,-3.1165d-001,+1.6439d-002,+2.1633d-004,
     .+1.7485d-004,+2.1587d-005,+6.1064d-006,-1.3755d-008,-7.8748d-008,
     .-5.9152d-001,-1.7676d-001,+8.1807d-003,+1.0445d-003,+2.3432d-004,
     .+9.3421d-006,+2.8104d-006,-1.5788d-007,-3.0648d-008,+2.6421d-010/ 

      data (bp_amp(i),i=1,55)/ 
     .+0.0000d+000,+0.0000d+000,+9.3340d-001,+0.0000d+000,+8.2346d-001,
     .+2.2082d-001,+0.0000d+000,+9.6177d-001,-1.5650d-002,+1.2708d-003,
     .+0.0000d+000,-3.9913d-001,+2.8020d-002,+2.8334d-002,+8.5980d-004,
     .+0.0000d+000,+3.0545d-001,-2.1691d-002,+6.4067d-004,-3.6528d-005,
     .-1.1166d-004,+0.0000d+000,-7.6974d-002,-1.8986d-002,+5.6896d-003,
     .-2.4159d-004,-2.3033d-004,-9.6783d-006,+0.0000d+000,-1.0218d-001,
     .-1.3916d-002,-4.1025d-003,-5.1340d-005,-7.0114d-005,-3.3152d-007,
     .+1.6901d-006,+0.0000d+000,-1.2422d-002,+2.5072d-003,+1.1205d-003,
     .-1.3034d-004,-2.3971d-005,-2.6622d-006,+5.7852d-007,+4.5847d-008,
     .+0.0000d+000,+4.4777d-002,-3.0421d-003,+2.6062d-005,-7.2421d-005,
     .+1.9119d-006,+3.9236d-007,+2.2390d-007,+2.9765d-009,-4.6452d-009/ 

      data (at_mean(i),i=1,55)/ 
     .+1.6257e+001,+2.1224e+000,+9.2569e-001,-2.5974e+001,+1.4510e+000,
     .+9.2468e-002,-5.3192e-001,+2.1094e-001,-6.9210e-002,-3.4060e-002,
     .-4.6569e+000,+2.6385e-001,-3.6093e-002,+1.0198e-002,-1.8783e-003,
     .+7.4983e-001,+1.1741e-001,+3.9940e-002,+5.1348e-003,+5.9111e-003,
     .+8.6133e-006,+6.3057e-001,+1.5203e-001,+3.9702e-002,+4.6334e-003,
     .+2.4406e-004,+1.5189e-004,+1.9581e-007,+5.4414e-001,+3.5722e-001,
     .+5.2763e-002,+4.1147e-003,-2.7239e-004,-5.9957e-005,+1.6394e-006,
     .-7.3045e-007,-2.9394e+000,+5.5579e-002,+1.8852e-002,+3.4272e-003,
     .-2.3193e-005,-2.9349e-005,+3.6397e-007,+2.0490e-006,-6.4719e-008,
     .-5.2225e-001,+2.0799e-001,+1.3477e-003,+3.1613e-004,-2.2285e-004,
     .-1.8137e-005,-1.5177e-007,+6.1343e-007,+7.8566e-008,+1.0749e-009/ 

      data (bt_mean(i),i=1,55)/ 
     .+0.0000e+000,+0.0000e+000,+1.0210e+000,+0.0000e+000,+6.0194e-001,
     .+1.2292e-001,+0.0000e+000,-4.2184e-001,+1.8230e-001,+4.2329e-002,
     .+0.0000e+000,+9.3312e-002,+9.5346e-002,-1.9724e-003,+5.8776e-003,
     .+0.0000e+000,-2.0940e-001,+3.4199e-002,-5.7672e-003,-2.1590e-003,
     .+5.6815e-004,+0.0000e+000,+2.2858e-001,+1.2283e-002,-9.3679e-003,
     .-1.4233e-003,-1.5962e-004,+4.0160e-005,+0.0000e+000,+3.6353e-002,
     .-9.4263e-004,-3.6762e-003,+5.8608e-005,-2.6391e-005,+3.2095e-006,
     .-1.1605e-006,+0.0000e+000,+1.6306e-001,+1.3293e-002,-1.1395e-003,
     .+5.1097e-005,+3.3977e-005,+7.6449e-006,-1.7602e-007,-7.6558e-008,
     .+0.0000e+000,-4.5415e-002,-1.8027e-002,+3.6561e-004,-1.1274e-004,
     .+1.3047e-005,+2.0001e-006,-1.5152e-007,-2.7807e-008,+7.7491e-009/ 

      data (at_amp(i),i=1,55)/ 
     .-1.8654e+000,-9.0041e+000,-1.2974e-001,-3.6053e+000,+2.0284e-002,
     .+2.1872e-001,-1.3015e+000,+4.0355e-001,+2.2216e-001,-4.0605e-003,
     .+1.9623e+000,+4.2887e-001,+2.1437e-001,-1.0061e-002,-1.1368e-003,
     .-6.9235e-002,+5.6758e-001,+1.1917e-001,-7.0765e-003,+3.0017e-004,
     .+3.0601e-004,+1.6559e+000,+2.0722e-001,+6.0013e-002,+1.7023e-004,
     .-9.2424e-004,+1.1269e-005,-6.9911e-006,-2.0886e+000,-6.7879e-002,
     .-8.5922e-004,-1.6087e-003,-4.5549e-005,+3.3178e-005,-6.1715e-006,
     .-1.4446e-006,-3.7210e-001,+1.5775e-001,-1.7827e-003,-4.4396e-004,
     .+2.2844e-004,-1.1215e-005,-2.1120e-006,-9.6421e-007,-1.4170e-008,
     .+7.8720e-001,-4.4238e-002,-1.5120e-003,-9.4119e-004,+4.0645e-006,
     .-4.9253e-006,-1.8656e-006,-4.0736e-007,-4.9594e-008,+1.6134e-009/ 

      data (bt_amp(i),i=1,55)/ 
     .+0.0000e+000,+0.0000e+000,-8.9895e-001,+0.0000e+000,-1.0790e+000,
     .-1.2699e-001,+0.0000e+000,-5.9033e-001,+3.4865e-002,-3.2614e-002,
     .+0.0000e+000,-2.4310e-002,+1.5607e-002,-2.9833e-002,-5.9048e-003,
     .+0.0000e+000,+2.8383e-001,+4.0509e-002,-1.8834e-002,-1.2654e-003,
     .-1.3794e-004,+0.0000e+000,+1.3306e-001,+3.4960e-002,-3.6799e-003,
     .-3.5626e-004,+1.4814e-004,+3.7932e-006,+0.0000e+000,+2.0801e-001,
     .+6.5640e-003,-3.4893e-003,-2.7395e-004,+7.4296e-005,-7.9927e-006,
     .-1.0277e-006,+0.0000e+000,+3.6515e-002,-7.4319e-003,-6.2873e-004,
     .-8.2461e-005,+3.1095e-005,-5.3860e-007,-1.2055e-007,-1.1517e-007,
     .+0.0000e+000,+3.1404e-002,+1.5580e-002,-1.1428e-003,+3.3529e-005,
     .+1.0387e-005,-1.9378e-006,-2.7327e-007,+7.5833e-009,-9.2323e-009/ 

C     parameter t
      t = dsin(dlat)

C     degree n and order m
      n = 9
      m = 9

C determine n!  (faktorielle)  moved by 1
      dfac(1) = 1
      do i = 1,(2*n + 1)
        dfac(i+1) = dfac(i)*i
      end do

C     determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
      do i = 0,n
        do j = 0,min(i,m)
          ir = int((i - j)/2)
          sum = 0
          do k = 0,ir
            sum = sum + (-1)**k*dfac(2*i - 2*k + 1)/dfac(k + 1)/
     .       dfac(i - k + 1)/dfac(i - j - 2*k + 1)*t**(i - j - 2*k)
          end do
C         Legendre functions moved by 1
          P(i + 1,j + 1) = 1.d0/2**i*dsqrt((1 - t**2)**(j))*sum
        end do
      end do

C     spherical harmonics
      i = 0
      do n = 0,9
        do m = 0,n
          i = i + 1
          aP(i) = P(n+1,m+1)*dcos(m*dlon)
          bP(i) = P(n+1,m+1)*dsin(m*dlon)
        end do
      end do

C     Geoidal height
      undu = 0.d0
      do i = 1,55
        undu = undu + (a_geoid(i)*aP(i) + b_geoid(i)*bP(i))
      end do

C     orthometric height
      hort = dhgt - undu

C     Surface pressure
      apm = 0.d0
      apa = 0.d0
      do i = 1,55
        apm = apm + (ap_mean(i)*aP(i) + bp_mean(i)*bP(i))
        apa = apa + (ap_amp(i) *aP(i) + bp_amp(i) *bP(i))
      end do
      pres0  = apm + apa*dcos(doy/365.25d0*2.d0*pi)

C     height correction for pressure
      pres = pres0*(1.d0-0.0000226d0*hort)**5.225d0

C     Surface temperature
      atm = 0.d0
      ata = 0.d0
      do i = 1,55
        atm = atm + (at_mean(i)*aP(i) + bt_mean(i)*bP(i))
        ata = ata + (at_amp(i) *aP(i) + bt_amp(i) *bP(i))
      end do
      temp0 =  atm + ata*dcos(doy/365.25d0*2*pi)

C     height correction for temperature
      temp = temp0 - 0.0065d0*hort
          
      return
      end

      subroutine gmf (dmjd,dlat,dlon,dhgt,zd,gmfh,gmfw)

      implicit none

C     This subroutine determines the Global Mapping Functions GMF
C     Reference: Boehm, J., A.E. Niell, H. Schuh: "Global Mapping Functions (GMF):
c     A new empirical mapping function based on numerical weather model data"
C     Geophys. Res. Lett., 33, L07304, doi:10.1029/2005GL025546, 2006.
C
C     input data
C     ----------
C     dmjd: modified julian date
C     dlat: latitude in radians
C     dlon: longitude in radians
C     dhgt: height in m
C     zd:   zenith distance in radians
C
C     output data
C     -----------
C     gmfh: hydrostatic mapping function
C     gmfw: wet mapping function
C
C     Johannes Boehm, 2005 August 30
C
c PT050831: remove implicit declarations and declare all the variables

      integer i,j,k,m,n,ir
      real*8 dfac(20),P(10,10),aP(55),bP(55)
     .      ,ah_mean(55),bh_mean(55),ah_amp(55),bh_amp(55)
     .      ,aw_mean(55),bw_mean(55),aw_amp(55),bw_amp(55)
     .      ,doy,dmjd,dlat,dlon,dhgt,zd,gmfh,gmfw,pi,t,sum
     .      ,bh,c0h,phh,c11h,c10h,ch,ahm,aha,ah,sine,beta
     .      ,gamma,topcon,a_ht,b_ht,c_ht,hs_km,ht_corr_coef
     .      ,ht_corr,aw,bw,cw,awm,awa
  
      pi = 3.1415926359d0

C     reference day is 28 January
C     this is taken from Niell (1996) to be consistent
      doy = dmjd  - 44239.d0 + 1 - 28

      data (ah_mean(i),i=1,55)/
     .+1.2517d+02, +8.503d-01, +6.936d-02, -6.760d+00, +1.771d-01,
     . +1.130d-02, +5.963d-01, +1.808d-02, +2.801d-03, -1.414d-03,
     . -1.212d+00, +9.300d-02, +3.683d-03, +1.095d-03, +4.671d-05,
     . +3.959d-01, -3.867d-02, +5.413d-03, -5.289d-04, +3.229d-04,
     . +2.067d-05, +3.000d-01, +2.031d-02, +5.900d-03, +4.573d-04,
     . -7.619d-05, +2.327d-06, +3.845d-06, +1.182d-01, +1.158d-02,
     . +5.445d-03, +6.219d-05, +4.204d-06, -2.093d-06, +1.540d-07,
     . -4.280d-08, -4.751d-01, -3.490d-02, +1.758d-03, +4.019d-04,
     . -2.799d-06, -1.287d-06, +5.468d-07, +7.580d-08, -6.300d-09,
     . -1.160d-01, +8.301d-03, +8.771d-04, +9.955d-05, -1.718d-06,
     . -2.012d-06, +1.170d-08, +1.790d-08, -1.300d-09, +1.000d-10/

      data (bh_mean(i),i=1,55)/
     . +0.000d+00, +0.000d+00, +3.249d-02, +0.000d+00, +3.324d-02,
     . +1.850d-02, +0.000d+00, -1.115d-01, +2.519d-02, +4.923d-03,
     . +0.000d+00, +2.737d-02, +1.595d-02, -7.332d-04, +1.933d-04,
     . +0.000d+00, -4.796d-02, +6.381d-03, -1.599d-04, -3.685d-04,
     . +1.815d-05, +0.000d+00, +7.033d-02, +2.426d-03, -1.111d-03,
     . -1.357d-04, -7.828d-06, +2.547d-06, +0.000d+00, +5.779d-03,
     . +3.133d-03, -5.312d-04, -2.028d-05, +2.323d-07, -9.100d-08,
     . -1.650d-08, +0.000d+00, +3.688d-02, -8.638d-04, -8.514d-05,
     . -2.828d-05, +5.403d-07, +4.390d-07, +1.350d-08, +1.800d-09,
     . +0.000d+00, -2.736d-02, -2.977d-04, +8.113d-05, +2.329d-07,
     . +8.451d-07, +4.490d-08, -8.100d-09, -1.500d-09, +2.000d-10/
       
      data (ah_amp(i),i=1,55)/
     . -2.738d-01, -2.837d+00, +1.298d-02, -3.588d-01, +2.413d-02,
     . +3.427d-02, -7.624d-01, +7.272d-02, +2.160d-02, -3.385d-03,
     . +4.424d-01, +3.722d-02, +2.195d-02, -1.503d-03, +2.426d-04,
     . +3.013d-01, +5.762d-02, +1.019d-02, -4.476d-04, +6.790d-05,
     . +3.227d-05, +3.123d-01, -3.535d-02, +4.840d-03, +3.025d-06,
     . -4.363d-05, +2.854d-07, -1.286d-06, -6.725d-01, -3.730d-02,
     . +8.964d-04, +1.399d-04, -3.990d-06, +7.431d-06, -2.796d-07,
     . -1.601d-07, +4.068d-02, -1.352d-02, +7.282d-04, +9.594d-05,
     . +2.070d-06, -9.620d-08, -2.742d-07, -6.370d-08, -6.300d-09,
     . +8.625d-02, -5.971d-03, +4.705d-04, +2.335d-05, +4.226d-06,
     . +2.475d-07, -8.850d-08, -3.600d-08, -2.900d-09, +0.000d+00/
       
      data (bh_amp(i),i=1,55)/
     . +0.000d+00, +0.000d+00, -1.136d-01, +0.000d+00, -1.868d-01,
     . -1.399d-02, +0.000d+00, -1.043d-01, +1.175d-02, -2.240d-03,
     . +0.000d+00, -3.222d-02, +1.333d-02, -2.647d-03, -2.316d-05,
     . +0.000d+00, +5.339d-02, +1.107d-02, -3.116d-03, -1.079d-04,
     . -1.299d-05, +0.000d+00, +4.861d-03, +8.891d-03, -6.448d-04,
     . -1.279d-05, +6.358d-06, -1.417d-07, +0.000d+00, +3.041d-02,
     . +1.150d-03, -8.743d-04, -2.781d-05, +6.367d-07, -1.140d-08,
     . -4.200d-08, +0.000d+00, -2.982d-02, -3.000d-03, +1.394d-05,
     . -3.290d-05, -1.705d-07, +7.440d-08, +2.720d-08, -6.600d-09,
     . +0.000d+00, +1.236d-02, -9.981d-04, -3.792d-05, -1.355d-05,
     . +1.162d-06, -1.789d-07, +1.470d-08, -2.400d-09, -4.000d-10/
       
      data (aw_mean(i),i=1,55)/
     . +5.640d+01, +1.555d+00, -1.011d+00, -3.975d+00, +3.171d-02,
     . +1.065d-01, +6.175d-01, +1.376d-01, +4.229d-02, +3.028d-03,
     . +1.688d+00, -1.692d-01, +5.478d-02, +2.473d-02, +6.059d-04,
     . +2.278d+00, +6.614d-03, -3.505d-04, -6.697d-03, +8.402d-04,
     . +7.033d-04, -3.236d+00, +2.184d-01, -4.611d-02, -1.613d-02,
     . -1.604d-03, +5.420d-05, +7.922d-05, -2.711d-01, -4.406d-01,
     . -3.376d-02, -2.801d-03, -4.090d-04, -2.056d-05, +6.894d-06,
     . +2.317d-06, +1.941d+00, -2.562d-01, +1.598d-02, +5.449d-03,
     . +3.544d-04, +1.148d-05, +7.503d-06, -5.667d-07, -3.660d-08,
     . +8.683d-01, -5.931d-02, -1.864d-03, -1.277d-04, +2.029d-04,
     . +1.269d-05, +1.629d-06, +9.660d-08, -1.015d-07, -5.000d-10/
       
      data (bw_mean(i),i=1,55)/
     . +0.000d+00, +0.000d+00, +2.592d-01, +0.000d+00, +2.974d-02,
     . -5.471d-01, +0.000d+00, -5.926d-01, -1.030d-01, -1.567d-02,
     . +0.000d+00, +1.710d-01, +9.025d-02, +2.689d-02, +2.243d-03,
     . +0.000d+00, +3.439d-01, +2.402d-02, +5.410d-03, +1.601d-03,
     . +9.669d-05, +0.000d+00, +9.502d-02, -3.063d-02, -1.055d-03,
     . -1.067d-04, -1.130d-04, +2.124d-05, +0.000d+00, -3.129d-01,
     . +8.463d-03, +2.253d-04, +7.413d-05, -9.376d-05, -1.606d-06,
     . +2.060d-06, +0.000d+00, +2.739d-01, +1.167d-03, -2.246d-05,
     . -1.287d-04, -2.438d-05, -7.561d-07, +1.158d-06, +4.950d-08,
     . +0.000d+00, -1.344d-01, +5.342d-03, +3.775d-04, -6.756d-05,
     . -1.686d-06, -1.184d-06, +2.768d-07, +2.730d-08, +5.700d-09/
       
      data (aw_amp(i),i=1,55)/
     . +1.023d-01, -2.695d+00, +3.417d-01, -1.405d-01, +3.175d-01,
     . +2.116d-01, +3.536d+00, -1.505d-01, -1.660d-02, +2.967d-02,
     . +3.819d-01, -1.695d-01, -7.444d-02, +7.409d-03, -6.262d-03,
     . -1.836d+00, -1.759d-02, -6.256d-02, -2.371d-03, +7.947d-04,
     . +1.501d-04, -8.603d-01, -1.360d-01, -3.629d-02, -3.706d-03,
     . -2.976d-04, +1.857d-05, +3.021d-05, +2.248d+00, -1.178d-01,
     . +1.255d-02, +1.134d-03, -2.161d-04, -5.817d-06, +8.836d-07,
     . -1.769d-07, +7.313d-01, -1.188d-01, +1.145d-02, +1.011d-03,
     . +1.083d-04, +2.570d-06, -2.140d-06, -5.710d-08, +2.000d-08,
     . -1.632d+00, -6.948d-03, -3.893d-03, +8.592d-04, +7.577d-05,
     . +4.539d-06, -3.852d-07, -2.213d-07, -1.370d-08, +5.800d-09/
       
      data (bw_amp(i),i=1,55)/
     . +0.000d+00, +0.000d+00, -8.865d-02, +0.000d+00, -4.309d-01,
     . +6.340d-02, +0.000d+00, +1.162d-01, +6.176d-02, -4.234d-03,
     . +0.000d+00, +2.530d-01, +4.017d-02, -6.204d-03, +4.977d-03,
     . +0.000d+00, -1.737d-01, -5.638d-03, +1.488d-04, +4.857d-04,
     . -1.809d-04, +0.000d+00, -1.514d-01, -1.685d-02, +5.333d-03,
     . -7.611d-05, +2.394d-05, +8.195d-06, +0.000d+00, +9.326d-02,
     . -1.275d-02, -3.071d-04, +5.374d-05, -3.391d-05, -7.436d-06,
     . +6.747d-07, +0.000d+00, -8.637d-02, -3.807d-03, -6.833d-04,
     . -3.861d-05, -2.268d-05, +1.454d-06, +3.860d-07, -1.068d-07,
     . +0.000d+00, -2.658d-02, -1.947d-03, +7.131d-04, -3.506d-05,
     . +1.885d-07, +5.792d-07, +3.990d-08, +2.000d-08, -5.700d-09/

C     parameter t
      t = dsin(dlat)

C     degree n and order m
      n = 9
      m = 9

C determine n!  (faktorielle)  moved by 1
      dfac(1) = 1
      do i = 1,(2*n + 1)
        dfac(i+1) = dfac(i)*i
      end do

C     determine Legendre functions (Heiskanen and Moritz, Physical Geodesy, 1967, eq. 1-62)
      do i = 0,n
        do j = 0,min(i,m)
          ir = int((i - j)/2)
          sum = 0
          do k = 0,ir
            sum = sum + (-1)**k*dfac(2*i - 2*k + 1)/dfac(k + 1)/
     .       dfac(i - k + 1)/dfac(i - j - 2*k + 1)*t**(i - j - 2*k)
          end do
C         Legendre functions moved by 1
          P(i + 1,j + 1) = 1.d0/2**i*dsqrt((1 - t**2)**(j))*sum
        end do
      end do

C     spherical harmonics
      i = 0
      do n = 0,9
        do m = 0,n
          i = i + 1
          aP(i) = P(n+1,m+1)*dcos(m*dlon)
          bP(i) = P(n+1,m+1)*dsin(m*dlon)
        end do
      end do

C     hydrostatic
      bh = 0.0029
      c0h = 0.062
      if (dlat.lt.0) then ! southern hemisphere
        phh  = pi
        c11h = 0.007
        c10h = 0.002
      else                ! northern hemisphere
        phh  = 0
        c11h = 0.005
        c10h = 0.001
      end if
      ch = c0h + ((dcos(doy/365.25d0*2*pi + phh)+1)*c11h/2 + c10h)*
     .           (1-dcos(dlat))

      ahm = 0.d0
      aha = 0.d0
      do i = 1,55
        ahm = ahm + (ah_mean(i)*aP(i) + bh_mean(i)*bP(i))*1d-5
        aha = aha + (ah_amp(i) *aP(i) + bh_amp(i) *bP(i))*1d-5
      end do
      ah  = ahm + aha*dcos(doy/365.25d0*2.d0*pi)

      sine   = dsin(pi/2 - zd)
      beta   = bh/( sine + ch  )
      gamma  = ah/( sine + beta)
      topcon = (1.d0 + ah/(1.d0 + bh/(1.d0 + ch)))
      gmfh   = topcon/(sine+gamma)

C     height correction for hydrostatic mapping function from Niell (1996)
      a_ht = 2.53d-5
      b_ht = 5.49d-3
      c_ht = 1.14d-3
      hs_km  = dhgt/1000.d0
C
      beta   = b_ht/( sine + c_ht )
      gamma  = a_ht/( sine + beta)
      topcon = (1.d0 + a_ht/(1.d0 + b_ht/(1.d0 + c_ht)))
      ht_corr_coef = 1/sine - topcon/(sine + gamma)
      ht_corr      = ht_corr_coef * hs_km
      gmfh         = gmfh + ht_corr

C     wet
      bw = 0.00146
      cw = 0.04391

      awm = 0.d0
      awa = 0.d0
      do i = 1,55
        awm = awm + (aw_mean(i)*aP(i) + bw_mean(i)*bP(i))*1d-5
        awa = awa + (aw_amp(i) *aP(i) + bw_amp(i) *bP(i))*1d-5
      end do
      aw =  awm + awa*dcos(doy/365.25d0*2*pi)

      beta   = bw/( sine + cw )
      gamma  = aw/( sine + beta)
      topcon = (1.d0 + aw/(1.d0 + bw/(1.d0 + cw)))
      gmfw   = topcon/(sine+gamma)
      
      end 
