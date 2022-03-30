CTITLE READ_OCEAN
 
      subroutine read_ocean ( ocean_file, oceamp, ocephs,
     .                    site_names, num_sites)

      implicit none 
 
*     This routine will open the ocean file and read in the
*     coeffiecients for the current sites in the solutions.  If the
*     file can be found then an error message is printed and the
*     amplitudes are set to zero.
 
 
* PASSED VARIABLES
 
*   num_sites   - Number of sites
 
      integer*4 num_sites
 
*   oceamp(11,3,num_sites)  - Amplitudes in mm (order of
*               - components is NEU (different from CALC but the
*               - other defintions are the same) (m)
*   ocephs(11,3,num_sites)  - Phases as defined in CALC (deg)
 
      real*8 oceamp(11,3,num_sites), ocephs(11,3,num_sites)
 
*   site_names(num_sites)       - Names of the sites
*   ocean_file  - Name of the file with coefficients (parameter
*               - kalman_param.h)
 
      character*(*) site_names(num_sites), ocean_file
 
* LOCAL VARIABLES
 
*   ierr, jerr  - IOSTAT error on file open and reads, and for
*               - decoding values
*   i,j,k       - Loop counters
*   indx        - Pointer in string
*   iel         - Station number for line just read.
*   found       - Bit mapped word with bit set as site is found
*               - (NOTE: Limited to 32 sites in this version)
 
      integer*4 ierr, jerr, i,j,k, indx, iel, found
 
*   vals(11)        - The 11 values read from data file line.  We
*               - then decide what to do with them.
 
      real*8 vals(11)
 
*   finished        - Indicates that we have found all of the
*               - coefficients
 
      logical finished
 
*   line        - Line read from input file.
 
 
      character*132 line
 
***** CLEAR all of the ampltitudes first
 
      do i = 1, num_sites
          do j = 1, 11
              do k = 1,3
                  oceamp(j,k,i) = 0.d0
                  ocephs(j,k,i) = 0.d0
              end do
          end do
      end do
 
****  Try to open the input file
      open(102, file=ocean_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',ocean_file,
     .                   0, 'READ_OCEAN')
 
****  There is not much can do, get out
      if( ierr.ne.0 ) RETURN
 
****  Start reading the file and try to find stations
 
      found = 0
      finished = .false.
 
      do while ( .not.finished )
 
          read(102,'(a)', iostat=ierr ) line
          if( ierr.eq.0 .and. line(1:2).eq.'  ' ) then
 
*             See if station is one of ours (Name is assummed to start
*             in column 3.)
              indx = 3
              call get_cmd( line, site_names, num_sites, iel, indx)
 
*****         See if we found station
              if( iel.gt.0 ) then
*                             ! Counter for number of lines read
                  i = 0
                  do while (i.le.5 .and. ierr.eq.0 )
                      read(102,'(a)', iostat=ierr) line
                      if( line(1:1).eq.' ' .and. ierr.eq.0 ) then
                          i = i + 1
                          read(line,*,iostat=jerr) vals
                          call report_error('IOSTAT',jerr,'read',
     .                        line, 0, 'READ_OCEAN')
 
*                         Clear the values if there is an error
                          if( jerr.ne.0 ) then
                              do j = 1,11
                                  vals(j) = 0.d0
                              end do
                          end if
 
*                         Now assign into the correct positions
                          if( i.eq.1 ) then
                              call cpoc(oceamp(1,3,iel), vals)
                          else if ( i.eq.2 ) then
                              call cpoc(oceamp(1,2,iel), vals)
                          else if ( i.eq.3 ) then
                              call cpoc(oceamp(1,1,iel), vals)
                          else if ( i.eq.4 ) then
                              call cpoc(ocephs(1,3,iel), vals)
                          else if ( i.eq.5 ) then
                              call cpoc(ocephs(1,2,iel), vals)
                          else if ( i.eq.6 ) then
                              call cpoc(ocephs(1,1,iel), vals)
                          end if
                          call sbit(found,iel,1)
*                                 ! Line not comment
                      end if
*                                 ! Finding all the lines
                  end do
*                                 ! This is not one our sites
              else
 
*                 Loop and pull of the data lines
                  i = 0
                  do while ( i.le.5 .and. ierr.eq.0 )
                      read(102, '(a)', iostat=ierr) line
                      if( line(1:1).eq.' ') i = i + 1
                  end do
              end if
*                         ! Not EOF and looking for station name
          end if
 
****      See if we have found all of the stations
          if( found.eq.2**num_sites-1 .or. ierr.ne.0 ) then
              finished = .true.
          end if
 
*                         ! Looping over data file
      end do
 
****  Close the input file and return
      close(102)
      return
      end
 
CTITLE OCEAN_LOAD
 
      subroutine ocean_load( site, epoch, site_pos, site_part, oceamp,
     .                    ocephs, ocean_rad, ocean_horz )
 
*     This routine will compute the ocean loading contributions to
*     to delay and rate for each of the sites in site.
 
      include '../includes/const_param.h'
 
* PASSED VARIABLES
 
*   site(2)     - Numbers of the sites in the baseline
 
      integer*4 site(2)
 
*   site_pos(3,*)   - Site positions XYZ (m)
*   oceamp(11,3,*)  - Amplitudes of the ocean loading
*   ocephs(11,3,*)  - Phases of the ocean loading
*   epoch           - JD of current measurement.
 
      real*8 site_pos(3,*), oceamp(11,3,*), ocephs(11,3,*), epoch
 
*   site_part(3,2,2)    - Site partial derivatives by XYZ, delay
*                   - and rate, and by site 1/2.
*   ocean_rad(2,2)  - Contribution of radial load to site 1/2
*                   - for delay and rate (ps, fs/s). Signs such
*                   - that contribution is subtracted for first
*                   - site and added for second
*   ocean_horz(2,2) - Contribution of horizontal load to
*                   - site 1/2 for delay and rate (ps, fs/s)
 
      real*4 site_part(3,2,2), ocean_rad(2,2), ocean_horz(2,2)
 
* LOCAL VARIABLES
 
*   i,j,k,l     - Loop counters
*   isgn        - Sign for contribution.
 
      integer*4 i,j,k,l, isgn
 
*   xjd             - Julian date at 0:00 UT
*   ut1             - UT of measurement (seconds)
*   oceneu(3,2)     - NEU values and rate of the ocean load
*   ocexyz(3,2)     - XYZ values and rate of the ocean load
*   angle(11)       - Angular arguments of the 11 tides
*   dangl(11)       - Rates of change of the tide arguments
*   loc_coord(3)        - Local coordinates (used by rotate_geod)
*   rot_matrix(3,3) - Rotation matrix from NEU to XYZ
*   radial(3,2)     - Radial tide correction
*   horz(3,2)       - Horizontal tide correction
 
      real*8 xjd, ut1, oceneu(3,2), ocexyz(3,2), angle(11), dangl(11),
     .    loc_coord(3), rot_matrix(3,3), radial(3,2), horz(3,2)
 
      DATA dangl  / 1.3557D-04, 1.4544D-04, 1.3026D-04, 1.4624D-04,
     1              7.312D-05,  6.245D-05,  7.232D-05,  5.714D-05,
     2              1.067D-05,  5.289D-06,  7.964D-07 /
      
      isgn = 0
 
****  Start by splitting the ut1 from the JD
 
      xjd = int( epoch - 0.5d0) + 0.5d0
      ut1 = (epoch - xjd ) * 86400.d0
 
*     Get the ocean loading arguments
      call ocearg( ut1, xjd, angle )
 
*     Now loop over each site
      do l = 1, 2
 
          if( l.eq.1 ) isgn = -1
          if( l.eq.2 ) isgn = +1
 
*         Loop over the components
          do k = 1,3
*                                 ! value
              oceneu(k,1) = 0.d0
*                                 ! Rate of change
              oceneu(k,2) = 0.d0
 
*             Now sum contributions
              do j = 1,11
                  oceneu(k,1) = oceneu(k,1) +
     .                    oceamp(j,k,site(l)) *
     .                    cos(angle(j)-ocephs(j,k,site(l))*pi/180)
                  oceneu(k,2) = oceneu(k,2) -
     .                    oceamp(j,k,site(l)) * dangl(j) *
     .                    sin(angle(j)-ocephs(j,k,site(l))*pi/180)
              end do
          end do
 
****      Now change the sign of the horizontal terms due to BUG in
*         original data files
          oceneu(1,1) = -oceneu(1,1)
          oceneu(2,1) = -oceneu(2,1)
          oceneu(1,2) = -oceneu(1,2)
          oceneu(2,2) = -oceneu(2,2)
 
*         Now save the radial and horizontal separately
          do i = 1,2
              radial(1,i) = 0.d0
              radial(2,i) = 0.d0
              radial(3,i) = oceneu(3,i)
 
              horz(1,i)   = oceneu(1,i)
              horz(2,i)   = oceneu(2,i)
               horz(3,i)   = 0.d0
           end do
 
*****     Now convert these values to XYZ so that we can use site
*         partials to get contribution: Do Radial first
 
          call rotate_geod(radial(1,1), ocexyz(1,1), 'NEU', 'XYZ',
     .            site_pos(1,site(l)), loc_coord, rot_matrix)
          call rotate_geod(radial(1,2), ocexyz(1,2), 'NEU', 'XYZ',
     .            site_pos(1,site(l)), loc_coord, rot_matrix)
 
****      Now compute the radial tide contribution: Delay
          ocean_rad(l,1) = (ocexyz(1,1)*site_part(1,l,1) +
     .                    ocexyz(2,1)*site_part(2,l,1) +
     .                    ocexyz(3,1)*site_part(3,l,1) ) * isgn
*                                               Rate
          ocean_rad(l,2) = (ocexyz(1,1)*site_part(1,l,2) +
     .                    ocexyz(2,1)*site_part(2,l,2) +
     .                    ocexyz(3,1)*site_part(3,l,2) +
     .                   (ocexyz(1,2)*site_part(1,l,1) +
     .                    ocexyz(2,2)*site_part(2,l,1) +
     .                    ocexyz(3,2)*site_part(3,l,1)) * 1000.d0)*
     .                    isgn
 
*****     Now convert these values to XYZ so that we can use site
*         partials to get contribution: Do Horizontal next
 
          call rotate_geod(horz(1,1), ocexyz(1,1), 'NEU', 'XYZ',
     .            site_pos(1,site(l)), loc_coord, rot_matrix)
          call rotate_geod(horz(1,2), ocexyz(1,2), 'NEU', 'XYZ',
     .            site_pos(1,site(l)), loc_coord, rot_matrix)
 
****      Now compute the radial tide contribution: Delay
          ocean_horz(l,1) = (ocexyz(1,1)*site_part(1,l,1) +
     .                    ocexyz(2,1)*site_part(2,l,1) +
     .                    ocexyz(3,1)*site_part(3,l,1)) * isgn
*                                               Rate
          ocean_horz(l,2) = (ocexyz(1,1)*site_part(1,l,2) +
     .                    ocexyz(2,1)*site_part(2,l,2) +
     .                    ocexyz(3,1)*site_part(3,l,2) +
     .                   (ocexyz(1,2)*site_part(1,l,1) +
     .                    ocexyz(2,2)*site_part(2,l,1) +
     .                    ocexyz(3,2)*site_part(3,l,1)) * 1000.d0)*
     .                    isgn
 
*     End looping over the sites
      end do
 
****  Thats all.  Contributions have been computed
      return
      end
 
CTITLE CPOC
 
      subroutine cpoc( ocent, vals )
 
*     Routine to copy the vals from the line in the ocean file into
*     the correct place in the ocean amp or phase arrays
 
*   ocent(11)   - The correct entry for the copy
*   vals(11)        - Values read from data file
 
      real*8 ocent(11), vals(11)
 
*   i           - Loop counter
 
      integer*4 i
 
****  Copy the values accross
      do i = 1,11
          ocent(i) = vals(i)
      end do
 
***** Thats all
      return
      end
 
      SUBROUTINE ocearg (UT1, XJD, ANGLE)
C
C   1      ocearg
C
C   1.1    ocearg PROGRAM SPECIFICATIONS
C
C   1.1.1  THIS SUBROUTINE COMPUTES THE ANGULAR ARGUMENTS
C          FOR SCHWIDERSKI COMPUTATION OF 11 OCEAN TIDES.
C
C          C A U T I O N
C          = = = = = = =
C
C          SCHWIDERSKI MODIFIES THE ANGULAR ARGUMENTS OF THE DIURNAL
C          TERMS BY +/- 90 DEGREES. THEREFORE HIS DIURNAL PHASES
C          CANNOT BE USED WITH THE STANDARD DOODSEN OR CARTWRIGHT
C          CONVENTIONS.
C
C
C   1.1.2  RESTRICTIONS  -  NONE
C
C   1.1.3  REFERENCES  -  MERIT STANDARS, APRIL 81, SECOND DRAFT,
C                         APPENDICES 7,11
C
C   1.2.   ocearg PROGRAM INTERFACE
C
C   1.2.1  CALLING SEQUENCE
C
C          INPUT VARIABLES  -
C
C          1.  UT1     -  THE UT1 TIME OF THE DAY. (SEC)
C
C          2.  XJD     -  THE JULIAN DATE AT ZERO HOURS UTC OF THE
C                         DATE IN QUESTION. (DAYS)
C
C          OUTPUT VARIABLES  -
C
C          1.  ANGLE(11) -THE ANGULAR ARGUMENTS FOR SCHWIDERSKI
C                         COMPUTATION OF THE OCEAN TIDES IN THE ORDER:
C                         M2, S2, N2, K2, K1, O1, P1, Q1, MF, MM, SSA
C
C   1.2.2. COMMON BLOCKS USED
C
C
      include '../includes/const_param.h'
CC
C
C  1.2.3   PROGRAM SPECIFICATIONS
C
      REAL*8 ANGLE(11), SPEED(11), convd, ut1, xjd 
 
      REAL*4 ANGFAC(4,11)
 
C
C  1.2.4   DATA BASE ACCESS  -  NONE
C
C  1.2.5   EXTERNAL INPUT/OUTPUT  -  POSSIBLE DEBUG
C
C  1.2.6   SUBROUTINE INTERFACE
C
C          CALLER SUBROUTINE  -  OCEG
C
C          CALLED SUBROUTINE  -  DMOD
C
C  1.2.7   CONSTANTS USED  -
C
C          1.  CONVD
C
C          2.  TWOPI
C
C          3.  ANGFAC(4,11) -  TABLE OF MULTIPLES OF ARGUMENTS.
C                              (UNITLESS)
C
C          4.  CENTJ        -  THE NUMBER OF JULIAN DAYS PER JULIAN
C                              CENTURY. (DAYS/CENT.) (CENTJ=36525.D0)
C
C          5.  SPEED(11)    -  COEFFICIENTS FOR THE COMPUTATION OF
C                              ARGUMENTS IN THE FOLLOWING ORDER OF TIDES
C                              M2, S2, N2, K2, K1, O1, P1, Q1, MF, MM,
C                              SSA. (RAD/SEC)
C
C          6.  XJD75        -  THE JULIAN DATE OF JAN 0.0 1975. (DAYS)
 
      real*8 centj, xjd75
      integer*4 k
 
C
      DATA ANGFAC / 2.D0, -2.D0, 0.D0, 0.D0,
     2              0.D0,  0.D0, 0.D0, 0.D0,
     3              2.D0, -3.D0, 1.D0, 0.D0,
     4              2.D0,  0.D0, 0.D0, 0.D0,
     5              1.D0,  0.D0, 0.D0, 0.25D0,
     6              1.D0, -2.D0, 0.D0,-0.25D0,
     7             -1.D0,  0.D0, 0.D0,-0.25D0,
     8              1.D0, -3.D0, 0.D0,-0.25D0,
     9              0.D0,  2.D0, 1.D0, 0.D0,
     A              0.D0,  1.D0,-1.D0, 0.D0,
     1              2.D0,  0.D0, 0.D0, 0.D0 /
C
      DATA CENTJ / 36525.D0 /
C
      DATA SPEED / 1.40519D-4, 1.45444D-4, 1.37880D-4,
     1             1.45842D-4, 0.72921D-4, 0.67598D-4,
     2             0.72523D-4, 0.64959D-4, 0.053234D-4,
     3             0.026392D-4, 0.003982D-4 /
C
      DATA XJD75 / 2442412.5D0 /
C
C
C
C  1.2.8   PROGRAM VARIABLES
C
C          1.  CAPT  -  THE NUMBER OF JULIAN CENTURIES BETWEEN JAN 0.5, 1900
C                       AND THE OBSERVATION
C
C          2.  FDAY  -  FRACTIONAL PART OF UNIVERSAL TIME (UT1) DAY IN
C                       SECONDS. (SEC)
C
C          3.  H0    -  MEAN LONGITUDE OF SUN AT BEGINNING OF DAY. (RAD)
C
C          4.  ICAPD -  JULIAN DAYS SINCE JANUARY 0.0 UT 1975. (DAYS)
C
C          5.  P0    -  MEAN LONGITUDE OF LUNAR PERIGEE AT BEGINNING
C                       OF DAY. (RAD)
C
C          6.  S0    -  MEAN LONGITUDE OF MOON AT BEGINNING OF DAY.
C                       (RAD)
 
      real*8 capt, fday, h0, icapd, p0, s0
 
C
C  1.2.9   PROGRAMMER  -  CLYDE GOAD
C          83:10:08 HAROLD M. SCHUH
C          89:10:08 Jim Ryan ANGFAC and SPEED moved to ema.
C
C  1.3     ocearg PROGRAM STRUCTURE
C
C     Compute FDAY  -  fractional part of UT1 day in seconds.
      FDAY = UT1
      convd = pi/180.d0
C
C     Compute ICAPD  -  Days since JAN 0, 0.00 UT, 1975
C     and  CAPT   -  Julian centuries since JAN 0.5, 1900.
      ICAPD = XJD - XJD75
      CAPT = ( 27392.500528D0 + 1.000000035D0 * ICAPD ) / CENTJ
C
C     Compute mean longitude of sun at beginning of day.
      H0 = (279.69668D0 + (36000.768930485D0 + 3.03D-4 * CAPT) * CAPT)
     1     * CONVD
C
C     Compute mean longitude of moon at beginning of day.
      S0 = ((( 1.9D-6 * CAPT - 0.001133D0 ) * CAPT + 481267.88314137D0)
     1     * CAPT + 270.434358D0 ) * CONVD
C
C     Compute mean longitude of lunar perigee at beginning of day.
      P0 = ((( -1.2D-5 * CAPT - .010325D0 ) * CAPT + 4069.0340329577D0)
     1     * CAPT + 334.329653D0 ) * CONVD
C
C     Calculate the angular arguments.
C     Run a loop over the 11 main tides.
C
      DO K = 1,11
        ANGLE(K) = SPEED(K)*FDAY + ANGFAC(1,K)*H0 + ANGFAC(2,K)*S0
     1             + ANGFAC(3,K)*P0 + ANGFAC(4,K)*2*pi
C
        ANGLE(K) = MOD( ANGLE(K), 2*pi )
        IF( ANGLE(K) .LT. 0.D0 ) ANGLE(K) = ANGLE(K) + 2*pi
      Enddo
C
C     Normal termination.
      RETURN
      END
 
