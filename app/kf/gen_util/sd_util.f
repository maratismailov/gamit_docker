 
CTITLE INIT_SD
 
      subroutine init_sd ( sd_file, apply_sd, print_mess ,
     .                     site_names, num_sites, mid_epoch )

      implicit none 
 
*     Routine to initialize the applying of diurnal and semidiurnal
*     UT1, pole and tide corrections.  This initial versions tries
*     to open the SD_DATA file.  If it can then its reads the values
*     and saves in common for use by the add_sd routine
 
      include '../includes/kalman_param.h'
c     include '../includes/obs_header.h'
c     include '../includes/obs_data.h'
      include '../includes/sd_common.h'
 
* PASSED VARIABLES

*   num_sites   - Number of sites

      integer*4 num_sites

*   mid_epoch   - Mid epoch of experiment

      real*8 mid_epoch
 
*   apply_sd    - Indicates we should apply sd data
*   print_mess  - Indicates that a message to used should be printed
 
      logical apply_sd, print_mess

*   sd_file     - Name of the file with the semi and diurnal
*                 corrections
*   site_names  - Name of the sites

      character*(*) sd_file, site_names(num_sites)
 
* LOCAL VARIBALES
 
*   ierr        - IOSTAT error
*   indx        - pointer to position in string
*   idum        - dummy for multiread and read_line
*   iel         - Pointer for getting station name
*   trimlen     - Length of string
 
      integer*4 ierr, indx, idum, iel, trimlen
 
*   eof         - Indicates eof of file reached
 
      logical eof
 
*   line        - Line read from input file
 
      character*128 line
 
*   sd_site_name    - Name of site read from file for the tide
*               - application
*   type        - Type of correction given. OPtions are UT1,
*               - XY and TIDES.
*   cdum        - Dummy character for multiread and read_line
 
      character*8 sd_site_name, type
      character*1 cdum

****  Try to open file
      apply_sd = .false.
 
      sd_ut1_num = 0
      sd_xy_num  = 0
      sd_tides_num = 0
 
      if( trimlen(sd_file).gt.0 ) then
          open(150, file = sd_file, iostat=ierr, status='old')
*                                 ! Open OK, Proceed
      else
          ierr = -99
      end if
      if( ierr.eq.0 ) then

          eof = .false.
          do while ( .not.eof )
              read(150,'(a)', iostat=ierr) line

*             See if none comment line
*                                                               ! Decode
              if( line(1:1).eq.' ' .and. ierr.eq.0 ) then
                  indx = 1
                  call GetWord(line, type, indx)
                  call casefold(type)
*                 Check type
                  if( type(1:3).eq.'UT1' ) then
                      sd_ut1_num = sd_ut1_num + 1
                      call check_sd_max(sd_ut1_num, max_sd_ut1,
     .                     'UT1')
 
                      call multiread(line,indx,'I4',ierr,
     .                    sd_ut1_arg(1,sd_ut1_num), cdum, 6)
                      call multiread(line,indx,'R8',ierr,
     .                    sd_ut1_val(1,sd_ut1_num), cdum, 2)
                  end if
 
                  if( type(1:2).eq.'XY' ) then
                      sd_xy_num = sd_xy_num + 1
                      call check_sd_max(sd_xy_num, max_sd_xy,
     .                     'xy pole postions')
 
                      call multiread(line,indx,'I4',ierr,
     .                    sd_xy_arg(1,sd_xy_num), cdum, 6)
                      call multiread(line,indx,'R8',ierr,
     .                    sd_xy_val(1,sd_xy_num), cdum, 2)
                  end if

                  if( type(1:2).eq.'NU' ) then
                      sd_nut_num = sd_nut_num + 1
                      call check_sd_max(sd_nut_num, max_sd_nut,
     .                     'nutation')

                      call multiread(line,indx,'I4',ierr,
     .                    sd_nut_arg(1,sd_nut_num), cdum, 6)
                      call multiread(line,indx,'R8',ierr,
     .                    sd_nut_val(1,sd_nut_num), cdum, 4)
                  end if

 
                  if( type(1:2).eq.'TI' ) then
                      sd_tides_num = sd_tides_num + 1
                      call check_sd_max(sd_tides_num, max_sd_tides,
     .                     'tides')
 
                      call multiread(line,indx,'I4',ierr,
     .                    sd_tides_arg(1,sd_tides_num), cdum, 6)
                      call multiread(line,indx,'R8',ierr,
     .                    sd_tides_val(1,sd_tides_num), cdum, 6)
                      call read_line(line,indx,'CH',ierr, idum,
     .                    sd_site_name)
*                     Get the number of the site corresponding
*                     to the site name.  If name is all use 99.
                      call casefold(sd_site_name)
                      if( sd_site_name(1:3).eq.'ALL' ) then
                          sd_tides_site(sd_tides_num) = 99
*                                  ! See if we can find the site
                      else
                          indx = 1
                          call get_cmd(sd_site_name, site_names,
     .                         num_sites, iel, indx)
*                                  ! If site found, then save
                          if( iel.gt.0 ) then
                             sd_tides_site(sd_tides_num) = iel
*                                  ! forget about this tide, does not
*                                  ! apply to any of our sites.
                          else
                             sd_tides_num = sd_tides_num-1
                          end if
                      end if
                  end if
*                             ! See if EOF
              else
                  if ( ierr.ne.0 ) eof = .true.
              end if
*                             ! Looping over file
          end do
          close(150)
 
*         See if we have values
          if( sd_ut1_num+sd_xy_num+sd_nut_num+
     .        sd_tides_num.gt.0                ) then
              apply_sd = .true.
          end if
*                             ! File open OK
      end if

****  Print message
      if( apply_sd .and. print_mess ) then
          write(*,200) sd_ut1_num, sd_xy_num, sd_nut_num,
     .                 sd_tides_num
 200      format(' Apply corrections from SD_FILE: ',I3,
     .           ' UT1 values,  ',I3,' XY pole, ',/,
     .    33x,I3,' Nutation and ',i3,' Tides')
      end if

*     Now evalute the values at mid_epoch (if no file then values are
*     zeroed.

      call get_sd_mid ( mid_epoch, apply_sd, num_sites )
 
****  Thats all
      return
      end
 
CTITLE ADD_SD
 
      subroutine add_sd
 
*     This routine will apply the semi-diurnal and diurnal corrections
*     for UT1, pole and tides to the theoretical model.  APPY_SD should
*     be checked before this routine is called.
*     ****Currently only UT1 and XY implemented
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
      include '../includes/sd_common.h'
 
*   i,j,k       - Loop counters
 
      integer*4 i,j,k
 
*   arg         - Angulatr argument for the correction (rads)
*   dt(2)       - Delay and rate correction (ps and fs/s)
*   fund_arg(6) - Values of the 5 Brown's arguments (l,l',F,D,
*               - Omega) and gst+pi (rads)
*   xy_part(2)  - Polar motion partial derivatives
*   dt2000      - Time since dj2000 
 
      real*8 arg, dt(2), fund_arg(6), xy_part(2), dt2000

*   tides_val(12)   - Values for the spherical harmonic
*                   - coefficients for the 12 extended
*                   - tidal displacements.

      real*8 tides_val(12)
 
***** Get the fundamental arguments at this epoch
 
      call fund_angles( epoch, fund_arg)
 
*     Clear the corrections
*                     ! Delay and rate values
      do j = 1,2
          dt(j) = 0.d0
      end do
 
*     Now loop over the UT1 contributions
      do i = 1, sd_ut1_num
 
*         Get the argument
          call sd_arg( sd_ut1_arg(1,i), fund_arg, arg, 6)
 
*         Compute delay and rate contriubutioin
*                         ! Loop on delay and rate
          do j = 1, 2
              dt(j) = dt(j) + ( sd_ut1_val(1,i)*cos(arg) +
     .               sd_ut1_val(2,i)*sin(arg) )*pmu_part(3,j)
          end do
      end do
 
****  Now do polar motion
      do i = 1, sd_xy_num
 
*         Get the argument
          call sd_arg( sd_xy_arg(1,i), fund_arg, arg, 6)
 
*         Compute delay and rate contriubutioin
*                         ! Loop on delay and rate
          do j = 1, 2
              xy_part(1) = -pmu_part(1,j)* cos(arg) +
     .                      pmu_part(2,j)* sin(arg)
              xy_part(2) =  pmu_part(2,j)* cos(arg) +
     .                      pmu_part(1,j)* sin(arg)
 
              dt(j) = dt(j) + ( sd_xy_val(1,i)*xy_part(1) +
     .               sd_xy_val(2,i)*xy_part(2) )
 
          end do
      end do

*     Now loop over the nut contributions
      do i = 1, sd_nut_num

*         Get the argument
          call sd_arg( sd_nut_arg(1,i), fund_arg, arg, 5)

*         Compute delay and rate contriubutioin
*                         ! Loop on delay and rate
*         check for FCN mode (6 th argument value non-zero)
          dt2000 = epoch - dj2000
          if( sd_nut_arg(6,i).ne.0 ) then
              arg = (dt2000/429.8)*2.d0*pi
          end if

*         check for prcession constant. (arguments zero)
          if( sd_nut_arg(1,i) + sd_nut_arg(2,i) + sd_nut_arg(3,i) +
     .        sd_nut_arg(3,i) + sd_nut_arg(4,i) + sd_nut_arg(5,i) +
     .        sd_nut_arg(6,i).eq.0 ) then
              arg = dt2000/365.25d0

*             Evaluate precession constant contribution (NOTE: we
*             also allow an offset in the positions so that we dont
*             move the coordinate system too much.)
              do j = 1, 2
                  dt(j) = dt(j) + ( sd_nut_val(1,i)*arg    + 
     .                   sd_nut_val(3,i) )* nut_part(1,j)  +
     .                            ( sd_nut_val(2,i)*arg    +
     .                   sd_nut_val(4,i) )* nut_part(2,j) 
              end do
*                      ! Normal nutation compute in regular fashion
          else 
              do j = 1, 2
                  dt(j) = dt(j) + ( sd_nut_val(1,i)*sin(arg)       +
     .                   sd_nut_val(3,i)*cos(arg) )* nut_part(1,j) +
     .                            ( sd_nut_val(2,i)*cos(arg)       +
     .                   sd_nut_val(4,i)*sin(arg) )* nut_part(2,j) 
              end do
          end if
      end do

***** Compute the tidal contribtions.  Add contributions for both sites.
*                     ! Loop over the two sites
      do i = 1,2

*         Get the values of the tides for this site
          call sd_get_tides(site(i), fund_arg, tides_val)

*         Now compute contribution
*                         ! Loop over 12 components
          do j = 1,12
*                         ! Loop over delay and rate
              do k = 1,2
                  dt(k) = dt(k) + tides_val(j)*etd_ext_part(j,i,k)
              end do
          end do
      end do

 
****  Now add the corrections to the theoretical model
      do j = 1, 3
          theoretical(j) = theoretical(j) + dt(1)
      end do

      theoretical(4) = theoretical(4) + dt(2)
 
***** Thats all
      return
      end
 
CTITLE SD_ARG
 
      subroutine sd_arg( sd_mult, fund_arg, arg, num_arg )
 
*     This routine computes the argumant of the tide given the
*     the fundamental argument values and the multipliers to
*     generate the tide
 
      include '../includes/const_param.h'
 
* PASSED VARIABLES
 
*   sd_mult(6)  - multipliers for Browns arguments and gst+pi
*   num_arg     - Number of arguments to be summed. (Allows the
*                 gst+pi argument to be skipped.)
 
      integer*4 sd_mult(6), num_arg
 
*   fund_arg(6) - Values for the fundamental artguments
*   arg         - Argumnent at this time
 
      real*8 fund_arg(6), arg
 
* LOCAL VRAIABLES
 
*   i           - Loop counter
 
      integer*4 i
 
****  Initialize argument and combine the fundamental angels with the
*     values passed
 
      arg = 0.d0
      do i = 1, num_arg
          arg = arg + sd_mult(i)*fund_arg(i)
      end do
 
      arg = mod(arg, 2*pi)
 
***** Thats all
      return
      end
 
 
CTITLE 'fund_angles'
 
      subroutine fund_angles( epoch, fund_arg )
 
 
*     Routine to compute the value of the fundamental argument
*     for Brown's arguments.  The sixth entry is returned as GST
*     push pi.  The additional pi is needed for compatability with
*     Doodson's Tide argument.
 
      include '../includes/const_param.h'
 
* epoch  - Julian date for arguments
* fund_arg(6) -  Brown's arguments plus GST+pi
 
 
      real*8 epoch, fund_arg(6)
 
*      cent             - Julian centuries to DJ2000.
*      el,eld           - Mean longitude of moon minus mean
*                       - longitude of moon's perigee (arcsec)
*      elc(5)           - Coefficients for computing el
*      elp,elpd         - Mean longitude of the sun minus mean
*                       - longitude of sun perigee (arcsec)
*      elpc(5)          - Coeffiecents for computing elp
*      f,fd             - Moon's mean longitude minus omega (sec)
*      fc(5)            - Coefficients for computing f
*      d,dd             - Mean elongation of the moon from the
*                       - sun (arcsec)
*      dc(5)            - coefficients for computing d
*      om,omd           - longitude of the ascending node of the
*                       - moon's mean orbit on the elliptic
*                       - measured from the mean equinox of date
*      omc(5)           - Coefficients for computing om.
*      gst              - Greenwich mean sidereal time (rad)
 
 
 
      real*8 cent, el,eld, elc(5), elp, elpd, elpc(5),
     .    f,fd, fc(5), d,dd, dc(5), om,omd, omc(5), gst
 
      data elc    /     0.064d0,    31.310d0,    715922.633d0,
     .             485866.733d0,    1325.0d0 /
      data elpc   /    -0.012d0,    -0.577d0,   1292581.224d0,
     .            1287099.804d0,      99.0d0 /
      data fc     /     0.011d0,   -13.257d0,    295263.137d0,
     .             335778.877d0,    1342.0d0/
      data dc     /     0.019d0,    -6.891d0,    1105601.328d0,
     .            1072261.307d0,    1236.0d0/
      data omc    /     0.008d0,     7.455d0,    -482890.539d0,
     .             450160.280d0,      -5.0d0/
 
****  Get the number of centuries to current time
 
      cent = (epoch-dj2000) / 36525.d0
 
****  Compute angular arguments
      el = elc(1) * cent**3 + elc(2) * cent**2 + elc(3) * cent
     .          + elc(4) + mod( elc(5) * cent, 1.d0 ) * sec360
      el = mod( el, sec360 )
      eld = 3.d0 * elc(1) * cent**2 + 2.d0 * elc(2) * cent + elc(3)
     .      + elc(5) * sec360
c
      elp = elpc(1) * cent**3 + elpc(2) * cent**2 + elpc(3) * cent
     .     + elpc(4) + mod( elpc(5) * cent, 1.d0 ) * sec360
      elp = mod( elp, sec360 )
      elpd = 3.d0 * elpc(1) * cent**2 + 2.d0 * elpc(2) * cent + elpc(3)
     .       + elpc(5) * sec360
c
      f = fc(1) * cent**3 + fc(2) * cent**2 + fc(3) * cent
     .     + fc(4) + mod( fc(5) * cent, 1.d0 ) * sec360
      f = mod( f, sec360 )
      fd = 3.d0 * fc(1) * cent**2 + 2.d0 * fc(2) * cent + fc(3)
     .     + fc(5) * sec360
c
      d = dc(1) * cent**3 + dc(2) * cent**2 + dc(3) * cent
     .     + dc(4) + mod( dc(5) * cent, 1.d0 ) * sec360
      d = mod( d, sec360 )
      dd = 3.d0 * dc(1) * cent**2 + 2.d0 * dc(2) * cent + dc(3)
     .     + dc(5) * sec360
c
      om = omc(1) * cent**3 + omc(2) * cent**2 + omc(3) * cent
     .     + omc(4) + mod( omc(5) * cent, 1.d0 ) * sec360
      om = mod( om, sec360 )
      omd = 3.d0 * omc(1) * cent**2 + 2.d0 * omc(2) * cent + omc(3)
     .      + omc(5) * sec360
c
*     Get GST
      call gst_jd( epoch, gst)
 
*     Now save the values.  Convert values from arcseconds to radians
      fund_arg(1) = el / (3600.d0*rad_to_deg)
      fund_arg(2) = elp/ (3600.d0*rad_to_deg)
      fund_arg(3) = f  / (3600.d0*rad_to_deg)
      fund_arg(4) = d  / (3600.d0*rad_to_deg)
      fund_arg(5) = om / (3600.d0*rad_to_deg)
      fund_arg(6) = gst + pi
 
***** Thats all
      return
      end

CTITLE CHECK_SD_MAX

      subroutine check_sd_max(sd_num, sd_max, type)

*     Routine to check the limits on the semidiurnal and diurnal
*     signals anwd warn user if two large.  If two large the last
*     value is over written.

      integer*4 sd_num, sd_max

      character*(*) type

****  Write meassage to std out
      if( sd_num.gt.sd_max ) then
          write(*,100) type, sd_max
 100      format('INIT_SD ERROR: Too many ',a,' values entered.',
     .           ' Max is ', I3,/,
     .            14x,' Last value being overwritten')
          sd_num = sd_max 
      end if

****  Thats all
      return
      end
 
CTITLE GET_SD_MID
 
      subroutine get_sd_mid ( epoch, apply_sd, num_sites )
 
*     Routine to return the accumulated values of the short period
*     ut1, xy and tidal displacements in each of the bands.  There are
*     some special treatments made for nutations.
*     This routine returns the value of the aliased signal at the epoch
*     and therefore the GST argument is ingored.
*     Band 2 -- semidiurnal (GST argument 2)
*     Band 1 -- DIurnal (GST argument 1)
*
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
c     include '../includes/obs_header.h'
      include '../includes/sd_common.h'
 
* PASSED VARIABLE

*   apply_sd - Indicates we should apply the semi and diurnal corrections

      logical apply_sd

*   num_sites - Number of sites

      integer*4 num_sites
 
*   epoch   - Mid epoch of experiment
 
      real*8 epoch
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   ja      - Argument index for GST to define the band (used
*           - for xy pole.  The order here is such that the
*           - array can be directly acessed in OUTSL
*   in      - Index for getting array positions.
 
      integer*4 i,j, ja, in
 
*   arg     - angular argument of each of the terms at eppch
*   fund_arg(6) - Values of the five fundamental arguments plus
*           - GST+pi as the sixth
*   dt2000  - Time from j2000
 
      real*8 arg, fund_arg(6), dt2000 

      ja = 0
 
****  See if we are appying semi-diurnal and diurnal corrections, if
*     we compute in the correct band.  Clear all the values first
      do i = 1,4
          sd_ut1_mid(i) = 0.d0
      end do
 
      do i = 1, 6
          sd_xy_mid(i) = 0.d0
      end do
 
      do i = 1,2
          sd_nut_mid(i) = 0.d0
      end do
 
      do i = 1, 12
          do j = 1, max_sites
              sd_tides_mid(i,j) = 0.d0
          end do
      end do
 
*     Now see if we are applyong
      if( apply_sd ) then
 
*         Compute the fundamental arguments
          call fund_angles( epoch, fund_arg)
 
*         Now do each of the signals in each of the bands.  Start with
*         UT1
 
*                     ! Loop over bands, diurnal and semi
          do i = 1,2
              do j = 1, sd_ut1_num
 
*                 See if in correct band
*                                                     ! Yes it is.
                  if( sd_ut1_arg(6,j).eq.-i ) then
                      call sd_arg( sd_ut1_arg(1,j), fund_arg, arg,5)

*                     compensate for the +pi in the GST argument
                      arg = arg + sd_ut1_arg(6,j)*pi

                      in = (i-1)*2
                      sd_ut1_mid(in+1) = sd_ut1_mid(in+1) +
     .                    sd_ut1_val(1,j)*cos(arg) + 
     .                    sd_ut1_val(2,j)*sin(arg) 
                      sd_ut1_mid(in+2) = sd_ut1_mid(in+2) +
     .                    sd_ut1_val(1,j)*sin(arg) -
     .                    sd_ut1_val(2,j)*cos(arg) 
*                             ! In correct band
                  end if
*                             ! Looping over UT1 signals
              end do
*                             ! Looping over the bands
          end do
 
*****     Now do xy pole position.
*                     ! Loop over bands.  Prograde diurnal, retro
          do i = 1,3
*                     ! semi, and prograde semi
              if( i.eq.1 ) ja =  1
              if( i.eq.2 ) ja = -2
              if( i.eq.3 ) ja =  2
 
              do j = 1, sd_xy_num
 
*                 See if in correct band
*                                                 ! Yes it is.
                  if( sd_xy_arg(6,j).eq.ja) then
                      call sd_arg( sd_xy_arg(1,j), fund_arg, arg,5)
                      arg = arg + sd_xy_arg(6,j)*pi
                      in = (i-1)*2
                      sd_xy_mid(in+1) = sd_xy_mid(in+1) +
     .                    sd_xy_val(1,j)*cos(arg) -
     .                    sd_xy_val(2,j)*sin(arg)
                      sd_xy_mid(in+2) = sd_xy_mid(in+2) +
     .                    sd_xy_val(1,j)*sin(arg) +
     .                    sd_xy_val(2,j)*cos(arg)
                  end if
*                             ! Looping over xy signals
              end do
*                             ! Looping over the bands
          end do
 
*****     Now do nutations.
          do i = 1, sd_nut_num
 
*              Get the argument
               call sd_arg( sd_nut_arg(1,i), fund_arg, arg,5)
 
*              Compute delay and rate contriubutioin
*              check for FCN mode (6 th argument value non-zero)
               dt2000 = epoch - dj2000
               if( sd_nut_arg(6,i).ne.0 ) then
                   arg = (dt2000/429.8)*2.d0*pi
               end if
 
*              check for prcession constant. (arguments zero)
               if( sd_nut_arg(1,i) + sd_nut_arg(2,i) + sd_nut_arg(3,i) +
     .             sd_nut_arg(3,i) + sd_nut_arg(4,i) + sd_nut_arg(5,i) +
     .             sd_nut_arg(6,i).eq.0 ) then
                   arg = dt2000/365.25d0
 
*                  Evaluate precession constant contribution (NOTE: we
*                  also allow an offset in the positions so that we dont
*                  move the coordinate system too much.)
                   sd_nut_mid(1) = sd_nut_mid(1) +
     .                    sd_nut_val(1,i)*arg + sd_nut_val(3,i)
                   sd_nut_mid(2) = sd_nut_mid(2) +
     .                      sd_nut_val(2,i)*arg + sd_nut_val(4,i)
*                    ! Normal nutation compute in regular fashion
              else
                   sd_nut_mid(1) = sd_nut_mid(1)      +
     .                    sd_nut_val(1,i)*sin(arg) +
     .                    sd_nut_val(3,i)*cos(arg)
                   sd_nut_mid(2) = sd_nut_mid(2)      +
     .                    sd_nut_val(2,i)*cos(arg)  +
     .                    sd_nut_val(4,i)*sin(arg)
              end if
          end do

****      Get the tidal values at each of the sites
          do i = 1, num_sites
              call sd_get_tides( i, fund_arg, sd_tides_mid(1,i))
          end do

*****     Thats all of the bands and signals
*                 ! SD corrections were applied
      end if
 
***** Thats all
      return
      end
 
 
CTITLE SD_GET_TIDES
 
      subroutine sd_get_tides( site_num, fund_arg, tides_val )
 
*     This routine scans the list of tide corrections to be
*     applied and computes the contributions at site 'site_num'
*     and saves the values in tides_val.  (The diurnal tides
*     are in the first six elements and the semidiurnals in the
*     last six elements.
*     Here we call sd_arg with only 5 arguments because GST is
*     already built into the tidal partials.
 
      include '../includes/const_param.h'
      include '../includes/kalman_param.h'
      include '../includes/sd_common.h'
 
* PASSED VARIABLES
 
*   site_num        - Number of current site
 
      integer*4 site_num
 
*   fund_arg(6)     - Values of the fundamental arguments
*   tides_val(12)   - Values of the diurnal and semidiurnal
*                   - tide Spherical harmonic coefficients
*                   - for radial, north and East (all m)
 
      real*8 fund_arg(6), tides_val(12)
 
* LOCAL VARIABLES
 
*   i,j     - Loop counters
*   offset  - Change in index between diurnal and semidiurnal
 
      integer*4 i,j, offset
 
*   arg     - Argument of the tide (computed ignoring GST)
 
      real*8 arg
 
*   reported(max_sd_tides)  - Set true when an error in the
*           - GST argumented has been reported.
 
      logical reported(max_sd_tides)
 
      data reported / max_sd_tides*.false. /

* MOD TAH 910130: Quick fix to not having enough tidal sites

      if( site_num.gt.max_sites ) RETURN
 
 
***** Clear all of the values first
 
      do i = 1,12
          tides_val(i) = 0.0d0
      end do
 
***** Scan the list of tides and evaluate those which apply to this
*     site.
 
      do i = 1, sd_tides_num
          if( sd_tides_site(i).eq. site_num .or.
     .        sd_tides_site(i).eq. 99           ) then
              call sd_arg( sd_tides_arg(1,i), fund_arg, arg, 5)
 
*             Now sum up the values
*                                                     ! Diurnal
              if( sd_tides_arg(6,i).eq.-1 ) then
                  offset = 0
*                                                     ! Semidiurnal
              else if( sd_tides_arg(6,i).eq.-2 ) then
                  offset = 6
*                                                 ! Invalid
              else
*                                                 ! argument
                  if( .not.reported(i) ) then
                      write(*,100) i, (sd_tides_arg(j,i),j=1,6)
 100                  format('SD_GET_TIDES ERROR: Invalid GST ',
     .                    ' argument for tide # ',i3,'.',/,
     .                20x,'Arguments are ',6I3,/,
     .                20x,'Remaining tide corrections ignored')
                      reported(i) = .true.
                  end if
 
                  RETURN
              end if
 
*                             ! Loop in pairs for cos and sin
* MOD TAH 930401 Added the sin term into the expression for the
*             tide values.
              do j = 1, 6, 2
                  tides_val(j+offset) = tides_val(j+offset) +
     .                sd_tides_val(j,i)   * cos(arg) -
     .                sd_tides_val(j+1,i) * sin(arg)
                  tides_val(j+offset+1) = tides_val(j+offset+1) +
     .                sd_tides_val(j,i)   * sin(arg) +
     .                sd_tides_val(j+1,i) * cos(arg)

              end do
*                         ! Site matches
          end if
*                         ! Looping over all tides input.
      end do
 
****  Thats all
      return
      end
 
CTITLE EVAL_SD
 
      subroutine eval_sd( epoch, pmu )

*     This routine will compute the instantaneous values of x and y pole
*     and ut1 at a given time.   
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include '../includes/sd_common.h'

* PASSED VARIABLES

* epoch  - Time at which value needed. (jd + fraction of day)
* pmu(3) - Pole and ut1 (all mas)
*

      real*8 epoch, pmu(3)

* LOCAL VARIABLES
 
*   i,j,k       - Loop counters
 
      integer*4 i,j
 
*   arg         - Angulatr argument for the correction (rads)
*   fund_arg(6) - Values of the 5 Brown's arguments (l,l',F,D,
*               - Omega) and gst+pi (rads)
*   xy_part(2)  - Polar motion partial derivatives
*   dt2000      - Time since dj2000 
 
      real*8 arg, fund_arg(6) 

***** Get the fundamental arguments at this epoch
 
      call fund_angles( epoch, fund_arg)
 
*     Clear the corrections
*                     ! Delay and rate values
      do j = 1,3
          pmu(j) = 0.d0
      end do
 
*     Now loop over the UT1 contributions
      do i = 1, sd_ut1_num
 
*         Get the argument
          call sd_arg( sd_ut1_arg(1,i), fund_arg, arg, 6)
 
*         Compute total UT1 change
          pmu(3) = pmu(3) + sd_ut1_val(1,i)*cos(arg) +
     .               sd_ut1_val(2,i)*sin(arg) 
      end do
 
****  Now do polar motion
      do i = 1, sd_xy_num
 
*         Get the argument
          call sd_arg( sd_xy_arg(1,i), fund_arg, arg, 6)
 
*         Compute delay and rate contriubutioin
*                         ! Loop on delay and rate

          pmu(1) = pmu(1) - sd_xy_val(1,i)* cos(arg) +
     .                      sd_xy_val(2,i)* sin(arg)
          pmu(2) = pmu(2) + sd_xy_val(1,i)* sin(arg) +
     .                      sd_xy_val(2,i)* cos(arg)
 
      end do

***** For the moment that is all.
      return
      end


