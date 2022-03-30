CTITLE PREREAD_CF
 
      subroutine preread_cf

      implicit none 
 
*     This routine will preread the cfiles and compute the clocks
*              good data and bias flag.  May add low elevation
*     after each epoch of data has been read.  These clocks are saved
*     and will be used when the cfile is updated.
 
      include '../includes/kalman_param.h'
      include '../includes/cfile_def.h'
 
      include 'modear.h'
 
* LOCAL variables
 
*   ierr    - IOSTAT error
*   i,j,k   - Loop counters
*   pn      - PRN number of current satellite
 
 
      integer*4 ierr, i,j,k, pn
 
*   mjd  - Julian date
*   lat, lng, rad   - Lat, long and radius of site
 
 
 
      real*8 mjd, lat, lng, rad
 
****  Open the cfile
      call open_cf( 100, in_cf, ierr )
      call report_error('IOSTAT',ierr,'open', in_cf,1,'preread_cf')
 
**** Now read the header records of the cfile
      call read_cf1(100, 'ALL', ierr)
      call report_error('IOSTAT',ierr,'cf1 read', in_cf,
     .                 1,'preread_cf')
 
      call read_cf2(100, 'ALL', ierr)
      call report_error('IOSTAT',ierr,'cf2 read', in_cf,
     .                 1,'preread_cf')
 
      call read_cf3(100, 'ALL', ierr)
      call report_error('IOSTAT',ierr,'cf3 read', in_cf,
     .                 1,'preread_cf')
 
**** Now save off some the information that we need.
      lat = cf_preval(1)
      lng = cf_preval(2)
      rad = cf_preval(3)*1.d3
 
*       Convert to XYZ
      site_xyz(1) = rad*cos(lat)*cos(lng)
      site_xyz(2) = rad*cos(lat)*sin(lng)
      site_xyz(3) = rad*sin(lat)
 
****  Now loop over the cfile
      do i = 1, cf_nepoch
          call read_cf4(100,'ALL', ierr )
          call report_error('IOSTAT',ierr,'cf4 read', in_cf,
     .                1,'preread_cf')
 
****      Now loop over the data at this epoch
          do j = 1, num_sat
              P1o(j) = 0.d0
              P2o(j) = 0.d0
              omc_OK(j) = .false.
          end do
 
****      Get the nominal time for this epoch
          call ydsd_to_mjd( cf_iyr, cf_idoy, cf_sod, mjd )
          data_epoch = mjd
 
          do j = 1, cf_msat
              call read_cf5(100,'ALL', ierr )
 
*             Save the observed pseudorange values (Save for
*             good data and bias flag.  May add low elevation
*             later).
              if( cf_ierfl.eq.0 .or. cf_ierfl.eq.10 ) then
                  pn = cf_iprn
                  omc_OK(pn) = .true.
                  do k = 1, cf_ndat
                      if( cf_dattyp(k).eq.3 .or.
     .                    cf_dattyp(k).eq.5     ) then
                          P1o(pn) = cf_obsv(k)
                      else if( cf_dattyp(k).eq.4 ) then
                          P2o(pn) = cf_obsv(k)
                      end if
                  end do
              end if
          end do
 
****      Now we have the observed Psuedo ranges at L1 and L2 for
*         all oberved satellites.
 
          call comp_clock( site_clock(i), i )
      end do
 
***** Thats all
      close(100)
      return
      end
 
CTITLE YDSD_TO_MJD
 
      subroutine ydsd_to_mjd( yr, doy, sec, mjd )

      implicit none 
 
*     Routine to compute Julian date from Year, day-of-year and
*     seconds of day when seconds is a double precision value
 
* PASSED VARIABLES
 
*   yr, doy     - Year and day of year
 
      integer*4 yr, doy
 
*   sec     - Seconds of day
*   mjd     - Modified Julian date
 
      real*8 sec, mjd
 
* LOCAL VARIABLES
 
*   date(5)     - Calender date with hrs and min
 
      integer*4 date(5)
 
*   sectag      - Seconds tag for ymdhms_to_mjd
 
      real*8 sectag
 
****  Save values for call to ymdhms_to_mjd
 
      date(1) = yr
      date(2) = 1
      date(3) = doy
      date(4) = 0
      date(5) = 0
      sectag = 0
 
      call ymdhms_to_mjd(date, sectag, mjd )
      mjd = mjd + sec/86400.d0
 
***** Thats all
      return
      end
 
CTITLE 'YMDHMS_TO_MJD'
 
      SUBROUTINE YMDHMS_to_MJD ( date, seconds, epoch )

      implicit none 

*     -------------------------
*    .,      26JUL85                <871210.1317>
 
*     Author: T.Herring            3:53 PM  THU., 10  JULY, 1986
*
*$S "YMDHMS_to_MJD"
*-----------------------------------------------------------------------
*     Routine to convert a calender date with hours, minutes and
*     seconds to a Modified Julian date. The calender date is
*     ordered as year, month, day of month, hours, and
*     minutes.  These values are stored in a single I*2 array.
*
*     If the year is greater than 1000 then the it is assumed to
*     contain the centuries.
*
*     This routine is only valid for dates after 1600 Jan 0.
*
*     CALLING SEQUENCE:
*     =================
*     CALL ymdhms_to_JD ( date, seconds, epoch)
*
*     WHERE:
*     date    is a 5 element array containing year, month, day of month,
*             hour, and minutes.
*             (I*4 5 element array INPUT)
*     seconds is the seconds part of the epoch
*             (REAL*8 INPUT)
*     epoch   is the JD with fractional days corresponding to date and
*             seconds.
*             (REAL*8 OUTPUT)
*
*----------------------------------------------------------------------
*$E
 
 
*         day           - day of month
*         date(1)     - 5 element date array with year, month, day,
*               - hours and minutes.
 
*   days_to_month(12)   - number of days from start of year to
*               - each month
 
*   leap_days   - the number of leap days we need to include
 
*   month       - month of year
*   year        - the years since 1900.0
 
*   years_from_1600 - number of years since 1600 Jan 0.
 
*   days_from_1600  - number of days since 1600 Jan 0.
 
      integer*4 day, date(5), days_to_month(12), leap_days, month,
     .    year, years_from_1600, days_from_1600
 
*       epoch   - The JD corresponding to date and seconds
*   fraction    - Fraction of a day.  We compute this separately
*               - to avoiod some rounding error when added to MJD
*   mjd         - the computes MJD with fractional days
*   seconds     - the seconds part of the epoch.
 
      real*8 epoch, fraction, mjd, seconds
 
*       leap_year   - Indicates that this is a leap year
 
      logical leap_year
 
      data  days_to_month /   0,  31,  59,  90, 120, 151,
     .                      181, 212, 243, 273, 304, 334 /
 
***** START, Make sure year is from 1900
 
      year      = date(1)
      month     = date(2)
      day       = date(3)
 
      if( year.lt.1000 ) year = year + 1900
 
***** Compute number of years from 1600
 
      years_from_1600 = year - 1600
 
***** Now compute number of leap days upto the start of the year
*     we are in (i.e., if the year we are in is a leap year do not
*     add the leap day yet)
 
      leap_days =   (years_from_1600 -   1)/  4
     .            - (years_from_1600 +  99)/100
     .            + (years_from_1600 + 399)/400  + 1
 
      if( years_from_1600.eq.0 ) leap_days = leap_days - 1
 
***** Now see if we are in leap year
 
      leap_year = .false.
      if(   mod(years_from_1600,  4).eq.0     .and.
     .     (mod(years_from_1600,100).ne.0.or.
     .      mod(years_from_1600,400).eq.0)       ) leap_year = .true.
 
***** Now compute number of days sinec 1600
 
      days_from_1600 = years_from_1600*365  + leap_days +
     .                 days_to_month(month) + day - 142350
 
***** Add extra day if we are after Februrary and this is a leap year
      if( month.gt.2 .and. leap_year ) then
          days_from_1600 = days_from_1600 + 1
      end if
 
***** Compute the mjd and add in the fraction of a day part
*     The MJD of 1600 Jan 0 is -94554
 
      fraction  = seconds/86400.d0  + date(5)/1440.d0 + date(4)/24.d0
 
C     mjd   = -94554.d0 + days_from_1600 + fraction
      mjd   =  days_from_1600 + fraction
      epoch = mjd 
 
***** THATS ALL
      RETURN
      end
 
