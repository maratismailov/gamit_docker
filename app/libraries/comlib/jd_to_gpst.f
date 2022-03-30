CTITLE JD_TO_GPST

      subroutine jd_to_gpst( jd, gpsw, gpsd, gpss)

      implicit none

*     Routine to convert JD or MJD to GPS Week, GPS doy of week
*     and GPS seconds of week.

* PASSED variables
      real*8 jd   ! Either julian date or modified julian date
      real*8 gpss ! GPS seconds of week (returned)

      integer*4 gpsw  ! GPS Week
      integer*4 gpsd  ! GPS day of week.

* LOCAL
      real*8 MJD_gpsst ! Week 0 MJD
      real*8 mjd       ! Input JD converted to MJD

      data MJD_gpsst / 44244.0d0 /

***   See if have JD or MJD
      if( jd.gt.2000000.d0 ) then
          mjd = jd - 2 400 000.5d0
      else
          mjd = jd
      end if

      gpsw = (mjd - MJD_gpsst)/7.d0
      gpsd =  mjd - gpsw*7 - MJD_gpsst
      gpss = (mjd - gpsw*7 - MJD_gpsst)*86400.d0

      return
      end
