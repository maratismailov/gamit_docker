      program mit_to_iers

      implicit none 
 
*     This program will convert a pmu file obtained from a
*     globk solution using the iers.ext extract command file
*     into the standard IERS format.
*
*     The iers.ext extract file looks like:
*   *    Global solution version of PMU extract. (For GLOBK output)
*   |  field 1 "Solution refers to     :"  5 I4 format   0 "(1x,i4,4(1x,i2))"
*   |  field 2 "X-pole position        (mas)"        2 R8 readline 0 1 3
*   |  field 3 "Y-pole position        (mas)"        2 R8 readline 0 1 3
*   |  field 4 "UT1-AT                 (mts)"        2 R8 readline 0 1 3
*   |  field 5 "prefit chi**2 for"                   1 R8 readline 0 5
*   |  field 6 "Pole/UT1 correlations: XY, XU, YU"   3 R8 readline 0 1 2 3
*   |  field 7 "radio sources, and"                  2 I4 readline 1 3 9
 
*   |  outform 2 "(2(f8.2,1x))"
*   |  outform 3 "(2(f8.2,1x))"
*   |  outform 4 "(2(f12.3,1x))"
*   |  outform 5 "(2x,f6.2)"
*   |  outform 6 "(1x,3(f6.3,1x))"
*   |  outform 7 "(1x,2i4)"
 
*   |  run
*
 
* VARIBLE DECLARATIONS.
 
*   ierr    - IOSTAT error
*   jerr    - IOSTAT error on decoding input line
*   rcpar   - read runstring
*   trimlen - Length of string
*   iout    - Output unit number (200 if file, 6 if stdout)
*   len_run - Length of runstring
*   date(5) - Date array
*   nstat   - Number of stations in solution
*   nsvs    - NUmber of satellites
 
      integer*4 ierr, jerr, rcpar, trimlen, iout, len_run, date(5),
     .    nstat, nsvs
 
*   sec_tag - Seconds tag
*   jd_pmu  - Julian date of pmu measurement
*   mjd_pmu - Modified julians date of pmu measurement
*   x,sx,y,sy   - Pole postion and sigma
*   utmat, sut, utmutc  - UT1-AT, sigma and UT1-UTC
*   dut     - Short period UT1 correction (sec)
*   chi, rms    - Prefit chi**2 of solution and sqrt.
*   rxy, rxu, ryu   - Correlations
*   jd_leap_792 - Julians date of July 1, 1992 leap second.
 
      real*8 sec_tag, jd_pmu, mjd_pmu, x,sx,y,sy, utmat, sut, utmutc,
     .    dut, chi, rms, rxy, rxu, ryu, jd_leap_792
 
*   line    - Line read from input (ignored in col 1 non-blank)
*   in_file - Input file
*   out_file    - Output file (overwritten by this program)
 
      character*256 line, in_file, out_file
 
*     Get the runstring print help if there is problem
      len_run = rcpar(1,in_file)
      if( len_run.le.0 ) call proper_runstring('mit_to_iers.hlp',
     .                        'mit_to_iers',1)
 
*     Open the input
      open(100, file=in_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',in_file,1,'mit_to_iers')
 
****  Get the output.  If not given write to the screen
      len_run = rcpar(2,out_file)
      if( len_run.eq.0 ) out_file = '6'
      call open_lu(iout, out_file, ierr, 'unknown')
      call report_error('IOSTAT',ierr,'open',out_file,1,'mit_to_iers')
 
***** Get the julian date of the 1992 July 1 leap second
      date(1) = 1992
      date(2) = 7
      date(3) = 1
      date(4) = 0
      date(5) = 0
      sec_tag = 0.d0
      call ymdhms_to_jd( date, sec_tag, jd_leap_792)
 
****  Write out some header lines
      call systime(date, sec_tag)
      write(iout,100) date, in_file(1:trimlen(in_file))
 100  format('* MIT_TO_IERS: Run on ',i4,'/',i2,'/',i2,1x,i2,':',i2,/,
     .       '* Input file : ',a)
 
*     Now loop over input
      do while ( ierr.eq.0 )
          read(100,'(a)',iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0 ) then
 
*             Process the line
              read(line,*,iostat=jerr) date, x,sx,y,sy, utmat, sut,
     .                    chi, rxy, rxu, ryu, nstat, nsvs
 
*             report error if any, and process if no error:
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                    'mit_to_iers')
 
*             Start the conversions:
*             Get MJD
              sec_tag = 0.0
              call ymdhms_to_jd( date,sec_tag, jd_pmu)
              mjd_pmu = jd_pmu - 2400000.5d0
 
*             Convert pole position and UT1
              x  =  x /1000.d0
              sx = sx /1000.d0
              y  =  y /1000.d0
              sy = sy /1000.d0
              utmat = utmat /1000.d0
              sut = sut / 1000.d0
 
*****         Convert UT1-AT to UT1-UTC (Add back the short period
*             UT1 terms)
              call short_period_ut1( jd_pmu, dut )
              utmat = utmat + dut
 
*             Now convert to UT1-UTC.
*             WARNING: This is not general code, valid only for
*             1991-1993.
              utmutc = utmat + 26.d0
              if( jd_pmu.gt. jd_leap_792) utmutc = utmutc + 1.d0
 
*             Covert the chi**2 of solution to rms
              rms = sqrt(chi)
 
*             Now write out the output line
              write(iout, 200) mjd_pmu, x,y,utmutc, sx, sy, sut,
     .                        rms, rxy, rxu, ryu, nstat, nsvs
 200          format(1x,f7.1,1x,2F9.5,1x,F9.6,1x,2f8.5,1x,f8.6,
     .               f5.2,3f7.3,1x,2i4,' 0')
          end if
      end do
 
****  Thats all
      close(100)
      end
 
 
 
