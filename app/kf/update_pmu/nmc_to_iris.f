      program nmc_to_iris
 
*     Program to read the new format nmc AAM data into an
*     iris style format.

      include '../includes/const_param.h'

*     The runstring for the program is:
*     % nmc_to_iris [input file] [output file] w u
*     where w is add the wind contribution,
*           u is use the UK met office format (only 6 values
*                per line instead of 8)
*
 
*   ierr, jerr  - IOSTAT errors
*   rcpar   - Reads the runstring
*   len_run - Length of runstring
*   trimlen - Length of string
*   i,j     - Loop counters
*   date(4) - Date of determination (assummed 0 hrs)
*   date_full(5) - Full date 
 
      integer*4 ierr, jerr, rcpar, len_run, trimlen, date(4), j,
     .          date_full(5), i

*   chi(3)  - Three Chis (no inverted barometer)
*   chii(3) - Three Chis (with inverted barometer correction)
*   cwn, cws, cpn, cps, cpin, cpis - Components of the excitation.
*   sectag  - Seconds tag for mjd conversion
*   cn100, cs100 - Values at 100 mbar (not used)
*   xp, yp  - Estimates of the X-pole and Y-pole obtained from the
*             excitation.
*   xpi, ypi  - Estimates of the X-pole and Y-pole obtained from the
*             excitation. iWIth inverted barometer.
*   ut1     - Estimate of UT1 from the excitation.
*   ut1i     - Estimate of UT1 from the excitation.
*   dxp, dyp, dut1 - Changes in values between data points
*   dt      - Change in time between data point (days)
*   cwf     - Chandler wobble frequency (rads/day)
*   t1, t2  - Epochs of measuremenst
*   arg     - argument of exp(-icwf*(t2-t1))
 
 
      real*8 chi(3), chii(3), cwn, cws, cpn, cps, cpin, cpis,
     .       cn100, cs100, xp, yp, xpi, ypi, ut1, ut1i, 
     .       dxp, dyp, dut1, dt, cwf, t1, t2, sectag, arg
 
*   input_file  - Read from runstring
*   output_file - Read from runstring
*   pmu_file    - Output integrated polar motion/UT1 file
*   header      - Header line from input.
*   wind        - if = w then wind added
*   ukmo        - if = u then UKMO data file
 
      character*128 input_file, output_file, pmu_file, header
      character*4   wind, ukmo
 
*   line        - Line read from input.
 
      character*256 line

      cwf = (1/433.d0)*2*pi
 
****  Decode the runstring
      len_run = rcpar(1, input_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'nmc_to_iris.hlp', 'nmc_to_iris', 1)
      end if
      open(100, file=input_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',input_file,1,
     .                'nmc_to_iris')
 
****  Get the output and open
      len_run = rcpar(2,output_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'nmc_to_iris.hlp', 'mit_to_iris', 1)
      end if
      open(200, file=output_file, iostat=ierr, status='new')
      call report_error('IOSTAT',ierr,'creat',output_file,1,
     .                'nmc_to_iris')

****  Get the output pmu file and open
      len_run = rcpar(3,pmu_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'nmc_to_iris.hlp', 'mit_to_iris', 1)
      end if
      open(201, file=pmu_file, iostat=ierr, status='new')
      call report_error('IOSTAT',ierr,'creat',pmu_file,1,
     .                'nmc_to_iris')

****  See if wind to be added
      len_run = rcpar(4, wind)
      call casefold(wind)
      len_run = rcpar(5, ukmo)
      call casefold( ukmo )
 
***** Get the header of the GSFC file and write first records of iris
*     file.
 
      read(100,'(a)', iostat=ierr) header
      call report_error('IOSTAT',ierr,'read',input_file,1,
     .                'nmc_to_iris')
      write(200,100) header(1:max(1,trimlen(header))), wind, ukmo
  100 format('* AAM Excitations from: ',a,/,
     .      '* Wind option ',a4,' UKMO option ',a4,/,
     .      '*  Yr Mo Da  H  M    Chi(1)   Chi(1)I  Chi(2)  Chi(2)I',
     .      '  Chi(3)    Chi(3)I',/,
     .      '*                      arcsec            arcsec    ',
     .      '    arcsec       ')
      write(201,110) header(1:max(1,trimlen(header))), wind, ukmo
  110 format('* AAM Pole Position and UT1 from: ',a,/,
     .      '* Wind option ',a4,' UKMO option ',a4,/,
     .      '*  Yr Mo Da  H  M    XP (1)   XP (1)I  YP (2)  YP  2)I',
     .      '  UT1(3)    UT1(3)I',/,
     .      '*                      arcsec            arcsec    ',
     .      '    arcsec       ')
 
****  Now start the conversion
      t1   = 0.d0
C Old values, replaced for Jan 1983
C     xp   = -.13051d0
C     yp   = .16516d0
* Y_pole has 0.326 mean removed.
      xp   = -.186d0
      yp   = -.113d0
      ut1  = 0.d0
      xpi  = -.186d0
      ypi  = -.113d0
      ut1i = 0.d0
      ierr = 0
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr ) line
*                                     ! Decode
          if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0                 ) then
 
*             Read the date and time
              read(line,150, iostat=jerr) date
 150          format(1x,4i2)
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'nmc_to_iris')
              if( date(1).lt.500) date(1) = date(1) + 1900

*             Now read the three Chi values
              do j = 1,3
                 if( ukmo(1:1).eq.'U' ) then
                     read(100,160)  cwn, cpn,cpin,
     .                              cws, cps, cpis
  160                format(8f10.5)
                 else
                     read(100,160)  cn100, cwn, cpn,cpin,
     .                              cs100, cws, cps, cpis
                  end if
                  if( wind(1:1).ne.'W' ) then
                      cwn = 0
                      cws = 0
                  end if
                  if( j.le.2 ) then
*                    Polar motion excitation (asecs)
                     chi(j) = (cwn + cpn + cws + cps)*0.020626d0
                     chii(j) = (cwn + cpin + cws + cpis)*0.020626d0
                  else
                     chi(j) = - (cwn + cpn + cws + cps)*86400.d-7 
     .                         +0.121941d0
                     chii(j) = -(cwn + cpin + cws + cpis)*86400.d-7 
     .                         +0.121941d0
                  end if
              end do 
 
*             Write out file
              if( cpn.ne.-99.0 .and. cps.ne.-99.0 ) then
                 write(200,200, iostat=ierr) date, (chi(j),
     .                           chii(j), j=1,3)
  200            format(1x,i4,3(i3),' 0', 1x,6(f9.6,1x))
                 call report_error('IOSTAT',ierr,'writ',line,1,
     .                'nmc_to_iris')

                 do i = 1,4
                    date_full(i) = date(i)
                 end do
                 date_full(5) = 0
                 sectag = 0.d0

*                Now start integrating the excitations.
                 if( t1.eq.0.d0 ) then
                     call ymdhms_to_jd( date_full, sectag, t1) 
                 else
                     call ymdhms_to_jd( date_full, sectag, t2)
                     dt = t2 - t1
                  
*                    Do the non-inverted baramoter first
                     arg = cwf*(t2-t1)
                     dxp = xp*cos(arg) - yp*sin(arg) + chi(2)*cwf*dt
                     dyp = yp*cos(arg) + xp*sin(arg) - chi(1)*cwf*dt
                     ut1 = ut1 + chi(3)*dt
                     xp = dxp
                     yp = dyp
   
*                    Now do the inverted barometer term.
                     dxp = xpi*cos(arg) - ypi*sin(arg) + chii(2)*cwf*dt
                     dyp = ypi*cos(arg) + xpi*sin(arg) - chii(1)*cwf*dt
                     ut1i = ut1i + chii(3)*dt
                     xpi = dxp
                     ypi = dyp
                     write(201,210, iostat=ierr) date, xp, xpi,
     .                              yp, ypi, ut1, ut1i
  210                format(1x,i4,3(i3),' 0', 1x,6(f10.6,1x))
                     t1 = t2
                end if
              end if

*                     ! No error on read
          end if
*                     ! Looping until end of input.
      end do
 
***** Thats all
      close(100)
      close(200)
      end
 
