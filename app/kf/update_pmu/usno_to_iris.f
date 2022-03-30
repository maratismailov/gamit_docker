
      program usno_to_iris
 
*     This program will read a IERS eop format file and convert
*     it to NGS format iris format for use in update_pmu.
*     The runstring for the program is:
*     % usno_to_iris [input file] [output file] [units]
*
 
*   ierr, jerr  - IOSTAT errors
*   rcpar   - Reads the runstring
*   len_run - Length of runstring
*   trimlen - Length of string
*   i,j     - Loop counters
*   date(5) - Date of determination (assummed 0 hrs)
 
      integer*4 ierr, jerr, rcpar, len_run, trimlen, date(5)
 
*   x,y     - X and Y pole position
*   sx,sy   - Sigmas on x and y
*   mjd     - Modified julian date of dertermination
*   ut1     - ut1 - utc (seconds)
*   su      - Sigma for Ut1
*   psi, eps    - Nutations in psi and epsilon
*   sp, se  - Sigmas for psi and eps
*   sectag  - Seconds tag for mjd conversion
*   units   - Scale factor for the program.
*   mjd     - Mjd of measurements
 
 
      real*8 x,y, sx,sy, ut1, su, units, sectag, mjd, dpsi, deps
 
*   input_file  - Read from runstring
*   output_file - Read from runstring
*   header      - Header line from input.
 
      character*128 input_file, output_file, header
 
*   line        - Line read from input.
 
      character*256 line
 
****  Decode the runstring
      len_run = rcpar(1, input_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'usno_to_iris.hlp', 'usno_to_iris', 1)
      end if
      open(100, file=input_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',input_file,1,
     .                'usno_to_iris')
 
****  Get the output and open
      len_run = rcpar(2,output_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'usno_to_iris', 'usno_to_iris', 1)
      end if

*     Get optional units:
      len_run = rcpar(3, header)
      if( len_run.gt.0 ) then
          read(header,*, iostat=ierr) units
          call report_error('IOSTAT',ierr,'decod', header,1,
     .                      'Getting units')
      else
          units = 1.0d0
      end if
      open(200, file=output_file, iostat=ierr, status='new')
      call report_error('IOSTAT',ierr,'creat',output_file,1,
     .                'usno_to_iris')
 
***** Get the header of the GSFC file and write first records of iris
*     file.
 
      read(100,'(a)', iostat=ierr) header
      call report_error('IOSTAT',ierr,'read',input_file,1,
     .                'usno_to_iris')
      write(200,100) header(1:max(1,trimlen(header)))
  100 format('*EOR Determinations: ',a,/,
     .      '*  Yr Mo Da  H  M      X      SE        Y      SE  ',
     .      '  UT1-UTC    SE  ',/,
     .      '*                      arcsec            arcsec    ',
     .      '    Seconds      ')
 
****  Now start the conversion
      ierr = 0
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr ) line
*                                     ! Decode
          if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0 ) then
 
*             Read from position 7 on (skip date at front, use mjd for
*             epoch)
              read(line,*, iostat=jerr) mjd, x, y, ut1, dpsi, deps,
     .                                  sx, sy, su
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'usno_to_iris')
 
*             convert the sigmas
              x   = x   / units  
              y   = y   / units  
              ut1 = ut1 / units   
              sx  = sx  / units   
              sy  = sy  / units   
              su  = su  / units   
              call jd_to_ymdhms(mjd, date, sectag)
 
*             Write out file
              if( ierr.eq.0 ) 
     .        write(200,200, iostat=ierr) date, x, sx, y, sy, ut1, su
  200         format(1x,i4,4(i3),1x,4(f10.6),1x,2(f11.6,1x))
              call report_error('IOSTAT',ierr,'writ',line,1,
     .                'usno_to_iris')
*                     ! No error on read
          end if
*                     ! Looping until end of input.
      end do
 
***** Thats all
      close(100)
      close(200)
      end
 
