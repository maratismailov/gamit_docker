      program gsfc_to_iris
 
*     This program will read a GSFC eop format file and convert
*     it to NGS format iris format for use in update_pmu.
*     The runstring for the program is:
*     % gsfc_to_iris [input file] [output file]
*
 
*   ierr, jerr  - IOSTAT errors
*   rcpar   - Reads the runstring
*   len_run - Length of runstring
*   trimlen - Length of string
*   i,j     - Loop counters
*   date(5) - Date of determination (assummed 0 hrs)
 
      integer*4 ierr, jerr, rcpar, len_run, trimlen, i, date(5)
 
*   x,y     - X and Y pole position
*   sx,sy   - Sigmas on x and y
*   mjd     - Modified julian date of dertermination
*   ut1     - ut1 - utc (seconds)
*   su      - Sigma for Ut1
*   psi, eps    - Nutations in psi and epsilon
*   sp, se  - Sigmas for psi and eps
*   sectag  - Seconds tag for mjd conversion
 
 
      real*8 x,y, sx,sy, mjd, ut1, su, psi, eps, sp, se, sectag
 
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
          call proper_runstring( 'gsfc_to_iris', 'gsfc_to_iris', 1)
      end if
      open(100, file=input_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',input_file,1,
     .                'gsfc_to_iris')
 
****  Get the output and open
      len_run = rcpar(2,output_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'gsfc_to_iris', 'gsfc_to_iris', 1)
      end if
      open(200, file=output_file, iostat=ierr, status='new')
      call report_error('IOSTAT',ierr,'creat',output_file,1,
     .                'gsfc_to_iris')
 
***** Get the header of the GSFC file and write first records of iris
*     file.
 
      read(100,'(a)', iostat=ierr) header
      call report_error('IOSTAT',ierr,'read',input_file,1,
     .                'gsfc_to_iris')
      write(200,100) header(1:max(1,trimlen(header)))
  100 format(' EOR Determinations: ',a,/,
     .      '   Yr Mo Da  H  M      X      SE        Y      SE  ',
     .      '  UT1-UTC    SE    NUT PSI  SE   NUT EPS SE',/,
     .      '                       arcsec            arcsec    ',
     .      '    Seconds           mas           mas')
 
****  Skip next 5 records from input
      do i = 1,5
          read(100,'(a)', iostat=ierr) line
          call report_error('IOSTAT',ierr,'read',input_file,1,
     .                    'gsfc_to_iris')
 
      end do
 
****  Now start the conversion
      ierr = 0
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr ) line
*                                     ! Decode
          if( ierr.eq.0 ) then
 
*             Read from position 7 on (skip date at front, use mjd for
*             epoch)
              read(line(7:),*, iostat=jerr) mjd, x, y, ut1, psi, eps,
     .                sx, sy, su, sp, se
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'gsfc_to_iris')
 
*             convert the sigmas
              sx = sx / 1000.d0
              sy = sy / 1000.d0
              su = su / 10000.d0
 
*             Now convert date and write out line
              call jd_to_ymdhms( mjd, date, sectag)
 
*             Write out file
              write(200,200, iostat=ierr) date, x, sx, y, sy, ut1, su,
     .            psi, sp, eps, se
  200         format(1x,i4,4(i3),1x,4(f8.5),1x,2(f9.6),1x,4(f7.2))
              call report_error('IOSTAT',ierr,'writ',line,1,
     .                'gsfc_to_iris')
*                     ! No error on read
          end if
*                     ! Looping until end of input.
      end do
 
***** Thats all
      close(100)
      close(200)
      end
 
