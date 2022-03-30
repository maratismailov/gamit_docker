
      program mit_to_iris
 
*     This program will read a MIT  eop format file and convert
*     it to NGS format iris format for use in update_pmu.
*     The runstring for the program is:
*     % mit_to_iris [input file] [output file]
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
 
 
      real*8 x,y, sx,sy, ut1, su
 
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
          call proper_runstring( 'mit_to_iris.hlp', 'mit_to_iris', 1)
      end if
      open(100, file=input_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',input_file,1,
     .                'mit_to_iris')
 
****  Get the output and open
      len_run = rcpar(2,output_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'mit_to_iris', 'mit_to_iris', 1)
      end if
      open(200, file=output_file, iostat=ierr, status='new')
      call report_error('IOSTAT',ierr,'creat',output_file,1,
     .                'mit_to_iris')
 
***** Get the header of the GSFC file and write first records of iris
*     file.
 
      read(100,'(a)', iostat=ierr) header
      call report_error('IOSTAT',ierr,'read',input_file,1,
     .                'mit_to_iris')
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
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
 
*             Read from position 7 on (skip date at front, use mjd for
*             epoch)
              read(line,*, iostat=jerr) date, x, sx, y, sy, ut1, su
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'mit_to_iris')
 
*             convert the sigmas
              x   = x   / 1000.d0
              y   = y   / 1000.d0
              ut1 = ut1 / 1000.d0
              sx  = sx  / 1000.d0
              sy  = sy  / 1000.d0
              su  = su  / 1000.d0
 
*             Write out file
              write(200,200, iostat=ierr) date, x, sx, y, sy, ut1, su
  200         format(1x,i4,4(i3),1x,4(f10.6),1x,2(f11.6,1x))
              call report_error('IOSTAT',ierr,'writ',line,1,
     .                'mit_to_iris')
*                     ! No error on read
          end if
*                     ! Looping until end of input.
      end do
 
***** Thats all
      close(100)
      close(200)
      end
 
