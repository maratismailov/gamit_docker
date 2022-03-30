
      program mjd_to_date

*     Program to convert either mjd to date at the beginning of
*     of a line or date to mjd.  The direction is set by the
*     input type (mjd or date) 
*     % mjd_to_date [input file] [output file] [input type]
*     Where input type is MJD to convert mjd to date or
*                         DATE to convert date to MJD
*     Non-blank first characters are treated as comments
*
 
*   ierr, jerr  - IOSTAT errors
*   rcpar   - Reads the runstring
*   len_run - Length of runstring
*   trimlen - Length of string
*   i,j     - Loop counters
*   date(5) - Date of determination (assummed 0 hrs)
*   indx    - Position in string
 
      integer*4 ierr, jerr, rcpar, len_run, trimlen, date(5),
     .          indx
 
*   mjd     - Mjd of measurements
*   sectag  - Seconds tag
 
      real*8 sectag, mjd 
 
*   input_file  - Read from runstring
*   output_file - Read from runstring
 
      character*128 input_file, output_file
 
*   line        - Line read from input.
 
      character*256 line

*   intype      - Input date type

      character*8 intype, cdum
 
****  Decode the runstring
      len_run = rcpar(1, input_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'mjd_to_date.hlp', 'mjd_to_date', 1)
      end if
      open(100, file=input_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',input_file,1,
     .                'mjd_to_date')
 
****  Get the output and open
      len_run = rcpar(2,output_file)
*                                 ! Print out the help file, and
      if( len_run.le.0 ) then
*                                 ! stop
          call proper_runstring( 'mjd_to_date', 'mjd_to_date', 1)
      end if

*     Get optional input type (default is mjd)
      len_run = rcpar(3, intype)
      if( len_run.gt.0 ) then
          call casefold( intype )
      else
          intype = 'MJD'
      end if
      open(200, file=output_file, iostat=ierr, status='new')
      call report_error('IOSTAT',ierr,'creat',output_file,1,
     .                'mjd_to_date')
 
****  Now start the conversion
      ierr = 0
      sectag = 0
      do while ( ierr.eq.0 )
          read(100,'(a)', iostat=ierr ) line
*                                     ! Decode
          if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0 ) then

              indx = 1
              if( intype(1:2).ne.'DA' ) then
                  call read_line(line,indx,'R8',jerr, mjd, cdum )  
                  call jd_to_ymdhms(mjd,date,sectag)
                  write(200,210) date, line(indx:trimlen(line))
 210              format(1x,i4,4i3,1x,a)
              else
                  call multiread(line,indx,'I4',jerr, date, cdum, 5) 
                  call ymdhms_to_mjd(date, sectag, mjd)
                  write(200,220) mjd, line(indx:trimlen(line))  
 220              format(1x,f14.5,1x,a)
             end if
          else if( ierr.eq.0 ) then           
             write(200,'(a)') line(1:max(trimlen(line),1))
          end if
 
      end do
 
***** Thats all
      close(100)
      close(200)
      end
 
