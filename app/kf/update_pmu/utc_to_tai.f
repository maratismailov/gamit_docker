      program utc_to_tai

*     program to convert a standard globk pmu file with UT1-UTC
*     entries into UT1-TAI entries.

*     Program uses the GAMIT leap.sec file which it first tries
*     to open directly.  If this fails it tries to construct the
*     full path from the HELP_DIR name.

*     Runstring is:
*     % utc_to_tai <in file> <out file>

      integer*4 max_nl     ! Maximum number of leap seconds
      parameter ( max_nl = 1000 )

      integer*4  nl  	! NUmber of leap seconds
     .,          rcpar  ! Reads runstring
     .,          len_run  ! Length of runstring
     .,          trimlen  ! Length of string
     .,          ierr, jerr   ! IOSTAT error flags
     .,          i        ! Loop counter
     .,          date(5)  ! yr, mon, day, hr, min
     .,          indx     ! Position in string

      real*8 xjd, jd      ! PEP JD and standard JD
     .,      val(6)       ! Values of xp, +-, yp +-, ut1-utc, +-
                          ! in arcsec and time sec
     .,      sectag       ! seconds tag
     .,      leap_jd(max_nl)   ! Julian dates of leap seconds
     .,      tai_utc(max_nl)   ! Value of TAI-UTC after the date of
                               ! each leap second.

      character*4 cdum        ! Dummy string for multiread

      character*256 infile, outfile  ! In and output file names
     .,             help_dir         ! name of help directory
     .,             leap_file        ! name of leap second file
     .,             line             ! Line read from file.


***** Read the runstring to get the file names
      write(*,120)
 120  format(' UTC_TO_TAI conversion program')
      len_run = rcpar(1, infile)
      if( len_run.le.0 ) then
          call proper_runstring('utc_to_tai.hlp','UTC_TO_TAI',1)
      end if
      open(100,file=infile, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',infile,1,'utc_to_tai')

      len_run = rcpar(2, outfile)
      if( len_run.le.0 ) then
          call proper_runstring('utc_to_tai.hlp','UTC_TO_TAI',1)
      end if
      open(200,file=outfile, status='unknown', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',outfile,1,'utc_to_tai')

****  Now try to open leap.sec file
      open(101,file='leap.sec', status='old', iostat=ierr)
      if( ierr.ne.0 ) then

*         Could not open file.  Try to generate name
          call getenv('HELP_DIR', help_dir)
          leap_file = help_dir(1:trimlen(help_dir)) //
     .                '../gamit/tables/leap.sec'
          open(101,file=leap_file, status='old',iostat=ierr)
          if( ierr.ne.0 ) then
              write(*,140) leap_file(1:trimlen(leap_file)), ierr
 140          format('**ERROR** Could not open leap.sec, also tried ',
     .               a,' IOSTAT error ',i5,/,
     .               ' Cannot convert UTC to TAI')
              stop 'UTC_TO_TAI: Could not open leap.sec file'
          end if
      end if

****  Now read the leap.sec file and make tai minus utc
      nl = 1
      leap_jd(1) = 2444786.50d0
      tai_utc(1) = -20.d0

      read(101,'(a)' )  line
      read(101,'(a)' )  line
      ierr = 0 
      do while ( ierr.eq.0 ) 
         read(101,'(a)', iostat=ierr) line
         if( ierr.eq.0 .and. trimlen(line).gt.0 ) then
             read(line,*,iostat=jerr) xjd
             if( jerr.eq.0 ) then
                 nl = nl + 1
                 leap_jd(nl) = xjd + 0.5
                 tai_utc(nl) = tai_utc(nl-1) - 1.d0
             end if
         end if
      end do

***** Now loop over the input flle
      ierr = 0
      write(200,200) infile(1:trimlen(infile))
 200  format('* File:',a,', converted from UT1-UTC to UT1-TAI by',
     .       ' program utc_to_tai')
      do while ( ierr.eq.0 ) 
         read(100,'(a)', iostat=ierr) line
         if(  ierr.eq.0 ) then
             if( line(1:1).eq. ' ' .and. trimlen(line).gt.0 ) then
C                read(line,*)  date,val
                 indx = 1
                 call multiread(line, indx,'I4', jerr,date, cdum, 5)
                 call multiread(line, indx,'R8', jerr,val , cdum, 6)

*                convert date to jd
                 sectag = 0.d0
                 call ymdhms_to_jd( date, sectag, jd )

*                Find leap seconds
                 do i = 1, nl-1
                    if( jd.ge.leap_jd(i) .and.
     .                  jd.lt.leap_jd(i+1) ) then
                        if( abs(val(5)).lt.1.d0 ) then
                            val(5) = val(5) + tai_utc(i)
                        end if
                    end if
                 end do
                 if( abs(val(5)).lt.1.d0 ) then
                     val(5) = val(5) + tai_utc(nl)
                 end if

                 write(200,210) date, val,
     .                    line(indx:max(indx,trimlen(line)))
 210             format(i5,4i3,4F9.6,1x,2F12.7,1x,a)
             else
                 write(200,'(a)') line(1:trimlen(line))
             end if
          end if
      end do

****  Thats all
      close(100)
      close(101)
      close(200)
      end 
                    



