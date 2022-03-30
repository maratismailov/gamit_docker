      program tgrep 

      implicit none 
*
*     Program to extract lines from a globk output file 
*     in a way similar to grep but with the time added to 
*     beginnings of the lines.
*
      character*100 infile   ! Name of input file
      character*256 line     ! Line read from file
      character*10 stime     ! Date in form of string
      character*128 search   ! String to be search for (^ for
*                            ! forced blank)
      integer*4 ierr         ! IOSTAT error
      integer*4 trimlen      ! Length of string
      integer*4 rcpar        ! Reads runstring
      integer*4 lenrun       ! Length of runstring parameter
      integer*4 unit         ! Input unit number (50 if file given,
*                            ! 5 if stdin)
      integer*4 lensearch    ! Length of search string
      integer*4 nr           ! Counter

*     Get the string passed by used
      lenrun = rcpar(1,search)
      if( lenrun.le.0 ) then
          call proper_runstring('tgrep.hlp','tgrep',1)
      endif
      lensearch = lenrun
*     Replace any ^ symbols with white space
      call sub_char(search,'^',' ')

*     Get the file name
      nr = 1
      do while ( lenrun.gt.0 ) 
         nr = nr + 1
         lenrun = rcpar(nr,infile)
         if( lenrun.gt.0 ) then
            unit = 50
            open(unit,file=infile,status='old',iostat=ierr)
            call report_error('IOSTAT',ierr,'open',infile,
     .                       1,'tgrep')
*           Set initial time to TIMESWRONG
            stime = 'TIMESWRONG'

            do while ( ierr.eq.0 )
               read(unit,'(a)',iostat=ierr) line
               if( index(line,'Solution refers').gt.0 ) then
                   stime = line(48:56)
               elseif( index(line,'EXPERIMENT date').gt.0 ) then
cEXPERIMENT date : 2007/10/ 7 11:59    (2007.7658) [Seconds tag  45.000]
                   stime = line(41:49)
               else if( index(line,
     .                search(1:lensearch)).gt.0 ) then
                   write(*,220) stime, line(1:trimlen(line))
 220               format(a10,1x,a)
               endif
            end do
         end if
      enddo
      end
         
* Solution refers to     : 1999/ 1/ 7 11:59    (1999.0178) [Seconds tag  45.000]
* EXPERIMENT date : 2007/10/ 7 11:59    (2007.7658) [Seconds tag  45.000]
*         1         2         3         4         5
*123456789012345678901234567890123456789012345678901234567890

