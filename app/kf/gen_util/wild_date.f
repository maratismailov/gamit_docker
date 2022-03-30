CTITLE WILD_DATE

      subroutine wild_date( infile, jd, option )

      implicit none

*     Subroutine to map wild card string in the
*     infile name

      integer*4 max_save   ! Maximum number of template names that can
                           ! saved.  Templates are based on changes in
                           ! file extent
      integer*4 num_keys   ! Number of keys available

      parameter ( max_save = 200 ) 
      parameter ( num_keys = 11 )

* INPUT variables
      real*8 jd   ! JD of experiment (may be MJD as well).  Used to get date

      character*(*) infile  ! Name of input file. On first call this
                            ! will have template which will be saved
                            ! for later calls
      character*(*) option  ! "I" as first character will increment file
                            !     name if no template strings are found
                            !     based on extent on file name (used when
                            !     infile will be used multiple times).
                            ! "N" will ignore the ext and simply substitute
                            !     into the name.

* LOCAL variables

      integer*4 j              ! Loop counter
     .,         nt             ! Number of template match
     .,         indx           ! Index to extent location in infile

      integer*4 gpsw, gpsd     ! GPS week and doy of week
      integer*4 date(5), doy   ! Date Y,M,D,hr,mn,sc and day of year
      real*8    sec            ! Seconds of day
      real*8    gpss           ! GPS seconds of week.

      logical   found_key      ! Set true if key wild found.  If not found
                               ! and option is I number added to file name.

      character*6 keys(num_keys)   ! Wild-card keys that can be used
      character*6 lckey,uckey      ! lower and upper case version of key
      character*4 code(num_keys)   ! JD coded into the key types

      character*16 inext, savext   ! File name extents (after last .) for 
                                   ! input and saved templates
* SAVE variables
      integer*4 num_save       ! Number of saved templates
     .,         count(max_save)  ! Counter with number of times called with given
                              ! templated (used in case no wild cards are
                              ! included and optiom set to I(ncrement)

      character*128 template(max_save)  ! Saved names of templates with wild cards

      data num_save  / 0 /    ! initial number of templates 
      data keys / '<WWWW>',   ! GPS week           1
     .            '<D>   ',   ! GPS day of week    2
     .            '<Y4>  ',   ! 4-digit year       3
     .            '<Y2>  ',   ! 2-digit year       4
     .            '<Y1>  ',   ! 1-digit year (arc) 5
     .            '<DOY> ',   ! Day of year        6
     .            '<MM>  ',   ! Month of year      7
     .            '<DD>  ',   ! Day of month       8
     .            '<HR>  ',   ! Hour of day        9
     .            '<MIN> ',   ! Minute of hour    10
     .            '<SEC> '  / ! Seconds of minute 11

      save num_save, count, template

****  Start: see if we have a template match for the extent in the file
      if( option.ne.'N' ) then   ! Option for when infilename may be  
                ! overwriten with update name.

         if( num_save.gt.0 ) then
*            We have saved templates.  See if match to infile
             call getextent(infile, inext, option)
             nt = 0   ! Template number
             do j = 1, num_save
                call getextent(template(j),savext, option)
                if( savext.eq.inext ) then
                    nt = j    ! Save number of template to use
                    exit      ! Exit from loop
                endif
             enddo
*            If nt is still zero, then there was no match so save new template
             if( nt.eq.0 ) then
                 num_save = num_save + 1
                 if( num_save.gt.max_save ) then
                     write(*,110) max_save, trim(infile)
 110                 format('**FATAL** Too many WILD_DATE templates. ',
     .                      'Exceeds max ',I3,' Template ',a)
                     stop 'WILD_DATE: Too many templates'
                 endif
                 template(num_save) = infile
                 nt = num_save
             endif
         else
*            Save infile as new template
             num_save = 1
             nt = 1
             template(nt) = infile
             call getextent(infile, inext, option)
             count(nt) = 0
         endif
      endif

****  Create the strings with the date encoded into the key formats
      call jd_to_gpst( jd, gpsw, gpsd, gpss )
      call jd_to_ymdhms( jd, date, sec )
      call ymd_to_doy( date, doy) 

      write(code( 1),'(I4.4)') gpsw
      write(code( 2),'(I1)')   gpsd
      write(code( 3),'(I4)')   date(1)
      write(code( 4),'(I2.2)') mod(date(1),100)
      write(code( 5),'(I1)')   mod(date(1),10)
      write(code( 6),'(I3.3)') doy
      write(code( 7),'(I2.2)') date(2)
      write(code( 8),'(I2.2)') date(3)
      write(code( 9),'(I2.2)') date(4)
      write(code(10),'(I2.2)') date(5)
      write(code(11),'(I2.2)') nint(sec)

****  Now take the template and replace the keys with date information
      found_key = .false.
      if( option.ne.'N' ) then   ! Only replace infile name if N option not
                                 ! selected.
         infile = template(nt)
      endif
      do j = 1, num_keys
          uckey = keys(j)
          indx = index(infile,trim(uckey))
          if( indx.eq.0 ) then
*            Try lower case version
             lckey = keys(j)
             call caseunfold(lckey)
             indx = index(infile,trim(lckey))
             if( indx.gt.0 ) then
                call sub_char(infile,lckey,code(j))
                found_key = .true.
             endif
          else

*             Sub_char will replace all occurences
              call sub_char(infile,uckey,code(j))
              found_key = .true.
          endif
      end do

****  See if key was not found. If not found and I used as option
*     add count(nt) to file name
      if( .not.found_key .and. option(1:1).eq.'I' ) then
          count(nt) = count(nt)+1
*         Only update the name if we have used the name more than
*         once.
          if( count(nt).gt.1 ) then 
             indx = index(infile,inext)
             write(infile(indx-1:),220) count(nt),inext
 220         format('_',I5.5,'.',a)
          endif
      endif
      return
      end
   
CTITLE GETEXTENT

      subroutine getextent(file, ext, option)  

      implicit none

*     Function to return extent string after last period.

* PASSED 
      character*(*) file    ! File name
      character*(*) ext     ! Extent string
      character*(*) option  ! Option: If N then missing extent is ignored.

* LOCAL
      integer*4 j   ! Counter
     .,         trimlen  ! Function to get length of string

***** Find the '.' in file name
      j = trimlen(file)
      ext = ' '
      do while ( j.gt.0 )
         if( file(j:j).eq.'.' ) then
*           Extent found
            ext = file(j+1:)
            j = -1
         else
            j = j - 1
         endif
      enddo
      if( j.eq.0 .and. option(1:1).eq.'I' ) then
*        No period found
         write(*,120) trim(file)
 120     format('**FATAL** no extent found in ',a,/,
     .          'Use of wild_date requires file with extent')
         stop 'WILD_DATE/GETEXTENT: No extent found'
      endif

****  Thats all
      return
      end



  


