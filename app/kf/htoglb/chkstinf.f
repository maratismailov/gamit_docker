
      Program CHKSTINF

      implicit none

*     This program will read the headers from a set of rinex files and 
*     check the match between the rinex information and station,info
*     Any differences will be listed as a warning.  If the rinex file
*     does nor exist, the site name and date will be check to see if
*     there is an entry.

*     Program only works with new format station.info in the standard format.
           
*     Calling arguments:
*     chkstinf <station.info> < ... List of rinex files ...>
*



* LOCAL VARIABLES     
* ierr -- IOSTAT error
* i    -- Loop counter

      integer*4 ierr, i 
     .,  rcpar   ! Reads runstring
     .,  len_run ! Length of runstring

      integer*4 yr_in, doy_in  ! Year and doy of data
       
      real*8 anth, antn, ante  ! Antenna Height, North and East

      character*16 anttyp  ! Type of antenna
      character*4  radome  ! Type of radome
     .,            site    ! Site 4-char code
      character*20 rctype  ! Receiver type
      character*256 ref_stinf ! Name of station.info file
     .,            rx_file  ! Name of RX file

****  Get the name of station.info
      len_run = rcpar(1,ref_stinf)
      if ( len_run.eq.0 ) then
          call proper_runstring('chkstinf.hlp','chkstinf',1)
      endif

*     Now open he reference station.info file
      open(100, file=ref_stinf, status='old', iostat=ierr)
      if( ierr.ne.0 )  call report_stat('FATAL','CHKSTINF'
     .   ,'chksting',ref_stinf 
     .   ,'Error opening reference station.info file', ierr)  

****  OK now loop over rinex files
      i = 0
      do while ( len_run.gt.0 )
         i = i + 1
         len_run = rcpar(i+1,rx_file)
         if( len_run.gt.0 ) then

*            OK, found a name see if we can open RX file
             open(101,file=rx_file, status='old',iostat=ierr)
             call get_rxinf(101,ierr, rx_file,site, yr_in, doy_in, 
     .                 anth, antn, ante, rctype, anttyp, radome)
             close(101)

*             Now see if match found is station.info
              call check_stinf(100, rx_file,site, yr_in, doy_in, 
     .                 anth, antn, ante, rctype, anttyp, radome)

         endif
      end do

***** Thats all 
      close(100)
      end

CTITLE GET_RXINF

      subroutine get_rxinf(unr ,ierr, rx_file,site, yr_in, doy_in, 
     .                 anth, antn, ante, rctype, anttyp, radome)

      implicit none


      integer*4 unr    ! Unit number for RINEX file
     .,         ierr   ! IOSTAT error: Passed in and if error only
                       ! site name and doy_in set to valid values
     .,         yr_in, doy_in  ! Year and doy of data
     .,         trimlen ! Length of string 

      real*8 anth, antn, ante  ! Antenna Height, North and East

      character*(*) anttyp  ! Type of antenna
      character*(*) radome  ! Type of radome
      character*(*) rctype  ! Receiver type
      character*(*) site    ! 4-char site name (returned upper case)
      character*(*) rx_file  ! Name of rinex file

* LOCAL VARIABLES
      integer*4 ymd(3)   ! Yr, month day
     .,   rxe   ! index to end of rinex file name 
     .,   ind   ! index for finding strings
     .,   jerr  ! IOSTAT error

      character*80 line  ! Line read from rx file
      character*4 doy_ch, yr_ch  ! Character version of doy and yr



****  Get the site name from the rx_file name.  Remove any compresssion
*     extents
      call sub_char(rx_file,'.Z',' ' )
      call sub_char(rx_file,'.gz',' ' )
      rxe = trimlen(rx_file)
      site = rx_file(rxe-11:rxe-8)
      call casefold(site)
      doy_ch = rx_file(rxe-7:rxe-5)
      read(doy_ch,*,iostat=jerr) doy_in
      yr_ch =  rx_file(rxe-2:rxe-1)
      read(yr_ch,*,iostat=jerr) yr_in
      if( yr_in .lt. 50 ) then
          yr_in = yr_in + 2000
      else
          yr_in = yr_in + 1900
      end if

***** Now if there was no error opening file start readind rx_file
      rctype  = 'UNK'
      anttyp = 'UNK'
      radome  = 'UNK'
      do while ( ierr.eq.0 ) 
         read(unr,'(a)',iostat=ierr) line

*        Look for data records we need
         ind = index(line,'END OF HEADER')
         if( ind.gt.0 ) ierr = -1
         ind = index(line,'REC # / TYPE / VERS')
         if( ind.gt.0 ) then
             rctype = line(21:40)
         endif
         ind = index(line,'ANT # / TYPE')
         if( ind.gt.0 ) then
             anttyp = line(21:36)
             radome  = line(37:40)
         endif
         ind = index(line,'TIME OF FIRST OBS')
         if( ind.gt.0 ) then
             read(line,*,iostat=jerr) ymd
             call ymd_to_doy( ymd, doy_in )
             yr_in = ymd(1)
         endif
         ind = index(line,'ANTENNA: DELTA H/E/N')
         if( ind.gt.0 ) then
             read(line,*,iostat=jerr) anth, antn, ante
         endif

      end do

***** Thats all
      return
      end
   
CTITLE CHECK_STINF

      subroutine check_stinf(uns, rx_file,site, yr_in, doy_in, 
     .                 anth, antn, ante, rctype, anttyp, radome)

      implicit none

      integer*4 uns    ! Unit number for RINEX file
     .,         yr_in, doy_in  ! Year and doy of data

      real*8 anth, antn, ante  ! Antenna Height, North and East

      character*(*) anttyp  ! Type of antenna
      character*(*) radome  ! Type of radome
      character*(*) rctype  ! Receiver type
      character*(*) site    ! 4-char site name (returned upper case)
      character*(*) rx_file  ! Name of rinex file

* LOCAL VARIABLES
      integer*4   jerr  ! IOSTAT error
     .,   ierr   ! IOSTAT error: Passed in and if error only
                       ! site name and doy_in set to valid values
     .,   yr_st, doy_st  ! Start year and doy for stinf entry
     .,   yr_en, doy_en  ! End yera and doy for stinf entry
     .,   pass   ! Pass number in reading stinf file.

      real*8 ht   ! Height in station.info

      character*256 line  ! Line read from station.info
      character*80  mess  ! Report line 

      character*4 site_in  ! site in stinf line

      logical found  ! Set true when entry found or eof
     .,       OK     ! True if line match


***** Continue reading done station.info (based on names being alphabetical)
      found = .false.
      pass = 0

      do while ( .not.found )
         read(uns,'(a)',iostat=ierr) line
         if( ierr.ne.0 ) then
*           We have reached EOF.  If this is pass 0, rewind file and
*           try again
            if ( pass.eq.0 ) then
                 pass = 1
                 rewind(uns)
            else
*                Site was not found.  Print warning and return
                 found = .true.
                 write(mess, 220 ) site, yr_in, doy_in
 220             format(A4,1x,I4,1x,I3,' not in station.info')
                 call report_stat('WARNING','CHKSTINF','check',rx_file,
     .                             mess, 0)
            endif
         else
            if( line(1:1).eq.' ' ) then
*               see if name matches
                site_in = line(2:5)
                call casefold(site_in)
                if( site.eq.site_in ) then
*                   Name match.  Get the time range for this entry
                    read(line(26:33),*,iostat=jerr) yr_st, doy_st
                    read(line(45:52),*,iostat=jerr) yr_en, doy_en
                    if( yr_st+doy_st/366.0 .le. 
     .                  yr_in+doy_in/366.0) then
*                       Data if after start; see if before end
                        if( yr_en+doy_en/366.0 .gt. 
     .                      yr_in+doy_in/366.0) then
*                           OK: Times are in range.  See if rcv match
                            OK = .true.
                            found = .true.

*                           Do someclean up of names for common RINEX
*                           errors
                            call fixnames(rctype, anttyp, radome)
 
                            if( rctype .ne.line( 98:117) ) OK = .false.
                            if( anttyp.ne.line(171:186) ) OK = .false.
                            if( radome .ne.line(188:191) ) OK = .false.
                            if( .not. OK ) then
                               write(mess, 320) site, yr_in, doy_in
  320                          format(A4,1x,I4,1x,I3,1x,
     .                               'mismatch in station.info')
                               call report_stat('WARNING','CHKSTINF',
     .                            'check',rx_file, mess, 0)
                               write(*,330) site, yr_in, doy_in, 
     .                             rctype, line( 98:117),
     .                             anttyp, line(171:186),
     .                             radome, line(188:191)
  330                          format(A4,1x,I4,1x,I3,1x,
     .                             'RCVR ',2(a,1x),'| ANT ',2(a,1x),
     .                             '| Radome ',2(a,1x))
 

                            endif
*                           Check the antenna height
                            read(line(62:70),*,iostat=jerr) ht
                            if( ht.ne.anth ) then
                               write(mess, 360) site, yr_in, doy_in
  360                          format(A4,1x,I4,1x,I3,1x,
     .                               'height error in station.info')
                               call report_stat('WARNING','CHKSTINF',
     .                            'check',rx_file, mess, 0)
                            end if    

                        end if
                    end if
                end if
            end if
         end if
      end do

***** Thats all
      return
      end

CTITLE FIXNAMES

      subroutine fixnames(rctype, anttyp, radome)

      implicit none

*     Routine to fix common errors in names

* PASSED
      character*(*) anttyp  ! Type of antenna
      character*(*) radome  ! Type of radome
      character*(*) rctype  ! Receiver type

****  Common errors
      call sub_char(rctype,'MICROZ',       'UZ-12')
      call sub_char(anttyp,'AOAD/M_TA_NGS','AOAD/M_T')
      call sub_char(anttyp,'ASH701945.02B','ASH701945B_M')
      call sub_char(radome,'    ',         'NONE' )

***** Thats all
      return
      end
