CTITLE READ_SP3
 
      subroutine read_sp3 ( sp3_file )
 
      implicit none

*     Thuis routine will read sp3 epehemeris file. All ephemeris
*     entries are read
 
      include '../includes/sp3_def.h'
      include '../includes/const_param.h'

* PASSED 
      character*(*) sp3_file
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   i           - Loop counter
*   lpn         - Read Prn number: Modified for constellation type
*   trimlen     - length of string
*   ll          - Length of line string  (not used)
*   date(5)     - ymdhm
*   rinex_version - Version of rinex file
*   pn          - Assigned PRN number (lpn+constellation offset)
*   in          - Index for running current GNSS PRN
 
 
      integer*4 ierr, i,j,  trimlen,  date(5), lpn, pn, in
      real*4 rinex_version
      integer*4 un   ! Unit number for Sp3 file
 
*   sectag      - Seconds tag in date
 
      real*8 sectag
 
*   still_header    - Indicates that we are sill in the header
*                   - part of file.
*   have            - Set true if we already have this satellite (not used)
 
      logical still_header
 
*   line            - line read from file
 
      character*256 line

* MOD TAH 180312: Added for GNSS
      integer*4 conadd  ! Constant added to account for constellation 
                        ! type
     .,         conoff  ! Function to return constellation offset

      integer*4 concount(7)  ! Counts of number of oonstellations 
                        ! found

*   cr - Carriage return (for handling dos files)
      character*1 cr

      cr = char(13)
 
****  Open the sp3_file 
      un = 97
      open(un, file=sp3_file, status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',sp3_file,1,
     .            'svsp3/read_nav')
 
*     Loop over the file reading the ephemeris entries.  First clear
*     all of the PRN numbers so we know when a PRN has been read
 
      do i = 1, max_sat
          prn_sp3(i) = 0
      end do
      rinex_version = 1
 
      still_header = .true.
      do while ( still_header )
          read(un,'(a)', iostat=ierr) line
          if( line(1:1).eq.'*' ) still_header = .false.
          if( ierr.ne.0 ) then
              write(*,120) sp3_file(1:trimlen(sp3_file))
 120          format(' Error reading ',a,' before end of ',
     .               'header found')
              stop 'svsp3: Reading sp3 file'
           end if
      end do
 
      num_sat = 0
      num_sp3 = 1
      read(line(2:),*) date, sectag
      call ymdhms_to_mjd( date, sectag, sp3_time(num_sp3))
      sp3_refmjd = int(sp3_time(num_sp3)+1.d-6)
      sp3_sec(num_sp3) = date(4)*3600.0d0 + date(5)*60.d0 +
     .                   sectag
 
*     Now start reading the entries
      do while ( ierr.eq.0 )

*         Read down the elements of the satellites
          read(un,'(a)',iostat=ierr) line
*         Skip the EOF test so we can read concatenated files
!          if( line(1:3).eq.'EOF' ) then
!              ierr = -1
*         Clear line if end of file so at it would be processed.
          if( ierr.ne.0 ) line = ' '
          if (line(1:1).eq.'P' ) then

*            See what the constellation is
             if( line(2:2).eq.' ' ) line(2:2) = 'G'
             conadd = conoff(line(2:2))
             read(line(3:4),*) lpn
             pn = lpn + conadd
*            See if we already have this one
             in = num_sat + 1
             do j = 1, num_sat
                if( prn_sp3(j).eq. pn ) then
                   in = j
                   exit
                endif
             end do
             if( in.gt.num_sat ) then
                num_sat = in
                prn_sp3(in) = pn
             endif 
             if( in.gt.max_sat ) then
                 call report_stat('FATAL','svpos','read_sp3',' ',
     .              'Too many satellites: MAX ',max_sat)
             endif

             read(line(5:),*) 
     .              sp3_xyz(1,in,num_sp3),
     .              sp3_xyz(2,in,num_sp3),
     .              sp3_xyz(3,in,num_sp3),
     .              sp3_clk(in,num_sp3)
              sp3_xyz(1,in,num_sp3) = sp3_xyz(1,in,num_sp3)*1000.d0
              sp3_xyz(2,in,num_sp3) = sp3_xyz(2,in,num_sp3)*1000.d0
              sp3_xyz(3,in,num_sp3) = sp3_xyz(3,in,num_sp3)*1000.d0

*             Make sure clock is valid
              if( abs(sp3_clk(in,num_sp3)-999999.999999d0).lt.1) then
*                 No clock, either set to zero or to previous value
                  if( num_sp3.gt.3 ) then
                     sp3_clk(in,num_sp3) = sp3_clk(in,num_sp3-1)*1.d6
                  else
                     sp3_clk(in,num_sp3) = 0.d0
                  endif
              endif
*             Convert clock to seconds
              sp3_clk(in,num_sp3)   = sp3_clk(in,num_sp3)*1.d-6
          else if (line(1:1).eq.'*' ) then
              num_sp3 = num_sp3 + 1
              read(line(2:),*) date, sectag
              call ymdhms_to_mjd( date, sectag, sp3_time(num_sp3))
              sp3_sec(num_sp3) = date(4)*3600.0d0 + date(5)*60.d0 +
     .                   sectag + 
     .                   int(sp3_time(num_sp3)-sp3_refmjd)*86400.d0
          end if
      end do
 
      write(*,300) num_sat, num_sp3, sp3_file(1:trimlen(sp3_file))
 300  format('* ', i5,' satellites, at ',i4,' epochs found in ',a)
      write(*,320) prn_sp3(1:num_sat)
 320  format('Mapped PRNS found: ',100(I4))
 

      if( num_sat.eq.0 ) stop ' SVPOS: No satellites found'

****  Assign the frequencies assoicated with each satellite
      call assign_freqs
 
      return
      end

CTITLE ASSIGN_FREQS


      subroutine assign_freqs

      implicit none

*     Routine to assign the two primary frequencies and a possible
*     third frequency to each satellite found in the SP3 file.

      include '../includes/sp3_def.h'
      include '../includes/const_param.h'
      include '../../libraries/includes/freq_def.h'   ! GAMIT Frequency table


* LOCAL VARAIBLES
      integer*4 num_glonass  ! Number of GLONASS satellites (if none no need to
                             ! svnav.dat)
     .,         i,j          ! Loop counters
     .,         ierr         ! IOSTAT error

* Variables read from svnav,dat
      integer*4 gpn     ! Glonass PRN 
     .,         gslot   ! Slot of Glonaaa
     .,         gstryd(2)   ! Start YEAR DOY
     .,         gendyd(2)   ! End YEAR DOY

      real*8 gsjd, gejd     ! Start and end JD
     .,      sec            ! Seconds of day (not used).


      character*128 svnav_file  ! Name of svnav.dat file in ~/gg/tables
      character*128 line    ! Line read from svnav.dat file (99 needed).

      character*1 offcon    ! Function to return constellation letter based on
                            ! offset PRN (original PRN for mod(prn,100))


****  Loop over the PRNS to see type.  For GLOSNASS we also need to read svnav.dat
*     to get the slot for each prn.
*     See any Glonass
      num_glonass = 0

      do i = 1, num_sat
         if( prn_sp3(i).gt.100 .and. prn_sp3(i).lt.200 ) then
             num_glonass = num_glonass + 1
         end if
      end do

***** See if glonass found
      sec = 0.0d0
      if( num_glonass.gt.0 ) then
          write(*,'("GLONASS FOUND: Getting slots")') 
* MOD TAH 210201: Point svnav to original format svnav.dat.all.gnss explictly
*         because new link to svnav.dat is to igs_metadata.snx which track
*         does not year read.
          svnav_file = '~/gg/tables/svnav.dat.allgnss'
          call subhome(svnav_file)
          open(103,file=svnav_file,iostat=ierr,status='old')
          call report_error('IOSTAT',ierr,'open',svnav_file,0,
     .                      'SVSP3/ASSIGN_FREQS')
          do while ( ierr.eq.0 )
             read(103,'(a)',iostat=ierr) line
             if( ierr.eq.0 ) then
*               See if GLONASS
                if( line(1:3).eq.' R ' ) then
*                   Decode line
                    read(line, 120) gpn, gslot, gstryd, gendyd
 120                format(10x,I2,1x,I3,53x,I4,1x,I3,8x,I4,1x,i3)
CR   802   9  -2  GLONASS-K1              935000.     U      0.2500  2016  47  0  0  2100   1  0  0
                    call yds_to_mjd( gstryd(1), gstryd(2), sec, gsjd)
                    call yds_to_mjd( gendyd(1), gendyd(2), sec, gejd)
                    if( sp3_time(1)      .ge. gsjd .and. 
     .                  sp3_time(num_sp3).le. gejd ) then
*                       in time range, see if we have in SP3
                        do j = 1, num_sat
                           if( gpn.eq.prn_sp3(j)-100) then
*                              Found match, so save frequencies
                               fL1(j) = glonass_f1 + glonass_df1*gslot
                               fL2(j) = glonass_f2 + glonass_df2*gslot
                               fL5(j) = glonass_f3 + glonass_df3*gslot
                               exit
                           end if
                        enddo
                    end if
                endif  
             end if
          end do
          close(103)
      endif

****  OK now assign everyone one
      do i = 1,num_sat
*        Check PRN numbers
         if( prn_sp3(i).gt.0 .and.prn_sp3(i).lt.100 ) then   ! GPS
             fL1(i) = gps_f1 
             fL2(i) = gps_f2
             fL5(i) = gps_f5
         elseif ( prn_sp3(i).gt.200 .and.prn_sp3(i).lt.300 ) then ! Galileo
*            Galileo has more frequencies so maybe expand later
             fL1(i) = galileo_f1   ! E1
             fL2(i) = galileo_f5   ! E5a
             fL5(i) = galileo_f7   ! E5b              
         elseif ( prn_sp3(i).gt.300 .and.prn_sp3(i).lt.400 ) then ! Beidou
* MOD TAH/RWK 200223: Updated frequencies for BDS-3 to reflect most common
*            observed frequencis
C            fL1(i) = beidou_f2    ! B1
C            fL2(i) = beidou_f7    ! B2
C            fL5(i) = beidou_f6    ! B3
             fL1(i) = beidou_f2    ! B1-2
             fL2(i) = beidou_f6    ! B3
             fL5(i) = beidou_f7    ! B2b  
         elseif ( prn_sp3(i).gt.400 .and.prn_sp3(i).lt.500 ) then ! QZSS
*            QZSS has a LEX 6 frequency as well/
             fL1(i) = gps_f1       ! Same as GPS
             fL2(i) = gps_f2
             fL5(i) = gps_f5
         elseif ( prn_sp3(i).gt.500 .and.prn_sp3(i).lt.600 ) then ! SBAS
             fL1(i) = gps_f1       ! L1
             fL2(i) = gps_f5       ! L5
             fL5(i) = 0.0d0        ! Not defined
         elseif ( prn_sp3(i).gt.600 .and.prn_sp3(i).lt.700 ) then ! IRNSS
             fL1(i) = irnss_f5     ! IRNSS I5
             fL2(i) = irnss_f5     ! IRNSS S-band  (L9)
             fL5(i) = 0.d0         ! not defined
         endif
      end do

****  Write out what we have for frequencies
      if( out_sp3freqs ) then
         write(*,220) 
 220     format('* CONSTELLATION FREQUENCY TABLE',/,
     .          '*SYS PRN   F1 (MHz)   F2 (MHz)   F3 (MHz)')
         do i = 1, num_sat
            write(*,240) offcon(prn_sp3(i)), mod(prn_sp3(i),100), 
     .                   fL1(i)/1.d6, fL2(i)/1.d6, fL5(i)/1.d6
 240        format('* ',a1,2x,I3,1x,3(F10.3,1x))
         end do
      end if

***** Thats all
      return
      end

CTITLE PtoL

      integer*4 function PtoL( prn )

*     Function to return the list entry for a given PRN number
*     If the PRN is not in the list -1 is return.  (To be called
*     after the SP3 file has been read).

      include '../includes/sp3_def.h'

* PASSED 
      integer*4 prn   ! PRN number offset for the GNSS constellation.

* LOCAL
      integer*4 i     ! Loop counter

      PtoL = -1
      do i = 1, num_sat
         if( prn.eq.prn_sp3(i) ) then
           PtoL = i
           exit
         end if
      end do

****  Thats all
      return
      end

CTITLE READ_ANTMOD_SVS

      subroutine read_antmod_svs

      implicit none

*     Routine to read and save satellite antenna offsets for the satelites in the
*     sp3 file being used.
*
      include '../includes/sp3_def.h'


* LOCAL VARIABLES 
      integer*4 i,j        ! Loop counters
     .,         ierr, jerr ! IOSTAT error
     .,         trimlen    ! Length of string
     .,         date(5)    ! YMDHM for calender date
     .,         antprn     ! PRN number for satellite, set 0 for site 
                           ! and -1 for entry not be processed.
     .,         lv         ! List vechicle number from (PtoL function)
     .,         ifrq       ! Frequency being read: 1 -- L1, 2 -- L2, 
                           ! 12 -- L1 and L2
     .,         nfrq       ! Number of frequencies to be read (ignored)
     .,         nnad       ! Number of Nadir  values to read
     .,         naz        ! Number of AZ values to read
     .,         inaz       ! Index for AZ lines.
     .,         un         ! Unit number

      integer*4 PtoL       ! Function to return satellite number from list
                           ! given PRN (offset by constellation).

      real*8    sec_tag    ! Seconds tag
     .,         mjd_vfrom, mjd_vuntil ! MJD for Valid from and Valid
                           ! until entries.
     .,         dPos(3)    ! Either NEU for site or XYZ for satellite (mm)
     .,         vals(max_dna)   !  Phase center values
     .,         azval      ! Value of azimuth

      character*1  svsys   ! Satellite systems in antex file
      character*8  noazi   ! String with NOAZI values

      character*20 label   ! Head label of file type
     .,            anttypx ! Antenna or Satellite type
     .,            antsnx  ! Satellite PRN number or serial number

      character*10 antid3,antid4  ! Satellite SV number and COSPAR number

      real*8 antex_vers    ! Version of antex file
     .,      dazi          ! Spacing in azimuth (0 for no depedence).
     .,      dzen(3)       ! Start, stop and dzen values

      logical eoh          ! Set true when end of header reached
     .,       eof          ! Set true when EOF found
     .,       site_found   ! Set true if antenna matches at least one site
     .,       read_dph     ! Set true when the next set of lines are the
                           ! phase center model

      character*1024 line  ! Lines from files (allows 91 values with f8.2 
                           ! format to be saved).
      character*256  file  ! Name of antmod.dat in ~/gg/tables

***** Try to open files
      un = 98
      file = '~/gg/tables/antmod.dat' 
      call subhome(file) 
      open(un,file=file,status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',file,0,
     .                  'read_antsvs')
      if( ierr.ne.0 ) RETURN

****  OK Start reading file
      read(un,'(a)',iostat=ierr) line
      call report_error('IOSTAT',ierr,'head read',line,0,
     .                  'read_antsvs/Header')
      if ( ierr.ne.0 ) RETURN
      read(line,'(f8.1,12x,a1,39x,a20)', iostat=ierr) 
     .    antex_vers, svsys, label
      call report_error('IOSTAT',ierr,'head decod',line,0,
     .                  'read_antsvs/Header')
      if ( ierr.ne.0 ) RETURN
      call casefold(line)
      if( label(1:5).ne.'ANTEX' ) then
          write(*,120) file(1:trimlen(file)), line(1:trimlen(line))
 120      format('ANTMOD_FILE ',a,' is not ANTEX.  First line is ',/,a)
          RETURN
      end if

****  Now skip down until we find end of header
      eoh = .false.
      do while ( .not. eoh )
          read(un,'(a)',iostat=ierr) line
          call casefold(line)
          if( line(61:73).eq.'END OF HEADER' ) eoh = .true.
          call report_error('IOSTAT',ierr,'read',file,0,'read_antsvs')
          if( ierr.ne.0 ) RETURN
      end do

****  OK: See what comes next
      eof = .false.
      antprn = -1 
      lv = -1 
      do while ( .not.eof )
         read(un,'(a)',iostat=ierr) line
         if ( ierr.eq. -1 ) eof = .true.
         if( ierr.ne.0 .and. ierr.ne.-1 ) then
            call report_error('IOSTAT',ierr,'read',file,0,'read_antsvs')
            close(un)
            RETURN
         end if
         call casefold(line)
         IF( ierr.eq.0 .and. trimlen(line).gt.0 .and. 
     .                                   line(61:67).ne.'COMMENT' ) THEN

***      See type of line
         if( line(61:64).eq. 'TYPE' ) then
            read(line,'(a20,a20,a10,a10)',iostat=ierr) anttypx,antsnx
     .         ,antid3,antid4
            if( anttypx(1:5).eq.'BLOCK' ) then
*               Set the PRN number to which the following entries
*               applies
                read(antsnx(2:3),'(I2)',iostat=jerr) antprn
                call report_error('IOSTAT',jerr,'decod',antsnx,0,
     .                'read_antsvs/TYPE')
                lv = PtoL(antprn) 
            else if( anttypx(1:7).eq.'GLONASS' ) then
                read(antsnx(2:3),'(I2)',iostat=jerr) antprn
                call report_error('IOSTAT',jerr,'decod',antsnx,0,
     .                'read_antsvs/TYPE')
                antprn = antprn + 100
                lv = PtoL(antprn) 
            else if( anttypx(1:7).eq.'GALILEO' ) then
                read(antsnx(2:3),'(I2)',iostat=jerr) antprn
                call report_error('IOSTAT',jerr,'decod',antsnx,0,
     .                'read_antsvs/TYPE')
                antprn = antprn + 200
                lv = PtoL(antprn) 
            else if( anttypx(1:7).eq.'BEIDOU' ) then
                read(antsnx(2:3),'(I2)',iostat=jerr) antprn
                call report_error('IOSTAT',jerr,'decod',antsnx,0,
     .                'read_antsvs/TYPE')
                antprn = antprn + 300
                lv = PtoL(antprn) 
            else
                antprn = -1
            endif
         endif

         if( line(61:64).eq. 'DAZI' .and. lv.gt.0 ) then
             read(line,'(2x,f6.1)',iostat=jerr) dazi 
             call report_error('IOSTAT',jerr,'decod',line,0,
     .               'read_antsvs/DAZI')
         endif
*
         if( line(61:64).eq. 'ZEN1' .and. lv.gt.0 ) then
             read(line,'(2x,3f6.1)',iostat=jerr) dzen 
             call report_error('IOSTAT',jerr,'decod',line,0,
     .               'read_antsvs/DZEN')
             do j = 1,3
                svs_dna(j,1,lv) = dzen(j)
             end do
             num_svs_dph(1,lv) = 
     .                   nint((dzen(2)-dzen(1))/dzen(3))+1
             num_svs_dph(2,lv) = 1
         endif
         if( line(61:70).eq. 'VALID FROM' .and. lv.gt.0 ) then
             read(line,'(5i6,f13.7)',iostat=jerr) date, sec_tag
             call report_error('IOSTAT',jerr,'decod',line,0,
     .               'read_antsvs/VALID FROM')
*            Check that our data is after this date
             call ymdhms_to_mjd(date, sec_tag, mjd_vfrom)
             if ( sp3_time(1).lt.mjd_vfrom ) then
*                Our data is before this time, so don't use
                 lv = -1
             end if
         endif
         if( line(61:71).eq. 'VALID UNTIL' .and. lv.gt.0 ) then
             read(line,'(5i6,f13.7)',iostat=jerr) date, sec_tag
             call report_error('IOSTAT',jerr,'decod',line,0,
     .               'read_antsvs/VALID UNTIL')
*            Check that our data is before this date
             call ymdhms_to_mjd(date, sec_tag, mjd_vuntil)
             if ( sp3_time(1).gt.mjd_vuntil ) then
*                Our data is before this time, so don't use
                 lv = -1
             end if
         endif

*
         if( line(61:73).eq. 'START OF FREQ' .and. lv.gt.0 ) then
             read(line,'(3x,1x,i2)',iostat=jerr) ifrq
             call report_error('IOSTAT',jerr,'decod',line,0,
     .               'read_antsvs/START OF FREQ')
             read_dph = .true.
         endif

         if( line(61:71).eq. 'END OF FREQ' .and. lv.gt.0 ) then
             read(line,'(3x,1x,i2)',iostat=jerr) ifrq
             call report_error('IOSTAT',jerr,'decod',line,0,
     .               'read_antsvs/END OF FREQ')
             read_dph = .false.
         endif

         if( line(61:77).eq. 'NORTH / EAST / UP' .and. 
     .                                           lv.gt.0 ) then
             read(line,'(3f10.2)',iostat=jerr) dPos
             call report_error('IOSTAT',jerr,'decod',line,0,
     .               'read_antsvs/NORTH / EAST / UP')
*            Save values in correct slots
*            This is satellite so save values in slot
             if( ifrq.le.2 ) then
                do j = 1, 3
                   svs_L12(j,ifrq,lv) = dPos(j)/1000.d0
                enddo
             else
                do j = 1, 3
                   svs_L12(j,1,lv) = dPos(j)/1000.d0
                   svs_L12(j,2,lv) = dPos(j)/1000.d0
                enddo
             end if
         endif

*****    See if we are in READ_DPH mode.  Just read lines
         if( read_dph .and. lv.gt.0 .and.
     .       line(61:73).ne. 'START OF FREQ' .and.
     .       line(61:77).ne. 'NORTH / EAST / UP' ) then
*            Get values from line
             nnad = num_svs_dph(1,lv)
             naz  = num_svs_dph(2,lv)

             if( nnad.gt.max_dna )
     .          call report_stat('FATAL','SP3_LIB','read_antmod_svs',
     .              file, 'Too Many Nadir values, Num ',nnad)
             

             read(line,'(a8,91f8.2)',iostat=jerr) 
     .                 noazi, (vals(i),i=1,nnad)
             call report_error('IOSTAT',jerr,'decod',line,0,
     .               'read_antsvs/PHASE Values')

*            See what we have
             if( noazi.eq.'   NOAZI' .and. naz.eq.1 .and. 
     .                                                 lv.gt.0 ) then
*                This is the azimuth averaged values.  If dazi is zero
*                use these values.
*                Satellite values
                 if ( ifrq.le.2 ) then
                    do j = 1, nnad
                        svs_dphs(j,1,ifrq,lv) = vals(j)/1000.d0
                    enddo
                 else
                    do j = 1, nnad
                        svs_dphs(j,1,1,lv) = vals(j)/1000.d0
                        svs_dphs(j,1,2,lv) = vals(j)/1000.d0
                    enddo
                 endif
                 if( out_sp3PCV ) 
     .           write(*,200) prn_sp3(lv), nnad, ifrq
     .,               (svs_L12(j,1,lv), j=1,3)
     .,               (svs_L12(j,2,lv), j=1,3)
     .,               (svs_dphs(j,1,ifrq,lv)*1000.d0,j=1,nnad)
 200             format('SVS PRN ',i3,' #NADIR ',i3, ' IFRQ ',i3,/
     .,                 'SVS L1 ',3F8.4, ' L2 ',3F8.4,/
     .,                 'DPH L1 ',15F8.2) 
             end if
         end if

*        Done reading file
         ENDIF  ! Not EOF
      end do

      close(un)

      return
      end

