CTITLE GET_UNR_ANT

      subroutine get_unr_ant( code )

      implicit none

*     Routine to open and read station.info to get the
*     times and antenna types for a given 4-char site name.
*     Once all entries are read, open and read antmod.dat
*     to get the G01 and G02 dNEU values.  
*     These values are later used to fix the East coordinate
*     error where the values for East have been applied with
*     the wrong sign.
*
      include 'tssum.h'
      include 'tscon.h'

* PASSED
      character*(*) code ! 4-character site code
      integer*4 nent_start  ! First entry number for this site
     .,         nent_end    ! End entry number for this site

* LOCAL

      integer*4 ierr, jerr  ! IOSTAT errors
     .,    yr, doy    ! Year and DOY 
     .,    nfnd       ! Number antenna NEU values founf
     .,    j          ! Loop counter.
     .,    freq       ! Gxx frequency number.  Need G01 and G03
     .,    trimlen    ! GG function to retrun length of string


      real*8 mjd, sec   ! MJD and Second of day
     .,      dNEU(3,2)  ! North/East/Up antenna offsets at G01/G02
     .,      dE_LC      ! Estimated offset in East LC

      logical done, found  ! Boolean values for ending file reads.
      logical vals_read    ! Set true when we know dNEU value (same 
                           ! antenna type multiple times)
      logical update       ! Set true if antenna or radome change and
                           ! we need to save new information

      character*16  ant     ! Antenna type
      character*8   radome  ! Radome type
      character*256 line    ! Line read from file. 

*     Open station info on unit 104.
      if( trimlen(stinf_file).eq.0 ) stinf_file = 'station.info'
*     Try local version firt
      open(104,file=stinf_file,iostat=ierr, status='old')
      if( ierr.ne.0 ) then   ! Try ~/gg/tables/station.info
          stinf_file = '~/gg/tables/station.info'
          call subhome( stinf_file )
          open(104,file=stinf_file,iostat=ierr, status='old')
          if( ierr.ne.0 ) then
              call report_error('IOSTAT',ierr,'open',stinf_file,
     .                          1,'GET_UNR_ANT')
          endif
      endif

*     Now start reading file for entries
      done = .false.
      num_stinf = 0
      sec = 0
      do while ( .not.done )
         read(104,'(a)',iostat=ierr) line
         if( code(1:4).eq.line(2:5) .and. line(1:1).eq.' ') then
*           Found site, get antenna info and time
            read(line(25:),*,iostat=jerr) yr,doy
            read(line(171:186),'(a)') ant
            read(line(188:193),'(a)') radome
*           Convert time to MJD
            call yds_to_mjd( yr, doy, sec, mjd)
            update = .true.
            if( num_stinf.gt.0 ) then
*               Only save if changed antenna type
                 if ( stinf_ant(num_stinf).eq. ant .and.
     .                stinf_rad(num_stinf).eq. radome ) then
*                  No change so no need to update
                   update = .false.
                endif
            endif
            if( update ) then

               num_stinf = num_stinf + 1
               if( num_stinf.gt. max_stinf ) then
                  call report_error('FATAL',max_stinf,'Too many stinf',
     .                 code,1,'get_unr_ant')
               endif
               stinf_mjd(num_stinf) = mjd
               stinf_ant(num_stinf) = ant
               stinf_rad(num_stinf) = radome
               stinf_sta(num_stinf) = 0       ! Set ant entry not set yet
            end if
         elseif( num_stinf.gt.0 ) then
            done = .true.
         endif
         if( ierr.ne.0 ) then
            write(*,120) code(1:4) 
 120        format('**WARNING** NO STATION.INFO entries for ',a)
            done = .true.
         endif
      end do

***** Clear any existing dNEU values (If antenna not found then zeros
*     will not change time series values.)
      do j = 1, num_stinf
          stinf_neu(:,:,j) = 0.0d0
      end do

***** Close station.info and open antmod.dat to see what the NEU offsets 
*     are for each antenna
      close(104)
      if( trimlen(antmod_file).eq.0 ) antmod_file = 'antmod.dat'
*     Try local version firt
      open(104,file=antmod_file,iostat=ierr, status='old')
      if( ierr.ne.0 ) then   ! Try ~/gg/tables/station.info
          antmod_file = '~/gg/tables/antmod.dat'
          call subhome( antmod_file )
          open(104,file=antmod_file,iostat=ierr, status='old')
          if( ierr.ne.0 ) then
              call report_error('IOSTAT',ierr,'open',antmod_file,
     .                          1,'GET_UNR_ANT')
          endif
      endif

****  Start reading the file and saving the G01 and G02 dNEU values
      done = .false.
      do while ( .not. done )
         read(104,'(a)',iostat=ierr ) line
         if( line(61:64).eq.'TYPE' ) then
*           OK check our list and see if one we need
            ant = line(1:16)
            radome = line(17:20)
*           Loop over our list to see if match
            nfnd = 0   ! Number of entries found
            vals_read = .false.
            do j = 1,num_stinf
               if( stinf_sta(j).le.0 ) then
*                  Check antenna
                   if( ant.eq.stinf_ant(j) ) then
*                      Antenna type matches.  If radome is NONE
*                      save NONE if radome has not been found
                       if( stinf_rad(j).eq.radome .or. 
     .                     radome.eq.'NONE' ) then
*                          Read and save values
                           found = vals_read
                           do while ( .not.found )
                              read(104,'(a)') line
                              if( line(61:65).eq.'START' .and.
     .                            line(4:4).eq.'G' ) then
                                  read(line(5:6),*) freq
*                                 read next line to get dNEU values
                                  read(104,*) dNEU(:,freq)
                                  if( freq.eq.2 ) found = .true.
                              endif
                              vals_read = .true.
                           enddo
*                          Save values
                           stinf_neu(:,1,j) = dNEU(:,1)
                           stinf_neu(:,2,j) = dNEU(:,2)
 
                           if( stinf_rad(j).eq.radome ) then
                              stinf_sta(j) = 1
                           else
                              stinf_sta(j) = -1   ! Show value is NONE radome
                           endif
                       endif
                   endif
               else
                   if( stinf_sta(j).gt.0 ) nfnd = nfnd + 1
               endif
            enddo
*           See if we have values for all entries
            if( nfnd.eq. num_stinf ) done = .true.
         elseif( ierr.ne.0 ) then
            done = .true.
            write(*,220) code
 220        format('**WARNING** Not all antenna entries found for ',a)
         endif
      end do
*     Close file
      close(104)

****  Report the values being used.
      do j = 1,num_stinf
         dE_LC = 2*( 2.54573*stinf_neu(2,1,j) - 
     .               1.54573*stinf_neu(2,2,j) )
         write(*,240) code, j, stinf_mjd(j), stinf_ant(j), stinf_rad(j),
     .      stinf_sta(j), stinf_neu(:,1,j), stinf_neu(:,2,j), dE_LC
 240     format('UNR ',a,1x,I3,1x,F13.4,1x,a,1x,a,1x,I2,
     .      ' G01 ',3(1x,F8.2),' G02 ',3(1x,F8.2),' dE LC ',F8.2)
      end do

***** Thats all
      return 
      end 


