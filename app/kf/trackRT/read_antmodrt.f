CTITLE READ_ANTMOD

      subroutine read_antmodRT(file)

      implicit none

*     Routine to read and save phase center models.
*
      include 'trackRT.h'

* PASSED  VARIABLES
      character*(*) file  ! Name of file. If open OK, added to
*           list of files

* LOCAL VARIABLES 
      integer*4 i,j        ! Loop counters
     .,         ierr, jerr ! IOSTAT error
     .,         trimlen    ! Length of string
     .,         date(5)    ! YMDHM for calender date
     .,         antprn     ! PRN number for satellite, set 0 for site 
                           ! and -1 for entry not be processed.
     .,         antsit(max_site)   ! Site numbers for antenna type
     .,         ifrq       ! Frequency being read: 1 -- L1, 2 -- L2, 
                           ! 12 -- L1 and L2
     .,         nfrq       ! Number of frequencies to be read (ignored)
     .,         nzen       ! Number of zenith values to read
     .,         naz        ! Number of AZ values to read
     .,         inaz       ! Index for AZ lines.
     .,         un         ! Unit number

      real*8    sec_tag    ! Seconds tag
     .,         mjd_vfrom, mjd_vuntil ! MJD for Valid from and Valid
                           ! until entries.
     .,         dPos(3)    ! Either NEU for site or XYZ for satellite (mm)
     .,         vals(max_zen)   !  Phase center values
     .,         azval      ! Value of azimuth

      character*1  svsys   ! Satellite systems in antex file
      character*8  noazi   ! String with NOAZI values

      character*20 label   ! Head label of file type
     .,            anttypx ! Antenna or Satellite type
     .,            antsnx  ! Satellite PRN number or serial number
     .,            antassigned(max_site)  ! Antenna assigned to site (radome may be missing).

      character*10 antid3,antid4  ! Satellite SV number and COSPAR number

      real*8 antex_vers    ! Version of antex file
     .,      dazi          ! Spacing in azimuth (0 for no depedence).
     .,      dzen(3)       ! Start, stop and dzen values

      logical eoh          ! Set true when end of header reached
     .,       eof          ! Set true when EOF found
     .,       site_found   ! Set true if antenna matches at least one site
     .,       read_dph     ! Set true when the next set of lines are the
                           ! phase center model
     .,       fullfnd(max_site), partfnd(max_site)

      character*1024 line  ! Lines from files (allows 91 values with f8.2 
                           ! format to be saved).

***** Try to open files
      un = 300
      open(un,file=file,status='old', iostat=ierr)
      call report_error('IOSTAT',ierr,'open',file,0,
     .                  'read_antmod')
      if( ierr.ne.0 ) RETURN

****  See if antenna information is known:
      if( .not. ante_off ) then
         write(*,100)
 100     format('***FATAL** antmod_file command issued before ',
     .          'ante_off information for sites')
         stop 'TRACK: Use ANTE_OFF command before ANTMOD_FILE'
      endif

****  OK Start reading file
      read(un,'(a)',iostat=ierr) line
      call report_error('IOSTAT',ierr,'head read',line,0,
     .                  'read_antmod/Header')
      if ( ierr.ne.0 ) RETURN
      read(line,'(f8.1,12x,a1,39x,a20)', iostat=ierr) 
     .    antex_vers, svsys, label
      call report_error('IOSTAT',ierr,'head decod',line,0,
     .                  'read_antmod/Header')
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
          call report_error('IOSTAT',ierr,'read',file,0,'read_antmod')
          if( ierr.ne.0 ) RETURN
      end do

****  Reached end of header.  Now read entries in file.
      num_antmods = num_antmods + 1
      write(*,140) num_antmods, antex_vers, file(1:trimlen(file))
 140  format('Reading ANTMOD File ',i2,' Ver ',F4.1,1x,a)
      antmod_file(num_antmods) = file

****  OK: See what comes next
      eof = .false.
      antprn = 0
      read_dph = .false.

      do i = 1, max_site
        antsit(i) = 0
        antassigned(i) = 'NONE'
        fullfnd(i) = .false.
        partfnd(i) = .false.
      end do

      do while ( .not.eof )
         read(un,'(a)',iostat=ierr) line
         if ( ierr.eq. -1 ) eof = .true.
         if( ierr.ne.0 .and. ierr.ne.-1 ) then
            call report_error('IOSTAT',ierr,'read',file,0,'READ_ANTMOD')
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
     .                'READ_ANTMOD/TYPE')
            else if( anttypx(1:7).eq.'GLONASS' ) then
                antprn = -1
            else
                antprn = 0
*               See if we can match to a site.  For each matching
*               site save site number, else set to 0.  If none
*               match, then set antprn = -1 so that entries are
*               ignored. Make radome type be none if not present. 
                if ( anttypx(17:20).eq.'    ' ) anttypx(17:20) = 'NONE'
                site_found = .false.
                do i = 1, num_site
                   if ( anttypx.eq. ant_name(i) ) then
                      antsit(i) = i
                      site_found = .true.
                      fullfnd(i) = .true.
                      antassigned(i) = anttypx
                   else
*                     See if we can match with no-radome type and the
*                     full name has not been matched.
                      if ( anttypx(1:16).eq.ant_name(i)(1:16) .and.
     .                     anttypx(17:20).eq.'NONE' .and.
     .                     .not.fullfnd(i)   ) then
                          antsit(i) = i
                          partfnd(i) = .true.
                          antassigned(i) = anttypx 
                      else 
*                         No match
                          antsit(i) = 0
                      endif
                   endif
                end do
                if( .not. site_found ) antprn = -1

            endif
            dazi = 0.0d0
          endif
*
          if( line(61:64).eq. 'DAZI' .and. antprn.ge.0 ) then
              read(line,'(2x,f6.1)',iostat=jerr) dazi 
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'READ_ANTMOD/DAZI')
              if( antprn.eq.0 ) then
*                Set range of values for site dependent azimuth
                 do i = 1, num_site
                    if( antsit(i).eq.i ) then
                       sit_dzn(1,2,i) = 0.0d0
                       sit_dzn(2,2,i) = 360.d0
                       sit_dzn(3,2,i) = dazi
                       if( dazi.gt.0 ) then
                           num_sit_dph(2,i) = nint(360.d0/dazi) + 1
                       else
                           num_sit_dph(2,i) = 1  ! Just one value, no
*                                   ! azimuth dependence
                       endif
                    end if
                 end do
              endif

          endif

          if( line(61:64).eq. 'ZEN1' .and. antprn.ge.0 ) then
              read(line,'(2x,3f6.1)',iostat=jerr) dzen 
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'READ_ANTMOD/DZEN')
              if( antprn.gt.0 ) then
                  do j = 1,3
                     svs_dna(j,1,antprn) = dzen(j)
                  end do
                  num_svs_dph(1,antprn) = 
     .                        nint((dzen(2)-dzen(1))/dzen(3))+1
                  num_svs_dph(2,antprn) = 1
              else
                  do i = 1, num_site
                     if( antsit(i).eq.i ) then
                        do j = 1,3
                          sit_dzn(j,1,i) = dzen(j)
                        end do
                        num_sit_dph(1,i) = 
     .                        nint((dzen(2)-dzen(1))/dzen(3))+1

                     end if
                  end do
              endif

          endif
          if( line(61:69).eq. '# OF FREQ' .and. antprn.ge.0 ) then
              read(line,'(i6)',iostat=jerr) nfrq 
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'READ_ANTMOD/#_OF_FREQ')
          endif
          if( line(61:70).eq. 'VALID FROM' .and. antprn.ge.0 ) then
              read(line,'(5i6,f13.7)',iostat=jerr) date, sec_tag
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'READ_ANTMOD/VALID FROM')
*             Check that our data is after this date
              call ymdhms_to_mjd(date, sec_tag, mjd_vfrom)
              if ( sp3_time(1).lt.mjd_vfrom ) then
*                 Our data is before this time, so don't use
                  antprn = -1
              end if
          endif
          if( line(61:71).eq. 'VALID UNTIL' .and. antprn.ge.0 ) then
              read(line,'(5i6,f13.7)',iostat=jerr) date, sec_tag
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'READ_ANTMOD/VALID UNTIL')
*             Check that our data is before this date
              call ymdhms_to_mjd(date, sec_tag, mjd_vuntil)
              if ( sp3_time(1).gt.mjd_vuntil ) then
*                 Our data is before this time, so don't use
                  antprn = -1
              end if
          endif

*
          if( line(61:73).eq. 'START OF FREQ' .and. antprn.ge.0 ) then
              read(line,'(3x,1x,i2)',iostat=jerr) ifrq
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'READ_ANTMOD/START OF FREQ')
              read_dph = .true.
          endif

          if( line(61:71).eq. 'END OF FREQ' .and. antprn.ge.0 ) then
              read(line,'(3x,1x,i2)',iostat=jerr) ifrq
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'READ_ANTMOD/END OF FREQ')
              read_dph = .false.
          endif

          if( line(61:77).eq. 'NORTH / EAST / UP' .and. 
     .                                            antprn.ge.0 ) then
              read(line,'(3f10.2)',iostat=jerr) dPos
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'READ_ANTMOD/NORTH / EAST / UP')
*             Save values in correct slots
              if ( antprn.gt.0 ) then
*                This is satellite so save values in slot
                 if( ifrq.le.2 ) then
                    do j = 1, 3
                       svs_L12(j,ifrq,antprn) = dPos(j)/1000.d0
                    enddo
                 else
                    do j = 1, 3
                       svs_L12(j,1,antprn) = dPos(j)/1000.d0
                       svs_L12(j,2,antprn) = dPos(j)/1000.d0
                    enddo
                 end if
              else  ! Save for sites
                 do i = 1, num_site
                    if( ifrq.le.2 .and. antsit(i).eq.i ) then
                       do j = 1, 3
                          sit_L12(j,ifrq,i) = dPos(j)/1000.d0
                       enddo
                    elseif ( antsit(i).eq.i ) then
                       do j = 1, 3
                          sit_L12(j,1,i) = dPos(j)/1000.d0
                          sit_L12(j,2,i) = dPos(j)/1000.d0
                       enddo
                    end if
                 end do
              end if
          endif

*****     See if we are in READ_DPH mode.  Just read lines
          if( read_dph .and. antprn.ge.0 .and.
     .        line(61:73).ne. 'START OF FREQ' .and.
     .        line(61:77).ne. 'NORTH / EAST / UP' ) then
*             Get values from line
              if( antprn.gt.0 ) then
                 nzen = num_svs_dph(1,antprn)
                 naz  = num_svs_dph(2,antprn)
              else
                 do i = 1,num_site
                    if( i.eq.antsit(i) ) then
                       nzen = num_sit_dph(1,i)
                       naz  = num_sit_dph(2,i)
                    endif
                 enddo
              endif
              read(line,'(a8,91f8.2)',iostat=jerr) 
     .                  noazi, (vals(i),i=1,nzen)
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                'READ_ANTMOD/PHASE Values')

*             See what we have
              if( noazi.eq.'   NOAZI' .and. naz.eq.1 ) then
*                 This is the azimuth averaged values.  If dazi is zero
*                 use these values.
                  if( antprn.gt.0 ) then
*                    Satellite values
                     if ( ifrq.le.2 ) then
                        do j = 1, nzen
                            svs_dphs(j,1,ifrq,antprn) = vals(j)/1000.d0
                        enddo
                     else
                        do j = 1, nzen
                            svs_dphs(j,1,1,antprn) = vals(j)/1000.d0
                            svs_dphs(j,1,2,antprn) = vals(j)/1000.d0
                        enddo
                     endif
                     if( debug(9).gt.0 ) 
     .               write(*,200) antprn, nzen, ifrq
     .,                   (svs_L12(j,1,antprn), j=1,3)
     .,                   (svs_L12(j,2,antprn), j=1,3)
     .,                   (svs_dphs(j,1,ifrq,antprn)*1000.d0,j=1,nzen)
 200                 format('SVS PRN ',i2.2,' NZEN ',i3, ' IFRQ ',i3,/
     .,                     'SVS L1 ',3F8.4, ' L2 ',3F8.4,/
     .,                     'DPH L1 ',15F8.2) 
                  else ! Do the site values
*                    Site value and no-azimuth dependence and so use
*                    these values
                     do i = 1, num_site
                        if( antsit(i).eq.i ) then 
                           if ( ifrq.le.2 ) then
                              do j = 1, nzen
                                  sit_dphs(j,1,ifrq,i) = vals(j)/1000.d0
                              enddo
                           else
                              do j = 1, nzen
                                  sit_dphs(j,1,1,i) = vals(j)/1000.d0
                                  sit_dphs(j,1,2,i) = vals(j)/1000.d0
                              enddo
                           endif
                           if( debug(9).gt.0 ) 
     .                     write(*,220) i, nzen, naz, ifrq
     .,                        (sit_L12(j,1,i), j=1,3)
     .,                        (sit_L12(j,2,i), j=1,3)
     .,                        (sit_dphs(j,1,ifrq,i)*1000.d0,j=1,nzen)
 220                       format('SITE NUM ',i2.2,' NZEN ',i3, 
     .                          ' NAZ ',i3,' IFRQ ',i3,/
     .,                         'SIT L1 ',3F8.4, ' L2 ',3F8.4,/
     .,                         'DPH L1 ',91F8.2) 

                        endif
                     enddo
                  end if
*             OK: Process azimuth dependent value (only for sites)
              elseif( noazi.ne.'   NOAZI' .and. naz.gt. 1 ) then
*                 Get the value of the azimuth and compute index
                  read(noazi,*,iostat=jerr) azval
                  call report_error('IOSTAT',jerr,'decod',noazi,0,
     .                  'read_antmod/azimuth')
                  if( jerr.eq.0 ) then
                      inaz = nint(azval/(360.d0/(naz-1)))+1
*                     Save the values
                      do i = 1, num_site
                         if( antsit(i).eq.i ) then 
                            if ( ifrq.le.2 ) then
                               do j = 1, nzen
                                   sit_dphs(j,inaz,ifrq,i) = 
     .                                                 vals(j)/1000.d0
                               enddo
                            else
                               do j = 1, nzen
                                   sit_dphs(j,inaz,1,i) = 
     .                                                 vals(j)/1000.d0
                                   sit_dphs(j,inaz,2,i) = 
     .                                                 vals(j)/1000.d0
                               enddo
                            endif
                            if( inaz.eq.1 ) then
                               if( debug(9).gt.0 ) 
     .                         write(*,240) i, nzen, naz, ifrq
     .,                             (sit_L12(j,1,i), j=1,3)
     .,                             (sit_L12(j,2,i), j=1,3)
                               if( debug(9).gt.0 ) 
     .                         write(*,260) inaz, azval
     .,                             (sit_dphs(j,1,ifrq,i)*1000.d0,
     .                                                      j=1,nzen)
 240                           format('SITE NUM ',i2.2,' NZEN ',i3, 
     .                               ' NAZ ',i3,' IFRQ ',i3,/
     .,                              'SIT L1 ',3F8.4, ' L2 ',3F8.4)
                            else
                                if( debug(9).gt.0 ) 
     .                          write(*,260) inaz, azval
     .,                           (sit_dphs(j,inaz,ifrq,i)*1000.d0,
     .                                                      j=1,nzen)
 260                            format(I3,F6.1,2x,91F8.2)
                            endif  

                         endif

                     enddo
                  endif
              end if
          end if

*     Done reading file
          ENDIF  ! Not EOF
      end do
      close(un)

*     Report on antennas found
      do i = 1, num_site
         if( fullfnd(i) ) then
            write(*,320) site_names(i), antassigned(i), 
     .                   file(1:trimlen(file))
 320        format('Matched antenna and radome at ',a4,' to ',a20,
     .             ' in antmod_file ',a)
         elseif( partfnd(i) ) then
            write(*,340) site_names(i), antassigned(i), 
     .                   file(1:trimlen(file))
 340        format('Matched antenna only       at ',a4,' to ',a20,
     .             ' in antmod_file ',a)
         endif
      enddo
 
      return
      end



CTITLE READ_DCBS

      subroutine read_dcbsRT(file, obs_mjd )

      implicit none

*     Routine to read and save the DCB values for current time range.
*
      include '../includes/const_param.h'
      include 'trackRT.h'

* PASSED  VARIABLES
      real*8 obs_mjd      ! Observation time
      character*(*) file  ! Name of file. If open OK, added to
*           list of files

* LOCAL VARIABLES 
      integer*4 i          ! Loop counters
     .,         ierr, jerr ! IOSTAT error
     .,         eindx      ! Index of Epoch on line
     .,         svn, pn       ! SVN and PRN number
     .,         un         ! Unit number for file 

      integer*4 yr, doy    ! Year and DOY
     .,         trimlen    ! Length of string
     .,         dats(5), date(5)   ! start and end date

      real*4    dcb_ver    ! Version of DCB file

      real*8 sec           ! Seconds tag
     .,      jd, mjd       ! JD and MJD of YR and DOY
     .,      dcb_ns        ! DCB in ns
     .,      mjds, mjde    ! Start and end MJD

      character*80 line    ! Line read from file (Epoch year day or
                           ! PRN dcb (ns) sig (ns)

***** Initialize the values
      do i = 1, num_prn
         dcbs(i) = 0
      end do
      dcb_mjd = -1 

****  see if name
      if ( trimlen(file).eq.0 ) RETURN

****  Try to open file
      un = 83
      open(un,file=file,iostat=ierr,status='old')
      call report_error('IOSTAT',ierr,'open',file,0,'read_dcbsRT')
      if ( ierr.ne.0 ) RETURN

****  See which version we have
      read(un,'(a)',iostat=ierr) line
      if( line(1:9).eq.'* dcb.dat' ) then
         read(line(18:22),*,iostat=jerr) dcb_ver
         write(*,120) trim(file), dcb_ver
 120     format('Found dcb.dat file ',a,' Verion ',F6.1)
      endif

****  Now read the file
      if( dcb_ver.eq. 0.0 ) then 
         do while ( ierr.eq.0 )
            read(un,'(a)',iostat=ierr) line
            if( ierr.eq.0 .and. line(1:1).ne.'*' .and.
     .          trimlen(line).gt.0 ) then
                call casefold(line)
                eindx = index(line,'EPOCH')
                if( eindx.ge.1 ) then   ! Epoch line
                    read(line(7:),*,iostat=jerr) yr,doy
                    call report_error('IOSTAT',jerr,'decode',line,0,
     .                   'read_dcbsRT/Epoch')
                    if( jerr.eq.0 ) then
                       call yds_to_jd(yr,doy,sec, jd)
                       mjd = jd - 2400000.5d0
                    else
                       mjd = 0.d0
                    end if
                    if( mjd.gt.obs_mjd ) then
                       exit  !  Exit loop
                    end if
                else     ! Read the DCB line
                    read(line,*,iostat=jerr) pn, dcb_ns
                    call report_error('IOSTAT',jerr,'decode',line,0,
     .                   'read_dcbsRT/dcb')
                    if( jerr.eq.0 ) then
                        if( pn.le.max_prn ) then
                            dcbs(pn) = dcb_ns*1.d-9*vel_light
                        else
                            call report_stat('warning','trackRT',
     .                         'read_dcbsRT',line,
     .                         'pn number too large',pn)
                        endif
                    end if
                    dcb_mjd = mjd   ! Actual MJD for entry

                 end if
             end if
          end do
       else
          dats(2) = 1
          date(2) = 1
          do while ( ierr.eq.0 )
            read(un,'(a)',iostat=ierr) line
            if( ierr.eq.0 .and. line(2:2).eq.'G' ) then
                read(line,220, iostat=jerr ) svn, pn, dats(1),
     .               dats(3:5), date(1), date(3:5),  dcb_ns
 220            format(5x,I3,1x,I3,2x,I4,1x,I3,1x,I2,1x,I2,
     .                             2x, I4,1x,I3,1x,I2,1x,I2,1x,F10.3)
                if( jerr.ne.0 )
     .          call report_stat('warning','trackRT','read_dcbsRT',line,
     .                 'IOSTAT error',jerr)

                call ymdhms_to_mjd(dats,sec, mjds)
                call ymdhms_to_mjd(date,sec, mjde)
                if( obs_mjd.ge.mjds .and.  obs_mjd.le.mjde ) then
                     dcbs(pn) = dcb_ns*1.d-9*vel_light 
                endif
            endif
          enddo
        endif



****   Thats all
       call report_dcbs(6)
       close(un)
       return
       end







                 
      


