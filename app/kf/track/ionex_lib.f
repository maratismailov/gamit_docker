CTITLE READ_IONEX

      subroutine read_ionex

*     Routine to read an IONEX file and save values in memory
*     for later evaluation by theory.  If there are problems 
*     with the file, use_ionex is set false.

      implicit none

      include 'track_com.h'

* LOCAL VARIABLES
      integer*4 ierr  ! IOSTAT error
     .,         i,j,k ! Loop and index counters
     .,         nline ! Number of lines to read (5)
     .,         irow  ! Row in matrix (row 1 = 87.5, row 71=-87.5)
     .,         ncol  ! Number of columns to read
      integer*4 date(5)  ! YMDHMS
      real*8 sectag      ! Seconds tag

      real*8 lat, lon1, lon2   ! Lat and long range of map line,
      logical OK_to_read       ! Set true when OK to read ION block
      character*128 line       ! Line read from file.


****  Try to open file
      open(99,file=ionex_file, status='old',iostat=ierr)
      call report_error('IOSTAT',ierr,'open',ionex_file,0,
     .     'read_ionex')
      if( ierr.ne.0 ) then
          use_ionex = .false.
          RETURN
      endif

***** OK, start reading the file looking for START OF TEC MAP
      k = 0  ! Number of times found
      ionex_dlat = 2.5d0  ! 2.5 deg spacing
      OK_to_read = .false.

      do while ( ierr.eq. 0 )
         read(99,'(a)',iostat=ierr ) line
*        See if start of tec map
         if( ierr.eq.0 ) then 
            call casefold(line)
            if( index(line,'START OF TEC MAP').eq.61 ) then
*               Found start of map, now start reading in the values
                read(line,*) k   ! Map number
                OK_to_read = .true.
            elseif ( index(line,'END OF TEC MAP').eq.61 ) then
*               End of map. so no more reading
                OK_to_read = .false.
            elseif ( index(line,'EPOCH OF CURRENT MAP').eq.61 ) then
*               Read the time of this map
                read(line,*) date, sectag
                call ymdhms_to_mjd(date, sectag,ionex_times(k))
            elseif( index(line,'LAT/LON1/LON2/DLON/H').eq.61 .and.
     .                                            OK_to_read ) then
*               Read the latiude, long and height
                read(line,220,iostat=ierr) lat, lon1, lon2, ionex_dlng, 
     .                                     ionex_ht
 220            format(2x,5f6.1)
                call report_error('IOSTAT',ierr,'decod',line,0,
     .               'read_ionex/lat/lng')
                if ( ierr.ne.0 ) use_ionex = .false.
*               Now read in the block of data
                nline = int(360/ionex_dlng/16) + 1
                irow = nint((87.5d0-lat)/ionex_dlat) + 1
                do i = 1, nline
                   if( i.eq.nline ) then
                       ncol = 9
                   else
                       ncol = 16
                   endif
                   read(99,240) (ionex_tec(irow,(i-1)*16+j,k),j=1,ncol)
 240               format(16I5)
                enddo
             end if
          end if
      end do

      ionex_num = k  ! Number of times read in

****  Thats all
      close(99)

****  Debug DUMP TEC for Matlab
C     do i = 1, ionex_num
C        do j = 1,71
C           write(500+i,320) (ionex_tec(j, k, i),k=1,73)
C320        format(73I6)
C        end do
C     end do

      return
      end

CTITLE INTERP_IONEX

      subroutine interp_ionex( ep, ns, ch, site_pos, azd, 
     .                         elevd, mjd, tec_los )

      implicit none

*     Routine to return the TEC value (line-of-site) for a given
*     lat, long (rads) with az elev (rad) at time mjd and returns
*     teqc

      include '../includes/const_param.h'
      include 'track_com.h'

* PASSED VARIABLES
      integer*4 ep    ! Epoch (for saving and getting values)
     .,         ns    ! Site number
     .,         ch    ! Channel number

      real*8 site_pos(3) ! XYZ coordinates (convert to geocentric
                         ! lat and long (rads)
      real*8 azd, elevd  !  Azimuth/Elev (deg)
      real*8 az , elev   !  Azimuth/Elev (rad)
      real*8 mjd       !  MJD of measurement
      real*8 tec_los   !  Line of sight TEC (TECU) 

      real*8 Re        ! Radius of earth (km)

      parameter ( Re = 6371.0 )  

*
* LOCAL VARIABLES
      real*8 latc, lngc  ! Geocenter lat/long of site (rads)
      real*8 lati, lngi, lngr ! Lat/long of sub-ion point and
                         ! longitude mapped to account for motion of sun (deg)
      real*8 beta, gamma ! ZD at ion-point and angle of center of earth
      real*8 th2         ! Co-latitude sub-ion point
      real*8 sdl, cdl, dlng    ! Sin and Cos of difference in longitude
      real*8 dt          ! Time spacing on ion maps (hrs)

      real*8 tecs, tece, tec  ! Spatially interpolate tec at start, end
                         ! and interpolated to correct time 
      real*8 p,q         ! Fractional long and latitude interpolation

      integer*4 it       ! Time index
      integer*4 irl, iru ! Row numbers of latitude at lower and upper points
                         ! (Since lat starts at 87.5 and decreased in table
                         ! iru=irl-1
      integer*4 icl, icu ! Columnes of longitude lower and upper
      integer*4 iep      ! Epoch number

****  See if we already have values
      if( ion_known ) then   ! Just get the values from tables
!        print *,'TEC_LOS SAVE ',ch,ns,ep, tec_los_cse(ch,ns,ep)
         tec_los = tec_los_cse(ch,ns,ep)
         RETURN
      endif

****  We need to calculate: Convert Az, elev
      az = azd * pi/ 180
      elev = elevd * pi/180


****  Get geocentric lat/long
      latc = atan(site_pos(3)/sqrt(site_pos(1)**2+site_pos(2)**2))
      lngc = atan2(site_pos(2),site_pos(1))
      if( lngc.gt.pi ) lngc = lngc - 2*pi

****  Now compute ion 450 km interestion lat/long.  These we convert
*     back to deg to allow easy calc of grid point.
      beta  = asin(Re*cos(elev)/(Re+ionex_ht))
      gamma = pi/2 - beta - elev
      th2   = acos(cos(pi/2-latc)*cos(gamma) + 
     .             sin(pi/2-latc)*sin(gamma)*cos(az))
      sdl   = sin(gamma)*sin(az)/sin(th2)
      cdl   = (cos(gamma)-cos(pi/2-latc)*cos(th2))/
     .                   (sin(pi/2-latc)*sin(th2))
      dlng  = atan2(sdl, cdl)
      lati  = (90 - th2*180/pi)
      lngi  = (lngc + dlng)*180/pi
      if( lngi.lt.-180  ) lngi = lngi + 360
      if( lngi.ge. 180  ) lngi = lngi - 360
    

****  Get time idices we need (will be used to adjust the longitude 
*     of the grid as well).
      dt = (ionex_times(2)-ionex_times(1))*24
      it = int((mjd-ionex_times(1))*24/dt) + 1
      if( it.lt.1 .or. it.ge.ionex_num ) then
          write(*,120) it,mjd, ionex_times(1), dt
 120      format('**ERROR** Time outside IONEX file. Index ',i4,
     .           'MJD ',F8.2,' Start IONEX ',F8.2,' DT ',F6.2,' hrs')
          tec_los = 0
          use_ionex = .false.
          stop 'IONEX file out of time range'
      end if

****  Get grid and spatial interpolation for first time
      lngr = lngi  + (mjd-ionex_times(it))*360
      if( lngr.lt.-180 ) lngr = lngr + 360
      if( lngr.ge. 180 ) lngr = lngr - 360
 
*     Latitude cells
      irl = int((90-lati)/ionex_dlat)+1
      iru = irl - 1
      if( iru.lt.1 ) iru = 1
      if( irl.gt.71 ) irl = 71
*     Longitude cells 
      icl = int((lngr+180)/ionex_dlng) + 1
      icu = icl + 1
*     Now do spatial interpolation
      q = (irl*ionex_dlat-(90-lati))/ionex_dlat
      p = ((lngr+180) - (icl-1)*ionex_dlng)/ionex_dlng
      tecs = (1-p)*(1-q)*ionex_tec(irl,icl,it) +
     .          p *(1-q)*ionex_tec(irl,icu,it) +
     .       (1-p)* q   *ionex_tec(iru,icl,it) +  
     .          p * q   *ionex_tec(iru,icu,it)  

****  Now do spatial interpolation at end time
      lngr = lngi  + (mjd-ionex_times(it+1))*360
      if( lngr.lt.-180 ) lngr = lngr + 360
      if( lngr.ge. 180 ) lngr = lngr - 360

*     Longitude cells (re-interpolate indices)
      icl = int((lngr+180)/ionex_dlng) + 1
      icu = icl + 1
*     Now do spatial interpolation
      q = (irl*ionex_dlat-(90-lati))/ionex_dlat
      p = ((lngr+180) - (icl-1)*ionex_dlng)/ionex_dlng
      tece = (1-p)*(1-q)*ionex_tec(irl,icl,it+1) +
     .          p *(1-q)*ionex_tec(irl,icu,it+1) +
     .       (1-p)* q   *ionex_tec(iru,icl,it+1) +  
     .          p * q   *ionex_tec(iru,icu,it+1)  

****  Now do spatial interpolation at end time
      tec = 24*((ionex_times(it+1)-mjd)*tecs + 
     .          (mjd-ionex_times(it))*tece)/dt

****  Convert to TECU (/10 from integer values read) and map to LOS
      tec_los = tec/cos(beta)/10

****  Save the values
      tec_los_cse(ch,ns,ep) = tec_los

****  Save sub-ion point
      tec_subI_cse(1,ch,ns,ep) = lngi
      tec_subI_cse(2,ch,ns,ep) = lati

***** Test output
      iep = nint((mjd - ref_start)*86400/usr_interval)+1
      if( iep - int(iep/360)*360.eq.1 .and. 1.eq.2 ) 
     .write(*,220) iep, mjd, it,irl, icl, tecs, tece, tec, tec_los,
     .          90-beta*180/pi, elev*180/pi, az*180/pi, p, q, 
     .          lati, lngi, ns, ctop_cse(ch,ns,ep)  
 220  format('IONEX EP ',i6,' MJD ',F14.6,' IT/IRL/ICL ',3I3,
     .       ' tecs, tece, TEC ',4F10.3,' co-Beta Elev Az ',3F10.4,
     .       ' p,q ',2F6.3,' Lat/Long ', 2F10.4,' S/PN ',2i3)

****  Thats all
      return
      end

CTITLE  read_ionlos
      
      subroutine read_ionlos(ns)

      implicit none

*     Routine to read all the TEC data into memory.  The indexing is
*     based on the epoch in the file ION_DELAY_FILE.

      include '../includes/const_param.h'
      include '../includes/xfile_def.h'
      include 'track_com.h'

* PASSED VARIABLES
* ns  -- Station number to read (used to connect to the logical unit also)

      integer*4 ns

* LOCAL VARIABLES
* j    -- Loop counters to find channel
* ch   -- Channel number for the PRN at each epoch
* ierr, jerr   -- status of read statement
* ep    -- Epoch number for the current estimates read from ion_delay_file
* prn -- PRN corresponding to estimate read from ion_delay_file

      integer*4 j, ch, ierr, jerr,trimlen, ep, prn, nsaved, nread

* ion_tec -- the estimated TEC read from ion_delay_file

      real*8 ion_tec

      character*128 line   ! Line read from file 

***** Open delay file: There needs to one file per station or
*     the site number needs to be in the records of the file.

      open(99,FILE=ionlos_file(ns),iostat=ierr,status='old')
      call report_error('IOSTAT',ierr,'open',ionlos_file(ns),
     .       0,'READ_ION_DEL')
      use_ionlos(ns) = .false.

      if( ierr.ne.0 ) RETURN

      nsaved = 0
      nread  = 0

      do while ( ierr.eq.0  )

*         Read only the entries we need here
          read(99,'(a)',iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' .and.
     .        trimlen(line).gt.0 ) then 
              read(line,*,IOSTAT=jerr) ep, prn, ion_tec
              nread = nread + 1
              if( jerr.eq.0 ) then   ! Process line

*                 Find the channel number for this PRN at this epoch and site
                  ch = -1
                  do j = 1,num_chan_se(ns,ep)
                      if( ctop_cse(j,ns,ep).eq.prn ) then
                          ch = j
                          exit
                      endif
                  end do

*                 Save the value of TEC
                  if( ch.gt.0 ) then
                      tec_los_cse(ch,ns,ep) = ion_tec
                      nsaved = nsaved + 1
                  endif
              else
                  call report_error('IOSTAT',jerr,'read',line,0,
     .                 'READ_IONLOS')
              end if
          end if
      end do
      write(*,220) nsaved, nread, site_names(ns), trim(ionlos_file(ns))
 220  format('IONLOS: ',i6,2x,I6,
     .       ' LOS values saved and read for site ',a,' from file ',a)
      close(99)

****  Thats all
      return
      end

