CITTLE READ_IN_PBO

      subroutine read_in_pbo( nf )

      implicit none

*     Routine to read in PBO files and part of the processing.  Existing pos files
*     will be read also but these are not transformed.
* MOD TAH 140903: Added passing file number to save as the name of the site
*     for use in tscomp.  Passing 0 has routine act as before.

      include 'tsfit.h'
      include 'tscon.h'
      include '../includes/const_param.h'

* PASSED 
      integer*4 nf   ! File number for saving by file


      character*16 gsite_full

      integer*4 date(5), ierr, jerr, j, ns,ne, nu
      integer*4 isec    ! Integer seconds in epoch
     .,         unit    ! Unit number for reading file.
      integer*4 trimlen, indx
* MOD TAH 191014: Added rouding of times to 5-minute to handle GYPSY-X
      integer*4 allsec  ! hr/min/sec converted to seconds so that
                        ! nearest 5-minutes can be computed).
      real*8  gmjd, pos_xyz_fin(3),xyz_std(6), unc_geod(3),
     .        pos_neu_fin(3), neu_std(6), sec

      real*8 unc_llu(3)  ! Computed Lat/long/height
      real*8 rot_mat(3,3) ! Rotation from XYZ to NEU (and transpose)

      logical line_read ! Indicates header line have been read (so that 
                        ! earilier formats can be read)
      
      character*5 ts_in_type   ! Solution type read from file.
 
      character*512 line

***   First open the file 
      open(100,file=in_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',in_file,0,'tscon')
      unit = 100

      if( ierr.ne.0 ) RETURN

****  Start with header records
      read(unit,'(a)',iostat=ierr) line
      indx = index(line,'Frame :')
      if( indx.gt.0 ) reference_frame = line(indx+8:)
* MOD TAH 051129: See if version line next
      call report_error('IOSTAT',ierr,'read',line,0,'read_ts')
      read(unit,'(a)',iostat=ierr) line
      if( index(line,'Version').gt.0 ) then
          read(line,105,iostat=ierr) ts_ver_read
 105      format(16x,a)
c105      format('Format Version: ',a)
          call report_error('IOSTAT',ierr,'read',line,0,'read_ts')
          read(unit,'(a)',iostat=ierr) line
      endif

      read(line,110,iostat=ierr) site_name
 110  format(16x,a)
      site_name(5:8) = '_GPS'
      call casefold(site_name)

c110  format('4-character ID: ',a4)
      call report_error('IOSTAT',ierr,'read',line,0,'read_ts')
      read(unit,120,iostat=ierr) gsite_full
 120  format(16x,a16)
c120  format('Station name  : ',a16)
      call report_error('IOSTAT',ierr,'read','Station',0,'read_ts')
      read(unit,'(a)',iostat=ierr) line
c     read(line,130,iostat=ierr) date, isec
c130  format(16x,i4,i2,i2,1x,i2,i2,i2)
c     sec = isec
c130  format('First Epoch   : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
c     call report_error('IOSTAT',ierr,'read',line,0,'read_ts')
c     call ymdhms_to_mjd(date,sec,ts_first)
cFirst Epoch   : 20060521 120000
      read(unit,'(a)',iostat=ierr) line
c     read(unit,140,iostat=ierr) date, isec
c     sec = isec
c140  format(16x, i4,i2,i2,1x,i2,i2,i2)
c     call report_error('IOSTAT',ierr,'read','Last Epoch',0,'read_ts')
c140  format('Last Epoch    : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
c     call ymdhms_to_mjd(date,sec,ts_last)

C**** Do not read release date, since we will generate a new one here
C     read(unit,150,iostat=ierr) date_rel, nint(sec_rel)
      read(unit,'(a)',iostat=ierr) line
C150  format(16x, i4,i2,i2,1x,i2,i2,i2)
c150  format('Release Date  : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
      call report_error('IOSTAT',ierr,'read',line,0,'tssum')

****  We don't need the reference values because we will use the 
*     XYZ coordinates to compute everything.  Replace read values
*     with read line calls.
* MOD TAH 121213: Read the reference coordinates so they can be saved
      read(unit,'(a)',iostat=ierr) line
      read(line,200,iostat=ierr) ref_xyz
 200  format(25x,3F15.5)
c200  format('XYZ Reference position : ',3F15.5)
c     call report_error('IOSTAT',ierr,'read',line,0,'tssum')

c     call XYZ_to_GEOD(rot,ref_xyz, ref_llu) 
c     call loc_to_geod(ref_llu, ref_neu)

C     Use compute version
C      read(unit,220,iostat=ierr) (pi/2-ref_llu(1))*180/pi,
C     .                            ref_llu(2)*180/pi,ref_llu(3)
      read(unit,'(a)',iostat=ierr) line
C220  format(25x,2F16.10,1x,F10.5)
      call report_error('IOSTAT',ierr,'read',line,0,'tssum')
c220  format('NEU Reference position : ',2F16.10,1x,F10.5)

* MOD TAH 080108: See if field descrription lines are present
      read(unit,'(a)',iostat=ierr) line
      line_read = .true.
      if( index(line,'Start Field Description').gt.0 ) then
          do while (index(line,'End Field Description').eq.0 )
             read(unit,'(a)',iostat=ierr) line
             call report_error('IOSTAT',ierr,'read',ts_file,0,
     .                         'tssum/Field Description')
          enddo
*         Now read header in 
          read(unit,'(a)') line
          line_read = .false.
      endif

* MOD TAH 121212: Moved site name match to here from inside data
*     read (not clear why is was embedded in reading entries)
*     See if we can match site name
      indx = 0
      call get_cmd(site_name,gsite_names,num_site,ns,indx)
      if( ns.le.0 ) then
          num_site = num_site + 1
          if( num_site.gt.max_site ) then 
              call report_stat('FATAL','TSCON','read_in_pbo',
     .             '','Too many sites',max_site)
          endif

          num_code = num_site
          ns = num_site
          gsite_names(ns) = site_name
          in_code(ns) = site_name(1:4)
          in_full(ns) = gsite_full
      end if
      save_ref_xyz(:,ns) = ref_xyz(:)
      if( nf.eq.0 ) then
         nu = ns    ! Old default method
      else
         nu = nf    ! Save by file number for tscomp
      endif

****  Finished reading the header records.  Now we start down the
*     data records.

****  Now start reading the full PBO lines
      do while ( ierr.eq.0 )
         if( .not.line_read ) then
             read(unit,'(a)',iostat=ierr) line
* MOD TAH 210309: Explicit test of blanke line added
             if( trimlen(line).lt.240 ) then
                call report_error('SHORTLine',trimlen(line),
     .                            'read',line,0,'read_ts')
                ierr = -1
             end if
         else
             line_read = .false.
         endif

         if( ierr.eq.0 ) then
             read(line,300,iostat=jerr) date, isec, gmjd,
     .          (pos_xyz_fin(j),j=1,3), (xyz_std(j),j=1,6),
     .          (unc_llu(j),j=1,3),
     .          (pos_neu_fin(j),j=1,3), (neu_std(j),j=1,6),
     .           ts_in_type
 300         format(1x,i4,i2,i2,1x,i2,i2,i2,1x,F10.4,
     .           3F15.5,3F9.5,3F7.3,3x,2F16.10,1x,F10.5,3x,
     .           3(F9.5,1x),1x,3F9.5,3F7.3,1x,a5)
             call report_error('IOSTAT',jerr,'read',line,0,'read_ts')

* MOD TAH 191007: Modulo result to nearest 600 seconds
             allsec = date(4)*3600+date(5)*60+isec
             allsec = nint(allsec/600.0)*600
             date(4) = int(allsec/3600)
             date(5) = int((allsec-date(4)*60)/60)
             isec = allsec-date(4)*3600-date(5)*60
             sec = isec+1.d-3
             call ymdhms_to_mjd(date,sec,gmjd)

             if( gmjd.lt.first_mjd ) first_mjd = gmjd
             if( gmjd.gt.last_mjd  ) last_mjd  = gmjd

****         Now get NEU values
             call XYZ_to_GEOD(rot_mat, pos_xyz_fin, unc_geod )
             call loc_to_geod(unc_geod, pos_neu_fin)

****         OK: If read was OK, add to data set.
             if( jerr.eq.0 ) then


*****            OK save this entry
                 num_ent = num_ent + 1
                 if( num_ent.gt.max_ent ) then
                     call report_stat('FATAL','tscon','read_in_pbo','',
     .                   'Too many enties Max ',max_ent)
                 endif
                 ne = num_ent
                 in_ns(ne)  = nu
                 in_cs(ne)  = ns
                 in_mjd(ne) = gmjd
                 in_type(ne) = ts_in_type
                 do j = 1,3
                    in_xyz(j,ne) = pos_xyz_fin(j)
                    in_neu(j,ne) = pos_neu_fin(j)
                    in_llu(j,ne) = unc_llu(j)
                 end do
* MOD TAH 140829: Scale the sigmas
              xyz_std(1:3) = xyz_std(1:3)*sigma_scale
              neu_std(1:3) = neu_std(1:3)*sigma_scale
                 do j = 1,6
                    in_xyz_std(j,ne) = xyz_std(j)
                    in_neu_std(j,ne) = neu_std(j)
                 end do
             end if
*        Now read error
         end if
      end do

****  Tell user were we are
      write(*,310) in_file(1:trimlen(in_file)), num_site, num_code, 
     .             num_ent, sigma_scale 
 310  format('File: ',a,' Sites ',i5,' Codes ',i5,' Entries ',i10, 
     .       ' SigScale ',F10.3)
      return
      end 

CTITLE WRITE_TS

      subroutine write_ts(unit,ns)

      implicit none

      include '../includes/const_param.h'
      include 'tssum.h'
      include 'tscon.h'

      integer*4 ns  ! code site number
      integer*4 unit  ! unit number

      integer*4 ierr, i, j, k, date(5)
      integer*4 trimlen, lr, lm
      real*8 sec, rot(3,3)
      real*8 ts_mag, Cts_mag  ! Magnitude of XYZ/NEU and function to
                                  ! compute it.

***** Scan the ts list and get start and stop times
      ts_first = 1.d6
      ts_last  = 0.d0

      call remove_ejmp(-ns) 
      do i = 1, num_ts
*        For tscon we can set start and end times to account for
*        these values.
         if( tsprog(1:5).eq.'tscon' ) then
            if( ts_mjd(i).lt.ts_first .and. 
     .          ts_mjd(i).ge.first_mjd ) ts_first = ts_mjd(i)
            if( ts_mjd(i).gt.ts_last .and.
     .          ts_mjd(i).le.last_mjd ) ts_last = ts_mjd(i)
         else
            if( ts_mjd(i).lt.ts_first ) ts_first = ts_mjd(i)
            if( ts_mjd(i).gt.ts_last  ) ts_last = ts_mjd(i)
         end if
      end do

      if( tsprog(1:5).ne.'tscon' ) then
          first_mjd = ts_first
          last_mjd  = ts_last
      endif

***** Write out the header lines
      if( trimlen(reference_frame).eq.0 ) 
     .    reference_frame = runstring_refframe

      lr = trimlen(reference_frame) 
      if( lr.gt.0 ) then
         write(unit,100,iostat=ierr) reference_frame(1:lr)
 100     format('PBO Station Position Time Series.',
     .          ' Reference Frame : ',a)
      else
         write(unit,102,iostat=ierr)
 102     format('PBO Station Position Time Series')
      endif
* MOD TAH 051129: Write out version number
      write(unit,105,iostat=ierr) ts_ver
 105  format('Format Version: ',a)

      write(unit,110,iostat=ierr) in_code(ns)
 110  format('4-character ID: ',a4)
      write(unit,120,iostat=ierr) in_full(ns)
 120  format('Station name  : ',a16)
      call jd_to_ymdhms(ts_first+1d-6, date,sec)
      write(unit,130,iostat=ierr) date, nint(sec)
 130  format('First Epoch   : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
      call jd_to_ymdhms(ts_last+1d-6, date,sec)
      write(unit,140,iostat=ierr) date, nint(sec)
 140  format('Last Epoch    : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
      write(unit,150,iostat=ierr) date_rel, nint(sec_rel)
 150  format('Release Date  : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)

      call XYZ_to_GEOD(rot,ref_xyz, ref_llu) 
      call loc_to_geod(ref_llu, ref_neu)

      if( lr.gt.0 ) then
         lm = 1
         do while ( lm.lt.lr .and. 
     .              reference_frame(lm:lm).ne.' ' )
            lm = lm + 1
         end do
         if( reference_frame(lm:lm).eq.' ' ) lm = lm - 1
         write(unit,200,iostat=ierr) ref_xyz, reference_frame(1:lm)
 200     format('XYZ Reference position : ',3F15.5,' (',a,')')
         write(unit,220,iostat=ierr) (pi/2-ref_llu(1))*180/pi,
     .         ref_llu(2)*180/pi,ref_llu(3),
     .         reference_frame(1:lm)
 220     format('NEU Reference position : ',2F16.10,1x,F10.5,
     .          ' (',a,'/WGS84)')
      else
         write(unit,240,iostat=ierr) ref_xyz
 240     format('XYZ Reference position : ',3F15.5,' (',a,')')
         write(unit,260,iostat=ierr) (pi/2-ref_llu(1))*180/pi,
     .                            ref_llu(2)*180/pi,ref_llu(3)
 260     format('NEU Reference position : ',2F16.10,1x,F10.5)
      endif

      if( 1.eq.1 ) call tsdesc( unit )   ! Write full headers


****  OK, Now write out the time series
      do i = 1, num_ts
          call jd_to_ymdhms(ts_mjd(i)+1.d-6,date,sec)
*
*         Check size
          ts_mag = Cts_mag(ts_neu(1,i),ref_neu)
 
*         Only out put values inside time range
          if( ts_mjd(i).ge. first_mjd .and. 
     .        ts_mjd(i).le. last_mjd+0.25d-5 ) then 
             if( ts_mag.lt.99.999d0 ) then 
*               Use Standard format                  
                write(unit,300) date, nint(sec), ts_mjd(i),
     .              (ts_xyz(j,i),j=1,3), (ts_xyz_std(j,i),j=1,6),
     .              (ts_llu(j,i),j=1,3),
     .              (ts_neu(j,i)-ref_neu(j),j=1,3), 
     .              (ts_neu_std(j,i),j=1,6),
     .               ts_type(i)
 300            format(1x,i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2,1x,F10.4,
     .              3F15.5,3F9.5,3F7.3,3x,2F16.10,1x,F10.5,3x,
     .              3(F9.5,1x),1x,3F9.5,3F7.3,1x,a) 
             else
                write(unit,310) date, nint(sec), ts_mjd(i),
     .              (ts_xyz(j,i),j=1,3), (ts_xyz_std(j,i),j=1,6),
     .              (ts_llu(j,i),j=1,3),
     .              (ts_neu(j,i)-ref_neu(j),j=1,3), 
     .              (ts_neu_std(j,i),j=1,6),
     .               ts_type(i)
 310            format(1x,i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2,1x,F10.4,
     .              3F15.5,3F9.5,3F7.3,3x,2F16.10,1x,F10.5,3x,
     .              3(F9.3,1x),1x,3F9.5,3F7.3,1x,a) 
             endif
          end if
 

      end do

****  Thats all 
      return

      end

CTITLE CTS_MAG

      real*8 function CTS_mag(NEU, Ref)

      implicit none

*     routine to compute magnitude of max element
      real*8 NEU(3), Ref(3)
      real*8 dNEU(3)
      integer*4 i

      CTS_mag = 0.d0
      do i = 1,3
         dneu(i) = NEU(i)-REF(i)
         if( abs(dneu(i)).gt.CTS_mag ) CTS_mag = abs(dneu(i))
      end do

      return
      end

   
