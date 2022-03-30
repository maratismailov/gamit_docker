      program hfupd

      implicit none

*     This program will check the contents of a hfile against the
*     receiver/antenna/eccentricity info in the igs.snx file and
*     adjust the hfile entries as needed (including changing the 
*     station position as need be)

      include '../includes/kalman_param.h' 
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
      include 'hfupd.h'

* LOCAL VARIABLES
* nr -- Number of runstring entries at end of reading base list
* len_run -- Length of runstring
* kerr, ierr    -- Error flag while reading hfiles.
* rcpar   -- Reads runstring 

      integer*4 nr, len_run, kerr, ierr, rcpar

* line -- line read from sinex/station.info file
      character*128 line 

* end  -- Logical to indicate that we have processed all entries in runstring
* old_stinf -- Set true if old station.info format

      logical end, old_stinf                    

****  Initialize the program variables
      call init_hfupd  

* MOD RWK 150213: We do not know what GNSS is used, and it must be only one per h-file,
*                 so restrict to GPS for now
      call report_stat('WARNING','HFUPD','hfupd ',' ' 
     .,'Frequencies for updating eccentries assumed to be GPS L1 L2 ',0)

****  Decode the runstring for the program
      call decode_snx_site_run(nr)

****  Now open and read the igs.snx file (name optionally passed by
*     user)
      
      if( check_snx ) then
          open(100,file=igs_snx_file, iostat=ierr,status='old')
          call report_error('IOSTAT',ierr,'open',igs_snx_file,0,'HFUPD')
          if( ierr.eq.0 ) then
*            Read the first line to see type of file
             read(100,'(a)',iostat=ierr) line
             if( line(1:5).eq.'%=SNX' ) then
                 call read_igs_snx
             else 
**                rwk 100219:  old-format station.info no longer supported
**                 call check_oldstnfo( 100, old_stinf )
**                 if( old_stinf ) then
**                    Read the first line again since this is asssumed
**                    to have been read by read_stinfo
**                     read(100,'(a)',iostat=ierr) line
**                     call read_stinfo
**                 else
                     call read_stinfo_newf
**                 end if
             end if
          end if
      end if
      
      if ( report_snx ) then
          call snx_rep
      end if

      if( edt_hf ) then
          call read_edit
      end if

*     Loop over the hfile names that have been passed
      end = .false.
      do while ( .not. end )

*        Get the next hfile name
         nr = nr + 1
         len_run = rcpar(nr, hfile )
         if( len_run.gt.0 ) then
             write(*,150) hfile(1:len_run)
 150         format(/,70('+'),/,'Processing ',a)
             call fmpopen(cglb_dcb,kerr,hfile,'rwo',0)  
             call report_error('FmpOpen',kerr,'open', hfile,
     .                         0,'SNX_SITE_INFO')

             if( kerr.eq.0 ) then
                 call rw_glb_header('R', kerr)

*                Check if header OK.
                 if( kerr.ne.0 ) then
                    if( kerr.eq.2001 ) then
                        write(*,*) '**WRONG BYTE ORDER** Run swaph'
                    end if
                    call report_error('RW_GLB_HEADER',kerr,'read',
     .                  hfile,0,'SNX_SITE_INFO')
                 end if
             end if

****         Contunue processing only of error OK
             if( kerr.eq.0 ) then

****              Set update needed 
                  upd_needed = .false.

****              Get the information blocks that we need
                  call read_hfinf

*                 See if we are editing the file
                  if( edt_hf ) then
                      upd_needed = .true.
                      call apply_edit
                  end if

*                 Check the file contents against the sinex/station.info
*                 values
                  if( check_snx ) then 
                      upd_needed = .true.
                      call check_sinf( kerr )
                  endif

*                 See if CWU time updates neede
                  if( cwu_upd ) then
                      call check_cwu( kerr )
                  endif  

* MOD TAH 1911226: Check the satellite Name/PRN/SVS values (upd_needed
*                 returned through common).
                  call check_svinf
                  
*                 See if we need to update the solution and headers as well
                  if( upd_hf .and. upd_needed ) then
                      call upd_hf_soln( kerr )
                      call upd_hf_hdr ( kerr )
                  endif
 
                  if( .not.upd_hf .and. upd_needed ) then
                      call report_stat('WARNING','HFUPD','Update Need',
     .                      hfile, '-u option not set',0)
                  endif

             endif
         else
*            Runstring empty so we are done
             end = .true.
         end if
      end do

****  Thats all
      end

CTITLE DECODE_SNX_SITE

      subroutine  decode_snx_site_run(nr)

      implicit none

*     Routine to decode the runstring for hfupd

      include 'hfupd.h'

* PASSED VARIABLES

* nr  -- Current number of runstring entry.  Returns pointing to
*        first hfile in list

      integer*4 nr

* LOCAL VARIABLES

* len_run  -- Length of runstring
* rcpar    -- Reads runstring

      integer*4 len_run, rcpar

* runstring -- Entry returned

      character*128 runstring

* done      -- Logical to indicate that we are still reading
*     options
      logical done

****  Get the first entry and see of -s option
      nr = 1
      done = .false.

*     Loop over the possible entries in the runstring
      do while ( .not.done )

         len_run = rcpar(nr, runstring)
         if( len_run.eq.0 .and. nr.eq. 1 ) then
              call proper_runstring('hfupd.hlp','hfupd',1)
          else

*             See if sinex header file name passed
              if( runstring(1:2).eq.'-s' .or.
     .            runstring(1:2).eq.'-S' ) then
                  nr = nr + 1
                  len_run = rcpar(nr, igs_snx_file)
                  nr = nr + 1
                  check_snx = .true.
              end if

*             See if we just want to report sinex file (-rsnx)
              if( runstring(1:3).eq.'-rs' .or.
     .            runstring(1:3).eq.'-RS'       ) then
                  nr = nr + 1
                  report_snx = .true.
              endif

*             See if we just want to report sinex file (-rec)
              if( runstring(1:3).eq.'-re' .or.
     .            runstring(1:3).eq.'-RE'       ) then
                  nr = nr + 1
                  report_rec = .true.
              endif

*             See if want to update hfiles (Force reporting of
*             differences as well).
              if( runstring(1:2).eq.'-u' .or.
     .            runstring(1:2).eq.'-U'       ) then
                  upd_hf = .true.
                  nr = nr + 1
                  report_diff = .true.
              endif

*             See if we just want to report just difference
              if( runstring(1:2).eq.'-d' .or.
     .            runstring(1:2).eq.'-D'       ) then
                  nr = nr + 1
                  report_diff = .true.
              endif

*             See if want to edit hfile (i.e. remove sites)
              if( runstring(1:2).eq.'-e' .or.
     .            runstring(1:2).eq.'-E' ) then
                  nr = nr + 1
                  len_run = rcpar(nr, edit_file)
                  nr = nr + 1
                  edt_hf = .true.
              end if 

*             See if we want to apply the pole tide
              if( runstring(1:2).eq.'-p' .or.
     .            runstring(1:2).eq.'-P'       ) then
                  nr = nr + 1
                  app_ptide = .true.
*                 Check the next argument to see if 
*                 pmu_file name given or HF used to denote 
*                 use hfile
                  len_run = rcpar(nr, runstring)
                  pmu_file = runstring
                  call casefold(runstring)
                  if( runstring(1:2).eq.'HF' ) then
                     pmu_file = ' '
                  end if
                  nr = nr + 1
              endif

*             See if -honly update option passed
              if( runstring(1:2).eq.'-h' .or.
     .            runstring(1:2).eq.'-H'      ) then
                  nr = nr + 1
                  len_run = rcpar(nr, runstring)
                  nr = nr + 1
                  call casefold(runstring)
                  if( index(runstring,'ANT').gt.0 ) 
     .                            call sbit(honly_opt,1,1)
                  if( index(runstring,'PTD').gt.0 ) 
     .                            call sbit(honly_opt,2,1)
              end if

*             See if we want to update CWU file times
*             Added 091009.
              if( runstring(1:2).eq.'-c' .or.
     .            runstring(1:2).eq.'-C'       ) then
                  nr = nr + 1
                  cwu_upd = .true.
              end if

* MOD TAH 130201: New options
              if( runstring(1:2).eq.'-f' .or.
     .            runstring(1:2).eq.'-F'       ) then
                  nr = nr + 1
                  force_upd = .true.
              end if
              if( runstring(1:2).eq.'-i' .or.
     .            runstring(1:2).eq.'-I'       ) then
                  nr = nr + 1
                  ps_ignore = .true.
              end if


*             See if first character is not an option:
              if( runstring(1:1).ne.'-' ) done = .true.
              
          end if
      end do

****  Set the number of runstring entries back by one so that we will
*     read the next hfile name
      nr = nr - 1

****  Report options selected
      call report_opts

****  Thats all
      return
      end

CTITLE INIT_HFUPD

      subroutine init_hfupd

      implicit none

*     Routine to initialize the logicals and other constants
*     in hfupd

      include 'hfupd.h'

*
*     Sequence over the logicals
      check_snx   = .false.
      report_snx  = .false.
      report_rec  = .false.
      report_diff = .false.
      upd_hf      = .false.
      honly_upd   = .false.
      cwu_upd     = .false.
      force_upd   = .false.
      ps_ignore   = .false.

      honly_opt   = 0
      num_hfu_edits = 0
      sx_num_upd  = 0                  
  
*     File names
      igs_snx_file = ' '
      edit_file    = ' '
      pmu_file     = ' ' 

****  Thats all
      return 
      end

CTITLE REPORT_OPTS

      subroutine report_opts

      implicit none

*     Routine to the options selected for this run of hfupd

      include 'hfupd.h'

* LOCAL VARIABLES
* ---------------
* trimlen -- length of string

      integer*4 trimlen

****  Report the options
      write(*,120)
 120  format(/,'HFUPD Selected options',/,
     .         '----------------------')
      if( check_snx ) then
         write(*,200) igs_snx_file(1:max(1,trimlen(igs_snx_file)))
 200     format('Selected SINEX FILE Header:',t40,a)
      endif
      if( report_snx ) then
         write(*,220)
 220     format('Contents of sinex header will be reported')
      end if
      if( report_diff ) then
          write(*,240)
 240      format('Differences between hfiles will be reported')
      endif
      if( edt_hf ) then
          write(*,260) edit_file(1:max(1,trimlen(edit_file)))
 260      format('H-files will be edited using :',t40,a)
      endif
      if( app_ptide ) then
          if( trimlen(pmu_file).eq.0 ) then
             write(*,280)
 280         format('H-files will have pole-tide applied using ',
     .              'hfile pmu values')
          else
             write(*,285) pmu_file(1:trimlen(pmu_file))
 285         format('H-files will have pole-tide applied using ',a)
          end if
      end if
      if( honly_upd ) then
          write(*,290)
 290      format('Only header information will be updated')
      endif


      if( upd_hf ) then
          write(*,400)
 400      format(/,'****H-files will be updated****',/,
     .             '-------------------------------')
      endif

****  Thats all
      return
      end


CTITLE READ_IGS_SNX

      subroutine read_igs_snx

      implicit none

*     Routine to read the IGS sinex file and extract the site information

      include '../includes/kalman_param.h'
      include 'hfupd.h'

* LOCAL VARIABLES
* --------------- 
* ierr    -- IOSTAT error
* trimlen -- Length of string
      integer*4 ierr, trimlen

* line -- Line read from file

      character*128 line


***** Open the igs sinex file
      ierr = 0
      write(*,120) igs_snx_file(1:trimlen(igs_snx_file))
 120  format(/,'Reading SINEX header file: ',a) 

****  Loop over file reading and saving entries
      do while ( ierr.eq. 0 )
         read(100,'(a)', iostat=ierr) line

*        Check for each of the major blocks
         if( line(1:8).eq. '+SITE/ID' ) then
             call snx_id(ierr)
         end if
         if( line(1:9).eq. '+SITE/REC' ) then
             call snx_rec(ierr)
         end if
         if( line(1:9).eq. '+SITE/ANT' ) then
             call snx_ant(ierr)
         end if
         if( line(1:9).eq. '+SITE/GPS' ) then
             call snx_gps_pc(ierr)
         end if
         if( line(1:9).eq. '+SITE/ECC' ) then
             call snx_ecc(ierr)
         end if
      end do

****  Thats all
      return
      end

CTITLE SNX_ID

      subroutine snx_id(ierr)

      implicit none

*     Routine to read site ID block

      include '../includes/kalman_param.h'
      include 'hfupd.h'

* PASSED VARIABLES
* ----------------
* ierr -- IOSTAT error incase we hit end of file
      integer*4 ierr
       

* LOCAL VARIABLES
* --------------- 
* jerr -- IOSTAT error decoding line
* lng_deg, lng_min -- Longitude deg and mins
* lat_deg, lat_min -- Latitude deg and mins
* ns   -- Local counter for number of sites 
      integer*4 jerr, lng_deg, lng_min, lat_deg, lat_min, ns


* lng_sec, lat_sec -- Seconds part of longitude and latitude
* ht -- Site height
      real*8 lng_sec, lat_sec, ht

* block_end -- Logical to indicate that end of block reached 
      logical block_end

* code -- Code for station
* pt   -- Point character for site
* domes -- Domes number
* name  -- Long name for sites
* type  -- P indicates GPS
* line  -- Line read from file

      character*4 code
      character*2 pt
      character*9 domes
      character*22 name
      character*1 type
      character*128 line


****  Start reading the lines until we hit end of file or end of
*     block
      ns = 0 
      block_end = .false.

      do while ( ierr.eq.0 .and. .not.block_end ) 
         read(100,'(a)', iostat=ierr) line
         if( line(1:8).eq.'-SITE/ID' .or.ierr.ne.0 ) then
             block_end = .true.
             call report_error('IOSTAT',ierr,'read','SITE/ID block',
     .                          1,'SNX_ID')
         else

*            Decode the line
             if( line(1:1).eq.' ' ) then
                 read(line,120,iostat=jerr) code, pt, domes, type,
     .                      name,  lng_deg, lng_min, lng_sec, 
     .                      lat_deg, lat_min, lat_sec, ht
  120            format(1X,A4, 1X,A2,1X,A9, 1X,A1, 1X,A22, 1X,
     .                  I3, 1X,I2, 1X,F4.1, 1X,
     .                  I3, 1X,I2, 1X,F4.1, 1X,F7.1 )
                 call report_error('IOSTAT',jerr,'decod',line,0,
     .                  'SNX_ID')
                 if( jerr.eq.0 ) then
                    ns = ns + 1
                    call ckmax(ns, max_snx_sites, 'MAX_SNX_SITES')
                    call casefold(code)

                    if( pt.eq.' A' ) then 
                        sx_site_code(ns) = code(1:4) // '_GPS'
                    else
                        sx_site_code(ns) = code(1:4) // '_G1' // 
     .                                     pt(2:2)
                    end if
                    sx_site_dome(ns) = domes
                    sx_site_long(ns) = name
                    sx_site_gpos(1,ns) = lng_deg + lng_min/60.d0 +
     .                                   lng_sec/3600.d0
                    sx_site_gpos(2,ns) = lat_deg + lat_min/60.d0 +
     .                                   lat_sec/3600.d0
                    sx_site_gpos(3,ns) =  ht 

*                   Initialise the counters for this sites
                    sx_num_rec(ns) = 0
                    sx_num_ant(ns) = 0
                    sx_num_ecc(ns) = 0
                 end if
             end if
         end if
      end do

****  Save the number of sites
      num_snx_sites = ns
      write(*,200) num_snx_sites
 200  format('Total of ',i4,' sites found in igs.snx file')

***** Thats all
      return
      end

CTITLE SNX_REC

      subroutine snx_rec(ierr)

      implicit none

*     Routine to read site RECIEVER block

      include '../includes/kalman_param.h'
      include 'hfupd.h'

* PASSED VARIABLES
* ----------------
* ierr -- IOSTAT error incase we hit end of file
      integer*4 ierr
       

* LOCAL VARIABLES
* --------------- 
* jerr -- IOSTAT error decoding line
* syr,sdoy,ssec -- Start yr, doy and sec
* eyr,edoy,esec -- End yr, doy and sec
* ke, ne  -- Site entry and total entry number
* indx    -- Pointer in string
* js      -- Site number for code 
 
      integer*4 jerr,  syr,sdoy,ssec, eyr,edoy,esec, ke, ne,
     .          indx, js


* block_end -- Logical to indicate that end of block reached 
      logical block_end

* code -- Code for station
* occ  -- Occuption entry
* pt   -- Point character for site
* type  -- P indicates GPS
* rec_type -- Type of receiver
* rec_sn   -- Serial number of receiver
* rec_fw   -- Firmware version of receiver
* line     -- Line read from file
* sx_code  -- Check variable name

      character*4 code, occ
      character*2 pt
      character*1 type
      character*20 rec_type
      character*5  rec_sn
      character*11 rec_fw
      character*8  sx_code 
      character*128 line


****  Start reading the lines until we hit end of file or end of
*     block
      block_end = .false.
      ne = 0

      do while ( ierr.eq.0 .and. .not.block_end ) 
         read(100,'(a)', iostat=ierr) line
         if( line(1:9).eq.'-SITE/REC' .or.ierr.ne.0 ) then
             block_end = .true.
             call report_error('IOSTAT',ierr,'read','SITE/REC block',
     .                          1,'SNX_REC')
         else

*            Decode the line
             if( line(1:1).eq.' ' ) then
                 read(line,120,iostat=jerr) code, pt, occ, type,
     .                 syr, sdoy, ssec, eyr, edoy, esec,
     .                 rec_type, rec_sn, rec_fw 
  120            format(1X,A4,  1X,A2,  1X,A4,  1X,A1,  
     .                  1X,I2.2,  1x,I3.3, 1x,I5.5, 
     .                  1X,I2.2,  1x,I3.3, 1x,I5.5, 
     .                  1X,A20,  1X,A5, 1X,A11 )   
                 call report_error('IOSTAT',jerr,'decod',line,0,
     .                  'SNX_REC')
                 if( jerr.eq.0 ) then

*                   Now find the station
                    if( pt.eq.' A' ) then 
                        sx_code = code(1:4) // '_GPS'
                    else
                        sx_code = code(1:4) // '_G1' // 
     .                                     pt(2:2)
                    end if
                    indx = 1
                    call get_cmd(sx_code, sx_site_code, 
     .                    num_snx_sites, js, indx)


*                   Warn user if site name not found
                    if( js.le.0 ) then
                       write(*,190)  js, line(1:80)
 190                   format('**WARNING** Bad site code match: ',
     .                        i4,/, a80)
                    else
*                      Increment the number of receiver records
*                      and save the pointer for this station to that
*                      record
                       ne = ne + 1
                       call ckmax(ne, max_snx_recs, 'MAX_SNX_RECS')
                       sx_num_rec(js) = sx_num_rec(js) + 1
                       ke = sx_num_rec(js)
                       call ckmax(ke, max_snx_ent_per_site,
     .                            'MAX_SNX_ENT_PER_SITE')
                       sx_rec_rec(ke,js) = ne
            
                       call yds_to_jd(syr,sdoy,ssec, sxrecv_st(ne))
                       call yds_to_jd(eyr,edoy,esec, sxrecv_en(ne))
*                      Check to see if end time OK and not set to
*                      zero if last entry
                       if( ke.gt.1 .and. 
     .                     sxrecv_en(ne-1).lt.sxrecv_st(ne) ) then
                           sxrecv_en(ne-1) = sxrecv_st(ne)
                       end if
                       sxrecv_sty(ne) = rec_type
                       sxrecv_sn(ne)  = rec_sn
                       sxrecv_fw(ne)  = rec_fw
                    end if

                 end if
             end if
         end if
      end do

***** Thats all
      sx_ent_rec = ne
      write(*,200) sx_ent_rec
 200  format('Total of ',i4,' Receiver entries found')
      return
      end

CTITLE SNX_ANT

      subroutine snx_ant(ierr)

      implicit none

*     Routine to read site ANTENNA block

      include '../includes/kalman_param.h'
      include 'hfupd.h'

* PASSED VARIABLES
* ----------------
* ierr -- IOSTAT error incase we hit end of file
      integer*4 ierr
       

* LOCAL VARIABLES
* --------------- 
* jerr -- IOSTAT error decoding line
* syr,sdoy,ssec -- Start yr, doy and sec
* eyr,edoy,esec -- End yr, doy and sec
* ke, ne  -- Site entry and total entry number
* indx    -- Pointer in string
* js      -- Site number for code 
 
      integer*4 jerr,  syr,sdoy,ssec, eyr,edoy,esec, ke, ne,
     .          indx, js, j


* block_end -- Logical to indicate that end of block reached 
      logical block_end

* code -- Code for station
* occ  -- Occuption entry
* pt   -- Point character for site
* type  -- P indicates GPS
* ant_type -- Type of antenna
* ant_sn   -- Serial number of antenna
* line     -- Line read from file

      character*4 code, occ
      character*2 pt
      character*1 type
      character*20 ant_type
      character*5  ant_sn
      character*8  sx_code
      character*128 line


****  Start reading the lines until we hit end of file or end of
*     block
      block_end = .false.
      ne = 0

      do while ( ierr.eq.0 .and. .not.block_end ) 
         read(100,'(a)', iostat=ierr) line
         if( line(1:9).eq.'-SITE/ANT' .or.ierr.ne.0 ) then
             block_end = .true.
             call report_error('IOSTAT',ierr,'read','SITE/ANT block',
     .                          1,'SNX_ANT')
         else

*            Decode the line
             if( line(1:1).eq.' ' ) then
                 read(line,120,iostat=jerr) code, pt, occ, type,
     .                 syr, sdoy, ssec, eyr, edoy, esec,
     .                 ant_type, ant_sn 
  120            format(1X,A4,  1X,A2,  1X,A4,  1X,A1,  
     .                  1X,I2.2,  1x,I3.3, 1x,I5.5, 
     .                  1X,I2.2,  1x,I3.3, 1x,I5.5, 
     .                  1X,A20,  1X,A5 )   
                 call report_error('IOSTAT',jerr,'decod',line,0,
     .                  'SNX_ANT')
                 if( jerr.eq.0 ) then

*                   Now find the station
                    if( pt.eq.' A' ) then 
                        sx_code = code(1:4) // '_GPS'
                    else
                        sx_code = code(1:4) // '_G1' // 
     .                                     pt(2:2)
                    end if
                    indx = 1
                    call get_cmd(sx_code, sx_site_code, 
     .                    num_snx_sites, js, indx)

*                   Warn user if site name not found
                    if( js.le.0 ) then
                       write(*,190)  js, line(1:80)
 190                   format('**WARNING** Bad site code match: ',
     .                        i4,/, a80)
                    else
*                      Increment the number of receiver records
*                      and save the pointer for this station to that
*                      record
                       ne = ne + 1
                       call ckmax(ne,max_snx_recs,'MAX_SNX_RECS')
                       sx_num_ant(js) = sx_num_ant(js) + 1
                       ke = sx_num_ant(js)
                       call ckmax(ke,max_snx_ent_per_site,
     .                            'MAX_SNX_ENT_PER_SITE')
                       sx_rec_ant(ke,js) = ne
            
                       call yds_to_jd(syr,sdoy,ssec, sxante_st(ne))
                       call yds_to_jd(eyr,edoy,esec, sxante_en(ne))
*                      Check to see if end time OK and not set to
*                      zero if last entry
                       if( ke.gt.1 .and.
     .                     sxante_en(ne-1).lt.sxante_st(ne) ) then
                           sxante_en(ne-1) = sxante_st(ne)
                       endif
                       sxante_sty(ne) = ant_type
                       sxante_sn(ne)  = ant_sn
                       do j = 1, 3
                          sxL1a_ecc(j,ne) = -99.0
                          sxL2a_ecc(j,ne) = -99.0
                       end do
                    end if

                 end if
             end if
         end if
      end do

***** Thats all
      sx_ent_ant = ne
      write(*,200) sx_ent_ant
 200  format('Total of ',i4,' antenna entries found')
      return
      end

CTITLE SNX_ECC

      subroutine snx_ecc(ierr)

      implicit none

*     Routine to read site ECCENTRICITY block

      include '../includes/kalman_param.h'
      include 'hfupd.h'

* PASSED VARIABLES
* ----------------
* ierr -- IOSTAT error incase we hit end of file
      integer*4 ierr
       

* LOCAL VARIABLES
* --------------- 
* jerr -- IOSTAT error decoding line
* syr,sdoy,ssec -- Start yr, doy and sec
* eyr,edoy,esec -- End yr, doy and sec
* ke, ne  -- Site entry and total entry number
* indx    -- Pointer in string
* js      -- Site number for code 
 
      integer*4 jerr,  syr,sdoy,ssec, eyr,edoy,esec, ke, ne,
     .          indx, js


* block_end -- Logical to indicate that end of block reached 
      logical block_end

* ecc(3) -- Values of the Eccentricity (UNE)

      real*4 ecc(3)

* code -- Code for station
* occ  -- Occuption entry
* pt   -- Point character for site
* type  -- P indicates GPS
* ecc_type -- Type of Eccentricity
* line     -- Line read from file

      character*4 code, occ
      character*2 pt
      character*1 type
      character*3 ecc_type
      character*8 sx_code
      character*128 line


****  Start reading the lines until we hit end of file or end of
*     block
      block_end = .false.
      ne = 0

      do while ( ierr.eq.0 .and. .not.block_end ) 
         read(100,'(a)', iostat=ierr) line
         if( line(1:9).eq.'-SITE/ECC' .or.ierr.ne.0 ) then
             block_end = .true.
             call report_error('IOSTAT',ierr,'read','SITE/ECC block',
     .                          1,'SNX_ECC')
         else

*            Decode the line
             if( line(1:1).eq.' ' ) then
                 read(line,120,iostat=jerr) code, pt, occ, type,
     .                 syr, sdoy, ssec, eyr, edoy, esec,
     .                 ecc_type, ecc
  120            format(1X,A4,  1X,A2,  1X,A4,  1X,A1,  
     .                  1X,I2.2,  1x,I3.3, 1x,I5.5, 
     .                  1X,I2.2,  1x,I3.3, 1x,I5.5, 
     .                  1X,A3, 1x,F8.4, 1x, F8.4, 1x, F8.4 )   
                 call report_error('IOSTAT',jerr,'decod',line,0,
     .                  'SNX_ECC')
                 if( jerr.eq.0 ) then

*                   Now find the station
                    if( pt.eq.' A' ) then 
                        sx_code = code(1:4) // '_GPS'
                    else
                        sx_code = code(1:4) // '_G1' // 
     .                                     pt(2:2)
                    end if
                    indx = 1
                    call get_cmd(sx_code, sx_site_code, 
     .                    num_snx_sites, js, indx)

*                   Warn user if site name not found
                    if( js.le.0 ) then
                       write(*,190)  js, line(1:80)
 190                   format('**WARNING** Bad site code match: ',
     .                        i4,/, a80)
                    else
*                      Increment the number of receiver records
*                      and save the pointer for this station to that
*                      record
                       ne = ne + 1
                       call ckmax(ne,max_snx_recs,'MAX_SNX_RECS')
                       sx_num_ecc(js) = sx_num_ecc(js) + 1
                       ke = sx_num_ecc(js)
                       call ckmax(ke,max_snx_ent_per_site,
     .                            'MAX_SNX_ENT_PER_SITE')
                       sx_rec_ecc(ke,js) = ne
            
                       call yds_to_jd(syr,sdoy,ssec, sxecce_st(ne))                       
                       call yds_to_jd(eyr,edoy,esec, sxecce_en(ne))

*                      Check to see if end time OK and not set to
*                      zero if last entry
                       if( ke.gt.1 .and.
     .                     sxecce_en(ne-1).lt.sxecce_st(ne) ) then
                           sxecce_en(ne-1) = sxecce_st(ne)
                       endif
                       sxarp_ecc(1,ne) = ecc(2)
                       sxarp_ecc(2,ne) = ecc(3)
                       sxarp_ecc(3,ne) = ecc(1)
                    end if

                 end if
             end if
         end if
      end do

***** Thats all
      sx_ent_ecc = ne
      write(*,200) sx_ent_ecc
 200  format('Total of ',i4,' eccentricity entries found')
      return
      end

CTITLE SNX_GPS_PC

      subroutine snx_gps_pc(ierr)

      implicit none

*     Routine to read site GPS_PHASE_CENTER block

      include '../includes/kalman_param.h'
      include 'hfupd.h'

* PASSED VARIABLES
* ----------------
* ierr -- IOSTAT error incase we hit end of file
      integer*4 ierr
       

* LOCAL VARIABLES
* --------------- 
* jerr -- IOSTAT error decoding line
* ne -- Number of entries found in block.  If none are found then
*        we read local copy.
 
      integer*4 jerr,   ne,   i, j, k

* block_end -- Logical to indicate that end of block reached
* upd       -- Indicates that we should update antenna phase
*              center information.
* notadded  -- Set true if an antenna has not been added to
*     list.  (Occurrs when antenna name does not match any stations
*     usually because of radome type.)  Here add an enrty. 
      logical block_end, upd, notadded

* L1a_ecc(3) -- L1 antenna offsets
* L2a_ecc(3) -- L2 antenna offsets

      real*4 L1a_ecc(3), L2a_ecc(3)

* ant_type -- Type of receiver
* ant_sn   -- Serial number of receiver
* line     -- Line read from file

      character*20 ant_type
      character*5  ant_sn
      character*128 line


****  Start reading the lines until we hit end of file or end of
*     block
      block_end = .false.
      ne = 0

      do while ( ierr.eq.0 .and. .not.block_end ) 
         read(100,'(a)', iostat=ierr) line
         if( line(1:9).eq.'-SITE/GPS' .or.ierr.ne.0 ) then
             block_end = .true.
             call report_error('IOSTAT',ierr,'read',
     .                         'SITE/GPS_PC block',  1,'SNX_GPS_PC')
         else

*            Decode the line
             if( line(1:1).eq.' ' ) then
                 call sub_char(line(27:),'------','  -10.')
                 read(line,120,iostat=jerr)   ant_type, ant_sn,
     .                  L1a_ecc, L2a_ecc
  120            format(1X,A20, 1X,A5, 6(1X,F6.4) )   
                 call report_error('IOSTAT',jerr,'decod',line,0,
     .                  'SNX_GPS_PC')
                 if( jerr.eq.0 ) then

*                   Now we need to match this antenna phase center
*                   values to all antenna entries with this type
*                   of antenna.
                    notadded = .true.
                    ne = ne + 1
                    do i = 1, sx_ent_ant
                       if( ant_type.eq.sxante_sty(i) .and.
     .                     (ant_sn.eq.sxante_sn(i) .or.
     .                      ant_sn.eq.'-----')          ) then
*                          Match found to both type and serial number
*                          or serial number is generic (only update
*                          generic if we do not already have values).
                           upd = .true.
                           if( ant_sn.eq.'-----' .and. 
     .                         sxL1a_ecc(1,i).ne. -99.0 ) upd = .false.
                           if( upd ) then
                               notadded = .false.
                               sxL1a_ecc(1,i) = L1a_ecc(2)
                               sxL1a_ecc(2,i) = L1a_ecc(3)
                               sxL1a_ecc(3,i) = L1a_ecc(1)
                               sxL2a_ecc(1,i) = L2a_ecc(2)
                               sxL2a_ecc(2,i) = L2a_ecc(3)
                               sxL2a_ecc(3,i) = L2a_ecc(1)
                           end if
                       end if
                    end do
* MOD TAH 061110:   See if not added 
                    if( notadded ) then
                        write(*,140) ant_type
 140                    format('No direct use of antenna ',a20,
     .                         ' adding')
                        i = sx_ent_ant + 1
                        sxante_sty(i)  = ant_type
                        sxante_sn(i)   = ant_sn
                        sxL1a_ecc(1,i) = L1a_ecc(2)
                        sxL1a_ecc(2,i) = L1a_ecc(3)
                        sxL1a_ecc(3,i) = L1a_ecc(1)
                        sxL2a_ecc(1,i) = L2a_ecc(2)
                        sxL2a_ecc(2,i) = L2a_ecc(3)
                        sxL2a_ecc(3,i) = L2a_ecc(1)
                        sx_ent_ant = i
                    end if

                end if
             end if
          end if
      end do

* MOD TAH 061109: Check to see if have some antennas that have not matched 
*     yet.  This is case for Radome in site entries that has no corresponding
*     value in the PCV entries.  In this case we match to NONE entry.
      do i = 1, sx_ent_ant  
         if( abs(sxL1a_ecc(1,i)).gt.10.0 ) then
*            This antenna has never been updated.  Getting the antenna
*            type and match to radome values
             ant_type = sxante_sty(i)
             upd = .false.
*            Now scan list find NONE match
             do j = 1,sx_ent_ant
                if( i.eq.2 ) then
                   write(*,210) j, ant_type(1:16),sxante_sty(j)(1:16),
     .                    sxante_sty(j)(17:20)
 210               format('M ',i4,1x,'|',a16,'|',a16,'|',a4,'|') 
                endif
                if( ant_type(1:16).eq.sxante_sty(j)(1:16) .and.
     .              sxante_sty(j)(17:20).eq.'NONE' .and. 
     .              .not.upd ) then
                    write(*,220) i,ant_type,j, sxante_sty(j)
 220                format('No PCV for antenna Entry ',i4,1x,a20,
     .                     ' Matching to ',i4,1x,a20)
                    do k = 1,3
                       sxL1a_ecc(k,i) =sxL1a_ecc(k,j)
                       sxL2a_ecc(k,i) =sxL2a_ecc(k,j)
                    end do
                    upd = .true.
                end if
             end do
*            See if we found match
             if( .not. upd ) then
                 write(*,240) i,ant_type
 240             format('**ERROR** No PCV offsets for Antenna entry ',
     .                  i4,1x,a20)
             end if
         end if
      end do

***   If there are no entries: Stop
      if( ne.eq.0 ) then
         call report_stat('FATAL','HFUPD','snx_gps_pc',igs_snx_file,
     .                    'No GPS_PHS_CENTER data',0)
      endif

          
***** Thats all
      return
      end

CTITLE SNX_REP

      subroutine snx_rep

      implicit none

*     Routine to report the results of reading the igs.snx file.

      include '../includes/kalman_param.h'
      include 'hfupd.h'


* LOCAL VARIABLES
* --------------- 
* syr,sdoy,ssec -- Start yr, doy and sec
* eyr,edoy,esec -- End yr, doy and sec
* i,j, k        -- Loop counters 
      integer*4  syr,sdoy,ssec, eyr,edoy,esec, i,j,k,l


***** Loop over the stations
      do i = 1, num_snx_sites
         write(*,120) i, sx_site_code(i), sx_site_dome(i),
     .                sx_site_long(i), sx_num_rec(i), 
     .                sx_num_ant(i), sx_num_ecc(i)
 120     format(i4,'. ',a8,1x,a9,1x,a22,' RECV, ANTE, ECC ',3i3)

*        Now output each type of record
         do j = 1, sx_num_rec(i)
            k = sx_rec_rec(j,i)
            call jd_to_yds(sxrecv_st(k), syr, sdoy, ssec)
            call jd_to_yds(sxrecv_en(k), eyr, edoy, esec)
            write(*,160) k, syr, sdoy, ssec, eyr, edoy, esec,
     .                   sxrecv_sty(k), sxrecv_sn(k), sxrecv_fw(k)
 160        format(4x,'REC ',i4,'. ',i2,':',i3.3,':',i5.5,1x,
     .             i2,':',i3.3,':',i5.5,1x,a20,1x,a5,1x,a11)
         end do

         do j = 1, sx_num_ant(i)
            k = sx_rec_ant(j,i)
            call jd_to_yds(sxante_st(k), syr, sdoy, ssec)
            call jd_to_yds(sxante_en(k), eyr, edoy, esec)
            write(*,220) k, syr, sdoy, ssec, eyr, edoy, esec,
     .                   sxante_sty(k), sxante_sn(k),
     .                   (sxL1a_ecc(l,k), l=1,3),
     .                   (sxL2a_ecc(l,k), l=1,3)
 220        format(4x,'ANT ',i4,'. ',i2,':',i3.3,':',i5.5,1x,
     .             i2,':',i3.3,':',i5.5,1x,a20,1x,a5,1x,
     .             3(F6.4,1x),1x,3(F6.4,1x))
         end do

         do j = 1, sx_num_ecc(i)
            k = sx_rec_ecc(j,i)
            call jd_to_yds(sxecce_st(k), syr, sdoy, ssec)
            call jd_to_yds(sxecce_en(k), eyr, edoy, esec)
            write(*,260) k, syr, sdoy, ssec, eyr, edoy, esec,
     .                   (sxarp_ecc(l,k), l=1,3)
 260        format(4x,'ECC ',i4,'. ',i2,':',i3.3,':',i5.5,1x,
     .             i2,':',i3.3,':',i5.5,1x,
     .             3(F8.4,1x))
         end do
      end do

****  Thats all
      return
      end

CTITLE CHECK_SINF
 
      subroutine check_sinf( kerr ) 

      implicit none

*     Routine to check the site information in the sinf_def block
*     against the values found in the igs.snx format

      include '../includes/kalman_param.h'
      include '../includes/const_param.h' 
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
      include 'hfupd.h'

* PASSED VARIABLES
* kerr -- IOSTAT file error

      integer*4 kerr


* LOCAL VARIABLES
* syr,sdoy,ssec -- Start yr, doy and sec
* eyr,edoy,esec -- End yr, doy and sec
* i,j, k        -- Loop counters
* js            -- Sinex site number
* indx          -- Pointer in string
* nu            -- Short version of number of site to be updated
* trimlen       -- Length of string
 
      integer*4  syr,sdoy,ssec, eyr,edoy,esec, i,j,k,l, js, indx,
     .           nu, trimlen

* values_OK     -- Logical to indicate if the antenna offset 
*                  information matches between sinex file and
*                  hfile.
* ant_ent_OK    -- Logical set false if all antenna entries are
*                  zero.
* all_ant_zero  -- Set true if all antennas have zero L1 and L2 offset
*                  indicating that the ANTENNA block is not in the sinex file.
* rec_ent_OK    -- Set true if the receivers entries are OK

* MOD RWK 150213: Frequencies for ecc update (tempoary for GPS only)
      real*8 fL1,fL2 

      logical values_OK, ant_ent_OK, all_ant_zero, rec_ent_OK

* dome_num     -- Domes number from long site name
      character*12 dome_num
* radome       -- Converted raydome type
* short_name   -- 4-char name of site
      character*4 radome
      character*8 short_name
      character*20 outant  ! Antenna radome name with no +- sign 

 
****  Start be reading the names from the binary files
      sx_num_upd = 0

****  Check the antenna information to see if not all values
*     zero suggesting no values were present
      ant_ent_OK = .false.
      rec_ent_OK = .false.
      i = 0
      do while ( .not.ant_ent_OK .and. i.lt. cnum_sites )
         i = i + 1
         do j = 1, 3
            if( qL1a_ecc(j,i).ne.0 ) ant_ent_OK = .true.
            if( qL2a_ecc(j,i).ne.0 ) ant_ent_OK = .true.
         end do
      end do

      all_ant_zero = .not.ant_ent_OK

      if( .not.ant_ent_OK ) then
          write(*,100)
 100      format('***WARNING*** All antenna offsets are zero.',
     .           'NOT UPDATING ENTRIES' )
          call report_stat('WARNING','HFUPD','check_sinf',hfile,
     .           'HFILE has no antenna information',0)
* MOD TAH 031112: Removed the return so that code will continue
*         and adopt standard values for the antenna offsets.
C         RETURN
      end if

****  Tell user what is happening
      do i = 1, cnum_sites

*        See if we can find matching igs.snx entry
         values_OK = .true.
         indx = 1
! MOD TAH 120823: Check just the 4-char code not the whole name
         short_name = qsite_names(i)(1:4)
         call get_cmd(short_name, sx_site_code,  
     .                    num_snx_sites, js, indx)
!        call get_cmd(qsite_names(i)(1:4), sx_site_code,  
!    .                    num_snx_sites, js, indx)
         if( js.le.0 ) then
             write(*,120) js, i, qsite_names(i), qfull_names(i) 
120          format('NoMatch ',i4,1x, i4,'. ',a8,1x,a32)

*            No match on site code.  See if we can match on 
*            domes number
             dome_num = qfull_names(i)(23:31)
             indx = 1
             call get_cmd(dome_num, sx_site_dome,
     .                num_snx_sites, js, indx)

!             if( js.gt.0 ) then 
!             Only match if 4-char ID also matches (avois problems
!             with multiple receivers on same mark).
              if( js.gt.0  ) then 
                 if( qsite_names(i)(1:4).eq.sx_site_code(js)(1:4) ) then
                    write(*,125) js, sx_site_code(js), dome_num
 125                format('Match found with site ',i4,1x,a8,
     .                     ' using dome number ',a9)
                    qsite_names(i) = sx_site_code(js)
                 else
                    write(*,127) js, sx_site_code(js), qsite_names(i)
 127                format('Match found with site ',i4,1x,a8,
     .                     ' using dome number but ID different ',a8)
                 endif
             end if
         end if

* MOD TAH 990328: Report warning if no match
        if( js.eq.0 ) then
             call report_stat('WARNING','HFUPD','check_sinf',
     .                sx_site_code(js),
     .                'NoMatch to site name ',0)
        endif
           

* MOD TAH 990328: Check to see if we have antenna offsets
         if( js.gt.0 ) then
            ant_ent_OK = .false.
            if( qL1a_ecc(3,i).ne.0 ) ant_ent_OK = .true.
            if( qL2a_ecc(3,i).ne.0 ) ant_ent_OK = .true.
* MOD TAH 990402: Also check the antenna description:
* MOD TAH 031112: Only check the dscription if all offsets are
*           not zero.
            if( qante_ty(i)(1:8).ne.'--------' .and.
     .          trimlen(qante_ty(i)).ne.0      .and.
     .          ichar(qante_ty(i)(1:1)).ne.0 .and.
     .          .not. all_ant_zero ) ant_ent_OK = .true.
 
            if( .not. ant_ent_OK ) then
                 write(*,130) sx_site_code(js)
 130             format('NoAntenna information for site ',a8)
                 call report_stat('WARNING','HFUPD','check_sinf',
     .                    sx_site_code(js),
     .                   'Site has not antenna information',0)
* MOD TAH 031112: Allow code to continue and use standard values
C                js = 0
             endif
         endif
                 
*        OK, Site seems to be OK
         if( js.ne.0 ) then

*           Now output each type of record
            call jd_to_yds(qrecv_st(i),syr, sdoy, ssec)
            call jd_to_yds(qrecv_en(i),eyr, edoy, esec)

*           Now find the corresponding entry and make sure that it
*           matches.
            rec_ent_OK = .true.
            do j = 1, sx_num_rec(js)
               k = sx_rec_rec(j,js)
               if( cepoch_expt.ge. sxrecv_st(k) .and.
     .             (cepoch_expt.lt. sxrecv_en(k) .or.
     .              sxrecv_en(k).eq.jd_zero ) ) then

*                 Get times 
                  call jd_to_yds(sxrecv_st(k), syr, sdoy, ssec)
                  call jd_to_yds(sxrecv_en(k), eyr, edoy, esec)

*                 We have found the correct record:  See if matches
                  if ( qrecv_ty(i).ne.sxrecv_sty(k) ) 
     .                                     rec_ent_OK = .false.
                  if ( qrecv_sn(i)(1:5).ne.sxrecv_sn(k)(1:5) )
     .                                     rec_ent_OK = .false.
                  if ( qrecv_fw(i)(12:15).ne. sxrecv_sty(k)(17:20) )
     .                                     rec_ent_OK = .false.
                  if( .not. rec_ent_OK .and. report_rec ) then
                        write(*,210) qsite_names(i), cowner,
     .                        k, syr, sdoy, ssec, eyr, edoy, esec,
     .                        sxrecv_sty(k), sxrecv_sn(k),
     .                        qrecv_ty(i)  , qrecv_sn(i)
 210                    format(a8,1x,a4,' REC ',i5,'. ',
     .                         i2,':',i3.3,':',i5.5,1x,
     .                         i2,':',i3.3,':',i5.5,1x,a20,1x,a5,1x,
     .                         ' Used ', a20,1x,a5)
                  end if


*                 Update informtation
                  if( sxrecv_sty(k)(1:4).ne.'----' )
     .                qrecv_ty(i) = sxrecv_sty(k)
                  if( sxrecv_sn(k)(1:4).ne.'----' )
     .                qrecv_sn(i) = sxrecv_sn(k)
*                 Save the last part of the recv type in the fw slot where
*                 have extra space for it. 
                  if( sxrecv_sty(k)(1:4).ne.'----' )
     .                qrecv_fw(i)(12:15) = sxrecv_sty(k)(17:20)


               end if
            end do 

*           Now check the anntena entries.  Here we may need to
*           update the actual solution information as well.
            call jd_to_yds(qante_st(i),syr, sdoy, ssec)
            call jd_to_yds(qante_en(i),eyr, edoy, esec)

            do j = 1, sx_num_ant(js)
               k = sx_rec_ant(j,js)

*              See if time range matches
               if( cepoch_expt .ge. sxante_st(k) .and.
     .             (cepoch_expt  .lt. sxante_en(k) .or.
     .              sxante_en(k) .eq. jd_zero)  ) then
*                  We have found a matching time.  Now check the values
                   if( ant_ent_OK ) then
                     values_OK = .true.
* MOD TAH 170815: See if antenna name matches
                     if( sxante_sty(k)(1:15).ne.qante_ty(i)(1:15) ) 
     .                   values_OK = .false.
                     do l = 1, 3
                        if( abs(qL1a_ecc(l,i)-sxL1a_ecc(l,k)).gt.2.d-4 ) 
     .                                            values_OK = .false.
                        if( abs(qL2a_ecc(l,i)-sxL2a_ecc(l,k)).gt.2.d-4 ) 
     .                                            values_OK = .false.
                     end do

* MOD TAH 170815: See if antenna name matches
* MOD TAH 061220: Temporary mod for transition to ABS phase center model
*                    values_OK = .true.
* MOD TAH 070123: Commented out transition logic.

                     if ( .not.values_OK ) then
*                      Convert the radome number
                       if( ichar(qradome_ty(i)(1:1)).eq.0 ) then
                           radome = qante_sn(i)(6:8)
                       else
                           radome = qradome_ty(i)
                       endif

                       outant = qante_ty(i)(1:15) 

                       write(*,220) qsite_names(i), cowner, 0,
     .                         syr, sdoy, ssec, eyr, edoy, esec,
     .                         outant(1:16) // radome, 
     .                         qante_sn(i)(1:5),
     .                      (qL1a_ecc(l,i), l=1,3),
     .                      (qL2a_ecc(l,i), l=1,3)


                        call jd_to_yds(sxante_st(k), syr, sdoy, ssec)
                        call jd_to_yds(sxante_en(k), eyr, edoy, esec)
                        write(*,220) qsite_names(i), cowner,
     .                        k, syr, sdoy, ssec, eyr, edoy, esec,
     .                        sxante_sty(k), sxante_sn(k),
     .                        (sxL1a_ecc(l,k), l=1,3),
     .                        (sxL2a_ecc(l,k), l=1,3)
                     
 220                    format(a8,1x,a4,' ANT ',i5,'. ',
     .                         i2,':',i3.3,':',i5.5,1x,
     .                         i2,':',i3.3,':',i5.5,1x,a20,':',a5,1x,
     .                         3(F6.4,1x),1x,3(F6.4,1x))
                        write(*,225) cowner, qsite_names(i)(1:4),
     .                         outant(1:16) // radome, 
     .                         (qL1a_ecc(l,i), l=1,3),
     .                         (qL2a_ecc(l,i), l=1,3),
     .                         hfile(1:trimlen(hfile)),  
     .                         cowner,qsite_names(i)(1:4),
     .                         sxante_sty(k),
     .                         (sxL1a_ecc(l,k), l=1,3), 
     .                         (sxL2a_ecc(l,k), l=1,3),
     .                          hfile(1:trimlen(hfile))
 225                    format('Err: ',a4,1x,a4,1x,'Used ',a20,' dL1 ',
     .                         3(F8.4,1x),' dL2 ',3(F8.4,1x),' m ',a,/,
     .                         'Err: ',a4,1x,a4,1x,'Real ',
     .                         a20,' dL1 ', 3(F8.4,1x),' dL2 ',
     .                         3(F8.4,1x), ' m ',a)

****                    Now save the adjustments to be updated
* MOD RWK 150213: Set frequencies to GPS values (see warning in main hfupd)
                        fL1 = 154*10.23d6     
                        fL2 = 120*10.23d6   
                        sx_num_upd = sx_num_upd + 1
                        nu = sx_num_upd
                        call ckmax(nu, max_snx_sites, 'MAX_SNX_SITES')
*                       Save the site number and NEU adjustments
                        dNEU_site(nu) = i
                        do l = 1,3
                           dNEU_ecc(l,nu) = (fL1**2*
     .                             (qL1a_ecc(l,i)-sxL1a_ecc(l,k)) -
     .                              fL2**2*
     .                             (qL2a_ecc(l,i)-sxL2a_ecc(l,k)))/
     .                             (fL1**2-fL2**2)
*                          Now save the correct values
                           qL1a_ecc(l,i) = sxL1a_ecc(l,k)
                           qL2a_ecc(l,i) = sxL2a_ecc(l,k)
                        end do
                     end if
                   else
* MOD TAH 031112:     Just save the values with out update
                      do l = 1,3
                         qL1a_ecc(l,i) = sxL1a_ecc(l,k)
                         qL2a_ecc(l,i) = sxL2a_ecc(l,k)
                      end do
                   end if 

*                  Update the antenna information to match sinex
*                  file. 
                   if( sxante_sty(k)(1:4).ne.'----' )
     .                 qante_ty(i) = sxante_sty(k)
                   if( sxante_sn(k)(1:4).ne.'----' )
     .                 qante_sn(i) = sxante_sn(k)
*                  Save the last part of the ante type in the fw slot 
*                  where have extra space for it. 
                   if( sxante_sty(k)(1:4).ne.'----' )
     .                qante_sn(i)(6:8) = sxante_sty(k)(17:19)

               end if
            end do

*****       Do the eccentricity
            call jd_to_yds(qdata_st(i),syr, sdoy, ssec)
            call jd_to_yds(qdata_en(i),eyr, edoy, esec)

             do j = 1, sx_num_ecc(js)
               k = sx_rec_ecc(j,js)
               if( cepoch_expt.ge. sxecce_st(k) .and.
     .            (cepoch_expt.lt. sxecce_en(k) .or.
     .             sxecce_en(k) .eq. jd_zero)  ) then 

*                  Found the corrent entry: Now see if OK
                   values_OK = .true.
                   do l = 1, 3
                      if( abs(qarp_ecc(l,i)-sxarp_ecc(l,k)).gt.2.d-4)
     .                    values_OK = .false.
                   end do

                   if( .not.values_OK ) then

                      call jd_to_yds(sxecce_st(k), syr, sdoy, ssec)
                      call jd_to_yds(sxecce_en(k), eyr, edoy, esec)
                      write(*,260) qsite_names(i), cowner,
     .                      0, syr, sdoy, ssec, eyr, edoy, esec,
     .                      (qarp_ecc(l,i), l=1,3)
                      write(*,260) qsite_names(i), cowner,
     .                       k, syr, sdoy, ssec, eyr, edoy, esec,
     .                      (sxarp_ecc(l,k), l=1,3)
 260                  format(a8,1x,a4,' ECC ',i5,'. ',
     .                      i2,':',i3.3,':',i5.5,1x,
     .                      i2,':',i3.3,':',i5.5,1x,
     .                       3(F8.4,1x))
                      write(*,265) cowner, qsite_names(i)(1:4),
     .                      (qarp_ecc(l,i), l=1,3),
     .                      (sxarp_ecc(l,k), l=1,3),
     .                       hfile(1:trimlen(hfile))
265                   format('Err: ',a4,1x,a4,' Eccentricity: Used',8x,
     .                       3(F8.4,1x),'Real ',3(F8.4,1x),' m ',a)


*                     Now save the updates.  See if we already have
*                     site.
                      if( sx_num_upd.gt.0 .and.
     .                    dNEU_site(sx_num_upd).eq.i ) then
*                         Already have update, so add new corrections
                          nu = sx_num_upd
                          do l = 1,3
                             dNEU_ecc(l,nu) = dNEU_ecc(l,nu) +
     .                           (qarp_ecc(l,i)-sxarp_ecc(l,k))
                             qarp_ecc(l,i) = sxarp_ecc(l,k)
                          end do
                      else        
                          sx_num_upd = sx_num_upd + 1
                          nu = sx_num_upd
                          call ckmax(nu,max_snx_sites,'MAX_SNX_SITES')
                          dNEU_site(nu) = i
                          do l = 1,3
                             dNEU_ecc(l,nu) = 
     .                           (qarp_ecc(l,i)-sxarp_ecc(l,k))
                             qarp_ecc(l,i) = sxarp_ecc(l,k)
                          end do
                      end if
                   end if 
                end if
            end do
         end if
     
      end do

****  Thats all
      return
      end

CTITLE READ_HFINF
 
      subroutine read_hfinf 

      implicit none

*     Routine to read the information blocks from the hfile.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h' 
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
      include 'hfupd.h'


****  Read the globk hfile blocks that we will need.

      call sr_names

***** Get the full names
      call sr_full

***** Get the solution description
      call sr_description

***** Get the codes for the estimated parameters
      call sr_codes

***** Get the apriori parameter codes
      call sr_aprioris( temp_r8 )

****  Thats all
      return
      end

CTITLE UPD_HF_SOLN

      subroutine upd_hf_soln( kerr )

      implicit none

*     Routine to to update the hfile and write out the updated hfile/
*     (The original Hfile is overwritten---strictly only the updated
*     records are re-written).

      include '../includes/kalman_param.h'
      include '../includes/const_param.h' 
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
      include 'hfupd.h'

* PASSED VARIABLES
* kerr    -- IOSTAT error reading files

      integer*4 kerr

* LOCAL VARIABLES
* nwrd_cov -- Number of I*4 words needed for the covariance matrix
* nrec_cov -- Number of 128 word records that this represents
* nrec_sol -- Number of records needed to read solution and the tail
*             of the covariance matrix (fractional 128 word record part)
* crec_sol -- Record number for start of end of covariance matrix and
*             start of solution vector
* start_sol -- last element number of covariance matrix before start
*             of solution (i.e., parameter 1 is at start_sol+1)
* len_read -- Length of record read from file
* i,j      -- Loop counter
* sn       -- Station number
* np       -- Parameter number to be updated
* npa(3)   -- XYZ parameter numbers
* ierr, jerr -- IOSTAT errors reading polar motion table
* date(5)    -- Date read from polar motion tables
* trimlen    -- Length of string
* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 nwrd_cov
      integer*4 nrec_cov, start_sol, nrec_sol, 
     .          crec_sol, len_read, i, j, sn, np, npa(3),
     .          ierr, jerr, date(5), trimlen

* loc_coord(3)  -- Local coordinates of site
* rot_matrix(3,3) -- Rotation matrix
* dXYZ(3)       -- Change in XYZ coordinates
* dXYZ_ptd(3), dNEU_ptd(3) -- Pole tide corrections in XYZ and
*                  NEU (m)
* pmx, pmy      -- Polar motion values (mas)
* mpx, mpy      -- Mean pole position to which correction needs
*                  to be applied (IERS2010 model used)
* pmvals(4)     -- Polar motion values read from pmu file
* dt, dtmin     -- Time difference between experiment epoch and
*                  polar motion table epoch and min value
* sectag        -- Seconds tag for ymdhms_to_jd
* jd_pmu        -- Julian date of pmu table entry

      real*8 loc_coord(3), rot_matrix(3,3), dXYZ(3),
     .       dXYZ_ptd(3), dNEU_ptd(3), pmx, pmy, pmvals(4),
     .       dt, dtmin, sectag, jd_pmu
      real*8 mpx, mpy

* done  -- Set true when polar motion entries found
* appOK -- Set true if OK to apply pole tide
* kbit  -- Checks bit status 
* write_soln -- Set true when we actually need to write 
*     solution vector

      logical done, appOK, kbit, write_soln

* line  -- Line read from polar motion file
      character*256 line

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      data I8 / 1 /
 
****  See if we have any update to actually make
      write_soln = .false.
      if( sx_num_upd.eq.0 .and. .not.app_ptide ) RETURN

*     OK, we have updates.  First need to read the solution vector
*     in so that we can make the updates.   Compute the starting block
*     number for the solution vector resides in the file.
      nwrd_cov = ((I8*cnum_parn)*cnum_parn)*2
      nrec_cov = (nwrd_cov-1)/128
      start_sol = (nwrd_cov - nrec_cov*128)/2 
      nrec_sol = cnum_par_vals - nrec_cov
      crec_sol = crec_par_vals + nrec_cov 

*     Read the last part of covariance matrix and solution vector
      call readd(cglb_dcb, kerr, temp_r8, 128*nrec_sol, len_read, 
     .           crec_sol)

****  Now start updating the entries
      do i = 1, sx_num_upd
         sn = dNEU_site(i)
         call rotate_geod(dNEU_ecc(1,i), dxyz, 'NEU',
     .                    'XYZ', site_pos(1,sn), loc_coord,
     .                     rot_matrix) 

*        Now get the parameter numbers
         do j = 1, 3
            np = qparn_sites(j,sn)
*           Update the solution unless we have been told not to
*           update.
            if( np.gt.0 .and. .not.kbit(honly_opt,1) ) then
                write_soln = .true.
                temp_r8(start_sol+np) = temp_r8(start_sol+np) +
     .                                  dxyz(j)
            end if
            npa(j) = np
         end do
         if( .not.kbit(honly_opt,1) ) then 
             write(*,200) cowner, i, qsite_names(sn), 
     .                   (dNEU_ecc(j,i), j =1,3)
 200         format('Update ',a4,1x,i5,1x,a8,' dNEU ',3f8.4,' m')
         end if
      end do

***** See if we should apply the ptide correction
*     Check for swap problems
! MOD TAH 110512: Removed all test for swaph problem
!     if( cgamit_mod.gt.2**30 )  call swap_bytes(4, cgamit_mod, 1)
!     if( app_ptide .and. .not.kbit(cgamit_mod,19) ) then
      if( app_ptide ) then
          appOK = .false.
          if( (cwob_apr(1).ne.0 .and.cwob_apr(2).ne.0) .and.
     .         trimlen(pmu_file).eq.0        ) then
              pmx = cwob_apr(1)
              pmy = cwob_apr(2)
              appOK = .true.
!              write(*,310) pmx, pmy
! 310          format('Applying pole tide with PM values ',
!     .               2F10.4,' mas from hfile')

*         See if polar motion file, if so use this.
          elseif( trimlen(pmu_file).gt.0 ) then
              open(110, file=pmu_file, iostat=ierr, status='old')
              call report_error('IOSTAT',ierr,'open',pmu_file,0,
     .                          'HFUPD')
              if( ierr.eq.0 ) then 
                 done = .false.
                 dtmin = 1.d10 
                 do while ( .not.done )
                    read(110,'(a)',iostat=ierr) line
                    if( ierr.eq.0 .and. line(1:1).eq.' ') then
                        read(line,*,iostat=jerr) date, pmvals
                        if( jerr.eq.0 ) then
                           sectag = 0
                           call ymdhms_to_jd(date, sectag, jd_pmu)
                           dt = cepoch_expt - jd_pmu
                           if( abs(dt).lt.dtmin ) then
                               dtmin = abs(dt)
                               pmx = pmvals(1)*1000
                               pmy = pmvals(3)*1000
                           end if
                           if( dt.lt.0 ) done = .true.
                        endif
                    end if
                    if( ierr.ne.0 ) done = .true.
                 end do
                 if( ierr.eq.0 .and. dtmin.lt.2 ) then
                     appOK = .true.
!                    write(*,320) pmx, pmy, 
!    .                        pmu_file(1:trimlen(pmu_file)),dtmin
!320                 format('Applying pole tide with PM values ',
!    .               2F10.4,' mas from ',a,' DT ',f5.3,' days')
                 else
                     write(*,325) ierr, dtmin
 325                 format('Unable to apply Pole tide ierr, dtmin ',
     .                      i5,F15.2)
                 end if
              end if
          else
              write(*, 330) 
 330          format('Pole can not be applied:  No pmu values ',
     .               'give pmu file with -p option')
          endif

****      See if can apply the pole correction
          if( appOK ) then
*            Compute the difference from the mean pole
*            that is needed
             if ( .not. kbit(cgamit_mod,19) ) then
!                Pole tide never applied; apply to IERS2010 mean
                 call mean_pole( cepoch_expt,'IERS10', mpx, mpy)
                 pmx = pmx - mpx
                 pmy = pmy - mpy
             elseif( kbit(cgamit_mod,21) ) then ! See if applied
                                                ! to old mean pole 
                 call mean_pole( cepoch_expt,'IERS96', pmx, pmy)
                 call mean_pole( cepoch_expt,'IERS10', mpx, mpy)
!                Adjust the position
                 pmx = pmx - mpx
                 pmy = pmy - mpy
             elseif( kbit(cgamit_mod,23) ) then
!                Correction already done, no update
                 pmx = 0.d0
                 pmy = 0.d0
                 appOK = .false.
             endif
                
             if( appOK ) then 
                 write(*,340) pmx, pmy 
340              format('Applying pole tide with offsets of ',
     .                   2F10.4,' mas ')
                 do i = 1, cnum_sites

                    call comp_ptide(site_pos(1,i), pmx, pmy,
     .                           dNEU_ptd, dXYZ_ptd)
                    do j = 1,3 
                       np = qparn_sites(j,i)
                       if( np.gt.0 .and. .not.kbit(honly_opt,2) ) then
                           write_soln = .true.
                           temp_R8(start_sol+np) =  
     .                            temp_R8(start_sol+np) - dXYZ_ptd(j)
                       endif
                    enddo
                 end do
                 call sbit(cgamit_mod,19,1)   ! Set pole tide applied
                 call sbit(cgamit_mod,21,0)   ! Remove IERS96 Mean Pole
                 call sbit(cgamit_mod,23,1)   ! Set IERS2010 Mean Pole.
             else 
*               Correction already applied
                write(*,370 ) 
 370            format('Pole tide already applied to this hfile')
             endif
          end if
      end if
                   
*     Now write out the global (Only if we have been told to do so): 
*     Try not to confuse the user by only printing message is we 
*     are actually updating and at least the ant or ptd values are
*     being updated.
      if( upd_hf .and. write_soln  ) then
          write(*,'(a)') 'Updating solution vector'
          call writd(cglb_dcb, kerr, temp_r8, 128*nrec_sol,  crec_sol)
          call report_error('IOSTAT',kerr,'writ','Binary Hfile',
     .                      1,'UPD_HF_SOLN')
      end if

      return
      end

 
CTITLE UPD_HF_HDR

      subroutine upd_hf_hdr( kerr )

      implicit none

*     Routine to to update the hfile and write out the updated hfile/
*     (The original Hfile is overwritten---strictly only the updated
*     records are re-written).

      include '../includes/kalman_param.h'
      include '../includes/const_param.h' 
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
      include 'hfupd.h'

* PASSED VARIABLES
* kerr  -- IOSTAT error
      integer*4 kerr


***** This routine will write back header information of the
*     global hfile.
      write(*,120) 
 120  format('Updating the hfile header information')
     
****  Do each of the major blocks
 
      call writd( cglb_dcb, kerr, glb_header, 128*cnum_header, 1)
      call report_error('FmpWrite',kerr,'writ','global header',1,
     .                  'qw_glb_header')           

      call qw_glb_names
      call qw_glb_full
      call qw_codes

*     set the number of stations and write the station information
      qnum_sites = cnum_sites
      call qw_description 

****  Thats all we need to do
      return
      end

CTITLE CKMAX

      subroutine ckmax(na, max, name )

      implicit none

*     Routine to check bounds on arrays
*
* na  -- Current number of values
* max -- Maxiumum number allowed
      integer*4 na, max

* name -- Name of parameter being exceeded

      character*(*) name

      if( na.le.max ) RETURN
      write(*,120) name, max
 120  format('***FATAL*** Parameter ',a,' has been exceeded.',
     .       ' Maximum allowed in hfupd.f ',i5)
      call report_stat('FATAL','HFUPD',' ',name,
     .                 ' dimension exceeded. Max ', max)

      return
      end

CTITLE CHECK_CWU
 
      subroutine check_cwu( kerr ) 

      implicit none

*     Routine to check the start and stop times of CWU files and update 
*     if needed.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h' 
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
      include 'hfupd.h'

* PASSED VARIABLES
* kerr -- IOSTAT file error

      integer*4 kerr


* LOCAL VARIABLES
      real*8 sectag
     .,      rt_jd    ! Creation date of file (JD)
     .,      gpss     ! Start GPS week number
     .,      dt       ! Time between end and start
     
      integer*4 gps_week, gps_dow  !  GPS week and day-of-week
     .,      trimlen
     .,      yr, doy, sec  ! Sinex format time

      data gpss / 2444244.50d0 / 


****  Check the start and end times
      dt = cepoch_end - cepoch_start
*     Normally this time is 1.25 days (+- 3hrs on day)
      if( dt.gt.1.30d0 ) then
*         We need to update
          upd_needed = .true.
          cepoch_start = cepoch_expt - 0.625d0  ! Make consistent
          cepoch_end   = cepoch_expt + 0.625d0  ! Make consistent
*         Now output line to say update made
          sectag = crun_time(6)+crun_time(7)/100.d0
          call ymdhms_to_jd(crun_time,sectag, rt_jd)
          call jd_to_yds( rt_jd, yr, doy, sec)
          gps_week = (cepoch_expt - gpss)/7.d0
          gps_dow  = cepoch_expt - gps_week*7 - gpss
*         Write out update line
          write(*,120) gps_week,gps_dow, yr, doy, sec, 
     .                 hfile(1:trimlen(hfile))
 120      format('CWUUDP cwu',i4.4,I1,' Runtime ',I4.4,':',I3.3,
     .           ':',I5.5,' Hfile ',a)
      end if

****  Thats all
      return
      end


CTITLE CHECK_SVINF

      subroutine check_svinf

      implicit none

*     Routine to fix incorrect assignment of SV number to PRN
*     for satellites that have been re-commissioned.


      include '../includes/kalman_param.h'
      include '../includes/const_param.h' 
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
      include 'hfupd.h'

* PASSED VARIABLES
*     NONE 

* LOCAL VARIABLES
      integer*4 i   ! Loop counter
      integer*4 old_svn(max_glb_svs)   ! Values for SV number in h-file
      logical OK    ! Set false if SVN numbers disagree

      upd_needed = .false.

      qnum_svs = cnum_svs
      do i = 1,cnum_svs
        old_svn(i) =  qsvi_svn(i)
      end do  

*     Now get updated information
      call Get_SVnum( cepoch_expt ) 

*     Test:
      OK = .true.
      do i = 1,cnum_svs
        if( old_svn(i).ne.qsvi_svn(i) ) then
           OK = .false.
           if( qsvs_names(i)(1:3).ne.'PRN' ) then
*              New format
               write(qsvs_names(i),120) qsvs_names(i)(1:1), 
     .                       qsvi_svn(i), mod(qsvi_prn(i),100)
 120           format(a1,I3.3,'_',i2.2)
           else     ! Old format
               write(qsvs_names(i),140) mod(qsvi_prn(i),100),
     .                       qsvi_svn(i)
 140           format('PRN_',i2.2,I2.2)
           endif 

           write(*,160) qsvs_names(i), qsvi_prn(i), old_svn(i), 
     .                  qsvi_svn(i)
 160       format('Updated ',a,' PRN ',i3,' Old SVN ',i3,
     .            ' Correct SVN ',I3)
        endif 
      end do  

      if ( .not.OK ) upd_needed = .true.

      return 
      end




