CTITLE READ_EDIT

      subroutine read_edit

      implicit none

*     This routine will read the hfupd edit/rename file and 
*     save the entries.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include 'hfupd.h'

* LOCAL VARIABLES
* ierr, jerr  -- IOSTAT errors
* j           -- Loop counter
* date(5)     -- Date for getting JDs
* indx, indx_save -- Pointer to position in string
* trimlen     -- Length of string
* ne, te, tr, tn -- Short version for number of edits, total edits,
*                total restores and total renames.

      integer*4  ierr, jerr, j, date(5), indx, indx_save, trimlen,
     .           ne, te, tr, tn, eol 

* sectag      -- Seconds tag for jd conversion
* start, end  -- Start and end JD for commands
      real*8 sectag, start, end

* line        -- Line read from file
* word        -- Command word read from line
* name, newname -- Name and new name for site
* hfcode      -- Hfile code
* cdum        -- Dummy string for rename
      character*256 line
      character*8   word, name, newname, cdum
      character*16  hfcode 

  

****  Open the edit file
      open(101,file=edit_file, iostat=ierr, status='old')
      call report_error('IOSTAT',ierr,'open',edit_file,0,'hfupd')
      if( ierr.ne.0 ) edt_hf = .false.
      if( ierr.ne.0 ) RETURN

* MOD TAH 061110: Initialize the counters
      te = 0
      tr = 0
      tn = 0

***** Now read the edit file
      do while ( ierr.eq.0 )
         read(101,'(a)', iostat=ierr) line
         if( trimlen(line).gt.0 .and. line(1:1).eq.' ' .and.
     .       ierr.eq.0                                     ) then
* MOD TAH 130201: Remove comment part of line.
             eol = index( line,'!' )
             if( eol.eq.0 ) eol = index(line,'#')
             if( eol.gt.0 ) then
                 line(eol:) = ' '
             end if

*            See what we have
             indx = 0
             call getword(line, word, indx)
             call casefold(word)

*            Get the components from the line.  All lines are the
*            same except RENAME has an additional site name
             call getword(line, name, indx)
             call casefold(name)

*            If we have RENAME then read the new name
             if( word(1:3).eq.'REN' ) then
                 call getword(line, newname, indx)
                 call casefold(newname)
             else
                 newname = ' '
             end if
*            Get the Hfile code.  This is optional, so read value
*            and see if numeric
             indx_save = indx
             call getword(line, hfcode, indx)
             if ( trimlen(hfcode).eq.0 ) hfcode = 'ALL'
             call check_num(hfcode,jerr)
             if( jerr.eq.0 ) then
*                Entry is numeric.  So must be date.  Reset the
*                indx and save hfcode as all
                 hfcode = 'ALL'
                 indx = indx_save
             end if

*            See if epoch range passed
             do j = 1, 5
                call read_line(line,indx,'I4',jerr,date(j),cdum)
                if( jerr.ne.0 ) date(j) = 1
             end do
*            If there is an error, default the start to 1900
             if( jerr.ne.0 ) date(1) = 1900
             sectag = 0.0d0
             call ymdhms_to_jd(date, sectag, start)

****         Get the end time:
             do j = 1,5
                 call read_line(line,indx,'I4',jerr,date(j),cdum)
                 if( jerr.ne.0 ) date(j) = 1
             end do
*            If there is an error, default the end to 2100
             if( jerr.ne.0 ) date(1) = 2100
 
*            Process the end data
             sectag = 0.0d0
             call ymdhms_to_jd(date, sectag, end)

****         Add entries to list
             num_hfu_edits = num_hfu_edits + 1
             ne = num_hfu_edits
             call ckmax(ne,max_hfu_edits,'MAX_HFU_EDITS')
             hfu_edit_names(ne) = name
             hfu_edit_rname(ne) = newname
             hfu_edit_code(ne)  = hfcode
             hfu_edit_start(ne) = start
             hfu_edit_end(ne)   = end

*            Now see if we have edit or restore
             if( word(1:2).eq.'ED' .or. word(1:2).eq.'UN') then
                call sbit(hfu_edit_type,ne,1)
             else if ( word(1:3).eq.'RES' ) then
                call sbit(hfu_edit_type,ne,0)
             else if ( word(1:3).ne.'REN' ) then 
                call report_stat('WARNING','HFUPD','read_edit',
     .                           'Unknown option in ',line,0)
             end if

*            Now local sum to user about entries
             if( word(1:2).eq.'ED' .or. word(1:2).eq.'UN' ) te = te + 1
             if( word(1:3).eq.'RES' ) tr = tr + 1
             if( word(1:3).eq.'REN' ) tn = tn + 1
         end if
      end do

****  Tell user current status
      write(*,200) num_hfu_edits, te, tr, tn,
     .             edit_file(1:trimlen(edit_file))
 200  format('There are ',i5,' edit/rename entries, ',I4,' edits,',
     .       I4,' restores, ',i4,' renames in ',a)

****  Thats all 
      close(101)
      return
      end

CTITLE APPLY_EDIT 

      subroutine apply_edit

      implicit none

*     This routine will check the times of this hfile and
*     apply the edit/restores/renames that are needed.

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
      include 'hfupd.h' 


* LOCAL VARIABLES
* i,j, k  -- Loop counters
* trimlen -- Length of string
* indx    -- Pointer in string
* js      -- site/prn number from get_cmd 
* num_codes, code(6)  -- Number of parameter codes to be
*            tested.  For Site site edit/restore this is 6
*            (3 position/3 velocity); for PRN it is only 1
* entry  -- Site or PRN number
* pcode  -- Parameter code to be tested against edit/restore code
* type   -- Type of paramter (51 for PRN, 7-9 and 42-44 for sites)
* ent    -- type of element for PRN (will be 0 for sites)
* npe    -- Counter for number of parameters entries actually edited

      integer*4 i,j, k, trimlen, indx, js, num_codes, code(6),
     .          entry, pcode, type, ent, npe
      integer*4 ndup  ! Number of duplicate site names for internal renames
     .,         lname ! Length of original site name in rename entry

* kbit   -- Checks bit status
      logical kbit

* name    -- Name of site/prn to be edited 
      character*8 name

****  Loop over the edits that we have and see which are applicable
      write(*,120) num_hfu_edits
 120  format(/,'Checking ',i4,' edit/restore/renames',/,
     .         '   # Option   Site     Action')

****  See if we really need any update
      upd_needed = .false.

      do i = 1, num_hfu_edits

*        See if time range and hfcode matchs 
         if( cepoch_start.ge. hfu_edit_start(i) .and.
     .       cepoch_end  .le. hfu_edit_end(i)   .and.
     .      (index(hfile,hfu_edit_code(i)
     .                  (1:trimlen(hfu_edit_code(i)))).gt.0 .or.
     .       hfu_edit_code(i)(1:3).eq.'ALL')        ) then

*            OK: Update is needed
             upd_needed = .true.

*            OK, we are in time range and hfcode is either all or
*            matches the name of the hfile
*            See if we have the station or PRN
             name = hfu_edit_names(i)
             entry = 0
             if( hfu_edit_names(i)(1:4).eq.'PRN_' ) then
*                Its a satellite, so see if we have
                 indx = 0
                 call get_cmd(name, qsvs_names, cnum_svs, 
     .                 js, indx)
                 if( js.gt.0 ) then
*                   We have the satellite.  Set up the parameter
*                   codes that have to be edited
                    num_codes = 1
                    code(1) = 51
                    entry   = js
                 else
                    num_codes = 0
                 end if
             else
*                This is a station name, so check to see if we have
                 indx = 0
                 call get_cmd(name, qsite_names, cnum_sites, 
     .                 js, indx) 
                 if( js.gt.0 ) then
*                   We have the satellite.  Set up the parameter
*                   codes that have to be edited (These are 7,8,9
*                   for positions, and 42,43,44 for velocities)
                    num_codes = 6
                    do j = 1, 3
                       code(j) = 6+j
                       code(j+3) = 41+j
                    end do
                    entry   = js
                 else
                    num_codes = 0
                 end if
             end if
             
*            See what we need to do:
             if( trimlen(hfu_edit_rname(i)).eq.0 ) then
*                We are editing or restoring.  To edit we add 128
*                to the code number for the paramter
                 npe = 0
                 if( num_codes.gt.0 ) then 

*                   Loop over the parameter codes and see which
*                   ones match 

                    do j = 1, cnum_parn
                       call decode_code(qglb_codes(j), pcode, indx)
*                      For satellites we need to pull the PRN number
*                      from the high-order bytes.
                       if( pcode.eq.51 ) then
                           call decode_code(indx, type, ent)
                           indx = ent
                       endif 

*                      If we are restoring then check the parameter
*                      code + 128
                       if( .not.kbit(hfu_edit_type,i) ) then
                           pcode = pcode + 128
                       end if

*                      See if the parameter code matches what we
*                      need
                       do k = 1, num_codes
                           if( pcode.eq. code(k) .and.
     .                         entry.eq.indx ) then
*                             Match found.  There this parameter
*                             code must be changed.
                              if( pcode.lt.128 ) then
                                  qglb_codes(j) = qglb_codes(j) + 128
                              else
                                  qglb_codes(j) = qglb_codes(j) - 128
                              end if
                              npe = npe + 1
                           end if
                       end do
                    end do
                 end if

*                Tell user what is happening
                 if( pcode.lt.128 ) then
                    write(*,420) i, hfu_edit_names(i), npe
 420                format(i4,' Editing  ',a8,1x,i2,
     .                        ' parameters removed')
                 else
                    write(*,425) i, hfu_edit_names(i), npe
 425                format(i4,' Restoring ',a8,1x,i2,
     .                        ' parameters added') 
                 endif

*            This is a rename, so simply change the name of site
             else if( entry.gt.0 ) then
                 lname = trimlen(hfu_edit_rname(i))
                 if( ps_ignore ) lname = 6
                 if( qsite_names(entry)(1:lname).ne.
     .               hfu_edit_rname(i)(1:lname) ) then
                     name = qsite_names(entry)
                     name(1:lname) = hfu_edit_rname(i)(1:lname)
                    write(*,440) i, qsite_names(entry), name
 440                format(i4,' Renaming ',a8,' to ', a8)
                    qsite_names(entry)(1:lname) = 
     .                       hfu_edit_rname(i)(1:lname) 
                 endif 
             end if
*
* MOD TAH 130201: See if there is partial overlap.  Only apply to standard
*        site renames
         elseif(  hfu_edit_start(i).lt. cepoch_end     .and.
     .            hfu_edit_start(i).gt. cepoch_start   .and.
     .            trimlen(hfu_edit_rname(i)).gt.0      .and.   ! Make sure rename (not delete)
     .           (index(hfile,hfu_edit_code(i)
     .                  (1:trimlen(hfu_edit_code(i)))).gt.0 .or.   ! H-file code matched
     .             hfu_edit_code(i)(1:3).eq.'ALL')        ) then
*            Find the station.  Normally there should be two sites with
*            second needing the update, multiple renames o
             name = hfu_edit_names(i)
             indx = 0
*            Try to find name.  If we find it then orginal rename not done can't apply
             call get_cmd(name, qsite_names, cnum_sites, 
     .                 js, indx) 
             if( js.gt.0 .and..not. force_upd ) then
                write(*,510)  i, qsite_names(js), hfu_edit_rname(i)
 510            format(I4,' Cannot rename ',a,' to ',a,' at overlap',
     .                    ' start time because no original rename',/,
     .                    ' Site name left unaltered')
             else
*               Find second site name
                entry = 0
                ndup = 0
                lname = trimlen(hfu_edit_names(i))
                do j = 1, cnum_sites
                   if( qsite_names(j)(1:lname).eq.
     .                 hfu_edit_names(i)(1:lname) ) then
                       ndup = ndup + 1
                       entry = j 
                   end if
                enddo
                if( ndup.eq.2 .or. 
     .              (force_upd .and. ndup.eq.1) ) then ! Correct number
                   lname = trimlen(hfu_edit_rname(i))
                   if( ps_ignore ) lname = 6
                   if( qsite_names(entry)(1:lname).ne.
     .                 hfu_edit_rname(i)(1:lname) ) then
                       name = qsite_names(entry)
                       name(1:lname) = hfu_edit_rname(i)(1:lname)
                       write(*,520) i, qsite_names(entry), name
 520                   format(i4,' Renaming ',a8,' to ', 
     .                                     a8,' : Overlap Start')
                   endif 
                   qsite_names(entry)(1:lname) = 
     .                       hfu_edit_rname(i)(1:lname) 
                elseif( ndup.ne.0 ) then
                   write(*,530)  i, qsite_names(entry), 
     .                           hfu_edit_rname(i), ndup
 530               format(I4,' Cannot rename ',a,' to ',a,' at overlap',
     .                    ' start time because ',i3,' site names')
                end if
             endif           
         elseif(  hfu_edit_end(i).lt. cepoch_end     .and.
     .            hfu_edit_end(i).gt. cepoch_start   .and.
     .            trimlen(hfu_edit_rname(i)).gt.0    .and.   ! Make sure rename (not delete)
     .           (index(hfile,hfu_edit_code(i)
     .                  (1:trimlen(hfu_edit_code(i)))).gt.0 .or.   ! H-file code matched
     .             hfu_edit_code(i)(1:3).eq.'ALL')        ) then
*            Find the station.  Normally there should be two sites with
*            second needing the update, multiple renames o
             name = hfu_edit_names(i)
             indx = 0
*            Try to find name.  If we find it then orginal rename not done can't apply
             call get_cmd(name, qsite_names, cnum_sites, 
     .                 js, indx) 
             if( js.gt.0 .and..not. force_upd ) then
                write(*,610)  i, qsite_names(js), hfu_edit_rname(i)
 610            format(I4,' Cannot rename ',a,' to ',a,' at overlap',
     .                 ' end time because no original rename',/,
     .                 ' Site name left unaltered')
             else
*               Find second site name
                entry = 0
                ndup = 0
                lname = trimlen(hfu_edit_names(i))
                do j = 1, cnum_sites
                   if( qsite_names(j)(1:lname).eq.
     .                 hfu_edit_names(i)(1:lname) ) then
                       ndup = ndup + 1
                       if( ndup.eq.1 ) entry = j ! Save first
                   end if
                enddo
                if( ndup.eq.2 .or. 
     .              (force_upd .and. ndup.eq.1) ) then ! Correct number
                   lname = trimlen(hfu_edit_rname(i))
                   if( ps_ignore ) lname = 6
                   if( qsite_names(entry)(1:lname).ne.
     .                 hfu_edit_rname(i)(1:lname) ) then
                       name = qsite_names(entry)
                       name(1:lname) = hfu_edit_rname(i)(1:lname)
                       write(*,620) i, qsite_names(entry), name
 620                   format(i4,' Renaming ',a8,' to ', 
     .                                     a8,' : Overlap End')
                   endif 
                   qsite_names(entry)(1:lname) = 
     .                       hfu_edit_rname(i)(1:lname) 
                elseif( ndup.ne.0 ) then
                   write(*,630)  i, qsite_names(entry), 
     .                           hfu_edit_rname(i), ndup
 630               format(I4,' Cannot rename ',a,' to ',a,' at overlap',
     .                    ' end time because ',i3,' site names')
                end if
             endif           
 
         end if    ! Experiment in time range


      end do       ! Looping over all of the entries in the edit_file
*            
*     Thats all
      return 
      end
                     

CTITLE RACODE_to_NAME

      subroutine racode_to_name(antcod, ant_name,
     .                          rcvcod, rcv_name, dir)

      implicit none

*     Routine to convert receiver and antenna codes to full     
*     names.  The rcvant.dat file is checked locally and in
*     $HOME/gs/tables or $HOME/gg/tables
*
      include 'hfupd.h'

* PASSED VARIABLES
* rcvcod, rcv_name -- Receiver code (in) and name (out)
* antcod, ant_name -- Antenna code (in) and name (out)
* dir  -- Direction of conversion: TONAME or TOCODE

      character*(*) rcvcod, antcod 
      character*(*) rcv_name, ant_name
      character*(*) dir         

* LOCAL VARIABLES
* ierr, jerr  -- IOSTAT errors
* i   -- Loop counter
* trimlen -- Length of string
* na, nr  -- short version of number of receivers
                                                  
      integer*4  ierr, jerr, i, trimlen, na, nr

* first_call  -- Saved logical to indicate if routine has been
*     called yet

      logical first_call

* file -- Name for rcvant.dat
* home_dir -- Users home directory
* mode     -- Mode (either ant or rcv)
* line     -- Line read from file
* antcod   -- Antenna code after checking aliases
* code, name  -- Code and name read from rcvant.dat

      character*128 file, home_dir, line
      character*6   code
      character*20  name
      character*3   mode

      data first_call  / .true. /


***** If this is first call then read the rcvant file
      if( first_call ) then
          num_antcode = 0
          num_rcvcode = 0
          first_call = .false.

          file = 'rcvant.dat'
          open(104,file=file, iostat=ierr, status='old')
          if( ierr.ne.0 ) then
*             No local copy.   Get the users home directory and try
*             gs/tables 
              call getenv('HOME', home_dir)
              file = home_dir(1:trimlen(home_dir)) // '/gs/tables/' //
     .               file
              open(104,file=file, iostat=ierr, status='old')

* MOD TAH 000901: If not found in gs/tables, try gg/tables.
              if( ierr.ne.0 ) then
                 file = 'rcvant.dat'
                 file = home_dir(1:trimlen(home_dir)) // 
     .                          '/gg/tables/' // file
                 open(104,file=file, iostat=ierr, status='old')
              endif
              call report_error('IOSTAT',ierr,'open',file,0,'HFUDP')
              if( ierr.ne.0 ) then
                  write(*,120) file(1:trimlen(file))
 120              format('Error opening ',a,': Cannot get full ',
     .                   'reciever/antenna names')
                  RETURN
              end if
              write(*,'(a,1x,a)') 'Using rcvant.dat file',file
          end if

*         Now read the file:
          do while ( ierr.eq.0 )
              read(104,'(a)', iostat=ierr) line
              if( line(1:1).eq.' ' .and. trimlen(line).gt.0 .and.
     .            ierr.eq.0 ) then 

*                 See if we have a mode change
                  if( index(line,'RECEIVERS').eq.2 ) then
                      mode = 'rcv'
                  else if( index(line,'ANTENNAS').eq.2 ) then
                      mode = 'ant'
                  else
*                     Decode the line:
                      read(line,210,iostat=jerr) code, name
 210                  format(1x,a6,8x,a20) 
                      if( mode.eq.'rcv' .and.jerr.eq.0 ) then
                          num_rcvcode = num_rcvcode + 1
                          nr = num_rcvcode
                          call ckmax(nr,max_rcvant,'MAX_RCVANT')
                          rcv_codes(nr) = code
                          rcv_names(nr) = name
                      else  if( mode.eq.'ant' .and.jerr.eq.0 ) then
                          num_antcode = num_antcode + 1
                          na = num_antcode
                          call ckmax(na,max_rcvant,'MAX_RCVANT')
                          ant_codes(na) = code
                          ant_names(na) = name
                      else
                          write(*,240) jerr, mode
 240                      format('*** WARNING *** IOSTAT error ',i4,
     .                           ' reading rcvant.dat while in ',a,
     .                           ' mode')
                      end if
                  end if
              end if
          end do
          write(*,300) num_rcvcode, num_antcode
 300      format('Found ',i4,' receiver codes, and ',i4,
     .           ' antenna codes in rcvant.dat')
      end if

****  Now see if we can find a match to the codes passed
*     by user.  See which direction
      if( dir.eq.'TON' ) then
         rcv_name = ' '
         do i = 1, num_rcvcode
            if( rcvcod.eq.rcv_codes(i) ) then
                rcv_name =  rcv_names(i)
            end if
         end do

*        Now do antenna.  Check for aliases first  
         call ant_alias(antcod,antcod) 
         ant_name = ' '
         do i = 1, num_antcode
            if( antcod.eq.ant_codes(i) ) then
                ant_name =  ant_names(i)
            end if
         end do
*      Else connvert name to code (needed for Hisub)
      else  ! connvert name to code (needed for Hisub)
*             (Actually not needed anymore).
         rcvcod = ' '
         do i = 1, num_rcvcode
            if( rcv_name.eq.rcv_names(i) ) then
                rcvcod =  rcv_codes(i)
            end if
         end do

*        Now do antenna.  Check for aliases first  
         antcod = ' '
         do i = 1, num_antcode
            if( ant_name.eq.ant_names(i) ) then
                antcod =  ant_codes(i)
            end if
         end do
       end if

****  Check the results
      if( trimlen(rcv_name).eq.0 ) then
          write(*,400) rcvcod
 400      format('*** WARNING *** No match to receiver code ',a)
          rcv_name = 'RCV' // rcvcod
      end if
      if( trimlen(ant_name).eq.0 ) then
          write(*,410) antcod
 410      format('*** WARNING *** No match to antenna code ',a)
          ant_name = 'ANT' // antcod
      end if
      if( trimlen(antcod).eq.0 ) then
          write(*,420) ant_name
 420      format('*** WARNING *** No match to antenna name ',a)
          antcod = 'NONE'
      end if
          

****  Thats all
      return 
      end

CTITLE CK_RCV_DUP

      subroutine ck_rcv_dup(nr, ke, ns)

      implicit none

*     This routine checks to see if the most recently added
*     receiver entry for a site is a duplicate.  If it is the
*     new entry is removed.

      include 'hfupd.h'

* PASSED VARIABLES
* nr  -- Number of receiver entries
* ke  -- Number of entry at this station
* ns  -- Station number

      integer*4 nr, ke, ns

* LOCAL VARIABLES
* dup  -- Logical set true if entry is a duplicate
      logical dup

****  If this the first entry for a station then it can't be
*     a duplicate so return
      if( ke.eq.1 ) RETURN

*     Now check receiver name
      dup = .true.

*     Ckeck the entries
      if( sxrecv_sty(nr-1).ne. sxrecv_sty(nr) ) dup = .false.
      if( sxrecv_fw(nr-1) .ne. sxrecv_fw(nr)  ) dup = .false.
      if( sxrecv_sn(nr-1) .ne. sxrecv_sn(nr)  ) dup = .false.

****  See if we still have a duplicate
      if( dup ) then
*         Remove the new entry
          sx_num_rec(ns) =  sx_num_rec(ns) - 1
          sxrecv_en(nr-1)  = sxrecv_en(nr)

          nr = nr - 1
      end if

****  Thats all
      return
      end

CTITLE CK_ANT_DUP

      subroutine ck_ant_dup(na, ka, ns)

      implicit none

*     This routine checks to see if the most recently added
*     antenna entry for a site is a duplicate.  If it is the
*     new entry is removed.

      include 'hfupd.h'

* PASSED VARIABLES
* na  -- Number of receiver entries
* ka  -- Number of entry at this station
* ns  -- Station number

      integer*4 na, ka, ns

* LOCAL VARIABLES
* dup  -- Logical set true if entry is a duplicate
      logical dup

* i    -- Loop counter
      integer*4 i

****  If this the first entry for a station then it can't be
*     a duplicate so return
      if( ka.eq.1 ) RETURN

*     Now check receiver name
      dup = .true.

*     Ckeck the entries
      if( sxante_sty(na-1).ne. sxante_sty(na) ) dup = .false.

*     Check the antenna eccentricity offsets
      do i = 1, 3
         if( sxL1a_ecc(i,na-1).ne.sxL1a_ecc(i,na) ) dup = .false.
         if( sxL2a_ecc(i,na-1).ne.sxL2a_ecc(i,na) ) dup = .false. 
      end do

****  See if we still have a duplicate
      if( dup ) then
*         Remove the new entry
          sx_num_ant(ns) =  sx_num_ant(ns) - 1
          sxante_en(na-1)  = sxante_en(na)

          na = na - 1
      end if

****  Thats all
      return
      end

CTITLE CK_ECC_DUP

      subroutine ck_ecc_dup(ne, ke, ns)

      implicit none

*     This routine checks to see if the most recently added
*     eccentricity entry for a site is a duplicate.  If it is the
*     new entry is removed.

      include 'hfupd.h'

* PASSED VARIABLES
* ne -- Number of eccentricity entries
* ke  -- Number of entry at this station
* ns  -- Station number

      integer*4 ne, ke, ns

* LOCAL VARIABLES
* dup  -- Logical set true if entry is a duplicate
      logical dup

* i    -- Loop counter
      integer*4 i

****  If this the first entry for a station then it can't be
*     a duplicate so return
      if( ke.eq.1 ) RETURN

*     Now check receiver name
      dup = .true.

*     Check the antenna eccentricity offsets
      do i = 1, 3
         if( sxarp_ecc(i,ne-1).ne.sxarp_ecc(i,ne) ) dup = .false.
      end do

****  See if we still have a duplicate
      if( dup ) then
*         Remove the new entry
          sx_num_ecc(ns) =  sx_num_ecc(ns) - 1
          sxecce_en(ne-1)  = sxecce_en(ne)

          ne = ne - 1
      end if

****  Thats all
      return
      end

CTITLE READ_STINFO_NEWF

      subroutine read_stinfo_newf

      implicit none

*     Routine to read a new format station.info file and convert the entries to
*     sinex type headers

      include '../includes/kalman_param.h'
      include 'hfupd.h'

* LOCAL VARIABLES
* --------------- 
* ierr    -- IOSTAT error
* trimlen -- Length of string
* ns, nr, ne, na -- Number of stations, receiver entries and antenna 
* ke, ka     -- Number of receiver and antenna entries at specific
*     stations
* isec       -- Seconds of day 

      integer*4 ierr, jerr, trimlen, 
     .          ns, nr, na, ne, kr, ke, ka, isec, i, j

* jd         -- Julian date
      real*8 jd

* line -- Line read from file
* file   -- Name for antmod.dat/rcvant.dat
* home_dir  -- Users home directory
* prev_site -- Previous site name

      character*256 line, file, home_dir
      character*4   prev_site 

* STATION.INFO Gamit variables
* ----------------------------
* trkcod, sitcod  -- Tracker and site codes
* stanam          -- Long station name
* anthgt          -- Antenna height
* offstn, offste  -- Offsets north and east
c antdaz          -- Antenna deviation from true North (Alignment from True N in log)
c                    Added TAH 200203 for repro3
* rcvcod, antcod, hgtcod -- Codes from receiver, antenna and height
*                    measurements 
* iyr, idoy -- YR and DOY for start of antenna information
* isessn   -- Session number (ignored)
* istart(5)  I*4   start time (yr doy hr min sec) for values (yr 9999 if not found)
* istop(5)   I*4   stop time  (yr doy hr min sec) for values 

      character*4  trkcod, sitcod
      character*16 stanam
      character*6  rcvcod, antcod
      character*5  hgtcod 
      real*4  swvers
      real*8 anthgt, offstn, offste
      real*8 antdaz  ! Antenna aligment from True N (deg).
      integer*4 iyr, idoy, isessn, istart(5), istop(5)

* Additional variables needed for hisub 
* ipcv_unit  -- Unit number for antmod.dat
* offset_L1(3), offset_L2(3)  -- UNE offsets for L1 and L2
* dhpab    -- Direct height to preamp base
* amodel   -- Antenna phase center model
* rcvnam, antnam -- Full names for receivers and antenna
* MOD RWK 051005: variables need for hisub2 and get_antpcv  (eventually remove the above)
* Additional variables needed for hisub 
* ipcv_unit  -- Unit number for antmod.dat
* offset_L1(3), offset_L2(3)  -- UNE offsets for L1 and L2
* dhpab    -- Direct height to preamp base
* amodel   -- Antenna phase center model
* warnings -- print (T) or not (F) warnings for missing antenna info
* MOD RWK 051005: variables need for hisub2 and get_antpcv  (eventually remove the above)
* ihi_unit  -- Unit number for hi.dat                     
* ipcv_unit -- Unit number for antmod.dat        
* jdi -- integer PEP JD  
* offarp(3) -- UNE offsets of ARP from monument   ! Removed TAH 111202
* offset_L1(3), offset_L2(3)  -- UNE offsets for L1 and L2
* antmod -- requested antenna model
* amodel -- Antenna phase center model
* newant -- initial call to read antmod.dat
* radome_sub -- T if radome='none' subsituted for mising antenna/radome combination antmod.dat
* warnings -- print (T) or not (F) warnings for missing antenna info
* minelev  -- Minumum elevation angle for antenna model tables (from get+antpcv). Not used here.

      integer*4 ihi_unit, ipcv_unit, jdi  
      real*8 offset_L1(3), offset_L2(3), dhpab, rdum, minelev
      character*16 amodel
      character*5 radome             
      character*4 antmod
      logical warnings, radome_sub, newant

****  Variables needed for new format station.info
      real*8 arp_off(2)   ! height and radius this htcod wiht specific antenna
      integer*4 nlist
      integer*4 doy_in    ! Input day number (set zero here since are not using)
      character*6 item_list(20)
      character*20 values(20)
* MOD TAH 200317: Made length consistent with GAMIT (change 24 to 32).
      character*32 comment  

      character*30 rcvnam 
      character*20 antnam
      character*20 rcvrsn, antsn, rcvers

      character*204 stdhead  ! Station.info standard header line

      logical stdstinfo   ! Set true if standard format that can 
                          ! be read directly.
      

      stdstinfo = .false.

****  To use station.info we need two ancillary files: antmod.dat
*     and rcvant.dat  Open and read these files.   
*MOD RWK 051005: Now need a third file: hi.dat
      file = 'antmod.dat'
      open(102,file=file, iostat=ierr, status='old')
      if( ierr.ne.0 ) then
*         No local copy.   Get the users home directory and try
*         gs/tables 
          call getenv('HOME', home_dir)
          file = home_dir(1:trimlen(home_dir)) // '/gs/tables/' //
     .           file
          open(102,file=file, iostat=ierr, status='old')
* MOD TAH 000901: If not found in gs/tables, try gg/tables
          file = 'antmod.dat'
          if( ierr.ne.0 ) then
             file = home_dir(1:trimlen(home_dir)) // '/gg/tables/' //
     .              file
             open(102,file=file, iostat=ierr, status='old')
          endif 
          call report_error('IOSTAT',ierr,'open',file,0,'HFUDP')
          if( ierr.ne.0 ) then
              write(*,120) file(1:trimlen(file))
 120          format('Error opening ',a,': Cannot use station.info')
              RETURN
          end if
      end if
      write(*,'(a,1x,a)') 'Using PCV model',trim(file)
      ipcv_unit = 102     

      file = 'hi.dat'
      ihi_unit = 103
      close(ihi_unit)
      open(ihi_unit,file=file, iostat=ierr, status='old')  
***   See if open OK
      if( ierr.ne.0 ) then
*       MOD TAH 070823: Use the gg/tables version
        call getenv('HOME', home_dir)
        file = home_dir(1:trimlen(home_dir)) // '/gg/tables/' //
     .         file
        open(ihi_unit,file=file, iostat=ierr, status='old')
*       See if we still have problem
        if( ierr.ne.0 ) then  
           ihi_unit = 0
           write(*,'(a)') 'No hi.dat file, revert to hisub'
           ierr = 0 
        else
           write(*,'(a,1x,a)') 'Using hi.dat from tables',file
        endif
      else
        write(*,'(a,1x,a)') 'Using hi.dat from local ',file
      endif 

***** Report using station.info
      write(*,210) igs_snx_file(1:trimlen(igs_snx_file))
 210  format('Reading station.info file (new format): ',a)

*     Get the file format 
      call read_stnfohead( 100,line ) 

***   See if standard format
      stdhead = '*SITE  Station Name      Session Start' // 
     .          '      Session Stop       Ant Ht   HtCod' // 
     .          '  Ant N    Ant E    Receiver Type         Vers' //
     .          '                  SwVer  Receiver SN' //
     .          '           Antenna Type     Dome   Antenna SN' // 
     .          '           AntDAZ'

      if( line(1:204).eq.stdhead(1:204) ) then              
          stdstinfo = .true.
          write(*,'(a)') 'Reading standard format file'
      else
         call get_tokens( line,nlist,item_list )  
      end if 
         
      prev_site = ' '
      ns = 0
      na = 0
      ne = 0
      nr = 0

*     Now loop over file
      do while ( ierr.eq.0 )
         read(100,'(a)',iostat=ierr) line  
         if( line(1:1).eq.' ' .and. trimlen(line).gt.0 .and.
     .       ierr.eq.0 ) then

             warnings = .false.
             warnings = .true.
             antnam = ' '
             rcvnam = ' '      
             if( stdstinfo ) then
                call decode_stdstnf(line,  doy_in  
     .                         , sitcod,  stanam
     .                         , anthgt, offstn, offste, antdaz, rcvcod
     .                         , antcod,  rcvnam, rcvrsn 
     .                         , rcvers, antnam, antsn
     .                         , hgtcod, radome, swvers
     .                         , istart, istop, warnings) 
              else
*                New format decoding of line
                 call read_values(line,nlist,item_list,values,comment)  
                 doy_in = 0    ! Value used only for reporting
*   
* MOD TAH 050709: Initialize names.  These are declared C*15 and C*20
*                in decode_values and therefore will not be fully set. 
* MOD RWK 050719: Add warnings variable to call and set false to avoid report_stat calls
*                 reading all of station.info   
* MOD TAH 100617: Updated call to decode_values in gamit/lib/rstnfo.f
                 call decode_values( nlist, item_list, values, doy_in  
     .                         , sitcod,  stanam
     .                         , anthgt, offstn, offste, antdaz, rcvcod
     .                         , antcod,  rcvnam, rcvrsn 
     .                         , rcvers, antnam, antsn
     .                         , hgtcod, radome, swvers
     .                         , istart, istop, warnings) 

             end if

****         See we need antenna code for hisub call
             if( trimlen(antcod).eq.0 ) then
                 call racode_to_name(antcod,antnam,rcvcod,rcvnam,'TOC')
             end if

             iyr  = istart(1)
             idoy = istart(2)
             isessn = 0  
*            Convert the date to julian date 
             isec = istart(3)*3600.d0 + istart(4)*60 +
     .              istart(5)
             call yds_to_jd( iyr, idoy, isec, jd )  
   
c            get_antpcv needs integer PEP JD, not true Julian date
             jdi = idint(jd+0.5)

****         Call hisub to get the antenna conversions
             if( ierr.eq.0 ) then
* MOD TAH 030731: Removed the 1 from call after ipcv_unit. 
* MOD TAH 030731: Changed DHPAB to DHARP
                 if( hgtcod(1:5).eq.'DHPAB' ) hgtcod = 'DHARP'


*                Only call hisub if know the antenna type
                 if( antcod.ne. '     ' .and. 
     .               antcod.ne. '-----' ) then  


* MOD RWK 050719:  Add warnings argument and set to false for this routine        
                   warnings = .false. 
* MOD RWK 051005:  Temporarily use successful open of hi.dat to determine hisub call      
* MOD RWK 080522:  Remove ability to use old-style hisub
* MOD TAH 111102:  Add call to specifically read hisub file using.  This routine
*                  will return the direct offset and radius for this antenna
                  call hisub_save(ihi_unit, antnam, hgtcod, arp_off)

                  call antpcv_save( ipcv_unit, antnam, radome, 
     .                              offset_L1,offset_L2, 
     .                              radome_sub)
                 else
                    do i = 1,3 
                       offset_L1(i) = 0
                       offset_L2(i) = 0
                    enddo
                 endif

*                See if have a new site
                 if( sitcod.ne.prev_site ) then
*                    New site.  Add site to list of sites
                     ns = ns + 1
                     call ckmax(ns, max_snx_sites,'MAX_SNX_SITES')
! MOD TAH 120823: Remove end of name to make just 4 characters
                     sx_site_code(ns) = sitcod  ! // '_GPS'
                     sx_site_long(ns) = stanam
                     sx_site_dome(ns) = '-----M---'
                     do j = 1,3
                        sx_site_gpos(j,ns) = 0.d0
                     end do
*                    Initialise the counters for this sites
                     sx_num_rec(ns) = 0
                     sx_num_ant(ns) = 0
                     sx_num_ecc(ns) = 0
                     prev_site = sitcod
                 end if

*                Still same site.  See if the receiver code has
*                code
                 kr = sx_num_rec(ns)

*                Add the receiver information
                 nr = nr + 1 
                 call ckmax(nr, max_snx_recs,'MAX_SNX_RECS')
                 kr = kr + 1 
                 sx_rec_rec(kr,ns) = nr
                 call ckmax(kr, max_snx_ent_per_site ,
     .                      'MAX_SNX_ENT_PER_SITE' )
                 sx_num_rec(ns) = kr

                 sxrecv_sty(nr) = rcvnam(1:20)
                 sxrecv_sn(nr)  = rcvrsn
                 write(sxrecv_fw(nr),310) rcvcod, swvers
 310             format(a3,1x,f6.2)
                 sxrecv_st(nr)  = jd
                 sxrecv_en(nr)  = jd_zero
                 if( kr.gt.1 ) sxrecv_en(nr-1) = jd 
                 call ck_rcv_dup(nr,kr,ns)

*****            Now check the antenna type and eccentricity.
                 ka = sx_num_ant(ns)
                 ke = sx_num_ecc(ns)
                 na = na + 1
                 ne = ne + 1

                 ka = ka + 1
                 call ckmax(ka, max_snx_ent_per_site ,
     .                      'MAX_SNX_ENT_PER_SITE' )
                 ke = ke + 1
                 call ckmax(ke, max_snx_ent_per_site ,
     .                      'MAX_SNX_ENT_PER_SITE' )
                                        
                 sx_rec_ant(ka,ns) = na
                 sx_rec_ecc(ke,ns) = ne
                 sx_num_ant(ns) = ka
                 sx_num_ecc(ns) = ke

                 sxante_sty(na) = antnam(1:16) // radome
                 sxante_sn(na)  = antsn 
                 sxante_st(na)  = jd
                 sxante_en(na)  = jd_zero
                 if( ka.gt.1 ) sxante_en(na-1) = jd
                 sxecce_st(ne)  = jd
                 sxecce_en(ne)  = jd_zero
                 if( ke.gt.1 ) sxecce_en(ne-1) = jd    

*                Now save the eccentricity and L1/L2 offsets   
*MOD RWK 051005: change to use hisub2 variables (remove hisub and its variables eventually)
*MOD TAH 111102: change to use with new hisub_save
                 sxL1a_ecc(1,na) = offset_L1(1) 
                 sxL1a_ecc(2,na) = offset_L1(2)                          
                 sxL1a_ecc(3,na) = offset_L1(3) 

                 sxL2a_ecc(1,na) = offset_L2(1) 
                 sxL2a_ecc(2,na) = offset_L2(2)                           
                 sxL2a_ecc(3,na) = offset_L2(3)  

                 sxarp_ecc(1,ne) = offstn   ! offarp(2)
                 sxarp_ecc(2,ne) = offste   ! offarp(3)
                 if ( arp_off(1).eq.0 .and. arp_off(2).eq.0 ) then
                     sxarp_ecc(3,ne) = anthgt  
                 else
                     sxarp_ecc(3,ne) = sqrt((anthgt**2)-(arp_off(2)**2))
     .                               - arp_off(1)
                 endif

                 call ck_ant_dup(na, ka, ns)
                 call ck_ecc_dup(ne, ke, ns)
             end if
         end if
      end do

****  Now save the numbers of entries
      num_snx_sites = ns
      sx_ent_rec = nr
      sx_ent_ant = na
      sx_ent_ecc = ne
      write(*,510) num_snx_sites, trim(igs_snx_file)
 510  format('Total of ',i6,' sites found in ',a) 
      write(*,520) sx_ent_rec
 520  format('Total of ',i6,' Receiver entries found')
      write(*,530) sx_ent_ant
 530  format('Total of ',i6,' antenna entries found')
      write(*,540) sx_ent_ecc
 540  format('Total of ',i6,' eccentricity entries found')

      close(102)

****  Thats all
      return
      end 


CTITLE decode_stdstnf

      subroutine decode_stdstnf(line,  doy_in  
     .                         , sitcod,  stanam
     .                         , anth, antn, ante, antdaz, rcvcod
     .                         , antcod,  rctype, rcvrsn 
     .                         , rcvers, anttyp, antsn
     .                         , htcod, radome, swver
     .                         , istart, istop, warnings) 

      implicit none

*     Routine to read standard station.info format rather than using
*     the general decodings. This method is being used to speed up
*     the reading of this fire, particularly when station.info has
*     a large number of lines.

c       station.info variables
                                
      character*4 sitcod
      character*5 htcod,radome,aswver
      character*6 rcvcod,antcod
      character*20 rcvrsn
      character*16 anttyp  ! Dimensioned it 1 larger to avoid \0 
      character*16 stanam
      character*20 rctype,rcvers,antsn,anttyp20
      real*4 swver
      real*8 anth,antn,ante 
      real*8 antdaz  ! Antenna aligment from True N (deg).
      integer*4 isessn,istart(5),istop(5)
      character*(*) line   ! Line read from station.info
               
c       other variables
         
      logical warnings                                                
      character*1 pcncod
      character*5 char5
      character*6 char6
      character*20 char20       
      character*80 prog_name
      character*256 message    
      integer*4 doy_in,rcpar,len,ioerr,i,j 

***** Initialize values 
      stanam = ' '  
      isessn = 1
      anth = 0.0d0
      antn = 0.0d0
      ante = 0.0d0 
      antdaz = 0.0d0
      rcvcod = ' ' 
      rctype = ' ' 
      rcvrsn = ' ' 
      swver = 0.0
      rcvers =  ' '
      antcod = ' ' 
      anttyp = ' '   
      htcod = ' ' 
      antsn = ' '  
      radome = 'UNKN '
      htcod = 'DHARP'
      do i=1,5 
        istart(i) = 0
        istop(i) = 0
      enddo


****  Directly read the line using a fixed format
*SITE  Station Name      Session Start      Session Stop       Ant Ht   HtCod  Ant N    Ant E    Receiver Type         Vers                  SwVer  Receiver SN           Antenna Type     Dome   Antenna SN           AntDAZ
*0001  GEONET0001        2011  60  0  0  0  9999 999  0  0  0   0.0000  DHARP   0.0000   0.0000  TRIMBLE NETR9         Nav 4.17 Sig 0.00      4.17  --------------------  TRM29659.00      GSI    --------------------    0.
      read(line,120,iostat=ioerr) sitcod, stanam, istart,istop, 
     .        anth, htcod, antn,ante, rctype, rcvers,  aswver, rcvrsn,
     .        anttyp, radome, antsn, antdaz
 120  format(1x,a4,2x,a16,2x,I4,1x,I3,3(1x,I2),2x,I4,1x,I3,3(1x,I2),
     .             2x,F7.4,2x,A5,2x,F7.4,2x,F7.4,2x,A20,2x,a20,  
     .             2x,a5,2x,a20,2x,a15,2x,a5,2x,a20,1x,F5.0)  

***** Now test and compare as is done in rstnfo.f/decode_values
      if( aswver(1:2).eq.'--' ) then 
          swver = 0.0
      else
          read(aswver,'(f5.2)',iostat=ioerr) swver 
      endif

      if( radome.eq.'     '.or.radome.eq.'-----') then
          radome = 'UNKN '      
      else      
c rwk 101111: Make sure 'NONE' is uppercase (not always true in SOPAC file)
          call uppers(radome)   
      endif 

c     assign corresponding RINEX/GAMIT receiver, antenna, and firmware tokens   
!     if((rcvcod.eq.'      '.or.rcvcod.eq.'------') .and.
!    . (rctype(1:5).eq.'     '.or.rctype(1:5).eq.'-----') ) then
!       if( warnings ) then
!         write(message,'(a,a4,i5,i4,3i3,a)') 
!    .      'Neither RcvCod nor Receiver Type found for ',sitcod
!    .      ,(istart(i),i=1,5),' in station.info'
!         call report_stat('WARNING',prog_name,'lib/rstnfo',' '
!    .                    ,message,0)  
!       endif
!     elseif (rctype(1:5).eq.'     '.or.rctype(1:5).eq.'-----') then
!       call read_rcvant(1,2,char6,char20,char5,rcvcod,rctype,pcncod)
!     elseif (rcvcod.eq.'      '.or.rcvcod.eq.'------') then   
!       call read_rcvant(2,2,char6,char20,char5,rcvcod,rctype,pcncod) 
!     else
!       call report_stat('FATAL',prog_name,'lib/rstnfo',' '
!    .    ,'Receiver entry inconsistency in station.info',0)
!     endif   

****  Convert the full names to codes
*     Antenna to CODE 
!     anttyp20 = anttyp
!     call read_rcvantN(2,1,antcod,anttyp20,radome,char6,char20
!    .                    ,pcncod) 
*     Receiver to CODE
!     call read_rcvantN(2,2,char6,char20,char5,rcvcod,rctype,pcncod) 
                 

***** Thats all
      return
      end  
    
      subroutine read_rcvantN( dirn,type,antcod_in,antnam,radome
     .                      , rcvcod,rcvnam,pcncod )
c
c PURPOSE: This subroutine will exchange the antenna code for a full antenna name
c          and vice versa. The definitive place of Gamit where all antennas and receivers
c           are listed is in a data file rcvant.dat in gamit/tables. The
c          intention is that all antennas/receivers ever used will be uniquely identified by
c          their 6 character ant/rec code. Thus, passing a 6 char  code
c          through all gamit routines will suffice, as it can be turned into a
c          complete name at any time by calling this routine.
c
c          This routine will open the file rcvant.dat, read it for the required information
c          then close it again.
c
c VARIABLES:
c     dirn     - direction for information transfer: 1: code name  -> full name
c                                                    2: full name -> code cod   I*4
c     type     - flag for receiver or antenna request 1: antenna; 2: receiver   I*4
c     antcod   - antenna code used in gamit                                     C*6
c     antnam   - full antenna name                                              C*20  
c     radome   - radome name (same for gamit and IGS)                           C*5
c     rcvcod   - antenna code used in gamit                                     C*6
c     rcvnam   - full antenna name                                              C*20 
c     pcncod   - character indicating need for C1 or P2' bias correction        C*1
c                ( see gamit/makex/set_dcb_flag.f )


c P. Tregoning   12th September, 1995
c S. McClusky    20th October, 1995  
c R. King         5th February 2003
c T. Herring     10th November 2003: Fixed up opening rcvant.dat

      implicit none

      integer*4 dirn,type,i,ioerr,irec,indx,len,rcpar
                                    
      character*1 pcncod    
      character*5 radome
      character*6 antcod_in,antcod,rcvcod
      character*20 antnam,rcvnam,antenna(2),receiver(2)
      character*80 line,answer,prog_name
      character*256 message
* NOD TAH 031110: Addeded home_dir so that <home_dir>/gg/tables/rcvant.dat can be opened
*     if their is not a local copy.  Also need trimlen for length of string
      integer*4 trimlen
      character*256 home_dir, home_rcvant

      logical found,stopflag
      parameter (irec = 90)
      logical first_call

      data first_call / .true. /


c     get the module name for report_stat calls
      len = rcpar(0,prog_name)   
c     if mstinf2, just issue warning
      if( prog_name(1:6).eq.'mstinf'.or.
     .    prog_name(1:5).eq.'hfupd' .or.
     .    prog_name(1:6).eq.'conver' ) then
        stopflag = .false.
      else
        stopflag = .true.
      endif                                     
               
c   initialize controls and arrays
      ioerr = 0   
      answer = ' '
      found = .false.
      i = 1

c      print *,'READ_RCVANT antcod antnam radome '
c     .       ,antcod_in,antnam,radome
c       write(*,10) dirn,type,antcod_in,antnam,radome,rcvcod,rcvnam
c10     format(i2,i2,1x,a6,2x,a20,1x,a5,1x,a6,1x,a20)

c  open the file rcvant.dat
* MOD TAH 031110: Close unit 90 first.  It is opened someplace else and
*     this cause error on HPUX system
c      close(unit=90,iostat=ioerr)
      if( first_call ) then 
         open(unit = irec,file='rcvant.dat',status='old',iostat=ioerr)     
         if(ioerr.ne.0) then
* MOD TAH 031110: See if we can open a copy through the <home_dir>/gg/tables
*          Get the user's home directory:
           call getenv('HOME',home_dir)
           home_rcvant = home_dir(1:max(1,trimlen(home_dir))) //
     .                '/gg/tables/rcvant.dat'
           open(unit=irec,file=home_rcvant,status='old',iostat=ioerr)
         end if
         if(ioerr.ne.0 ) then
           write(message,'(2a)')
     .        'Cannot open file rcvant.dat'
     .        ,'-- may need link to gamit/tables/rcvant.dat'
           call report_stat('FATAL',prog_name,'lib/read_rcvant',' '
     .                  ,message,ioerr)
         endif
      endif
     
c now match up the receiver/antenna information as requested by the arguments
c passed to the subroutine

      if ( type.eq.1 ) then
c        it is an antenna request

c get the full RINEX name from the GAMIT code
        if( dirn.eq.1 ) call ant_alias(antcod_in,antcod)

c move through the file until we find the 'ANTENNA' block
15      do while(answer(1:7).ne.'ANTENNA')
          read(irec,'(a80)',iostat=ioerr) line
          if(ioerr.ne.0)  call report_stat('FATAL',prog_name
     .      ,'lib/read_rcvant',' '
     .      ,'Error looking for ANTENNA in rcvant.dat',ioerr)
          if (line(1:1) .ne. ' ' ) goto 15
          indx = index(line,'ANTENNA')
          read(line(indx:indx+6),'(a7)',iostat=ioerr) answer   
          if(ioerr.ne.0)  call report_stat('FATAL',prog_name
     .      ,'lib/read_rcvant',' '
     .      ,'Error reading line for ANTENNA in rcvant.dat',ioerr)
        enddo

20      do while (.not.found)
          read(irec,'(a)',iostat=ioerr) line   
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/read_rcvant',' '
     .       ,'Error reading antenna lines in rcvant.dat',ioerr)
          if (line(1:1) .ne. ' ' ) goto 20
          if (line(2:4) .eq. 'END' ) then  
            if( dirn.eq.1 )  write(message,'(a,a6,a,a6,a)') 
     .         'Input antenna type ',antcod_in,' with alias '
     .          ,antcod,' not in rcvant.dat' 
            if( dirn.eq.2 )  write(message,'(a,a20,a)')
     .         'Antenna name ',antnam,' not found in rcvant.dat'  
            if( stopflag ) then 
               call report_stat('FATAL',prog_name,'lib/read_rcvant'
     .                       ,' ',message,0)     
             else
               call report_stat('WARNING',prog_name,'lib/read_rcvant'
     .                       ,' ',message,0)    
               if( dirn.eq.1 ) antnam = ' '
               if( dirn.eq.2 ) antcod_in = ' '     
!               close(irec)
               rewind(irec)
               return
             endif 
          endif

          read(line,'(1x,a6,8x,a20)',iostat=ioerr) antenna(1),antenna(2)
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .     ,'lib/read_rcvant',' ','Error in antenna table format',ioerr)


          if (dirn.eq.1)then

c  we want to convert codnam into full name
            if(antenna(1)(1:6).eq.antcod) then
              found = .true.
              antnam = antenna(2)               
              if( radome(1:4).eq.'UNKN' .or. radome(1:4).eq.'DOME') then
c               check for radome types built into GAMIT codes
                if( antcod.eq.'ASHDMR' ) radome(1:4)='DOME'
                if( antcod.eq.'ASHDMD' ) radome(1:4)='DOME'
                if( antcod.eq.'LC303R' ) radome(1:4)='DOME'
                if( antcod.eq.'NO503R' ) radome(1:4)='DOME'
                if( antcod.eq.'TPSC3R' ) radome(1:4)='SNOW'
                if( antcod.eq.'TPSCC4' ) radome(1:4)='SNOW'
                if( antcod.eq.'AERCHR' ) radome(1:4)='DOME'
              endif
              antnam(17:20) = radome(1:4) 
c              print *,'found antcod antnam ',antcod,antnam
            endif

          elseif(dirn.eq.2)then

c  want to convert full name into codnam
            if(antenna(2)(1:16).eq.antnam(1:16)) then
              found = .true.
              antcod_in = antenna(1)(1:6)       
c             for a few cases, assign the radome from the name
              if( antnam.eq.'TPSCR3_GGD SNOW ') radome='SNOW '
              if( antnam.eq.'TPSCR4 SNOW     ') radome='SNOW '
c put the radome name in into the end of the antenna name slot
              antnam(17:20) = radome(1:4)            
            endif
          endif
        enddo

!        close(irec) 
        rewind(irec) 
        return

      elseif( type .eq. 2 ) then

c  it is a receiver request

c move through the file until we find the 'RECEIVER' block
120     do while(answer(1:8).ne.'RECEIVER')
          read(irec,'(a)',iostat=ioerr) line   
          if( ioerr.ne.0 )  call report_stat('FATAL',prog_name
     .        ,'lib/read_rcvant',' '
     .         ,'Error looking for RECEIVER in rcvant.dat file ',ioerr)
          if (line(1:1) .ne. ' ' ) goto 120
          indx = index(line,'RECEIVER')
          read(line(indx:indx+7),'(a8)',iostat=ioerr) answer  
          if( ioerr.ne.0 )  call report_stat('FATAL',prog_name
     .       ,'lib/read_rcvant',' '
     .       ,'Error reading line for RECEIVER in rcvant.dat ',ioerr)
        enddo

130     do while (.not.found)
          read(irec,'(a)',iostat=ioerr) line     
          if( ioerr.ne.0 ) call report_stat('FATAL',prog_name
     .       ,'lib/read_rcvant',' '
     .         ,'Error reading receiver lines in rcvant.dat',ioerr)
          if (line(1:1) .ne. ' ' ) goto 130
          if (line(2:4) .eq. 'END' ) then             
            if(dirn.eq.1) write(message,'(a,a6,a)') 'Receiver code '
     .               ,rcvcod,' not found in rcvant.dat'
            if(dirn.eq.2) write(message,'(a,a20,a)') 'Receiver name '
     .               ,rcvnam,' not found in rcvant.dat'   
            if( stopflag ) then
              call report_stat('FATAL',prog_name,'lib/read_rcvant',' '
     .                      ,message,0)
            else  
              if( dirn.eq.2.and.rcvnam(1:6).ne.'------') 
     .        call report_stat('WARNING',prog_name,'lib/read_rcvant',' '
     .                      ,message,0) 
             if( dirn.eq.1 ) rcvnam = ' '
             if( dirn.eq.2 ) rcvcod = ' ' 
!             close(irec) 
             rewind(irec)
             return  
            endif
          endif                              

          read(line,'(1x,a6,8x,a20,1x,a1)',iostat=ioerr) 
     .             receiver(1), receiver(2),pcncod 
          if(ioerr.ne.0) call report_stat('FATAL',prog_name
     .       ,'lib/read_rcvant',' '
     .         ,'Error in receiver tables format',ioerr)

          if (dirn.eq.1)then

c  we want to convert codnam into full name
            if(receiver(1)(1:6).eq.rcvcod)then
              found = .true.
              rcvnam = receiver(2)
            endif

          elseif(dirn.eq.2)then

c  want to convert full name into codnam
            if(receiver(2).eq.rcvnam)then
              found = .true.
              rcvcod = receiver(1)(1:6)
            endif
          endif
        enddo

!        close(irec)
        rewind(irec)
        return

      else

        call report_stat('FATAL',prog_name,'lib/read_rcvant',' '
     .                  ,'Unknown request type',0)

      endif
      end



CTITLE HISUB_SAVE 

      subroutine hisub_save(ihi_unit, antnam, hgtcod, arp_off)

      implicit none 

*     This routine will read hisub.dat and save the horizontal 
*     and vertical offsets of the antenna. With each call it
*     will check whether the antenna  and height type  has been 
*     previously seenand use the value from the saved list.
*     This should substantially speed up these routine
*     over the use of the standard hi_sub routine

      include 'hfupd.h'

* PASSED VARIABLES
      integer*4 ihi_unit   ! Unit number for hisub.dat file.
      character*(*) antnam ! Full name of antenna  
      character*(*) hgtcod ! Height code.  If DHARP then zeros are
                           ! returned and values not saved. 
      real*8 arp_off(2)    ! Vertical offset and horizontal radius (m)

* LOCAL VARIABLES
      integer*4 i, n       ! Loop counters
      integer*4 ierr, jerr ! IOSTAT errors
      integer*4 trimlen    ! Function

      character*256 line   ! Line read from file
      character*16 anttyp  ! Antenna type full name
      character*6  antcod  ! GAMIT code (not used)
      character*5  htcod   ! Height code read.
      real*8 vert, horiz   ! Vertical and radius for ARP (m)

      logical first_call   ! Set false after first call
      logical found        ! Set true when match found

      data first_call  / .true. /

      save first_call

****  Start by seeing if this DHARP.  If so simply return
      arp_off(1) = 0
      arp_off(2) = 0
      if( hgtcod(1:5).eq.'DHARP' ) then
          RETURN
      end if

****  OK, now scan existing entries to see if have already
      n = 0
      if( first_call .and. ihi_unit.gt.0 ) then

*         Readin all the entries and save
          ierr = 0
          do while ( ierr.eq.0 )
             read(ihi_unit,'(a)',iostat=ierr) line
             if( ierr.eq.0 .and. trimlen(line).gt.0 .and.
     .           line(1:1).eq.' ' ) then
                 read(line,'(1x,a16,4x,1x,a6,1x,a5,2f7.2)',iostat=jerr)
     .               anttyp,antcod,htcod,vert,horiz 
                 call report_error('IOSTAT',jerr,'decod',line,
     .               1, 'hisub_save')

                 n = n + 1
                 call ckmax(n, max_rcvant,'MAX_APR_HTSAVES' )
                 arp_anthtcod(n) = anttyp(1:16) // htcod
                 arp_offsave(1,n) = horiz
                 arp_offsave(2,n) = vert
             end if
          end do
          first_call = .false.
          num_arpcode = n
          close(ihi_unit)
          write(*,120) num_arpcode
 120      format('Read hisub.dat: ',I4,' entries found')
      end if

****  Now try to match to values passed.
      found = .false.
      do i = 1, num_arpcode
         if( antnam(1:16).eq.arp_anthtcod(i)(1:16) .and.
     .       hgtcod(1:5) .eq.arp_anthtcod(i)(17:21) ) then
*            Match found
             arp_off(1) = arp_offsave(1,i)
             arp_off(2) = arp_offsave(2,i)
             found = .true.
             exit
         endif
      end do
* 
      if( .not.found ) then
         write(*,220)  antnam(1:16),  hgtcod(1:5)
 220     format('No match found in hisub.dat  for ',a,1x,a)
*        Add this entry to table so that only reported once
         num_arpcode = num_arpcode+1
         call ckmax(num_arpcode, max_rcvant,'MAX_APR_HTSAVES_missing')
         arp_anthtcod(num_arpcode)  = antnam(1:16) // hgtcod
         arp_offsave(1,num_arpcode) = 0
         arp_offsave(2,num_arpcode) = 0 
      end if

****  Thats all

      end

CTITLE ANTPCV_SAVE

      subroutine antpcv_save( ipcv_unit, antnam, radome, 
     .           offset_L1,offset_L2, radome_sub)

      implicit none 

*     This routine will read the phase center model file and save 
*     the L1 and L2 offsets. 
*     With each call, the antenna tape is checked against previous 
*     antennas found, in the previous values will be used if the antenna
*     has been seen  before. Only the phase center offsets are saved
*     because these are the only values that can be checked. This 
*     method should be much faster in reading the standard antenna file.

      include 'hfupd.h'


* PASSED VARIABLES
      integer*4 ipcv_unit  ! Unit number for antmod.dat
      character*(*) antnam ! Antenna name (16 characters used)
      character*(*) radome ! Radome (4-characters used).

      real*8 offset_L1(3),offset_L2(3)
      logical radome_sub   ! Set true if radome not found and NONE used

* LOCAL VARIABLES
      integer*4 i, n       ! Loop counters
      integer*4 ierr, jerr ! IOSTAT errors
      integer*4 trimlen    ! Function


      character*256 line   ! line read from file
      character*4 svstr    ! G/R number for satellite

      logical first_call   ! Set false after first call
      logical found        ! Set true when match found
      logical read_offs	   ! Set true when offset should be read
      logical L1_true, L2_true  ! Each true when L1 or L2 offet

      data first_call  / .true. /

      save first_call

****  If this is first call; read the L1/L2 offsets from the 
*     PCV model file
      if( first_call ) then
         n = 0
*        Read until end of header
         found = .false.
         do while ( .not. found )
            read(ipcv_unit,'(a)',iostat=ierr) line
            if( ierr.eq.0 .and. 
     .          index(line,'END OF HEADER').gt.0 ) then
                found = .true.
            end if
            if( ierr.ne.0 ) found = .true.
         end do
*        Now start getting antenna information.
*        Loop till we find start of antenna
         do while ( ierr.eq.0 )
            read(ipcv_unit,'(a)',iostat=ierr) line
            if( index(line,'COMMENT').lt.60 .and.
     .          ierr.eq.0 ) then
*               Not a comment line
*               See if antenna type
                if( index(line,'TYPE').eq.61 ) then
*                   Make sure not a satellite
                    svstr = line(41:44)
                    if( trimlen(svstr).eq.0 ) then
*                      This is not a satellite
                       n = n + 1
                       call ckmax(n, max_rcvant,'MAX_ANTSAVES')
                       name_offsave(n) = line(1:20)
                       read_offs = .true.
                    else
                       read_offs = .false.  ! don't read values for sats
                    end if
                elseif( index(line,'START OF FREQUENCY').eq.61 ) then
                    L1_true = .false.
                    L2_true = .false.
                    if( line(4:6).eq.'G01' ) L1_true = .true.
                    if( line(4:6).eq.'G02' ) L2_true = .true.
                elseif( index(line,'NORTH / EAST / UP').eq.61 ) then
*                   If we have L1_true read values
* MOD TAH 191029: Only read offsets if not a satellite (n=0 causes problems)
                    if( L1_true .and. read_offs ) then
                        read(line,'(3F10.2)',iostat=jerr) 
     .                      (L1_offSave(i,n),i=1,3)
                        L1_offSave(:,n) = L1_offSave(:,n)/1000
* MOD TAH 130130: Set L1_true to false so that NORTH / EAST / UP line from
*                       possible FREQ RMS lines are not read
                        L1_true = .false.
                    end if

*                   If we have L2_true read values
* MOD TAH 191029: Only read offsets if not a satellite (n=0 causes problems)
                    if( L2_true .and. read_offs) then
                        read(line,'(3F10.2)',iostat=jerr) 
     .                      (L2_offSave(i,n),i=1,3)
                        L2_offSave(:,n) = L2_offSave(:,n)/1000
* MOD TAH 130130: Set L2_true to false so that NORTH / EAST / UP line from
*                       possible FREQ RMS lines are not read
                        L2_true = .false.
                    end if
                    call report_error('IOSTAT',jerr,'decod',
     .                   line,1,'antpcv_save')
                endif
             end if   ! Not comment
          end do      ! Reading file
          first_call = .false.
          num_offSave = n
          close(ipcv_unit)
          write(*,120) num_offSave
 120      format('Read antmod.dat: ',I4,' entries saved')
  
      endif            ! First call

****  Now see if we can find our antenna
      found = .false.
      radome_sub = .true.
      do i = 1, num_offSave
         if( antnam(1:16).eq.name_offsave(i)(1:16) .and.
     .       radome(1:4).eq.name_offsave(i)(17:20) ) then
            radome_sub = .false.
            found = .true.
            offset_L1 = L1_offSave(:,i)
            offset_L2 = L2_offSave(:,i)
            exit
         end if
      end do

*     IF we did not find match, then test with NONE of radome
      if( .not.found  ) then
         do i = 1, num_offSave
            if( antnam(1:16).eq.name_offsave(i)(1:16) .and.
     .          'NONE'.eq.name_offsave(i)(17:20) ) then
               found = .true.
               offset_L1 = L1_offSave(:,i)
               offset_L2 = L2_offSave(:,i)
               exit
            end if
         end do
      endif

****  If we did not find, warn user
      if( .not.found ) then
         write(*,220) antnam(1:16), radome(1:4)
 220     format('No match found in antmod.dat for ',a,1x,a)
         offset_L1 = 0
         offset_L2 = 0
*        Add this entry to table so that only reported once
         num_offSave = num_offSave+1
         call ckmax(num_offSave, max_rcvant,'MAX_ANTSAVES_missing')
         name_offsave(num_offSave) = antnam(1:16) // radome(1:4)
         L1_offSave(:,num_offSave) = 0
         L2_offSave(:,num_offSave) = 0
      end if

****  Thats all
      return
      end 
