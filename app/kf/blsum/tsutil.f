CTITLE READ_TSPOS

      subroutine read_tspos(unit)

*     Routine to read the tssum file attached to unit unit

      implicit none

      include 'tsfit.h'


      integer*4 unit   ! Unit number for position file

      integer*4 ierr, jerr, trimlen, date(5), i, j, isec
      integer*4 ts_srt(max_ts)  ! indexs to sorted time series 
      integer*4 kn, cnt, in 
      real*8 sec, rot(3,3), unc_geod(3)
      real*8 max_mjd, min_mjd   ! Max and min m,jd used in sorting
      logical line_read  ! Set true when line has been read checking for
                         ! Field dscriptions.
      logical time_order_out  ! Set true to denote file out of time order
      logical set ! Set true when ranked value put in ts_srt



      character*512 line

      rewind(unit)

****  Read the header records
      read(unit,'(a)',iostat=ierr) line
*     Exract the reference frame:
      read(line,'(52x,a)') reference_frame
c 100  format('PBO Station Position Time Series.',
c     .          ' Reference Frame : ',a)
* MOD TAH 051129: See if version line next
      call report_error('IOSTAT',ierr,'read',line,0,'read_tspos')
      read(unit,'(a)',iostat=ierr) line
      if( index(line,'Version').gt.0 ) then
          read(line,105,iostat=ierr) ts_ver_read
 105      format(16x,a)
c105      format('Format Version: ',a)
          call report_error('IOSTAT',ierr,'read',line,0,'read_tspos')
          read(unit,'(a)',iostat=ierr) line
      endif

      read(line,110,iostat=ierr) ts_code
 110  format(16x,a)
c110  format('4-character ID: ',a4)
      call report_error('IOSTAT',ierr,'read',line,0,'read_tspos')
      read(unit,120,iostat=ierr) ts_full
 120  format(16x,a16)
c120  format('Station name  : ',a16)
      call report_error('IOSTAT',ierr,'read','Station',0,'read_tspos')
      read(unit,'(a)',iostat=ierr) line
      read(line,130,iostat=ierr) date, isec
 130  format(16x,i4,i2,i2,1x,i2,i2,i2)
      sec = isec+1.d-3
c130  format('First Epoch   : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
      call report_error('IOSTAT',ierr,'read',line,0,'read_tspos')
      call ymdhms_to_mjd(date,sec,ts_first)
cFirst Epoch   : 20060521 120000
      read(unit,140,iostat=ierr) date, isec
      sec = isec+1.d-3
 140  format(16x, i4,i2,i2,1x,i2,i2,i2)
      call report_error('IOSTAT',ierr,'read','Last Epoch',
     .                   0,'read_tspos')
c140  format('Last Epoch    : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
      call ymdhms_to_mjd(date,sec,ts_last)

C**** Do not read release date, since we will generate a new one here
C     read(unit,150,iostat=ierr) date_rel, nint(sec_rel)
      read(unit,'(a)',iostat=ierr) line
C150  format(16x, i4,i2,i2,1x,i2,i2,i2)
c150  format('Release Date  : ', i4,i2.2,i2.2,1x,i2.2,i2.2,i2.2)
      call report_error('IOSTAT',ierr,'read',line,0,'tssum')

      read(unit,200,iostat=ierr) ref_xyz
 200  format(25x,3F15.5)
c200  format('XYZ Reference position : ',3F15.5)
      call report_error('IOSTAT',ierr,'read',line,0,'tssum')

      call XYZ_to_GEOD(rot,ref_xyz, ref_llu) 
      call loc_to_geod(ref_llu, ref_neu)

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

****  OK, Now read in the time series
      i = 0
      time_order_out = .false.
      do while ( ierr.eq.0 )
          if( .not.line_read ) then
              read(unit,'(a)',iostat=ierr) line
          endif
          line_read = .false.

          if( ierr.eq.0 ) then
             i = i + 1
             read(line,300,iostat=jerr) date, isec, ts_mjd(i),
     .        (ts_xyz(j,i),j=1,3), (ts_xyz_std(j,i),j=1,6),
     .        (ts_llu(j,i),j=1,3),
     .        (ts_neu(j,i),j=1,3), (ts_neu_std(j,i),j=1,6),
     .         ts_type(i)
  
             sec = isec+1.d-3
 300         format(1x,i4,i2,i2,1x,i2,i2,i2,1x,F10.4,
     .           3F15.5,3F9.5,3F7.3,3x,2F16.10,1x,F10.5,3x,
     .           3(F9.5,1x),1x,3F9.5,3F7.3,1x,a5)
             call report_error('IOSTAT',jerr,'read',line,0,'read_tspos')

* MOD TAH 140808: If there is a error, skip the line inless ts_keep set.
             if( jerr.eq.0 .or. ts_keep ) then 
*                Process line
                call ymdhms_to_mjd(date,sec,ts_mjd(i))

*
C MOD TAH 160821: Use the ts_neu values read from pos file rather than
C MOD TAH 191217: Comment: We seem to do this because the offsets in the
C              East values due to offsets have been "fixed" in the pos files
C              but this then can cause issues with the NEU values being 
C              converted back to LLU and XYZ values.  Current code and
C              convert_ltog in gen_util/loc_to_geod.f seems to fix this.
               do j = 1,3
                    ts_neu(j,i) = ts_neu(j,i) + ref_neu(j)
               end do
****           Now get NEU values: By coonversion of XYZ values
C MOD TAH 160821: Use the ts_neu values read from pos file rather than
c               re-computes.
C               call XYZ_to_GEOD(rot, ts_xyz(1,i), unc_geod )
C               call loc_to_geod(unc_geod, ts_neu(1,i))

****            Make sure that time series is in time order.
                if( i.gt.1 .and. ts_mjd(i).lt.ts_mjd(i-1) ) then
                   write(*,320) ts_file(1:trimlen(ts_file)),i
 320               format('**WARNING** Time series in ',a,
     .                    ' out of time order at record ',i5)
                   time_order_out = .true.
                endif 
* MOD TAH 090507: Check and remove zero sigma poinnts
                if( ts_neu_std(1,i).eq.0.0 .or.  ts_neu_std(2,i).eq.0.0 
     .              .or. ts_neu_std(3,i).eq.0.0 ) i = i - 1
             else
*               Line is bad so skip
                i = i - 1
             end if
          end if
      end do

* MOD TAH 080104: If time seris out of order, resort now.
      num_ts = i
      if( time_order_out ) then
*        Buildup up list of indices in time order.  In this
*        approach dupliucate times are replaced with later 
*        values.
         do i = 1,num_ts
             ts_srt(i) = 0
         end do
         max_mjd = -1e20
         cnt = 0
         do i = 1, num_ts
             min_mjd = 1e20
             set = .false.
             do j = 1, num_ts
                if( ts_mjd(j) .le. min_mjd .and.
     .              ts_mjd(j) .gt. max_mjd ) then
                    ts_srt(i) = j
                    min_mjd  = ts_mjd(j)
                    set = .true.
                endif
             enddo
             if( set ) then
                 cnt = cnt + 1
                 max_mjd = ts_mjd(ts_srt(i))
             end if
         enddo 
*        ts_srt now is sorted list.  Sort the list of values
         do i = 1, num_ts
            in = ts_srt(i)
*           if in > 0 then swap the ith entry with the in'th
*           value
            call save_ts(i,0)  ! cp i'th entries to save values (0)
            call save_ts(in,i) ! cp in'th entries to ith values
            call save_ts(0,in) ! cp saved values to in'th entry
*           Now update ts_srt to show values have moved
            kn = 0
            do j = i+1,num_ts
               if( ts_srt(j).eq.i ) then
                   kn =j
               endif
            enddo
            if( kn.gt. 0 ) then
                ts_srt(kn) = in
            endif
         end do
*        Save number of values in time series
         num_ts = cnt
      endif

      write(*,400) ts_file(1:trimlen(ts_file)), num_ts
 400  format('TS_file ',a,' sorted entries ',i5)
 
      return
      end

CTTLE TSREAD_EQFILE

      subroutine tsread_eqfile


*     Routine to read the eq file

      implicit none

      include 'tsfit.h'
      include '../includes/const_param.h'

* LOCAL VARIABLES

* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   trimlen - Length of string
*   iel     - Command number from list of earthquake commands
*   indx    - Position in string that we have decoded to.
 
*   jel     - Earthquake number
*   jndx    - Index when looking up code
*   date(5) - YMDHM of earthquake
*   jerr    - Error to test if hfile name passed in rename
*   indx_curr - Current value of indx incase no hfile name passed
 
      integer*4 ierr, trimlen, indx, jel, jndx, date(5), j,
     .          jerr, indx_curr, eol
 
*   geod(3) - Geodetic postion of earthquake (co-lat, long,
*           - ellipsiodal height.
*   sectag  - Seconds tag
*   lat, long       - Decimal degrees for epicenter
*   rad_km, dep_km  - Radius of influence and depth in
*           - km (converted to meters for internal use)
 
      real*8 geod(3), sectag, lat, long, rad_km, dep_km
 
*   buffer  - Line read from input file
 
      character*(256) buffer
 
*   test_code   - TEst code for Earth quake to see if
*           - duplicate or new event.
*	cdum    - Dummy string for read_line
*   rn_test     - Test string to see if hfile name passed in rename
 
      character*8 test_code, cdum, cmd
      character*16 rn_test

****  open the file
      open(101,file=eqfile,iostat=ierr,status='old')
      call report_error('IOSTAT',ierr,'open',eqfile,0,
     .                  'tsfit/tsread_eqfile')
      if( ierr.ne. 0 ) RETURN

*     Start reading file

      do while ( ierr.eq.0 )
          read(101,'(a)',iostat=ierr) buffer
*         Clear any comments
          eol = index( buffer,'!' )
          if( eol.eq.0 ) eol = index(buffer,'#')
          if( eol.gt.0 ) then
              buffer(eol:) = ' '
          end if

          if( ierr.eq.0 .and. buffer(1:1).eq.' ' .and.
     .        trimlen(buffer).gt.0 )  then
 
*             OK, Valid line for command.  Process:
              indx = 0
              call GetWord(buffer, cmd, indx )
              call casefold(cmd)

*     EQ_DEF: Define the parameters of a specifiic Earth quake.
*         Command must be given before the other commands
*         referring to this earthquake
              if( cmd(1:4).eq.'EQ_D' ) then   ! Eq_def
*                Get the Earthquake code
                 call getword(buffer, test_code, indx)
                 call casefold(test_code)
 
*                See if we have used this code before
                 jndx = 1
                 call get_cmd(test_code, eq_codes, num_eq, jel, jndx)
                 if( jel.gt.0 ) then
                     write(*,120) buffer(1:trimlen(buffer))
 120                 format('**WARNING** Earthquake defined by: ',a,/,
     .                      14x,'has been previously defined.  ',
     .                       'Using new definition')
                 else if ( jel.lt.-1 ) then
                     write(*,140) buffer(1:trimlen(buffer))
 140                 format('**WARNING** Earthquake defined by: ',a,/,
     .                  14x,'is not UNIQUE code. Event being ignored')
                 else
 
*                    Add this event
                     num_eq = num_eq + 1
                     jel = num_eq
                     if( num_eq.gt.max_eq ) then
                         write(*,160) max_eq, buffer(1:trimlen(buffer))
 160                     format('**WARNING** Maxiumum number of EQS',
     .                        '(',i4,') Exceeded.  Ignoring ',/,a)
                         jel = 0
                     end if
                 end if
 
****             Now process the command if jel is greater than 0 indicating
*                that we have a valid earthquake
                 if( jel.gt.0 ) then

                     eq_codes(jel) = test_code
                     read(buffer(indx:), *, iostat=jerr) lat, long, 
     .                       rad_km,dep_km, date
                     geod(1) = pi/2 - lat*pi/180
                     geod(2) = long*pi/180
                     geod(3) = -dep_km*1000.d0
                     call GEOD_to_XYZ( geod, eq_pos(1,jel) )
 
                     eq_rad(jel) = rad_km*1000.d0
                     eq_depth(jel) = dep_km*1000.d0
 
                     sectag = 0.d0
                     call ymdhms_to_mjd( date, sectag, eq_epoch(jel))

*                    Initialize the values associated with the earthquake
                     do j = 1, 6
                        eq_apr_coseismic(j,jel) = 0.d0
                        eq_log_sig(j,jel) = 0.d0
                        eq_exp_sig(j,jel) = 0.d0
                     end do
                     log_eq(jel) = .false.
                     exp_eq(jel) = .false.
                     dtr_eq(jel) = .false.

                 endif
              
* EQ_Rename Not needed: 
 
****  EQ_COSEI: Set the coseimic displacement apriori sigmas
              elseif( cmd(1:4).eq.'EQ_C' ) then
*                Get the Earthquake code
                 call getword(buffer, test_code, indx)
                 call casefold(test_code)
                 jndx = 1
                 call get_cmd(test_code, eq_codes, num_eq, jel, jndx)
                 if( jel.gt.0 ) then
 
*                    Now read the rest of the line.  sigmas in NEU (m) and
*                    spatial sigmas (NEU at 1 earthquake depth).
                     read(buffer(indx:),*,iostat=jerr) 
     .                       (eq_apr_coseismic(j,jel), j = 1,6)
                     call report_error('IOSTAT',jerr,'read',buffer,0,
     .                       'PROC_EQ_CMD/COSEISMIC')
                 else
                     write(*,210) buffer(1:trimlen(buffer))
 210                 format('**WARNING** Earthquake code not found.',
     .                      ' Command:',a)
                 end if

* EQ_Post not needed
* EQ_Pre  not needed
            elseif( cmd(1:4).eq.'EQ_L' ) then  
***** EQ_LOG : Read the log estimation information
*                Get the Earthquake code
                 call getword(buffer, test_code, indx)
                 call casefold(test_code)
                 jndx = 1
                 call get_cmd(test_code, eq_codes, num_eq, jel, jndx)
                 if( jel.gt.0 ) then
 
*                    Now read the rest of the line. Markov process values
*                    (fixed and spatially dependent)
                     read(buffer(indx:),*,iostat=jerr) eq_log_tau(jel),
     .                            (eq_log_sig(j,jel), j = 1,6)
                     call report_error('IOSTAT',jerr,'read',buffer,0,
     .                       'PROC_EQ_CMD/EQ_LOG')
                     log_eq(jel) = .true.
                 else
                     write(*,210) buffer(1:trimlen(buffer))
                 end if

* ne XP: This command is not in globk but can be used here
             elseif( cmd(1:4).eq.'EQ_E' ) then  
***** EQ_EXP : Read the exponent estimation information
*                Get the Earthquake code
                 call getword(buffer, test_code, indx)
                 call casefold(test_code)
                 jndx = 1
                 call get_cmd(test_code, eq_codes, num_eq, jel, jndx)
                 if( jel.gt.0 ) then
 
*                    Now read the rest of the line. Markov process values
*                    (fixed and spatially dependent)
                     read(buffer(indx:),*,iostat=jerr) eq_exp_tau(jel),
     .                            (eq_exp_sig(j,jel), j = 1,6)
                     call report_error('IOSTAT',jerr,'read',buffer,0,
     .                       'PROC_EQ_CMD/EQ_EXP')
                     exp_eq(jel) = .true.
                 else
                     write(*,210) buffer(1:trimlen(buffer))
                 end if

***** RENAME : General feature for renaming sites
             elseif( cmd(1:2).eq.'RE' ) then
   
*                Make sure we do not have too many renames already.
                 if( num_rn+1.gt.max_rn ) then
                     write(*,610) max_rn, (buffer(1:trimlen(buffer)))
 610                 format('**WARNING** Too many renames specified.',
     .                   ' Maximum allowed ',i5,/,
     .                   ' Ignoring ',a)
                     RETURN
                 end if
 
                 num_rn = num_rn + 1
                 call getword(buffer, rn_codes(1,num_rn),indx)
                 call getword(buffer, rn_codes(2,num_rn),indx)
                 call casefold(rn_codes(1,num_rn))
                 call casefold(rn_codes(2,num_rn))

* MOD TAH 971112: Check to see if a name restriction has been included
*                in the line
                 indx_curr = indx
                 call getword(buffer, rn_test, indx )
*                See if number or string
                 call check_num( rn_test, jerr )
                 if( jerr.ne.0 ) then
                     rn_hfiles(num_rn) = rn_test
                     lrn_hf(num_rn) = trimlen(rn_test)
                 else
                     rn_hfiles(num_rn) = ' '
                     indx = indx_curr
                     lrn_hf(num_rn) = 0
                 end if

*                Now start reading the oprional parts of the command.
*                Start time for the rename (ymdhms)
                 do j = 1,5
                     call read_line(buffer,indx,'I4',jerr,date(j),cdum)
                     if( jerr.ne.0 ) date(j) = 1
                 end do
*                If there is an error, default the start to 1900
                 if( jerr.ne.0 ) date(1) = 1900
*                Process the starting data
                 sectag = 0.0d0
                 call ymdhms_to_mjd(date,sectag,rn_times(1,num_rn))

****             Get the end time:
                 do j = 1,5
                     call read_line(buffer,indx,'I4',jerr,date(j),cdum)
                     if( jerr.ne.0 ) date(j) = 1
                 end do
*                If there is an error, default the end to 2100
                 if( jerr.ne.0 ) date(1) = 2100
 
*                Process the starting data
                 sectag = 0.0d0
                 call ymdhms_to_mjd(date,sectag,rn_times(2,num_rn))
 
****             See if there is a displacement associated with this
*                rename:
                 do j = 1,3
                     call read_line(buffer,indx,'R8',jerr,
     .                   rn_dpos(j,num_rn),cdum)
                     if( jerr.ne.0 ) rn_dpos(j,num_rn) = 0.d0
cd		     print*,'num_rn, j, rn_dpos:',num_rn,j,rn_dpos(j,num_rn)
                 end do
 
****             Now see what frame this is (XYZ/NEU)
                 call read_line(buffer,indx,'CH',jerr,sectag,
     .                   rn_types(num_rn))
                 call casefold(rn_types(num_rn))
                 if( jerr.ne.0 ) rn_types(num_rn) = 'XYZ'
		 
***** OFFSET <site year mo dd hr mn s.s n_off e_off u_offset>
             elseif (cmd(1:3).eq.'OFF' ) then  
*                Make sure we do not have too many renames already.
                 if( num_off+1 .gt. max_off ) then
                     write(*,670) max_off, (buffer(1:trimlen(buffer)))
 670                 format('**WARNING** Too many offsets specified.',
     .                   ' Maximum allowed ',i5,/,
     .                   ' Ignoring ',a)
                     return
                 end if
 
                 num_off = num_off + 1
                 call getword(buffer, off_codes(num_off),indx)
                 call casefold(off_codes(num_off))

*                Now start reading the ather parts of the command.
*                Time for the offset (ymdhms)
                 do j = 1,5
                     call read_line(buffer,indx,'I4',jerr,date(j),cdum)
                     if( jerr.ne.0 ) date(j) = 1
                 end do
*                If there is an error, default the start to 1900
                 if( jerr.ne.0 ) date(1) = 1900
*                Process the starting data
                 sectag = 0.0d0		 
                 call read_line(buffer,indx,'R8',jerr,sectag,cdum)
		 
                 call ymdhms_to_mjd(date,sectag,off_times(num_off))

****             Get the valuse associated with this offset
                 do j = 1,3
                     call read_line(buffer,indx,'R8',jerr,
     .                   off_dpos(j,num_off),cdum)
                     if( jerr.ne.0 ) off_dpos(j,num_off) = 0.d0
                 end do
 
****             Now see what frame this is (XYZ/NEU)
                 call read_line(buffer,indx,'CH',jerr,sectag,
     .                   off_types(num_off))
                 call casefold(off_types(num_off))
                 if( jerr.ne.0 ) off_types(num_off) = 'NEU'	     	  

***** BREAK: Command for adding break at a site
             elseif( cmd(1:2).eq.'BR' ) then
   
*                Make sure we do not have too many renames already.
                 if( num_rn+1.gt.max_rn ) then
                     write(*,710) max_rn, (buffer(1:trimlen(buffer)))
 710                 format('**WARNING** Too many breaks specified.',
     .                   ' Maximum allowed ',i5,/,
     .                   ' Ignoring ',a)
                     RETURN
                 end if
 
                 num_rn = num_rn + 1
                 call getword(buffer, rn_codes(1,num_rn),indx)
                 call casefold(rn_codes(1,num_rn))

*                Now set the rn_codes(2,num_rn) to denote that
*                this is a break command
                 rn_codes(2,num_rn) = rn_codes(1,num_rn)(1:4) 
     .                                     //  '_?PS' 


* MOD TAH 971112: Check to see if a name restriction has been included
*                in the line
                 indx_curr = indx
                 call getword(buffer, rn_test, indx )
*                See if number or string
                 call check_num( rn_test, jerr )
                 if( jerr.ne.0 ) then
                     rn_hfiles(num_rn) = rn_test
                     lrn_hf(num_rn) = trimlen(rn_test)
                 else
                     rn_hfiles(num_rn) = ' '
                     indx = indx_curr
                     lrn_hf(num_rn) = 0
                 end if

*                Now get the date
                 call multiread(buffer,indx,'I4',ierr,date,cdum,5)
*                Process the starting data
                 sectag = 0.0d0
                 call ymdhms_to_mjd(date, sectag,rn_times(1,num_rn))
*                set an open ended end-date
                 date(1) = 2100
                 do j = 2,5
                    date(j) = 1
                 end do
                 call ymdhms_to_mjd(date, sectag,rn_times(2,num_rn))

****             Finally see if open re-name string has been passed
                 call getword(buffer, rn_test, indx )
                 call casefold(rn_test)
                 jel = trimlen(rn_test)
                 if( jel.gt.0 ) then
                      rn_codes(2,num_rn)(5:5+jel) = '_' 
     .                                              // rn_test(1:jel)
                 endif

****             Set the offsets to zero for this rename
                 do j = 1,3
                     rn_dpos(j,num_rn) = 0.d0
                 end do
                 rn_types(num_rn) = 'NEU'

             endif
          endif
      end do

      write(*,320) num_eq, num_rn, eqfile(1:trimlen(eqfile))
 320  format('EQ: Found ',I4,' earthquakes, ',I8,' renames in ',a)

      close(101)
  
***** Thats all
      return
      END

