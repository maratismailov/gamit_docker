 
CTITLE MAKE_GLB
 
      subroutine make_glb( unit, num_in_file, vma_data, len_vma )
 
      implicit none

*     This routine will actually read, and make the global file.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
*   num_in_file - NUmber of covariance matrix in file
*   len_vma     - LEngth of vma space allowwed
*   vma_data(len_vma)   - Area for storing covariance matrix and
*               - solution vector
 
      integer*4 unit, num_in_file, len_vma, vma_data(len_vma)
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   indx        - Pointer in string
*   date(5)     - Date of run
*   i           - Loop counter
*   icov_parm, isol_parm    - Pointers to starts of covariance
*               - matrix and solution vector.
*   icov_copy, ieigvec,  ieigval, ibwork, izwork -- Pointers for doing 
*                 eigenvalues
*                 computation
*   iend        - Last word needed in memory
*   np          - Number of parameters in JPL file
*   etide_mod, erot_mod -- Tide and rotation models used.
*   otide_mod, oatml_mod, atmtd_mod, hydro_mod -- Models used for
*     loads

      integer*4 ierr, trimlen, indx, date(5), i, np, 
     .          etide_mod, erot_mod, jerr
      integer*4 otide_mod, oatml_mod, atmtd_mod, hydro_mod

      integer*8 icov_parm, isol_parm,iend, icov_copy, ieigvec,  ieigval, 
     .          ibwork, izwork
 
*   sectag      - Seconds tag of runtime
*   real_ver  - Real*8 version of solve version
      real*8 sectag, real_ver
 
*   finished        - Indicates that this q file is completed
*   svant_est       - Set true if antenna offsets estimated

      logical finished, svant_est, kbit
 
*   line        - Line read from file
 
      character*256 line
* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      data I8 / 1 /

      icov_parm = 0
      isol_parm = 0
 
***   Read off the header at the start
      num_in_file = 0
      finished = .false.

*     Initialize the satellite information
      call init_svinf
 
*                                             ! Blank line
      read(unit,'(a)', iostat=ierr) line
*
* MOD TAH 941005: See if this is a stacov file.
***** See if this is a JPL stacov file:
      if( index(line,'parameters').gt.0 .or.
     .    index(line,'PARAMETERS').gt.0        ) then
          read(line,*) np
          icov_parm = istart_vma
          isol_parm = icov_parm + 2*np*np

          call make_jpl( unit, line, np, vma_data(icov_parm), 
     .                   vma_data(isol_parm) )
          num_in_file = 1
          RETURN
      end if

* MOD TAH 950619: See if this a SINEX file

      if( index(line,'%=SNX').gt.0 ) then
          call preread_snx( unit, line, np, ierr, sepoch ) 
          icov_parm = istart_vma
          isol_parm = icov_parm + 2*np*np
          if( snx_constraint.gt.0 .or. snx_constraint.le.0) then 
             icov_copy = isol_parm + 2*np
             ieigvec   = icov_copy + 2*(I8*np)*np
             ieigval   = ieigvec   + 2*(I8*np)*np
             ibwork    = ieigval   + 2*np
             izwork    = ibwork    + 2*np
          else
             icov_copy = isol_parm + 2*np
             ieigvec   = icov_copy 
             ieigval   = ieigvec  
             ibwork    = ieigval 
             izwork    = ibwork  
          end if  
          
          if( izwork+2*np-istart_vma.gt.len_vma ) then
              write(*,95) (izwork+2*np-istart_vma)*4/1024.**2,
     .                    (icov_copy-istart_vma)*4/1024.**2
 95           format('**DISASTER** ',f8.2,' Mbytes needed for run.',
     .               ' Use the -M option to increase program size',/,
     .               ' Trying running without -C= option (if used).',
     .             /,' Only ',F8.2,' needed in this case')
                 
              call report_stat('FATAL', 'htoglb','make_snx',
     .              ' ','Not enough memory', 
     .              nint((izwork+2*np-istart_vma)*4/1024.**2) )
          end if
          
          call make_snx( unit, line, np, vma_data(icov_parm), 
     .                   vma_data(isol_parm),
     .                   vma_data(icov_copy),
     .                   vma_data(ieigvec), 
     .                   vma_data(ieigval),
     .                   vma_data(ibwork),
     .                   vma_data(izwork) )
     
          num_in_file = 1
          RETURN
      end if

* MOD TAH 941212: See if a GSFC SLR solution file
      if( index(line,'Name of solution').gt.0 ) then
*         This is a GSFC solution
          read(line,'(18x,a)') sdata_base
          call get_slr_parn(unit, line, np )
          rewind(unit)
          icov_parm = istart_vma
          isol_parm = icov_parm + 2*np*np
          call make_slr( unit, line, np, vma_data(icov_parm),
     .                   vma_data(isol_parm) ) 
          num_in_file = 1
          RETURN
      end if

* MOD TAH 950112: See if this is GFSC VLBI FILE.  (We need to
*     add the header line to this files to denote type).
      if( index(line,'GSFC VLBI').gt.0 ) then

*         GSFC VLBI cvr* file.
*         Read the next line to get the number of parameters
          read(unit,'(a)') line
          read(line,'(7x,i6)', iostat=ierr) np
          call report_error('IOSTAT',ierr,'read',line,0,
     .                      'GSFC VLBI parameter numbers')
          if( ierr.ne.0 ) then
              num_in_file = 0
          else
              icov_parm = istart_vma
              isol_parm = icov_parm + 2*np*np

              call make_gsfc_vlbi(unit, line, np, vma_data(icov_parm), 
     .                       vma_data(isol_parm) )
              num_in_file = 1
          end if
          RETURN
      end if

* MOD TAH 950112: See if a UT SLR solution file
      if( index(line,'CSR position').gt.0 ) then
*         This is a UT SLR solution
          read(line,'(22x,a)') sdata_base
          call get_ut_slr_parn(unit, line, np )
          rewind(unit)
          icov_parm = istart_vma
          isol_parm = icov_parm + 2*np*np
          call make_ut_slr( unit, line, np, vma_data(icov_parm),
     .                   vma_data(isol_parm) ) 
          num_in_file = 1
          RETURN
      end if

* MOD TAH 950301: See if a MIT CTOGOBS H-file
      if( index(line,'CTOGOBS H-FILE').gt.0 ) then
          read(unit,'(a)') line
          read(line(24:),*) np
          icov_parm = istart_vma
          isol_parm = icov_parm + 2*np*np

          call make_ctog( unit, line, np, vma_data(icov_parm), 
     .                    vma_data(isol_parm) )
          num_in_file = 1
          RETURN
      end if


*                                             ! GLOBALQ file line
      read(unit,'(a)', iostat=ierr) line
      call report_error('IOSTAT',ierr,'read','Line 2',0,'make_glb')
      if( ierr.ne.0 ) finished = .true.
 
*     Check to see if correct line, Set the hfile_type to be blank
*     and see if we have the correct file types.
      hfile_type = ' '
      gamit_hf_ver = 0
      indx = index(line,'GLOBALQ')
      if( indx.ne.0 ) hfile_type = 'GPS'
      indx = index(line,'GLOBALT')
      if( indx.ne.0 ) hfile_type = 'TER'
      indx = index(line,'GAMIT H-file')
      if( indx.ne.0 ) then   ! New gamit h-file version
          hfile_type = 'GPS'
          read(line(42:),*,iostat=ierr) real_ver
          call report_error('IOSTAT',ierr,'decod',line,0,
     .         'GAMIT H-file version')
          if( ierr.ne.0 ) real_ver = 2.0
          gamit_hf_ver = nint(real_ver*100)
      end if

*     See if we a legitimate type of file
      if( trimlen(hfile_type).eq.0 ) then
          write(*,100) hfile(1:trimlen(hfile)),
     .                line(1:trimlen(line))
  100     format(1x,a,' does not appear to be a hfile. ',
     .        ' Second line is:',/,a)
          finished = .true.
      end if

***** Skip next line

      read(unit,'(a)', iostat=ierr ) line
         
***** Get some more header information
      if( .not.finished ) then
          read(unit,'(a)', iostat=ierr) line
 
*         Generate the Fake KalObs file name
          sKalobs_file = hfile(1:trimlen(hfile)) 
          qexpt_title  = line(2:)
          orb_est = .false.  ! Set true when orbits estimated (non 
                         ! BASELINE mode solution).

* MOD TAH 050214: Read Solve version as well
          read(line,190,iostat=ierr) real_ver
 190      format(16x,f5.2)
          solve_ver = nint(real_ver*100)

*         Get the runtime from next line
          read(unit,'(a)', iostat=ierr) line
          read(line,200,iostat=ierr)  date, sectag
  200     format(15x,i4,1x,i2,1x,i2,2x,2(i2,1x),f3.0)
          do i = 1,5
              qrun_time(i) = date(i)
          end do
          qrun_time(6) = nint(sectag)
          qrun_time(7) = 0
 
 
*         Get the file owner from the next line.  Get the second
*         word from the line
          read(unit,'(a)', iostat=ierr) line
          indx = 1
          call getword( line, cowner, indx)
          call getword( line, cowner, indx)
          
*         Get Mfile name (claim this is database)
          read(unit,'(a)', iostat=ierr) line
          sdata_base = line(15:24)
 
*         Get the datum.  Check to make sure geocentric
          read(unit,'(a)', iostat=ierr) line

*         Set the datum type to blank and see if we find
*         one which matches
          datum_type = ' '
          indx = index(line,'geocentric')
          if( indx.gt.0 ) datum_type = 'LLR'
          indx = index(line,'cartesian')
          if( indx.gt.0 ) datum_type = 'XYZ'

*                                     ! Not geocentric, get out
          if( trimlen(datum_type).eq.0 ) then
              write(*,250) line(1:max(1,trimlen(line)))
  250         format(' Datum NOT found. Datum line is:',/,
     .                a,/,' Stopping processing of this file')
              finished = .true.
          end if
      end if
 
****  We have now finished non-repeating part of header.  If still
*     OK get the rest. Now start looping over individual soltions
*     Skip next line
* MOD TAH 980514: Check the model types and set the sgamit_mod
*     variable.
      read(unit,'(a)', iostat=ierr) line
* MOD TAH 050214: Read on old format if solve_ver < 1017
      if( solve_ver .lt. 1017 ) then
* Old line format
* Models:  E-tides 15   SP E-rot  7
          read(line,'(17x,i3,11x,i3)', iostat=jerr) etide_mod, erot_mod
          otide_mod = 0
          oatml_mod = 0
          atmtd_mod = 0
          hydro_mod = 0
      elseif( solve_ver .lt. 1042 ) then
* New line format (Later we will extract more from line)
* Models:  SP E-rot IERS96    7  E-tides IERS96   15  O-load OSO       3  Atm-load           0 Atm-tide           0  Hydrol-load           0
* MOD TAH 200219: Changed 1x,i2 to I3 to accomdate larger values.
* Models:  SP E-rot IERS10   11  SE-tide IERS03   19  O-load FES2004E  3  Atm-load           0 Atm-tide           0  Hydrol-load           0
*                                                 ^ Bit mapping changed to match sestbl.

          read(line,260, iostat=jerr) cspeopmod, erot_mod, 
     .       cetidemod, etide_mod, cotidemod, otide_mod, 
     .       coatmlmod, oatml_mod, catmtdmod, atmtd_mod, 
     .       chydromod, hydro_mod
 260      format(7x,12x,a8,  I3,10x,a8,  I3, 9x,a8,  I3,
     .              11x,a8,  I3,10x,a8,  I3,14x,a8,  I3)  
      else
!         Same format as above but we with new model, bit 4 (8) in the E-tide value
!         will mean the new mean pole position (IERS2010) model is used.
          read(line,260, iostat=jerr) cspeopmod, erot_mod, 
     .       cetidemod, etide_mod, cotidemod, otide_mod, 
     .       coatmlmod, oatml_mod, catmtdmod, atmtd_mod, 
     .       chydromod, hydro_mod
* MOD TAH 200220: Update for the IERS2020 mean pole model.  See if new model
*         is on; If so set appropriate bits.  If older model, save as before
          if( kbit(etide_mod,5) )  then
!             IERS 2020 Mean pole used: Bit 10 in GLOBK word.  In new version
*             no other bits should be on.  Check and reort.
              call sbit(etide_mod,10,1)  ! Turn on new mean IERS2020
              call sbit(etide_mod, 5,0)   ! Turn off old mean
              call sbit(etide_mod, 3,1)  ! Set to show poletide on.
          endif  
* MOD TAH 200220: Code below left unchanged from earilier versions
          if( kbit(etide_mod,4) )  then
!             IERS 2010 Mean pole used
              call sbit(etide_mod, 7,1)   ! Turn on new mean
              call sbit(etide_mod, 5,0)   ! Turn off old mean
              call sbit(etide_mod, 3,1)   ! Set to show poletide on.
          endif  
      endif 
      
      if( jerr.eq.0 ) then
          sgamit_mod = etide_mod*65536 + erot_mod
* MOD TAH 120714: Updated bits to show status that loads are available
*         when they have been applied
          if( kbit(oatml_mod,1) ) call sbit(oatml_mod,2,1)
          if( kbit(hydro_mod,1) ) call sbit(hydro_mod,2,1)

          sload_mod = otide_mod + 256*oatml_mod + 
     .                65536*(atmtd_mod + 256*hydro_mod)
* MOD TAH 051117: Set ocean tidal loading bit in gamit_mod if
*         the otide_mod is set
          if ( otide_mod.ne.0 ) call sbit(sgamit_mod,20,1)
      else
          write(*,'(a,1x,i4,1x,a,f5.2,1x,a35)') 
     .             'Error decoding Models line',
     .              jerr, 'SOLVE VER ',solve_ver/100., line
      endif
*     Atm Model line now follow
      if( solve_ver .ge. 1017 ) then
          read(unit,'(a)', iostat=ierr) line
!Atm models:  DryZen UFL   WetZen GP25  DryMap VMF1  WetMap VMF1  IonSrc GMAP  MagFld IGRF11
* MOD TAH 140403: Decode this line so that models are saved (allows
*     checking for inconsistencies later).
          read(line,280) cdryzen, cwetzen, cdrymap, cwetmap,
     .                   cionsrc, cmagfield
 280      format(21x,a4,4(9x,a4),9x,a8) 
      end if

      do while ( .not. finished )
 
*                                     ! Probably finished with this
          if( ierr.ne.0 .or. trimlen(line).eq.0 ) then
*                                 ! file
              finished = .true.
*                                 ! GEt the position of bias
          else
              indx = index(line,'biases')

*             if biases line found, read the next line:
              if( indx.ne.0 ) then
                  read(unit,'(a)', iostat=ierr) line
                  indx = 0
              end if

*             if biases not found, see if we can find keys:
              if( indx.eq.0 ) indx = index(line,'keys')

*                                     ! then problem, not correct line
              if( indx.eq.0 ) then

*                 Skip over the headers written by htoh (echo the lines
*                 to user just to remind them what has happened).
                  call skip_htoh( unit, line , ierr )
                  if( ierr.ne.0 ) then 
                      write(*,300) num_in_file+1,
     .                    line(1:max(1,trimlen(line)))
  300                 format(' Bias/keys line not where expected in',
     .                       ' solution ',i4,'.  Line read is:',/,a)
                      finished = .true.
                  end if
              end if

* MOD TAH 930208: Check the keys line to see what type of hfile this 
*             is.  Do only for GPS.
              if( index(line,'keys').gt.0 .and. 
     .            hfile_type.eq.'GPS' ) then
                    call get_gps_hext( line, hf_ext )
              else
                  hf_ext = ' '
              end if
          end if
 
****      Now start reading main part of file
          if ( .not. finished ) then
 
              num_in_file = num_in_file + 1	

*             The remaining headers are differnt between the GPS and
*             TER versions of the hfiles. Use different routines to
*             read the remaining header records.
              if( hfile_type.eq.'GPS' ) then
                   call read_rem_gpshead( unit, finished )
              else if( hfile_type.eq.'TER' ) then
                   call read_rem_terhead( unit, finished )
              else
                   write(*,310)  hfile_type
 310               format(' *** DISASTER *** Type of hfile lost at',
     .                    ' read remaining headers stage',/,
     .                    ' Current type is ',a,/,
     .                    '     Quitting on this hfile')
                    finished = .true.
              end if

*****         Now that we know the number of satellites and sites,
*             initialize the parameter number pointer arrays

              call init_qparn( qnum_sites, qnum_svs )
              call init_site(site_pos, site_vel, qnum_sites)
              call init_svs (svs_pos, qnum_svs, max_svs_elem)

              rad_known = .false.
              svant_known = .false.
 
*****         We are now ready to get the solution vector as well as
*             apriori values.  This routine also gets the live and dead
*             parameter numbers and the total number of parameters
              do i = 1, max_glb_parn
                 atoo(i) = 0
              end do
               
              call read_aprioris( unit )
 
*****         Now map the vma_data array for the number of parameters
*             which must be saved
 
              icov_parm = istart_vma
              isol_parm = icov_parm + qnum_parn*qnum_parn*2
              iend      = isol_parm + qnum_parn*2 - istart_vma + 1

*****         See if EOR turned on
              canal_type = 'S'
              if( qparn_pmu(1,1).ne.0 .or. qparn_pmu(1,2).ne.0 .or.
     .            qparn_pmu(1,3).ne.0 ) then
                  canal_type(trimlen(canal_type)+2:) = 'E'
              end if
              svant_est = .false.
              do i = 1, qnum_svs
* MOD TAH 210506: Fixed repeated -2 index; -2, -1, 0 now.
                if( qparn_svs(max_svs_elem-2,i).ne.0 .or.
     .              qparn_svs(max_svs_elem-1,i).ne.0 .or.
     .              qparn_svs(max_svs_elem  ,i).ne.0 ) 
     .              svant_est = .true.
              end do
              if( svant_est ) then
                 canal_type(trimlen(canal_type)+2:) = 'A'
              endif
 
*****         See if we have enough space
              if( iend.gt. len_vma ) then
                  write(*,500) iend*4/1024.**2, len_vma*4/1024.**2
 500              format(' ** DISASTER ** ',f8.3,' Mbytes required ',
     .                ' for this hfile and ',/, ' only ',f8.3,
     .                ' Mb assigned.  Resize with -m= option')
                  finished = .true.
              end if
*                         ! First Not Finished.
          end if
 
****      See if we should continue
          if( .not.finished ) then
 
*             Get the covarinace matrix.  This routine will also
*             reset the units.  Rotating the lat, long, radius
*             values to XYZ cartesian coordinates.is saved until
*             the covariance ,matrix is read (both converted at
*             the same time)
 
              call read_cov( unit, vma_data(icov_parm),
     .            vma_data(isol_parm) )

* MOD TAH 180202: See if option to apply rotation/translation de-consraints
*             for GAMIT BASELINE solutions has been set.
              if( (auto .and. .not. orb_est) .or.
     .             index(decon_str,'F').gt.0 ) then
                 
                 call decon_strapp(gamit_str, qnum_parn, 
     .               vma_data(icov_parm), vma_data(isol_parm) )
              endif
     
****          Read the next line and see if A priori covariance
*             matrix present
              read(unit,'(a)', iostat=ierr ) line
              if( index(line,'A priori').gt.0 .and. ierr.eq.0 ) then
                  call read_apr_cov(unit, line)
              end if
              if( ierr.eq.-1 ) finished = .true.
              

****          Now generate the apriori and parameter codes
 
              call qgen_codes
 
*****         Now write out the global file in GLOBK format.
 
              call gen_gname( hfile, sepoch, num_in_file, glb_dir,
     .                        glb_file, hf_ext, name_format,expt_code)
              write(*,600) glb_file(1:trimlen(glb_file))
 600          Format(' Creating Binary file ',a)

*****         Finish creating the GAMIT sln_def information
              call mk_sln_gamit
 
*****         Write out the header (create file if need be)
              call qw_glb_header
 
*             Write out the names
              call qw_glb_names
 
*             Write out the full names
              call qw_glb_full
 
*             Write out the solution description
              call qw_description
 
*             Write out the codes
              call qw_codes
 
*             Write out the apriori values
              call qw_aprioris
 
*             Write out the solution and covariance matrix
              call qw_soln( vma_data(istart_vma) )
 
*****         Thats all for this solution.  Now see if there is another
*             one
          end if
*                     ! Looping until finished
      end do
 
***** Thats all
      return
      end
 
CTITLE READ_REM_GPSHEAD
 
      subroutine read_rem_gpshead (unit, finished)
 
      implicit none

*     This rouitne will read the remaining headers upto the start of
*     the solution vector.  Finished is returned TRUE if there are
*     any problems.  Currently finished is not changed by this routine.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
 
      integer*4 unit
 
*   finished        - Indicates that this q file is completed
 
      logical finished
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   indx        - Pointer in string
*   date(5)     - Date of run
*   i           - Loop counter
*   num_obs_svs(max_glb_svs) - Number of observations on each satellite
*   num_obs_site(max_glb_sites) - Number of observations on each site.
*                 (used to elminate site if no data)
*   start_ep_num, end_ep_num - Start and stop epoch numbers (used to
*                 compute the real start and stop times of the data
*                 actually used in the GAMIT run.
*   dur_epoch   - Durataion of epoch (secs).  Used to adjust the start and
*                 stop times.
*   num_zen     - Number of zenith delays
*   ns          - Site number

      integer*4 ierr, trimlen, indx, date(5), i,
     .          num_obs_svs(max_glb_svs), num_obs_site(max_glb_sites),
     .          start_ep_num, end_ep_num, dur_epoch, num_zen, ns
 
*   sectag      - Seconds tag of runtime
*   L1_arp(3), L2_arp(3), arp_ecc(3) - Positions of L1, L2 and ARP from
*     ground mark
*   elev_cut   - Elevation cutoff
 
      real*8 sectag, L1_arp(3), L2_arp(3), arp_ecc(3), elev_cut
* MOD TAH 050622: Added
*   atmload(3) -- NEU atmospheric load (mm)
*   hydload(3) -- NEU hydorgraphic load (mm)
* MOD TAH 200205: Added antdaz
      real*4 antdaz  ! Antenna aligment from True N (deg).

      real*4 atmload(3), hydload(3)
 
*   code        - Code for site
*   ant_mod     - Antenna phase center model

      character*8  code       ! Site code
      character*16 ant_mod    ! Full Antex model name.

      character*4 nut    ! Check that Nuta appears on IC line
      character*1 offcon    ! Function to return constellation letter based on
                            ! offset PRN (original PRN for mod(prn,100))
      character*1 ssys   ! Satellite system for GNSS
     
*   line        - Line read from file
 
      character*256 line
 
***** read next line
 
      read(unit,'(a)', iostat=ierr) line

*     Skip down to number of sessions.  We will call this data
*     base version.
* MOD TAH 910110: Looped on Assumed lines since in some cases there is
*     is one or two (Maybe three but i haven't seen this)
* MOD TAH 950920: Changed loop until Number is found.
      line(2:8) = 'Assumed'
      do while ( line(2:7).ne.'Number' )
          read(unit,'(a)', iostat=ierr ) line
      end do 
 
*     Get number of sessions as version. (Start of line has changed so
*     should be Number of sessions.)
      read(line,'(21x,I8)', iostat=ierr) sversion
      call report_error('IOSTAT',ierr,'num session gett',
     .        line,0, 'make_glb')
 
*     Now get the number of stations
      read(unit,'(a)', iostat=ierr) line
      read(line,'(21x,i8)',iostat=ierr) qnum_sites
      call report_error('IOSTAT',ierr,'num site gett',
     .            line,0, 'make_glb')
 
*     Now get the number of satellites
      read(unit,'(a)', iostat=ierr) line
      read(line,'(16x,i8)',iostat=ierr) qnum_svs
      call report_error('IOSTAT',ierr,'num svs gett',
     .            line,0, 'make_glb')
 
*     Now skip down to the name of the Name of track stations
*     line
      do while ( index(line,'Name of track station').eq.0 ) 
          read(unit,'(a)', iostat=ierr) line
      end do 
 
      do i = 1, qnum_sites
          read(unit,'(a)', iostat=ierr) line
          indx = 5
          call GetWord(line, qsite_names(i), indx)

          read(line(indx+2:),50) qfull_names(i), num_obs_site(i),
     .                qrecv_ty(i), qrecv_sn(i), qrecv_fw(i),
     .                qante_ty(i), qradome_ty(i), qante_sn(i)
  50      format(a12,3x,i6,3x,a20,1x,a8,12x,1x,a8,1x,
     .           2x,a16,a4,1x,a8)

*         Clear character 6:8 of sn because these are not used in
*         sinex
          call trimlead(qante_sn(i))
          call trimlead(qrecv_sn(i))
          qante_sn(i)(6:8) = '  '
          qrecv_sn(i)(6:8) = '  ' 
       
          call sub_char(qfull_names(i)(1:trimlen(qfull_names(i))),
     .                  ' ', '_')
          qfull_names(i)(32:32) = 'P'
      end do
     
***** Now remove those sites with no data

      call compress_site( num_obs_site )
     
****  Now read the antenna position information (check that site is
*     still in the list
*     Skip two lines (Data files and code Cfile lines)
      read(unit,'(a)', iostat=ierr) line
      read(unit,'(a)', iostat=ierr) line
* MOD TAH 200205: Initialize antdaz for pre-310 h-files that do
*     not contain this value.
      antdaz = 0.0    
      do while (  line(2:10).ne.'Satellite' )
         read(unit,'(a)', iostat=ierr) line
         if( line(1:2).ne. ' S' ) then
             if( gamit_hf_ver.lt.200 ) then   ! Old format
                read(line, 60) code, arp_ecc, L1_arp, L2_arp, 
     .                         ant_mod, elev_cut, num_zen,
     .                         atmload, hydload
  60            format(6x,a4,18x,3(3(2x,f7.4)),
     .                 3x,a4,3x,f4.1,4x,i3,2x,3f7.2,1x,3f7.2)
             elseif( gamit_hf_ver.le.300 ) then
                read(line, 65) code, arp_ecc, L1_arp, L2_arp, 
     .                         ant_mod, elev_cut, num_zen,
     .                         atmload, hydload
  65            format(6x,a4,18x,3(3(2x,f7.4)),
     .                 3x,a16,2x,f4.1,4x,i3,2x,3f7.2,1x,3f7.2)
             else    ! Latest version 310 (TAH 200205)
                read(line, 70) code, arp_ecc, L1_arp, L2_arp, 
     .                         antdaz, ant_mod, elev_cut, num_zen,
     .                         atmload, hydload
* MOD TAH 200205: Added 1x,F5.1,1x, for antdaz
  70            format(6x,a4,18x,3(3(2x,f7.4)),1x,F5.1,1x,
     .                 3x,a16,2x,f4.1,4x,i3,2x,3f7.2,1x,3f7.2)
             end if

*            Find the station name
             indx = 1
             call get_cmd(code, qsite_names, qnum_sites, ns, indx)
             if( ns.gt.0 ) then
                 
                 qarp_ecc(1,ns) = arp_ecc(2)
                 qarp_ecc(2,ns) = arp_ecc(3)
                 qarp_ecc(3,ns) = arp_ecc(1)
                 
                 qL1a_ecc(1,ns) = L1_arp(2) - arp_ecc(2)
                 qL1a_ecc(2,ns) = L1_arp(3) - arp_ecc(3)
                 qL1a_ecc(3,ns) = L1_arp(1) - arp_ecc(1)
                 
                 qL2a_ecc(1,ns) = L2_arp(2) - arp_ecc(2)
                 qL2a_ecc(2,ns) = L2_arp(3) - arp_ecc(3)
                 qL2a_ecc(3,ns) = L2_arp(1) - arp_ecc(1)

                 qantdaz(ns)    = antdaz
                 
                 qant_mod(ns)   = ant_mod
                 qelev_cut(ns)  = elev_cut*pi/180.d0
                 qnum_zen(ns)   = num_zen

* MOD TAH 050622: Copy the loading values across
                 do i = 1,3
                    qatmload(i,ns) = atmload(i)
                    qhydload(i,ns) = hydload(i)
                 end do

              else
                 write(*,100) code
 100             format(' No data for site ',a)
              end if
          end if
      end do              
       
****  Tell user a little of what is going on
      Write(*,110) qnum_sites, hfile(1:trimlen(hfile))
 110  format(' There are ',i3,' sites in ',a,/,
     .       5x,'Name',6x,'Full name')
      do i = 1, qnum_sites
          write(*,120) i, qsite_names(i), qfull_names(i),
     .                 num_obs_site(i), qrecv_fw(i), qarp_ecc(3,i)
 120      format(1x,i3,1x,a8,2x,a16,1x,i6,2x,a12, 2x,f9.4)
      end do
 
***** Now skip down until Satellite used:
      do while ( line(2:10).ne.'Satellite' )
          read(unit,'(a)', iostat=ierr) line
      end do
 
****  Now get the satellite information.
* MOD TAH 101015: Base the reading of satellite information on file version
      if( gamit_hf_ver.lt.200 ) then
          call get_svs_inf( unit, num_obs_svs )
      else
          call get_svs_new( unit, num_obs_svs )
      endif

****  Tell user what is happening
* MOD TAH 180401: Mover code down because we can't get SVN numbers
*     until we know the epcoh of the h-file which is read later.
!      Write(*,140) qnum_svs, hfile(1:trimlen(hfile))
 140  format(' There are ',i3,' satellites in ',a,/, 5x,'Name')
!      do i = 1, qnum_svs
!          write(*,150) i, qsvs_names(i), num_obs_svs(i)
 150      format(1x,i3,1x,a8,2x,I6)
!      end do

****  Now compress those satellites not be used from list
      call compress_svs(  num_obs_svs )    

*     Now get number of observations total.
      read(unit,'(a)', iostat=ierr) line
      read(line,'(33x,i8,15x,i5,3x,i5,13x,i6)', iostat=ierr) snum_obs, 
     .               start_ep_num, end_ep_num, dur_epoch
 
***** Get the mid-epoch of the session.  (We use this to create
*     the file name as well)
      read(unit,'(a)', iostat=ierr)  line
      read(line,170) date, sectag, sgtime, sgpst_utc
 170  format(13x, i4,2x,i2,2x,i2,2x,i2,2x,i3,2x,f7.3,2x,a4,14x,i2)    
      call ymdhms_to_jd( date, sectag, qstart_epoch)
      read(unit,'(a)', iostat=ierr) line
      read(line,'(13x,i4,2x,i2,2x,i2,2x,i2,3x,i2,3x,f5.3)') date,sectag
      call ymdhms_to_jd( date, sectag, qend_epoch)

* MOH TAH 180401: Now we know date, get SVNumbers for satellites
      call Get_SVnum( qstart_epoch )

* MOD TAH 180401: Generate satellites names with SVN numbers
      do i = 1, qnum_svs
         ssys = offcon(qsvi_prn(i))
         write(qsvs_names(i),320) ssys, qsvi_svn(i),
     .                         mod(qsvi_prn(i),100)
 320     format(a1,I3.3,'_',i2.2)
      end do

*     New tell user what we have
      Write(*,140) qnum_svs, hfile(1:trimlen(hfile))
      do i = 1, qnum_svs
          write(*,150) i, qsvs_names(i), num_obs_svs(i)
      end do

****  Now adjust the start and stop times for the duration of data used.
      qend_epoch = qstart_epoch + (end_ep_num-1)*dur_epoch/86400.d0
      qstart_epoch = qstart_epoch + (start_ep_num-1)*dur_epoch/86400.d0
 
      sepoch = ( qstart_epoch + qend_epoch ) / 2
      qgpst_utc = sgpst_utc

****  Last peice of information that we need is the SV number.  We need 
*     read svnav.dat to get this information
      if( solve_ver .ge. 1017 ) then
          call Get_SVnum(qstart_epoch)
      end if

*     Now get the IC time 
      read(unit,'(a)', iostat=ierr) line
* MOD TAH 080103: Added cgnut, cggrav to line
* MOD TAH 140327: Added ceradmod and cantradmod to values read
*     from IC line.
      read(line,190) date,sectag, sgframe, sgprec, sgsrpmod, 
     .               nut, cgnut, cggrav, ceradmod, cantradmod
 190  format(13x, 4i4,2x,i3,2x,f7.3,3x,8x,
     .       a5,15x, a5, 13x ,a5, 2x,a4, 6x, a5, 17x, a5,
     .          25x, a5, 27x, a5 )
* MOD TAH 070926: Changed format so that gfortran would work.
c190  format(' ICs time  : ', 4i4,2x,i3,2x,f7.3,3x,' Frame: ',
c    .       a5,'   Precession: ',a5, '  SRP model: ',a5)
c MOD TAH 080102: New line being read
c ICs time  : 2007  12  24  12    0    0.000    Frame: J2000   Precession: IAU76  SRP model: BERNE  Nutation: IAU00  Gravity model: EGM96
c    .  ,'  Gravity model: ',gravmod
c    .  ,'  Earth-radiation model: ',eradmod
c    .  ,'  Antenna-radiation model: ',antradmod
cc     write(*,'(a)') line(1:90)
c     write(*,190) date,sectag, sgframe, sgprec, sgsrpmod
*     Make sure Nuta read
      if( nut.ne.'Nuta' ) then
         cgnut = 'IAU80'
         cggrav = 'EGM96'
      endif
     
      call ymdhms_to_jd( date, sectag, ssvs_epoch)

      if( sgframe(1:1) .eq.' ' ) sgframe  = 'B1950'
      if( sgprec(1:1)  .eq.' ' ) sgprec   = 'IAU68'
      if( sgsrpmod(1:1).eq.' ' ) sgsrpmod = 'SPHRC'
      if( sgtime(1:1)  .eq.' ' ) sgtime   = 'UTC'
* MOD TAH 140327: See if new format with erad and antenna thrust model
      if( ceradmod(1:1).eq.' ' ) ceradmod = 'UNKNOWN'
      if( cantradmod(1:1).eq.' ' ) cantradmod = 'UNKNOWN'
* MOD TAH 140403: See if models were read
      if( cdryzen(1:1).eq.' ') cdryzen = 'UNKN'
      if( cwetzen(1:1).eq.' ') cwetzen = 'UNKN'
      if( cdrymap(1:1).eq.' ') cdrymap = 'UNKN'
      if( cwetmap(1:1).eq.' ') cwetmap = 'UNKN'
      if( cionsrc(1:1).eq.' ') cionsrc = 'UNKNOWN'
      if( cmagfield(1:1).eq.' ') cmagfield = 'UNKNOWN'

      write(*,200) sgframe, sgprec, sgsrpmod, sgtime, cgnut, cggrav,
     .             ceradmod, cantradmod
 200  format(' GAMIT Hfile using ',8a8)
      write(*,210) cdryzen, cwetzen, cdrymap, cwetmap, 
     .             cionsrc, cmagfield
 210  format(' GAMIT Atm models  ',4(a4,1x), 
     .       ' 2nd Order Ion ',2(a8,1x))
 
*     Read the next line.  If it is the ICEpcoh line then this
*     is new format file.  If not skip over the next sections.
      read(unit,'(a)', iostat=ierr) line
      if( line(1:9).eq.' IC epoch' .or. 
     .    line(1:6).eq.' Earth'          ) then
          call read_orient( unit )
      end if

****  Other quanitites to get are:
*     stai_utc,
*     sut1_apr(2)
*     swob_apr(2,2)
*     snut_ang_apr(2,2)
*     secc_change(3,max_sites)  -- real*4.  There is slight
*     problem here since there can be many more sites in GPS
*     experiments than VLBI we do not have enough space for all
*     of the eccentricities.  Could extend the length of
*     sln_def record which would not be too hard to implement.
 
*
***** Get the number of parameters in solution.  Actually ingore this
*     line since the parameters numbers do not actually reflect what
*     is being output,
      read(unit,'(a)', iostat=ierr) line

*     Get the NRMS for solution      
      read(unit,'(a)', iostat=ierr) line
      read(line,'(44x,f12.5)') cchisq
      cchisq = cchisq**2
      write(*,220) sqrt(cchisq)
 220  format(' GAMIT Nrms for solution ',f12.5)
      
****  Skip Until we find the label line
      do while (index(line,'label').eq.0 )
          read(unit,'(a)', iostat=ierr) line
      end do
c     read(unit,'(a)', iostat=ierr) line
c     read(unit,'(a)', iostat=ierr) line
 
***** Thats all.  Now we are ready to read soltion
 
      return
      end
 
CTITLE READ_REM_TERHEAD
 
      subroutine read_rem_terhead (unit, finished)
 
      implicit none

*     This rouitne will read the remaining headers upto the start of
*     the solution vector.  Finished is returned TRUE if there are
*     any problems.  Currently finished is not changed by this routine.
*     This routine is for reading TER hfile headers
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
 
      integer*4 unit
 
*   finished        - Indicates that this q file is completed
 
      logical finished
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   indx        - Pointer in string
*   date(5)     - Date of run
*   i           - Loop counter
*   num_obs_svs(max_glb_svs) - Number of observations on each satellite
*   num_obs_site(max_glb_sites) - Number of observations on each site.
*                 (used to elminate site if no data)
*   start_ep_num, end_ep_num - Start and stop epoch numbers (used to
*                 compute the real start and stop times of the data
*                 actually used in the GAMIT run.
*   dur_epoch   - Durataion of epoch (secs).  Used to adjust the start and
*                 stop times.
*   cnt         - Generic counting variable to be used when reading in loops
*                 to make sure that we do not go to far.
 
      integer*4 ierr, trimlen, indx, date(5), i,
     .          num_obs_site(max_glb_sites), cnt
 
*   sectag      - Seconds tag of runtime
 
      real*8 sectag
 
*   line        - Line read from file
 
      character*256 line
 
***** for Terrestrial data just skip the next  block of lines until we
*     get to Number
*     Skip down to number of sessions.  We will call this data
*     base version.
*     is one or two (Maybe three but i haven't seen this)
      cnt = 0
      read(unit,'(a)', iostat=ierr ) line
      do while ( line(2:7).ne.'Number' .and. cnt.lt.100 )
          read(unit,'(a)', iostat=ierr ) line
          cnt = cnt + 1
      end do 

****  See if have an over-run
      if( cnt.ge.100 ) then
          write(*, 80) cnt, line(1:max(1,trimlen(line)))
  80      format(' ** DISASTER ** Overrun on reading remaining ',
     .           '  terrestrial headers.  Last line read: ',/,a)
          finished = .true.
          RETURN
      end if
 
*     Get number of sessions as version. (Start of line has changed so
*     should be Number of sessions.)
      read(line,'(21x,I8)', iostat=ierr) sversion
      call report_error('IOSTAT',ierr,'num session gett',
     .        line,0, 'make_glb')
 
*     Now get the number of stations
      read(unit,'(a)', iostat=ierr) line
      read(line,'(21x,i8)',iostat=ierr) qnum_sites
      call report_error('IOSTAT',ierr,'num site gett',
     .            line,0, 'make_glb')
 
*     set the number of  satellites to zero.
      qnum_svs = 0
 
*     Skip next line to get the site names
      read(unit,'(a)', iostat=ierr) line
 
      do i = 1, qnum_sites
          read(unit,'(a)', iostat=ierr) line
          indx = 5
          call GetWord(line, qsite_names(i), indx)

*         Use GetWord to allow for more flexiable formats
C         read(line(indx+2:),90) qfull_names(i), num_obs_site(i)
C 90      format(A16,I8)
          call GetWord(line, qfull_names(i), indx )
          read(line(indx+2:),*) num_obs_site(i)
          call sub_char(qfull_names(i)(1:trimlen(qfull_names(i))),
     .                  ' ', '_')
      end do

***** Now remove those sites with no data

      call compress_site( num_obs_site )
 
****  Tell user a little of what is going on
      Write(*,100) qnum_sites, hfile(1:trimlen(hfile))
 100  format(' There are ',i3,' sites in ',a,/,
     .       5x,'Name',6x,'Full name')
      do i = 1, qnum_sites
          write(*,120) i, qsite_names(i), qfull_names(i),
     .                 num_obs_site(i)
 120      format(1x,i3,1x,a8,2x,a16,1x,i6)
      end do
 
***** Now skip down until number of observations
      cnt = 0
      do while ( line(2:7).ne.'Number' .and. cnt.le.10)
          read(unit,'(a)', iostat=ierr) line
          cnt = cnt + 1
      end do

****  See if have an over-run
      if( cnt.ge.10 ) then
          write(*,125) cnt, line(1:max(1,trimlen(line)))
 125      format(' ** DISASTER ** Overrun on reading remaining ',
     .           '  terrestrial headers.  Last line read: ',/,a)
          finished = .true.
          RETURN
      end if
 
*     Now get number of observations total.
      read(line,'(25x,i10)', iostat=ierr) snum_obs
 
***** Get the mid-epoch of the session.  (We use this to create
*     the file name as well)
      read(unit,'(a)', iostat=ierr)  line
      read(line,'(13x,i4,2x,i2,2x,i2,2x,i2,3x,i2,3x,f5.3)') date,sectag
      call ymdhms_to_jd( date, sectag, qstart_epoch)
      read(unit,'(a)', iostat=ierr) line
      read(line,'(13x,i4,2x,i2,2x,i2,2x,i2,3x,i2,3x,f5.3)') date,sectag
      call ymdhms_to_jd( date, sectag, qend_epoch)

*     Now get the reference time (for positions when velocities
*     estimated
      read(unit,'(a)', iostat=ierr) line
      read(line,'(13x,i4,2x,i2,2x,i2,2x,i2,3x,i2,3x,f5.3)') date,sectag
      call ymdhms_to_jd( date, sectag, sepoch)

***** Get the number of parameters in solution.  Actually ingore this
*     line since the parameters numbers do not actually reflect what
*     is being output,
      read(unit,'(a)', iostat=ierr) line
 
****  Skip Until we find the label line
      do while (index(line,'label').eq.0 )
          read(unit,'(a)', iostat=ierr) line
      end do
c     read(unit,'(a)', iostat=ierr) line
c     read(unit,'(a)', iostat=ierr) line
 
***** Thats all.  Now we are ready to read soltion
 
      return
      end

CTITLE READ_APRIORIS
 
      subroutine read_aprioris( unit )
 
      implicit none

*     This routine reads the apriori values and adjustments from
*     the hfile.  While doing this it builds up the list of apriori
*     and parameter codes, and the parn arrays.
 
      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
 
      integer*4 unit
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   indx        - Pointer in string
*   np          - Running number of parameters actually
*               - estimated
*   sn          - Site number
*   vn          - Satellite vechile number
*   i           - Loop counter
 
      integer*4 ierr, trimlen, np

*   line_preread  - Logical that indicates that line from the
*                   file has already been read  in read_orbit
*                   while it searches for radiation parameters

      logical line_preread
 
*   line        - Line read from file
      character*256 line
 
***** Start reading the paramters until a blank line is found
 
      ierr = 0
      np = 0
      line_preread = .false.
 
      do while ( ierr.eq.0 )
 
          if( .not. line_preread) read(unit,'(a)', iostat=ierr ) line
          call report_error('IOSTAT',ierr,'solution read',hfile,0,
     .                'read_aprioris')
          if( trimlen(line).eq.0 ) then
*                             ! Finished with parameter estimates
              ierr = -1
          end if
 
*         If no error and not finished then decode
          if( ierr.eq.0 ) then
 
*****         See type of parameter (Only check for start of GEOC and
*             orbit.  All other parameters are ignored, and when these
*             are found a set sequence of reading is started.)
 
*****         Now based on type read the line
*                                         ! Read site lat/long/rad
              if( index(line,'GEOC LAT').gt.0  ) then

                  call read_geoc(unit, line, np )
*                             ! Site position.
              end if

* MOD TAH 040703: Check to see if ATM delay entry present
              if( index(line,'ATMZEN').gt.0 ) then
                  call read_atm(unit, line, np)
              end if

              if( index(line,'geoc x').gt.0 ) then
                  call read_xyz(unit, line, np)
              end if

              if( index(line,'Vx').gt.0 ) then
                  call read_vel(unit, line, np)
              end if

*             See if polar motion or UT1
              if( index(line,'X POLE').gt.0 .or.
     .            index(line,'Y POLE').gt.0 .or. 
     .            index(line,'UT1-TAI').gt.0 ) then
                  call read_pmu(unit, line, np )
              end if
 
*****         See if orbital parameter
*                                             ! Read Orbital elements
              if( index(line,'ORBIT').gt.0 ) then

                  call read_orbit(unit, line, np, line_preread)

              end if

*****         See if antenna offset for satellite found
              if( index(line,'SVANT').gt.0 ) then
                  call read_svant(unit, line, np, line_preread)  
              end if
 
*                                 ! Not blank line or file error
          end if
*                                 ! looping finding all the
*                                 ! parmeters
      end do
 
****  Thats all
      qnum_parn = np
      write(*,300) qnum_parn
  300 format(' Found ', I4,' parmeters estimated in solution')
      return
      end
 
CTITLE READ_COV
 
      subroutine read_cov( unit, cov_parm, sol_parm )
 
      implicit none

*     This routine reads the covarinance matrix. Once it is read,
*     it and the solution vector are scaled and rotated (for the
*     site coordinates)
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
 
      integer*4 unit
 
*   cov_parm(qnum_parn, qnum_parn)    - Covariance matrix
*   sol_parm(qnum_parn)          - Solution vector
 
      real*8 cov_parm(qnum_parn, qnum_parn), sol_parm(qnum_parn)
 
* LOCAL VARIABLES
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   indx        - Pointer in string
*   i,j         - Loop counter
*   nq, np      - Parmeter number for rotating solution and
*               - covariance matrix
*   num_orb_blks  - Number of orbital element blocks that need
*                 to be rotated for LLR coordinates
*   sb          - Start of block for orbital element rotations.
*   num_ant_blks  - Number of antenna offset elements that need to
*                 be rotated.

 
      integer*4 ierr, trimlen, indx, i,j,k, nq, np, num_orb_blks,
     .          sb,  nb, num_ant_blks

*   umat        - 3*3 unit matrix used for rotating the covariance
*                 elements between the station coordinates and the
*                 satellite orbits.  (For convenience we do the
*                 sv elements 3 at a time.  There are 9 for sv)

       real*8 umat(3,3)
 
*   line        - Line read from file
 
      character*256 line,message

      data umat / 1.d0, 0.d0, 0.d0,  0.d0, 1.d0, 0.d0, 
     .            0.d0, 0.d0, 1.d0 /
 
****  Skip the next (check to see that it says covariance)
      read(unit,'(a)', iostat=ierr ) line
      indx = index(line,'Covariance')
*                             ! We have a problem
      if( indx.eq.0 ) then
          write(*,100) line(1:max(1,trimlen(line)))
 100      format(' ** DISASTER ** Covariance matrix not where',
     .        ' expected.  Line read is:',/,a)
          RETURN
      end if
 
***** Now loop over all of the parameters which have been estimated
*     Read the column of the matrx
      do i = 1, qnum_parn
          read(unit,200, iostat=ierr) (cov_parm(j,i), j= 1,i)
 200      format(4x,2x,5000(5(1x,d23.16),:,/,6x))
          nb = atoo(i)  
          if( ierr.ne.0 ) then
              cov_parm(i,i) = 1.d0  
              write(message,'(a,2i5,a)')
     .             'Error reading row ',i,nb,' of covariance matrix'    
              call report_stat('WARNING','HTOGLB','htog_ut',' '
     .                        ,message,ierr)
          end if
      end do

****  Now fill in the symetric part of the matrix (This code actually
*     copies the diagonal as well but this should be no problem)
      do i = 1, qnum_parn
          call dwmov(cov_parm(1,i),1, cov_parm(i,1), qnum_parn,i)
      end do

****  Now rescale the matrix (meter everywhere except solar radiation
*     which will be ratio to direct effect)
      do i = 1, qnum_parn
          call dwsmy(qscale(i), cov_parm(1,i),1, 
     .                          cov_parm(1,i),1,  qnum_parn)
          call dwsmy(qscale(i), cov_parm(i,1),qnum_parn,
     .                          cov_parm(i,1),qnum_parn, qnum_parn)
 
          sol_parm(i) = qsol(i)*qscale(i)
 
      end do

****  Now rotate the station coordinates from local frame to cartesian
*     XYZ.  If needed

      if( datum_type.eq.'LLR' ) then 

*         Get the numbers of orbital parameters that need rotating
          num_orb_blks = (rad_num+6)/3
          num_ant_blks =  svant_num/3

          do i = 1, qnum_sites
 
*             Rotate the site coordinates
              call llr_to_xyz( llr(1,i), site_pos(1,i), rot_mat(1,1,i))
 
          end do
 
*****     Now rotate the solution vector and the covariance matrix
          do i = 1, qnum_sites
 
*             Solution vector
              nq = qparn_sites(1,i)

              if( nq.ne.0 ) then
                  call rotate_sol( rot_mat(1,1,i), sol_parm(nq))
 
*****             Now rotate the covariance matrx
                  do j = 1, qnum_sites
                      np = qparn_sites(1,j)
                      if( np.gt.0 ) then
                          call rotate_cov(rot_mat(1,1,i),
     .                        cov_parm(nq,np), qnum_parn,
     .                        rot_mat(1,1,j), 3, 3 )
                      end if
                  end do

*****             Now do the ATM  estimates.  Do one at
*                 a time.
                  do j = 1, qnum_sites
                      np = qparn_atm(j)
                      if( np.gt.0 ) then
                            call rotate_cov(rot_mat(1,1,i),
     .                           cov_parm(nq,np), qnum_parn,
     .                           umat, 3, 1 )
                      end if
                  end do

*****             Now do the rotation of the satellite/site covariance 
*                 elements.  Do in groups of three moving accross the
*                 rows.  When we are finished we will do the move the
*                 values to the columns.
                  do j = 1, qnum_svs
*                    There could be upto 5 groups of parameters for
*                    maximum of 15 orbit parameters.
                     np = qparn_svs(1,j)
                     if( np.gt.0 ) then
                         do k = 1,rad_num + 6
                            sb = np + (k-1)
                            call rotate_cov(rot_mat(1,1,i),
     .                           cov_parm(nq,sb), qnum_parn,
     .                           umat, 3, 1 )
                         end do
                     end if
*                    Now see if we need to rotate the svant block
                     np = qparn_svs(max_svs_elem-2,j) 
                     if( np.gt.0 .and. num_ant_blks.eq.1 ) then
                         call rotate_cov(rot_mat(1,1,i),
     .                        cov_parm(nq,np), qnum_parn,
     .                        umat, 3, 3 )
                     end if
                  end do

*****             Now do the polar motion/Ut1 estimates.  Do one at
*                 a time since we don't know how many there will be
                  do j = 1, 3
                     do k = 1,2
                        np = qparn_pmu(k,j)
                        if( np.gt.0 ) then
                            call rotate_cov(rot_mat(1,1,i),
     .                           cov_parm(nq,np), qnum_parn,
     .                           umat, 3, 1 )
                        end if
                     end do
                  end do

              end if
*                     ! Looping over sites
          end do

****      Now fill in the symetric part of the matrix (This code actually
*         copies the diagonal as well but this should be no problem)
	  
          do i = 1, qnum_parn
              call dwmov(cov_parm(1,i),1, cov_parm(i,1), qnum_parn,i)
          end do
*                     ! coordinates need rotating
      end if
      
****  Now construct the total solution (add apriori and adjustment
 
      do i = 1,qnum_sites
          do j = 1,3
              nq = qparn_sites(j,i)
              if( nq.gt.0 ) then
                  sol_parm(nq) = site_pos(j,i) + sol_parm(nq)
              end if
              nq = qparn_vel(j,i)
              if( nq.gt.0 ) then
                  sol_parm(nq) = site_vel(j,i) + sol_parm(nq)
              end if
              
          end do
*         Generate total atmospheric delay
          nq = qparn_atm(i)
          if( nq.gt.0 ) then
              sol_parm(nq) = site_atm(i) + sol_parm(nq)
          endif
      end do
 
***** Repeat for satellites
      do i = 1,qnum_svs
          do j = 1,max_svs_elem
              nq = qparn_svs(j,i)
              if( nq.gt.0 ) then
                  sol_parm(nq) = svs_pos(j,i) + sol_parm(nq)
              end if
          end do
      end do

***** Now add adjusmtent to get total values
      do i = 1,3
         do j = 1,2
            nq = qparn_pmu(j,i)
            if( nq.gt.0 ) then
               sol_parm(nq) = pmu_pos(j,i) + sol_parm(nq)
            end if
         end do
      end do

****  Thats all.  We are now ready to generate the solution
      return
      end
 
CTITLE LLR_TO_XYZ
 
      subroutine llr_to_xyz ( llr, xyz, rot_mat )
 
      implicit none

*     Routine to convert from lat, long and radius to XYZ cartersian
*     and to returm the rotation matrix from llr to xyz.
 
* PASSED VARIABLES
 
*   llr(3)  - lat, long, radius (rad, rad, kms)
*   xyz(3)  - XYZ geocentric (m)
*   rot_mat(3,3)    - Rotation from NEU ro XYZ
 
 
      real*8 llr(3), xyz(3), rot_mat(3,3)
 
****  Get the XYZ coordinates
      xyz(1) = llr(3)*cos(llr(1))*cos(llr(2)) * 1000.d0
      xyz(2) = llr(3)*cos(llr(1))*sin(llr(2)) * 1000.d0
      xyz(3) = llr(3)*sin(llr(1)) * 1000.d0
 
****  Now compute the rotation matrix
      rot_mat(1,1) = -sin(llr(1)) * cos(llr(2))
      rot_mat(1,2) = -sin(llr(2))
      rot_mat(1,3) =  xyz(1)/llr(3) / 1000.d0
 
      rot_mat(2,1) = -sin(llr(1)) * sin(llr(2))
      rot_mat(2,2) =  cos(llr(2))
      rot_mat(2,3) =  xyz(2)/llr(3) / 1000.d0
 
      rot_mat(3,1) =  cos(llr(1))
      rot_mat(3,2) =  0.d0
      rot_mat(3,3) =  xyz(3)/llr(3) / 1000.d0
 
****  Thats all
      return
      end
 
CTITLE ROTATE_COV
 
      subroutine rotate_cov(R1, cov, dim, R2, dim2, ns )
 
      implicit none

*     Routine to rotate the covariance matrix.  Cross element between
*     site 1 with rotation R1 and site site 2 with rotation R2
 
* PASSED VARIABLES
 
*   dim     - The row dimension of covaraince matrix
*   dim2    - dimension of R2 (assumed square)
*   ns      - number of parmaters to rotate
 
      integer*4 dim, dim2, ns
 
*   R1(3,3), R2(dim2,dim2)    - The rotation matrices for each site
*   cov(dim,ns)      - The cross covariance matrix (upper
*                   - right hand cornore pointed to)
 
 
      real*8 R1(3,3), R2(dim2, dim2), cov(dim,ns)
 
* LOCAL VARIABLES
 
*   i,j,k       - Loop couners
 
      integer*4 i,j,k
 
*   cov_R2T(3,3)    - Product of covariance matrix by R2 transposed
 
      real*8 cov_R2T(3,3)
 
****  Start.  First do:
*                      T
*     COV_R2T = Cov* R2
*
      do i = 1,3
          do j = 1,ns
              cov_R2T(i,j) = 0.d0
              do k = 1,ns
                  cov_R2T(i,j) = cov_R2T(i,j) + cov(i,k)*R2(j,k)
              end do
          end do
      end do
 
*     Now complete with R1 * cov_R2T
 
      do i = 1,3
          do j = 1,ns
              cov(i,j) = 0.d0
              do k = 1,3
                  cov(i,j) = cov(i,j) + R1(i,k)*cov_R2T(k,j)
              end do
          end do
      end do
 
***** Thats all
      return
      end
 
CTITLE ROTATE_SOL
 
      subroutine rotate_sol( R1, sol )
 
      implicit none

*     Routine to apply the R1 solutiuon matrix to solution vector
*
* PASSED VARIABLES
 
*   R1(3,3)     - The rotation matrice for this site
*   sol(3)      - The solution vector (adjustment)
 
      real*8 R1(3,3), sol(3)
 
* LOCAL VARIABLES
 
*   i,j     - Loop couners
 
      integer*4 i,j
 
*   R1_sol(3)   - Product of R1*sol
 
 
      real*8 R1_sol(3)
 
****  Loop over each element
      do i = 1,3
          R1_sol(i) = 0.d0
          do j = 1,3
              R1_sol(i) = R1_sol(i) + R1(i,j)*sol(j)
          end do
      end do
 
****  Now return the new values
      do i = 1,3
          sol(i) = R1_sol(i)
      end do
 
***** Thats all
      return
      end
 
CTITLE GEN_GNAME
 
      subroutine gen_gname( hfile, sepoch, num_in_file, glb_dir,
     .                      glb_file, hf_ext, name_format,expt_code)

      implicit none
 
*     This routine will generate the name of the output global
*     file from the date of the experiment
 
* PASSED VARIABLES
 
*   num_in_file - Bumber of covariance matrix in file.
 
      integer*4 num_in_file
 
*   sepoch      - JD of the experiment (mid-point)
 
      real*8 sepoch
 
*   hfile       - NAme of the input hfile
*   glb_dir     - Directory for output
*   glb_file        - Full name of output file.
*   hf_ext      - OPtional extent to be used on hfile name.
*                 (Read from keys line of GPS h-file)
*   name_format - Optional name format specification
*   expt_code   - Optional experiment name to be used
 
      character*(*) hfile, glb_dir, glb_file, hf_ext, name_format,
     .              expt_code
 
* LOCAL VARIABLES
 
*   len_name        - Generic lenght of name
*   len_expt    - Length of expt_code and expt_string
*   trimlen     - returns length of string
*   indx_q      - Index for the first . in the hfile name.
*   date(5)     - ymdhm date form
*   gps_week, gps_dow  - gps week number and doy of week
*
*   
*   j           - loop counter
 
      integer*4 len_name, trimlen, date(5), j, indx_q, gps_week,
     .          gps_dow, len_expt
 
*   sectag      - Seconds tag
*   nepoch      - Epoch to be used in naming files (rounded to
*                 nearest minute)
*   gpss        - JD of start of GPS week
 
      real*8 sectag, nepoch, gpss

*   expt_string - String to be added to name.  If not passed with
*                 -F option, it is derived from the input file name.

      character*8 expt_string

      data gpss / 2444244.50d0 / 

 
*     Start generating name
      glb_file = glb_dir
      len_name = trimlen(glb_file)
*                                 ! See if / there
      if( len_name.gt.0 ) then
          if( glb_file(len_name:len_name).ne.'/' ) then
              glb_file(len_name+1:) = '/'
              len_name = len_name + 1
          end if
      end if
 
*     Add 'h' as first letter
      len_name = len_name + 1
      glb_file(len_name:) = 'h'
 
*     Now add date:  See what the user wants as a file name format
      nepoch = nint((sepoch-gpss)*144.d0)/144.d0 + 1.d0/86400.d0 +
     .         gpss
      call jd_to_ymdhms( nepoch, date, sectag )
      if( name_format(1:4).eq.'GPSW' ) then
          gps_week = (nepoch-gpss)/7.d0
          gps_dow  = nepoch - gps_week*7 - gpss
          write(glb_file(len_name+1:),110) gps_week, gps_dow,
     .                                     (date(j),j=4,5)
 110      format(I4.4,i1,'_',2i2.2,'_')
      else if( name_format(1:5).eq.'Y2DOY' ) then
          call jd_to_yds( nepoch, date(1), date(2), date(3))
          date(1) = date(1) + 1900
          if( date(1).lt.2000 ) then
              date(1) = date(1) - 1900
          else
              date(1) = date(1) - 2000
          end if
          write(glb_file(len_name+1:),120) (date(j),j=1,2),
     .                                     (date(j),j=4,5)
 120      format(I2.2,i3.3,'_',2i2.2,'_')

      else if( name_format(1:5).eq.'Y4DOY' ) then
          call jd_to_yds( nepoch, date(1), date(2), date(3))
* MOD TAH 100302: test for Y2000
          if( date(1).lt.50 ) then
              date(1) = date(1) + 2000
          else
              date(1) = date(1) + 1900
          end if

          write(glb_file(len_name+1:),130) (date(j),j=1,2),
     .                                     (date(j),j=4,5)
 130      format(I4.4,i3.3,'_',2i2.2,'_')
      else if( name_format(1:2).eq.'Y4' ) then
           write(glb_file(len_name+1:),140) (date(j),j=1,5)
 140       format(i4.4,4i2.2,'_') 
*          Default to original case if none specified
      else
          if( date(1).lt.2000 ) then
              date(1) = date(1) - 1900
          else
              date(1) = date(1) - 2000
          end if
          write(glb_file(len_name+1:),150) (date(j),j=1,5)
 150      format(5i2.2,'_')
      end if
 
      len_name = trimlen(glb_file)

*     Now get character before '.' in name
      call get_indx_q( hfile, indx_q)

* MOD TAH 930208 to use optional extent if available.
* MOD TAH 050515: See if expt_code was passed with -F option
*     and base the string name on this.
      expt_string = hfile(indx_q-5:indx_q-2)
      if ( hf_ext(1:3).eq.'glj' ) 
     .     expt_string = hfile(indx_q-4:indx_q-1)
*     See if we need to change because expt_code passed
      if( trimlen(expt_code).gt.0 ) then
          expt_string = expt_code
      end if
      len_expt = trimlen(expt_string)
         

      if( trimlen(hf_ext).eq.0 ) then
          glb_file(len_name+1:) = expt_string(1:len_expt)  // '.gl'
          len_name = trimlen(glb_file)
          glb_file(len_name+1:) = char(97+num_in_file)
      else 
          if( hf_ext(1:3).eq.'glj') then
             glb_file(len_name+1:) = expt_string(1:len_expt) // '.'
     .                            // hf_ext(1:3)
          else
             glb_file(len_name+1:) =  expt_string(1:len_expt) // '.'
     .                            // hf_ext(1:3)
          endif
      end if
 
****  That all
      return
      end
 
CTITLE QGEN_CODES
 
      subroutine qgen_codes
 
      implicit none

*     This routine will scan the parmater entries generated while
*     the hfile was being read, and generate the apriori and estimated
*     parameter codes.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* LOCAL VARIABLES
 
*   i,j,k   - Loop counter
*   na      - number of apriori codes (working value)
*   np      - count of the number of parameter codes
*   nq      - Actual parameter position in covariance matrix
*   ne      - Counter for number of epoch depedent parameters
*   nf      - Occ number for each multi-day pmu
 
      integer*4 i,j, k, na, np, nq, ne, nf
 
****  Loop over the stations first
      na = 0
      np = 0
      ne = 0

* MOD TAH 070110: Initialize array in case some parameters have
*     no codes.
      do i = 1, qnum_parn
         qglb_codes(i) = 0
         qapr_codes(i) = 0
      end do
 
      do i = 1, qnum_sites
          do j = 1,3
              if( qparn_sites(j,i).gt.0 ) then 
                 na = na + 1
                 qapr_list(na) = site_pos(j,i)
                 qapr_codes(na) = 6+j + i*256
                 np = np + 1
                 nq = qparn_sites(j,i)
                 qglb_codes(nq) = 6+j + i*256
              end if
          end do

*         Now check the velocities
          do j = 1,3
              if( qparn_vel(j,i).ne.0 ) then
                  na = na + 1
                  qapr_list(na) = site_vel(j,i)
                  qapr_codes(na) = 41+j + i*256
                  np = np + 1
                  nq = qparn_vel(j,i)
                  qglb_codes(nq) = 41+j + i*256
              end if
          end do

*         Now check for axis offsets
 	  if( qparn_axo(i).ne.0 ) then
              na = na + 1
              qapr_codes(na) = 10 + i*256
              qapr_list(na) = site_axo(i)
	           np = np + 1
	           nq = qparn_axo(i)
	           qglb_codes(nq) = 10 + i*256
          end if

      end do

* MOD TAH 040703: Add apriori and values for atmospheric
*     delay estimates.
      do i = 1,qnum_sites
*         Now check for atmospheric delay offsets
          if( qparn_atm(i).ne.0 ) then
              na = na + 1
              qapr_codes(na) = 61 + i*256
              qapr_list(na) = site_atm(i)
              np = np + 1
              nq = qparn_atm(i)
              qglb_codes(nq) = 61 + i*256
          end if
      end do
 
****  Now do the satellites
 
      do i = 1, qnum_svs
          do j = 1,max_svs_elem
              if( qparn_svs(j,i).ne.0 ) then
                  na = na + 1
                  qapr_list(na) = svs_pos(j,i)
                  qapr_codes(na) = 51 + j*256 + i*65536
                  np = np + 1
                  nq = qparn_svs(j,i)
                  qglb_codes(nq) = 51 + j*256 + i*65536
              end if
          end do
      end do

****  Now add in polar motion UT1
      do i = 1,3
         do j = 1,2

*           Since the polar position UT1 are explicitly saved we
*           do not need to save apriori value here.
            if( qparn_pmu(j,i).ne.0 ) then
                np = np + 1
                nq = qparn_pmu(j,i)
                if( i.lt.3 ) then
*                   Save polar motion (wob) codes
                    qglb_codes(nq) = 13 + (i + (j-1)*2)*256
                else
*                   Save UT1 codes
                    qglb_codes(nq) = 14 + j*256
                end if
            end if
         end do
      end do
      
****  Now add multiday polar motion UT1 values
      ne = 0
      do i = 1,3
         do j = 1,2

*           Loop over the epochs for this station
            nf = 0
            do k = 1, qnum_mul_pmu(j,i)
                if( qparn_mul_pmu(j,i,k).ne.0 ) then
                   np = np + 1
                   ne = ne + 1
                   nf = nf + 1
                   nq = qparn_mul_pmu(j,i,k)
                   qglb_codes(nq) = 55+i + j*256 + ne*65536
                   qmul_par_ep(ne) = qref_ep(nq)                   
                   na = na + 1
                   qapr_list(na) = qpmu_mul_apr(j,i,nf)
                   qapr_codes(na) = 55+i + j*256 + ne*65536
                   qmul_apr_ep(ne) = qref_ep(nq)
                end if
            end do
         end do
         
      end do

****  Now add translation parameters
      do i = 1, 3
         if ( qparn_tran(i,1).ne.0 ) then
            np = np + 1
            na = na + 1
            nq = qparn_tran(i,1)
            qglb_codes(nq) = 52 + i*256
            qapr_codes(na) = 52 + i*256
         end if
      end do

***   Do translation rate
       do i = 1, 3
         if ( qparn_tran(i,2).ne.0 ) then
            np = np + 1
            na = na + 1
            nq = qparn_tran(i,2)
            qglb_codes(nq) = 53 + i*256
            qapr_codes(na) = 53 + i*256
         end if
      end do

****  Do Scale
      do i = 1,2
         if( qparn_scale(i).ne.0 ) then
             np = np + 1
             na = na + 1
             nq = qparn_scale(i)
             qglb_codes(nq) = 53 + i
             qapr_codes(na) = 53 + i
         end if
      end do
  
     

      
****  Now save the values
      qnum_apr_codes = na
      qnum_par_codes = np
      qent_par_codes = ne
      qent_apr_codes = ne

      
***** Thats all
      return
      end

CTITLE MAKE_JPL

      subroutine make_jpl( unit, line, np, cov_parm, sol_parm )

      implicit none

*     Routine to read and make a binary hfile for JPL files:
*     Lots of entries need to be dummied because information is
*     not in these files.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'


* unit   - Unit for reading the file
* np     - Number of parameters in this solution

      integer*4 unit, np

* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector

      real*8 cov_parm(np,np), sol_parm(np)

* line   - Line read from file

      character*(*) line 

* LOCAL VARIABLES

* date(5)  - Calender date
* num_month - Function to return numeric month given ascii month
* nq       - Generic parameter number
* ns       - Simple form of number of sites
* i,j      - Loop counters
* num_in_file = number of solutions in file (1 only for JPL solns)
* ierr     - IOSTAT error
* trimlen  - Length of string
* cjs      - Current JPL height site number.  Used to see if height
*            missed.

      integer*4 date(5), num_month, nq, ns, i,j, num_in_file, ierr,
     .          trimlen, cjs, indx

* ltog_parn(max_glb_parn) -- Mapping for JPL to global parameter numbers
* ng, nh -- global paramater number
      integer*4 ltog_parn(max_glb_parn), ng, nh


* sectag   - Seconds part of time
* param_sig - Parmeter sigma read from file
* corr      - Parameter corelation read from file

      real*8 sectag, param_sig, corr, sol_in

      character*4 site_code

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters


* mon      - ASCII Month 

      character*3 mon

      data I8 / 1 /

*
****  This is JPL file:  Read the number of paramters and date from
*     the first line.
      num_in_file = 1
      read(line,100) qnum_parn, date(1), mon, date(3)
 100  format(I5,15x,i2,a3,i2)

* MOD TAH 0103005: Initialization added because the JPL stacov
*     do not always include all the elements in the matrix
***** Initialize the covariance matrix and solution vector
      call dwint8(0.d0, cov_parm, 1, (I8*qnum_parn)*qnum_parn)
      call dwint(0.d0, sol_parm, 1, qnum_parn)

      
*     Convert character month to numeric month
      date(2) = num_month(mon)
      date(4) = 12
      date(5) =  0
      sectag  = 0.0
      call ymdhms_to_jd( date, sectag, sepoch )

****  Now preread start of file to get the parameter list and remap
*     velocities to be next to positions (file is rewound and first
*     line read).
      call preread_jpl( unit, ltog_parn)

****  Now get the liast of stations and read the parameter estimates
      ns = 0
      qnum_sites = 0
      do i = 1, qnum_parn
         read(unit,'(a)', iostat=ierr) line
*        Now see what the line is:
         if( index(line,'STA X').gt.0 ) then
             qnum_sites = qnum_sites + 1
             ns = qnum_sites
             read(line,200) np, qsite_names(ns), sol_in, 
     .                          param_sig
             ng = ltog_parn(np)
             cov_parm(ng,ng) = param_sig**2
             sol_parm(ng) = sol_in
 200         format(i5,2x,A4,7x,F30.16,4x,F23.16)
             qparn_sites(1,ns) = ng 
             qscale(ng) = 1.d0
             site_pos(1,ns) = sol_parm(ng)
         else if( index(line,'STA Y').gt.0 ) then
             read(line,200) np, qsite_names(ns), sol_in, 
     .                          param_sig
             ng = ltog_parn(np)
             cov_parm(ng,ng) = param_sig**2
             sol_parm(ng) = sol_in
             qparn_sites(2,ns) = ng 
             qscale(ng) = 1.d0
             site_pos(2,ns) = sol_parm(ng)
         else if( index(line,'STA Z').gt.0 ) then
             read(line,200) np, qsite_names(ns), sol_in, 
     .                          param_sig
             ng = ltog_parn(np)
             cov_parm(ng,ng) = param_sig**2
             sol_parm(ng) = sol_in
             qparn_sites(3,ns) = ng 
             qscale(ng) = 1.d0
             site_pos(3,ns) = sol_parm(ng)
         else if( index(line,'VEL X').gt.0 ) then
             read(line,200) np, site_code, sol_in, 
     .                          param_sig
             indx = 0
             ng = ltog_parn(np)
             call get_cmd(site_code,qsite_names, qnum_sites, ns, indx)
             cov_parm(ng,ng) = param_sig**2
             sol_parm(ng) = sol_in
             qparn_vel(1,ns) = ng 
             qscale(ng) = 1.d-3
             site_vel(1,ns) = sol_parm(ng)
         else if( index(line,'VEL Y').gt.0 ) then
             read(line,200) np, site_code, sol_in, 
     .                          param_sig
             ng = ltog_parn(np)
             indx = 0
             call get_cmd(site_code,qsite_names, qnum_sites, ns, indx)
             cov_parm(ng,ng) = param_sig**2
             sol_parm(ng) = sol_in
             qparn_vel(2,ns) = ng 
             qscale(ng) = 1.d-3
             site_vel(2,ns) = sol_parm(ng)

         else if( index(line,'VEL Z').gt.0 ) then
             read(line,200) np, site_code, sol_in, 
     .                          param_sig
             ng = ltog_parn(np)
             indx = 0
             call get_cmd(site_code,qsite_names, qnum_sites, ns, indx)
             cov_parm(ng,ng) = param_sig**2
             sol_parm(ng) = sol_in
             qparn_vel(3,ns) = ng 
             qscale(ng) = 1.d-3
             site_vel(3,ns) = sol_parm(ng)

         else if( index(line,'UT1-UTC   ').gt.0 ) then
             read(line, 220) np, sol_in, param_sig
 220         format(i5,13x,F30.16,4x,F23.16)
             ng = ltog_parn(np)
             cov_parm(ng,ng) = param_sig**2
             sol_parm(ng) = sol_in
             qparn_pmu(1,3) = ng
             qscale(ng) = 15.d3
             qut1_apr(1) = sol_parm(ng)*15.d3
         else if( index(line,'UT1-UTCRATE').gt.0 ) then
             read(line, 220) np, sol_in, param_sig
             ng = ltog_parn(np)
             cov_parm(ng,ng) = param_sig**2
             sol_parm(ng) = sol_in
             qparn_pmu(2,3) = ng
             qscale(ng) = 15.d3*86400.d0
             qut1_apr(2) = sol_parm(ng)*15.d3*86400.d0
         else if( index(line,'XPOLEMOTION').gt.0 ) then
             read(line, 220) np, sol_in, param_sig
             ng = ltog_parn(np)
             cov_parm(ng,ng) = param_sig**2
             qparn_pmu(1,1) = ng
             qscale(ng) = 1.d3
             qwob_apr(1,1) = sol_parm(ng)*1.d3
         else if( index(line,'XPOLERATE').gt.0 ) then
             read(line, 220) np, sol_in, param_sig
             ng = ltog_parn(np)
             cov_parm(ng,ng) = param_sig**2
             sol_parm(ng) = sol_in
             qparn_pmu(2,1) = ng
             qscale(ng) = 1.d3*86400.d0
             qwob_apr(1,2) = sol_parm(ng)*1.d3*86400.d0
         else if( index(line,'YPOLEMOTION').gt.0 ) then
             read(line, 220) np, sol_in, param_sig
             ng = ltog_parn(np)
             cov_parm(ng,ng) = param_sig**2
             qparn_pmu(1,2) = ng
             sol_parm(ng) = sol_in
             qscale(ng) = 1.d3
             qwob_apr(2,1) = sol_parm(ng)*1.d3
         else if( index(line,'YPOLERATE').gt.0 ) then
             read(line, 220) np, sol_in, param_sig
             ng = ltog_parn(np)
             cov_parm(ng,ng) = param_sig**2
             sol_parm(ng) = sol_in
             qparn_pmu(2,2) = ng
             qscale(ng) = 1.d3*86400.d0
             qwob_apr(2,2) = sol_parm(ng)*1.d3*86400.d0
         else

             write(*,290) line(1:max(1,trimlen(line)))
 290         format('Unkown line in JPL file: ',a)
         end if
      end do

*     Write out some information
      write(*,300) qnum_sites, hfile(1:trimlen(hfile)) 
 300  format(i5,' sites found in ',a)
      do i = 1, qnum_sites
         qsite_names(i)(5:8) = '_GPS'
         qfull_names(i) = qsite_names(i)(1:4) // ' ex JPL STACOV'
         qfull_names(i)(23:32) = '---------P'
         write(*,310) i, qsite_names(i), qfull_names(i)
 310     format(i4,'. ',a8,1x,a)
      end do

*     Now read the correlation maxtrix.
      do i = 2, qnum_parn
         do j = 1,i-1
* MOD TAH 010313: Added iostat error check incase less than number
*           of parameters expected.
            read(unit,*,iostat=ierr) np, nq, corr
            ng = ltog_parn(np)
            nh = ltog_parn(nq)
            cov_parm(ng,nh) = corr*sqrt(cov_parm(ng,ng)*cov_parm(nh,nh)) 
            cov_parm(nh,ng) = cov_parm(ng,nh)
         end do
      end do

****  Now rescale the matrix (meter everywhere except solar radiation
*     which will be ratio to direct effect)
      do i = 1, qnum_parn
          call dwsmy(qscale(i), cov_parm(1,i),1, 
     .                          cov_parm(1,i),1,  qnum_parn)
          call dwsmy(qscale(i), cov_parm(i,1),qnum_parn,
     .                          cov_parm(i,1),qnum_parn, qnum_parn)
 
          sol_parm(i) = sol_parm(i)*qscale(i)
      end do
 
*     Generate the Fake KalObs file name
      sKalobs_file = hfile(1:trimlen(hfile))
 
*     Get Mfile name (claim this is database)
      sdata_base = hfile(trimlen(hfile)-13:trimlen(hfile)-4)

*     Save some of the epochs:
      qstart_epoch = sepoch - 0.5d0
      sepoch_start = qstart_epoch
      qend_epoch   = sepoch + 0.5d0
      sepoch_end   = qend_epoch
      ssvs_epoch   = sepoch 
      cnum_soln_recs = 1
      qsowner(1) = 'JPL'
      sgamit_mod = 0
      sload_mod = 0
      qsepoch_start(1) = qstart_epoch
      qsepoch_end(1)   = qend_epoch
      qskalobs_file(1) = hfile
      qsexpt_title(1)  = 'JPL STACOV file'



****  Now read the stacov antenna heights.  Here we need to assume a
*     type of antenna
      call stacov_heights(unit, cjs) 
     
****  Now generate the apriori and parameter codes
 
      call qgen_codes
 
***** Now write out the global file in GLOBK format.
      hf_ext = 'glj'
 
      call gen_gname( hfile, sepoch, num_in_file, glb_dir,
     .                glb_file, hf_ext, name_format,expt_code)
      write(*,600) glb_file(1:trimlen(glb_file))
 600  Format(' Creating Binary file ',a)
 
***** Write out the header (create file if need be)
      call qw_glb_header
 
*     Write out the names
      call qw_glb_names
 
*     Write out the full names
      call qw_glb_full
 
*     Write out the solution description
      call qw_description
 
*     Write out the codes
      call qw_codes
 
*     Write out the apriori values
      call qw_aprioris
 
*     Write out the solution and covariance matrix
      call qw_soln( cov_parm )
 
***** Thats all for this solution.  Now see if there is another
*     one
      return
      end

CTITLE NUM_MONTH

      integer*4 function num_month ( mon )

*     Routine to return a numeric value for a month

* mon  - Month as characters

      character*(*) mon

* LOCAL
      integer*4 i
      character*3 months(12)

      data months / 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
     .              'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'  /


      num_month = 0

      do i = 1, 12
         if( mon(1:3).eq.months(i) ) num_month = i
      end do
      if( num_month.eq.0 ) then
          write(*,100) mon
 100      format('**ERROR** in num_month: Could not find ',a)
      end if

      return
      end

CTITLE GET_SLR_PARN 

      subroutine get_slr_parn( unit, line, np)

      implicit none

*     This routine will read the SLR(GSFC) solution part
*     of the covariance matrix to find out how many parameters
*     have been estimated.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'


* unit   - Unit for reading the file
* np     - Number of parameters in this solution

      integer*4 unit, np

* line   - Line read from file

      character*(*) line 

* LOCAL VARIABLES

* i,j      - Loop counters
* ierr     - IOSTAT error
* trimlen  - Length of string
* ns       - temp site counter.

      integer*4 date(5), ns, i, ierr,
     .          trimlen, sn, indx

* sectag   - Seconds part of time
* param_val - Value of parameter

       real*8 sectag,  param_val

* temp_name  - Temprary names for stations when rate checked.

       character*8 temp_name

*
****  This is an SLR file. Skip two lines to get to epoch of station
*     coordinates
      read(unit,'(a)', iostat=ierr ) line
      read(unit,'(a)', iostat=ierr ) line

*     Read the next line to get the epoch
      read(unit,'(a)', iostat=ierr ) line
      read(line,'(22x,3i2)') date(1), date(2), date(3)
      date(4) = 0
      date(5) = 0
      sectag = 0.d0

C     read(line,100) qnum_parn, date(1), mon, date(3)
C100  format(I5,15x,i2,a3,i2)
      
      call ymdhms_to_jd( date, sectag, sepoch )

      write(*,100) sdata_base, date
 100  format(/' Found SLR solution ',a,' Reference Epoch ',
     .        i4,'/',i2,'/',i2,1x,i2,':',i2)

*
****  Skip four lines
      do i = 1,4
         read(unit,'(a)') line
      end do

****  Now get the liast of stations and read the parameter estimates
      np = 0
      ns = 0
      do while ( trimlen(line).gt.0 )
         read(unit,'(a)', iostat=ierr) line
         if( trimlen(line).gt.0 ) then
             read(line(26:),*) param_val

*            See type of parameter
             np = np + 1
             if( index(line,'X coord (m)').gt.0 ) then
                 ns = ns + 1
                 qnum_sites = ns
                 qsite_names(ns) = line(6:9) // '_SLR'
                 qparn_sites(1,ns) = np
                 qscale(np) = 1.d0
                 site_pos(1,ns) = param_val
             else if (index(line,'Y coord (m)').gt.0 ) then
                 qparn_sites(2,ns) = np
                 qscale(np) = 1.d0
                 site_pos(2,ns) = param_val
             else if (index(line,'Z coord (m)').gt.0 ) then
                 qparn_sites(3,ns) = np
                 qscale(np) = 1.d0
                 site_pos(3,ns) = param_val

*            Now test for velocity.  Here we look up the site
*            name to see which site this is.
             else if (index(line,'X vel (m/y)').gt.0 ) then
                 temp_name = line(6:9) // '_SLR'
                 indx = 1
                 call get_cmd(temp_name, qsite_names, qnum_sites, 
     .                            sn,indx )
                 qparn_vel(1,sn) = np
                 qscale(np) = 1.d0
                 site_vel(1,sn) = param_val
             else if (index(line,'Y vel (m/y)').gt.0 ) then
                 qparn_vel(2,sn) = np
                 qscale(np) = 1.d0
                 site_vel(2,sn) = param_val
             else if (index(line,'Z vel (m/y)').gt.0 ) then
                 qparn_vel(3,sn) = np
                 qscale(np) = 1.d0
                 site_vel(3,sn) = param_val
             end if
          end if
      end do

      qnum_parn = np
             

***** Thats all
      return
      end

CTITLE MAKE_SLR

      subroutine make_slr( unit, line, np, cov_parm, sol_parm )

      implicit none

*     Routine to read the slr file and get the solution vector
*     and covariance matrix.  The file is assumed to be have
*     been rewound back to the beginning

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'


* unit   - Unit for reading the file
* np     - Number of parameters in this solution

      integer*4 unit, np

* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector

      real*8 cov_parm(np,np), sol_parm(np)

* line   - Line read from file

      character*(*) line 

* LOCAL VARIABLES

* nq       - Generic parameter number
* i,j      - Loop counters
* num_in_file = number of solutions in file (1 only for JPL solns)
* ierr     - IOSTAT error
* trimlen  - Length of string

      integer*4 nq, i,j, num_in_file, ierr,
     .          trimlen

*
****  This is an SLR(GSFC) file: Skip down to the start of the 
*     parameter estimates
      num_in_file = 1
      do while ( index(line,'LABEL').eq.0 )
         read(unit,'(a)' ) line
      end do

*     Now read the solution vector
      do i = 1, qnum_parn
         read(unit,'(a)', iostat=ierr) line
         read(line(26:),* ) sol_parm(i)
      end do

*     Skip two lines
      read(unit,'(a)') line
      read(unit,'(a)') line

*     Now read the covariance matrix
      do i = 1, qnum_parn
         read(unit, 200 ) nq, (cov_parm(j,i),j=1,i)
 200     format(i4,2x,3d24.16,:,/,1000(6x,3d24.16,:,/))
      end do

*     Write out some information
      write(*,300) qnum_sites, hfile(1:trimlen(hfile)) 
 300  format(i5,' sites found in ',a)
      do i = 1, qnum_sites
         qfull_names(i) = qsite_names(i) // ' from SLR(GSFC) file'
         write(*,310) i, qsite_names(i)
 310     format(i4,'. ',a8)
      end do

****  Now fill in the symetric part of the matrix (This code actually
*     copies the diagonal as well but this should be no problem)
      do i = 1, qnum_parn
          call dwmov(cov_parm(1,i),1, cov_parm(i,1), qnum_parn,i)
      end do

****  Now rescale the matrix (meter everywhere except solar radiation
*     which will be ratio to direct effect)
      do i = 1, qnum_parn
          call dwsmy(qscale(i), cov_parm(1,i),1, 
     .                          cov_parm(1,i),1,  qnum_parn)
          call dwsmy(qscale(i), cov_parm(i,1),qnum_parn,
     .                          cov_parm(i,1),qnum_parn, qnum_parn)
 
          sol_parm(i) = sol_parm(i)*qscale(i)
      end do
 
*     Generate the Fake KalObs file name
      sKalobs_file = hfile(1:trimlen(hfile))
 
*     Save some of the epochs:
      qstart_epoch = sepoch - 0.5d0
      qend_epoch   = sepoch + 0.5d0
      ssvs_epoch   = sepoch
     
 
****  Now generate the apriori and parameter codes
 
      call qgen_codes
 
***** Now write out the global file in GLOBK format.
      hf_ext = 'gls'
 
      call gen_gname( hfile, sepoch, num_in_file, glb_dir,
     .                glb_file, hf_ext, name_format,expt_code)
      write(*,600) glb_file(1:trimlen(glb_file))
 600  Format(' Creating Binary file ',a)
 
***** Write out the header (create file if need be)
      call qw_glb_header
 
*     Write out the names
      call qw_glb_names
 
*     Write out the full names
      call qw_glb_full
 
*     Write out the solution description
      call qw_description
 
*     Write out the codes
      call qw_codes
 
*     Write out the apriori values
      call qw_aprioris
 
*     Write out the solution and covariance matrix
      call qw_soln( cov_parm )
 
***** Thats all for this solution.  Now see if there is another
*     one
      return
      end

CTITLE INIT_HFR

      subroutine init_hfr

      implicit none

*     Routine to initialize the reading the hfile in terms of the
*     number of solutions etc.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'


*     Clear the initial values
      qnum_soln_recs = 1
      qnum_comb      = 0
      qnum_apr_cov   = 0
      qnum_apr_svs   = 0
      qnum_apr_diag  = 0

      return
      end

CTITLE STACOV_HEIGHTS

      subroutine stacov_heights(unit, cjs)

      implicit none

*     Routine to read and interpret the heights from a stacov file.

      include '../includes/kalman_param.h'
      include '../includes/const_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'

* PASSED VARIABLES
      integer*4 unit  ! Unit number for the stacov file
      integer*4 cjs   ! Current JPL height site number

* LOCAL VARIABLES
      integer*4 ierr  ! IOSTAT file read
      integer*4 jerr  ! Antenna height IOSTAT
      integer*4 js    ! Site number
      integer*4 j     ! Loop counter
      integer*4 indx  ! Position in string

      integer*4 havail(max_glb_site_wrds)  ! Bit set when height in stacov file
      logical kbit

      real*8    lc_h  ! LC height read from file

      character*80 line ! Line read from file
      character*8  tsite ! test site name


***** Read down until end of file
      ierr = 0 
      cjs = 0
      do j = 1, max_glb_site_wrds
         havail(j) = 0
      end do

      do while ( ierr.eq.0 )
         read(unit,'(a)',iostat=ierr) line
         if( ierr.eq.0 ) then
             tsite = line(2:5) // '_GPS'
*            see if we can file this site
             indx = 1
             call get_cmd(tsite, qsite_names, qnum_sites, js, indx)
             if( js.gt.0 ) then
*               OK, found the site
                read(line(18:),*,iostat=jerr) lc_h
                call report_error('IOSTAT',jerr,'decod',line,0,
     .                            'STACOV_HEIGHTS')
                
                qdata_st(js) = sepoch_start
                qdata_en(js) = sepoch_end
*
*               Now assume things about the antenna
                do j = 1,2
                   qarp_ecc(j,js) = 0.0d0
                   qL1a_ecc(j,js) = 0.0d0
                   qL2a_ecc(j,js) = 0.0d0
                end do
* MOD TAH 050513: Changed output to match CWU satsov
c               qL1a_ecc(3,js) = 0.1100d0
c               qL2a_ecc(3,js) = 0.1280d0
c               qarp_ecc(3,js) = lc_h - 
c     .               (pcf1*qL1a_ecc(3,js)+pcf2*qL2a_ecc(3,js))
                qL1a_ecc(3,js) = 0.d0
                qL2a_ecc(3,js) = 0.d0
                qarp_ecc(3,js) = lc_h
                qelev_cut(js) = 0.d0
                qnum_zen(js) = 0.d0
                qant_mod(js)  = '----'
                qante_ty(js)  = '----'
                qante_sn(js)  = '----'
                qrecv_ty(js)  = '----'
                qrecv_sn(js)  = '----'
                qrecv_fw(js)  = '----'
                write(*,220) js, qsite_names(js), qarp_ecc(3,js)
  220           format('Site ',i4,1x,a8,' Height ',F7.4,' m')
                call sbit(havail,js,1)

             end if
         end if
      end do

*     Report any missing sites at end
      do j = 1, qnum_sites
         if( .not.kbit(havail,j) ) then    
            write(*,230) j, qsite_names(j)
  230       format('Site ',i4,1x,a8,' Apparently no height')
         endif
      enddo 

*     Thats all
      return
      end


       
CTITLE PREREAD_JPL

      subroutine preread_jpl( unit, ltog_parn)

      implicit none

*     Routine to read and make a binary hfile for JPL files:
*     Lots of entries need to be dummied because information is
*     not in these files.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'


* unit   - Unit for reading the file
* ltog_parn(max_glb_parn) -- Mapping from JPL parameter number to
*     global number (used to put postions and velocities back 
*     together).

      integer*4 unit, ltog_parn(max_glb_parn)

* LOCAL
* line   - Line read from file

      character*(160) line 

* LOCAL VARIABLES

* i,j      - Loop counters
* num_in_file = number of solutions in file (1 only for JPL solns)
* ierr     - IOSTAT error
* trimlen  - Length of string
* cjs      - Current JPL height site number.  Used to see if height
*            missed.

      integer*4 nq, ns, i,indx, ierr, j

* velfound -- logical set true when velocity found
      logical velfound


      character*4 site_code

****  This is JPL file:  Read the number of paramters and date from
*     the first line.

****  Now get the liast of stations and read the parameter estimates
      ns = 0
      qnum_sites = 0
      velfound = .false.
      do i = 1, qnum_parn
         read(unit,'(a)', iostat=ierr) line
*        Now see what the line is:
         if( index(line,'STA X').gt.0 ) then
             qnum_sites = qnum_sites + 1
             ns = qnum_sites
         else if( index(line,'VEL X').gt.0 ) then
             velfound = .true.
         endif

      end do
      rewind(unit)
      read(unit,'(a)') line

****  Now see how we should order parameters
      if( velfound ) then
          write(*,220) ns, qnum_parn
 220      format('Velocities found for ',i4,' stations. Reordering ',
     .            i5,' parameters')
          do i = 1, ns
             do j = 1, 3
                ltog_parn(3*(i-1)+j) = 6*(i-1)+j 
                ltog_parn(3*ns+3*(i-1)+j) = 6*(i-1)+3+j 
             enddo
          enddo
          do i = 6*ns,qnum_parn
             ltog_parn(i) = i
          enddo
       else
           do i = 1,qnum_parn
             ltog_parn(i) = i
          enddo
       endif

*****  List values
       write(*,320) (i,ltog_parn(i),i=1,qnum_parn)
 320   format('PARAMETER MAPPING',/,300(5(I5,'->',I5,2x),/))
 
       return
       end

CTITLE DECON_STRAPP
 
      subroutine decon_strapp(str, np, cov_parm, sol_parm)
 
      implicit none

*     Routine to apply the Rotation and translaion (-d=TR) deconstraints
*     to the covariance matix.  (No change to sol_parm although it is 
*     passed.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES

*   str     - Deconstring passed from decon_str or gamit_str
      character*(*) str

*   np      - number of parameters expected
 
      integer*4 np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)

* LOCAL VARIABLES
* i,j, k, l, m    - Loop counters
* p1, p2  - Parameter numbers  

      integer*4 i,j,k,l,m 
      integer*4 p1, p2

* rot_var -- Variance of rotational noise added 
* tran_var -- Translational variatiance

      real*8 rot_var, tran_var

* pmu_rbrot(3,3,max_glb_sites) -- Partial for rotation.
      real*8 pmu_rbrot(3,3,max_glb_sites)

****  See if will apply Rotation de-constraint.
      if( index(str,'R').gt.0   .and.
     .    snx_constraint.eq.0               ) then 
         rot_var = 100.0d0**2
         write(*,'(a,F6.2,a)') ' Applying ',sqrt(rot_var),
     .                     ' mas Rotational Deconstraint'
         do i = 1, qnum_sites
*                     ! XYZ
             do j = 1,3
* MOD TAH 131231: Changed this to be a rigid body rotation and not
*                a coordinate system rotation.  Rigid body allows 
*                network to rotate correctly
                 call pmu_rigidbody( pmu_rbrot(1,j,i), site_pos(1,i),
     .               j, qut1_apr(1))
             end do
         end do

         do i = 1, qnum_sites
            do j = 1,3
               do k = 1, qnum_sites
                  do l = 1,3
                     p1 = qparn_sites(j,i)
                     p2 = qparn_sites(l,k)
                     do m = 1,3
                        cov_parm(p1,p2) = cov_parm(p1,p2) +
     .                       pmu_rbrot(m,j,i)*pmu_rbrot(m,l,k)*
     .                            (rot_var)
                     end do
                  end do
               end do
            end do
         end do

****     Now apply the covariance to the rotation portion of the
*        matrix
         do j = 1, 3
            do k = 1, qnum_mul_pmu(1,j)
               p1 = qparn_mul_pmu(1,j,k)
               do l = 1, qnum_mul_pmu(1,j)
                  p2 = qparn_mul_pmu(1,j,l)
                  cov_parm(p1,p2) = cov_parm(p1,p2) + rot_var
               end do
            end do
         end do

****     Now apply the cross terms between station position and 
*        EOP.  Only do this from X and Y pole because UT1 should
*        no be constrainable
         do j = 1, 3                 ! Loop over XY pole and UT1
            do k = 1, qnum_mul_pmu(1,j)
               p1 = qparn_mul_pmu(1,j,k)
               do i = 1, qnum_sites  ! Loop over Sites
                  do l = 1,3         ! Loop over XYZ at site
                     p2 = qparn_sites(l,i)
* MOD TAH 171005: Set negative correlation between site positions and
*                 rotation parameters. 
* MOD TAH 171005: Keep sign the same (related to rigid body rotation above?)
                     cov_parm(p1,p2) = cov_parm(p1,p2) +
     .                       pmu_rbrot(j,l,i)*(rot_var) 
                     cov_parm(p2,p1) = cov_parm(p1,p2)
                  end do
               end do
            end do
         end do
* MOD TAH 180205: Set bit 14 of sgamit_mod to show applied
         call sbit(sgamit_mod,14,1)

      end if


* MOD TAH 010308: Add translation noise to the GFZ files
      if( index(hfile,'gfz').gt.0 .or.
     .    index(str,'T').gt.0 ) then
         tran_var = 1.0d0
         write(*,'(a,F6.2,a)') ' Applying ',sqrt(tran_var),
     .                     ' m   Translation Deconstraint'
         do i = 1, qnum_sites
            do k = 1, qnum_sites
               do j = 1,3
                  p1 = qparn_sites(j,i)
                  p2 = qparn_sites(j,k)
                  cov_parm(p1,p2) = cov_parm(p1,p2) + tran_var
               end do
            end do
         end do
* MOD TAH 171005: See if translation parameters are estimated.  If
*        so, apply the negative covariance to the station positions
         do j = 1,3 
            p1 = qparn_tran(j,1)
            if( p1.gt.0 ) then 
                cov_parm(p1,p1) = tran_var
                do k = 1, qnum_sites
                   p2 = qparn_sites(j,k)
                   cov_parm(p1,p2) = cov_parm(p1,p2) - tran_var
                enddo
            endif
         enddo
* MOD TAH 180205: Set bit 15 of sgamit_mod to show applied
         call sbit(sgamit_mod,15,1)
      end if 

****  Thats all
      return
      end







