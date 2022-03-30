 
      subroutine get_svs_inf( unit, num_obs_svs)
 
      implicit none

*     This routine will read the satellite information after
*     the "Satellite used:" line is found.  The reason for
*     a subroutine is handle the complications of multiple
*     satellite lines in solve ver 9.28 when more than 24
*     satellites were used.

 
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
*   num_obs_svs(max_glb_svs) - Number of data for each satellite
  
      integer*4 unit, num_obs_svs(max_glb_svs)
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   indx    - Pointer in string
*   trimlen - Length of string
*   i   - Loop counter
*   nsvs_first  - Number of satellites on the
*       - first line
 
      integer*4 ierr, jerr, indx, trimlen, i, nsvs_first
 
*   split_svs_line  - Set true if the satellites lines
*           - have been split into two (solve
*           - version 9.28 with >24 satellites)
 
      logical split_svs_line
 
*   prn_lead        - Lead name for satellites (PRN typically)
*   prn_num     - PRN number as ascii
 
 
      character*4 prn_lead, prn_num, cdum
      character*4 sv_antMod   ! Model code from GAMIT
 
*   line        - Line read from file
 
 
 
      character*256 line
 
****  Now get the satellite information.
      read(unit,'(a)', iostat=ierr) line
*     This is line with PRN numbers
      read(unit,'(a)', iostat=ierr) line
 
* MOD TAH 941005: Check more carefully what the next line
*     contains.  Solve ver 9.28 could write another Satellite
*     used line if there were more than 24 satellites
      if( line(2:10).eq.'Satellite' ) then
          split_svs_line = .true.
*
*         Skip next line with channel numbers and then
*         Read the next line to get the first 'n' satellites
          read(unit,'(a)', iostat=ierr) line
          read(unit,'(a)', iostat=ierr) line
          indx = 1
          call getword(line, prn_lead, indx)
 
          i = 0
          do while ( indx.lt.trimlen(line))
             i = i + 1
             call getword(line, prn_num, indx)
             if( prn_num(2:2).eq.' ' ) then
                 prn_num(2:2) = prn_num(1:1)
                 prn_num(1:1) = '0'
             end if
 
             qsvs_names(i) = prn_lead(1:trimlen(prn_lead)) //
     .              '_' // prn_num(1:max(1,trimlen(prn_num)))

          end do
          nsvs_first = i

*         Read the next line so that we have the second line
*         of PRN numbers.
          read(unit,'(a)', iostat=ierr) line
 
      else
          split_svs_line = .false.
          nsvs_first = 0
      end if
 
*     Now pull off the PRN numbers
      indx = 1
      call getword(line, prn_lead, indx)
 
* MOD TAH 941005: run counter skipping the satellites on the first
*     line (if there is only one line then nsvs_first = 0
      do i = nsvs_first+1 , qnum_svs
          call GetWord(line, prn_num, indx )
 
*         Construct name
          if( prn_num(2:2).eq.' ' ) then
              prn_num(2:2) = prn_num(1:1)
              prn_num(1:1) = '0'
          end if
 
          qsvs_names(i) = prn_lead(1:trimlen(prn_lead)) //
     .            '_' // prn_num(1:max(1,trimlen(prn_num)))
      end do
 
***** Then use next line to get number of observations
 
*     If there are multiple lines then read first line
      if( split_svs_line ) then
          read(unit,'(a)', iostat=ierr) line
          indx = 1
          call multiread(line, indx, 'I4', ierr, num_obs_svs, cdum,
     .                   nsvs_first)
      end if
 
*     Get the number of observations on each satellite.
      read(unit,'(a)', iostat=ierr) line
* MOD TAH 050214: Start in column 5 since 10.17 solve put Obs at front of line
      indx = 5
      call multiread(line, indx, 'I4', ierr, num_obs_svs(nsvs_first+1),
     .               cdum, qnum_svs - nsvs_first)

***** For Post 10.17 solve there are more lines
      if( solve_ver .ge. 1017 ) then
          read(unit,'(a)', iostat=ierr) line
          indx = 8
          call multiread(line, indx, 'I4', ierr, qsvi_block,
     .               cdum, qnum_svs)
*         Read the line woth AntMod
          read(unit,'(a)', iostat=ierr) line
*         Loop over the line and save model information
          indx = 8
          do i = 1, qnum_svs
             call GetWord(line, sv_antMod, indx)
*            Convert the code  (ELI5 into a model name)
             call SMod_to_full( sv_antMod, qsvi_antmod(i), 
     .                          qsvi_ocode(i))

*            GAMIT Values are always LC so save (and only one)
* MOD TAH 200131: New version now L1/L2 so make '12' instead of '50'
             qsvi_ocode(i)(1:2) = '12'

*            Extract the PRN number from the name
             read(qsvs_names(i)(5:),*,iostat=jerr)  qsvi_prn(i)
          end do

*         Read and deocode DX, DY, DZ offsets
*  DY  0.0000-0.0012 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000-0.0012-0.0012-0.0012 0.0000-0.0012-0.0012-0.0012-0.0012-0.0012-0.0012-0.0012-0.0012 0.0000 0.0000 0.0000 0.0000-0.0012 0.0000 0.0000
          read(unit,'(a)', iostat=ierr) line
          indx = 6
          read(line,210,iostat=ierr) (qsvi_antpos(1,1,i),i=1,qnum_svs)
          call report_error('IOSTAT',ierr,'decod',line,1,'GET_SVS_INF')
 210      format(5x,32F7.4)
          read(unit,'(a)', iostat=ierr) line
          read(line,210,iostat=ierr) (qsvi_antpos(2,1,i),i=1,qnum_svs)
          read(unit,'(a)', iostat=ierr) line
          read(line,210,iostat=ierr) (qsvi_antpos(3,1,i),i=1,qnum_svs)

          if( index(line,'DZ').eq.0 ) then
             write(*,220) ierr,line
 220         format('DZ line not where expected: Ierr ',i4,1x,a30)
          endif
****      Last peice of information that we need is the SV number.  We need 
*         read svnav.dat to get this information.  We can't do this yet because
*         we need the epoch of these data and that has not yer been read
      endif

 
***** Thats all
      return
      end

CTITLE GET_SVS_NEW
 
      subroutine get_svs_new( unit, num_obs_svs)
 
      implicit none

*     This routine will read the satellite information after
*     the "Satellite used:" line is found.  This version reads
*     the 2.00 and greater version
 
* MOD TAH 150501: Block modified to handle 3.0 version hfiles (gamit_hf_ver>=300)
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - INput hfile unit
*   num_obs_svs(max_glb_svs) - Number of data for each satellite
  
      integer*4 unit, num_obs_svs(max_glb_svs)
 
* LOCAL VARIABLES
 
*   ierr    - IOSTAT error
*   trimlen - Length of string
*   i, j   - Loop counter
 
      integer*4 ierr, trimlen, i, j
 
*   prn_lead        - Lead name for satellites (PRN typically)
 
 
      character*1 offcon    ! Function to return constellation letter based on
                            ! offset PRN (original PRN for mod(prn,100))
      character*1 ssys       ! GNSS system 

      character*22 antbody   ! Antenna/Satellite Body type
 
*   line        - Line read from file
      character*256 line

****  Read the next line (Chan  PRN  Blk ...)
      read(unit,'(a)',iostat=ierr) line
* MOD TAH 180401: Removed getting LEAD since all GNSS use PRN in 
*     these columns
!     prn_lead = line(8:10)   ! Normally this is PRN

****  Now loop over each satellites
      do i = 1, qnum_svs
         read(unit,'(a)') line
         if( gamit_hf_ver.ge.310 ) then   ! Newest version (200127)
             read(line,120,iostat=ierr) qsvi_prn(i), antbody, 
     .           num_obs_svs(i), qsvi_antMod(i), 
     .          (qsvi_antpos(j,1,i),j=1,3), 
     .          (qsvi_antpos(j,2,i),j=1,3)
 120         format(5x,I4,2x,a22,I7,1x,A15,2(1x,3f8.4)) 
             call name_to_blk(+1, antbody, qsvi_block(i))
             qsvi_prn(i) = qsvi_prn(i) + int(qsvi_block(i)/10)*100
         elseif( gamit_hf_ver.eq.300 ) then   ! Previous
             read(line,130,iostat=ierr) qsvi_prn(i), antbody, 
     .           num_obs_svs(i), qsvi_antMod(i), 
     .          (qsvi_antpos(j,1,i),j=1,3)
 130         format(5x,I4,2x,a22,I7,1x,A15,1x,3F7.4) 
!  1     1  BLOCK IIF               1955  IGS08_1840 ELEV  0.3940 0.0000 1.5613
!  2     2  GALILEO-2               2578  IGS08_1930 ELEV  0.0000 0.0000 0.0000
!  3     3  GLONASS-M               2702  IGS08_1930 ELEV  0.0000 0.0000 0.0000
             call name_to_blk(+1, antbody, qsvi_block(i))
             qsvi_prn(i) = qsvi_prn(i) + int(qsvi_block(i)/10)*100

         else    !   (200-series files) 
            read(line, 220, iostat=ierr) qsvi_prn(i), qsvi_block(i),
     .           num_obs_svs(i), qsvi_antMod(i), 
     .          (qsvi_antpos(j,1,i),j=1,3)
 220        format(5x,I4,2x,I4,2x,I6,2x,a15,1x,3F7.4)
         endif

*        GAMIT Values are always LC so save (and only one)
* MOD TAH 200131: Make L12 '12' code rather than '50'
         qsvi_ocode(i)(1:2) = '12'
         qsvi_ocode(i)(3:3) = 'A'
         qsvi_ocode(i)(4:4) = qsvi_antMod(i)(12:12)
      end do

* MOD TAH 180401:Updated for GNSS with the addition of SVN numbers
***** Now get the SVS numbers from svnav.dat
* MOD TAH 180607: Only call if we know date (will be call latter as well)
      if( sepoch.gt.0.d0 ) call Get_SVnum( sepoch )

* MOD TAH 180401: Generate satellites names with SVN numbers
      do i = 1, qnum_svs
         ssys = offcon(qsvi_prn(i))
         write(qsvs_names(i),320) ssys, qsvi_svn(i),
     .                         mod(qsvi_prn(i),100)
 320     format(a1,I3.3,'_',i2.2)
      end do

***** Thats all
      return
      end


