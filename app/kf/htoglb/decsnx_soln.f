CTITLE DECSNX_SOLN
 
      subroutine decsnx_soln(unit, line, np, cov_parm, sol_parm)
 
      implicit none

*     Routine to decode the solution  blocks from SINEX:
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector

 
      real*8 cov_parm(np,np), sol_parm(np)
 
*   line        - Line read from file file
 
      character*(*) line

* MOD TAH 070704: Moved to common and intialized in make_snx so that
*  multiple sinex files can be processed at at same time.
C      logical estread, aprread, epread

* estread, aprread -- Set true once the matrix est and apriori
*     have been read.  We then can remove the constraints.

C      data  estread / .false. /, aprread / .false. /, epread / .false. /
      
* LOCAL VARIABLES
 
*   indx        - Position in string
*   ierr        - IOSTAT error
*   trimlen     - Length of string
*   ipivot(max_glb_parn)  - Pivot elements for matrix inverse
 
      integer*4 indx, trimlen, i, j, ipivot(max_glb_parn)
     
* scale(max_glb_parn) -- Scale for matrix inversion

      real*8 scale(max_glb_parn)

*       block_found - True if we know how to decode
*           - the block
 
 
      logical block_found

 
****  Start looking for blocks
      block_found = .false.
 
      indx = index(line,'SOLUTION/EPOCHS')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_soln_epochs(unit,line)
      end if
      
      indx = index(line,'SOLUTION/STATISTICS')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_soln_stats(unit,line)
      end if
 
      indx = index(line,'SOLUTION/ESTIMATE')
      if( indx.gt.0 ) then 
          block_found = .true.
          if ( .not. epread ) then
             call dec_soln_est(unit,line, np, cov_parm, sol_parm)
*            Only copy if these are not normal equations
             if( .not. neq_sys ) then
                do i = 1, np
                   sol_parm(i) = qsol(i)
                end do
             end if
             epread = .true.
          else
             call snx_finbl(unit)
          endif
  
      end if
 
      indx = index(line,'SOLUTION/NORMAL_EQUATION_VECTOR')
      if( indx.gt.0 ) then 
          block_found = .true.
          neq_sys = .true.
          vec_read = .true.
          call dec_soln_vec(unit,line, np, cov_parm, sol_parm)
          do i = 1, np
             sol_parm(i) = qsol(i)
          end do  
      end if

      indx = index(line,'SOLUTION/APRIORI')
      if( indx.gt.0 ) then 
          block_found = .true.
          if( epread ) then
             call dec_soln_apr(unit,line, np, cov_parm, sol_parm)
* MOD TAH 131231: Save the apriori values values qwob_apr and qut1_apr
*            if there was only 1 multiday PMU entry (1-day solution).
             do i = 1,2       ! Loop X/Y wob
                do j = 1,2    ! Value and rate
                   if( qnum_mul_pmu(j,i).eq.1 ) then
                       qwob_apr(i,j) = qpmu_mul_apr(j,i,1)
                   end if
                enddo
             end do
*            Now do UT1 and flip LOD
             if( qnum_mul_pmu(1,3).eq.1 ) 
     .           qut1_apr(1) =  qpmu_mul_apr(1,3,1)
*                          Flip LOD to dUT1/dt
             if( qnum_mul_pmu(2,3).eq.1 ) 
     .           qut1_apr(2) =  qpmu_mul_apr(2,3,1)

             apr_missed = .false.
          else
              call snx_finbl(unit)
          endif
          
      end if
 
      indx = index(line,'SOLUTION/MATRIX_ESTIMATE')
      if( indx.gt.0 ) then
          block_found = .true.
          estread = .true.
          cov_parm(1,1) = 1
          call dec_soln_mat(unit,line, np, cov_parm, sol_parm)
      end if
 
      indx = index(line,'SOLUTION/NORMAL_EQUATION_MATRIX')
      if( indx.gt.0 ) then
          block_found = .true.
          neq_sys = .true.
          neq_read = .true.
          cov_parm(1,1) = 1
          call dec_soln_neq(unit,line, np, cov_parm, sol_parm)
      end if
 
      indx = index(line,'SOLUTION/MATRIX_APRIORI')
      if( indx.gt.0 ) then
          block_found = .true.
          if( epread ) then
             aprread = .true.
             call dec_soln_con(unit,line, np, cov_parm, sol_parm)
          else
             call snx_finbl(unit)
          endif
      end if

      if ( .not.block_found ) then
          write(*,500) line(1:trimlen(line))
500       format('DECSNX_SOLN: Unknown block type',/,a)
      end if
     
****  Thats all
      return
      end
 
 
CTITLE DEC_SOLN_epochs
 
      subroutine dec_soln_epochs(unit,line )
 
      implicit none

*     Routine to decode the solution epochs to get the start and
*     stop times for each station
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
 
      integer*4 unit
 
*   line        - Line read from file file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   sn      - Site number from list of sites
*   jp      - Parameter number read from the input list
*   syr, sdoy, ssec - start yr, doy and ses
*   eyr, edoy, esec - stop yr, doy and see
*   indx        - Position in string
*   occ         - Occupation number
 
      integer*4 ierr, sn, syr, sdoy, ssec, indx, occ,
     .          eyr, edoy, esec, jerr
 
*   end_block   - Indicates end of block found
 
      logical end_block
 
*   code    - Site code read from Sinex file
*   pt      - Point characters
*   full_name  -  Full site name with occ and pt added
 
      character*8 code, full_name
      character*4 pt, chocc
 
****  Start reading the block
      end_block = .false.
      ierr = 0
 
      do while ( ierr.eq.0 .and. .not.end_block)
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then
  
              read(line,210,iostat=jerr) code,pt,chocc, syr, sdoy, ssec,
     .                       eyr, edoy, esec
 210          format(1x,a4,1x,a2,1x,a4,3x,i2,1x,i3,1x,i5,
     .               1x,i2,1x,i3,1x,i5)

*             Check the character occupancy
              if( chocc.eq.'----' ) chocc = '   1'
              read(chocc,*,iostat=jerr) occ

* MOD TAH 970716: Check for blank pt names, replace with ' A'
              if( pt(2:2).eq.' ' .or. pt(2:2).eq.'-' ) pt = ' A'

****          OK, construct the full name for the site and
*             see if we find it.
              write(full_name, 240) code, occ, pt(2:2)
              call casefold(full_name)
 240          format(a4,i3.3,a1)
              indx = 1
              call get_cmd(full_name, qsite_names, qnum_sites,
     .                        sn, indx )

              if( sn.gt.0 ) then  
                  call yds_to_jd( syr, sdoy, ssec, qdata_st(sn))
                  call yds_to_jd( eyr, edoy, esec, qdata_en(sn))
              end if
          end if
          if( line(1:1).eq.'-' ) then
              end_block = .true.
          end if
      end do
 
****  Thats all
      return
      end

CTITLE DEC_SOLN_stats
 
      subroutine dec_soln_stats(unit,line )
 
      implicit none

*     Routine to decode the solution stats block to get number of
*     data and chiqsq 
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
 
      integer*4 unit
 
*   line        - Line read from file file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   sn      - Site number from list of sites
*   jp      - Parameter number read from the input list
*   syr, sdoy, ssec - start yr, doy and ses
*   eyr, edoy, esec - stop yr, doy and see
*   indx        - Position in string
*   occ         - Occupation number
 
      integer*4 ierr, jerr

      real*8 rnum_obs
 
*   end_block   - Indicates end of block found
 
      logical end_block
 
****  Start reading the block
      end_block = .false.
      ierr = 0
 
      do while ( ierr.eq.0 .and. .not.end_block)
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then

*             See what entries we have:
              if( index(line,'NUMBER OF OBSERVATION').gt.0 ) then
                  read(line(32:), * , iostat=jerr) rnum_obs
                  snum_obs = rnum_obs
cc 210              format(32x,i22)
                  cnum_obs = snum_obs
              else if( index(line,'VARIANCE FACTOR').gt.0 ) then
                  read(line, * , iostat=jerr) cchisq
cc 220               format(32x,f22.15)
              endif
              
          end if
          if( line(1:1).eq.'-' ) then
              end_block = .true.
              write(*,240) cnum_obs, cchisq
 240          format(' Analysis with ',i12,' Obs, Chi**2 ',f22.6)
          end if
      end do
 
****  Thats all
      return
      end

CTITLE DEC_SOLN_EST
 
      subroutine dec_soln_est(unit,line, np, cov_parm, sol_parm)
 
      implicit none

*     Routine to decode the solution estimates part of the sinex
*     files.  These values are written directly into the sol_parm
*     vector.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
*   line        - Line read from file file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   sn      - Site number from list of sites
*   jp      - Parameter number read from the input list
*   cyr, cdoy, csec - Reference yr, doy and sec
*   indx        - Position in string
*   nr          - "real" parameter number in final covariance matrix
*   occ         - Occupation number
*   cons_type   - Constraint type on station (0,1 or 2)
 
      integer*4 ierr, sn, jp, cyr, cdoy, csec, indx, nr, occ, i,
     .          cons_type, jerr, ne, kerr, svn
 
*   est     - Estimate from file
*   est_sig - Sigma of estimate.  Needed when correlations are given
*   ref_jd  - JD of reference epoch for value
 
      real*8 est, est_sig, ref_jd
      
* dut1, dlod -- Short period corrections to UT1 amd LOD

      real*8 dut1, dlod 
      
*   end_block   - Indicates end of block found
 
      logical end_block
      logical elim   ! Set true if normal equations entry
                     ! to be removed due to being on wrong day (COD)
 
*   code    - Site code read from Sinex file
*   stype   - Parameter type read from sinex file
*   estu        - Units for the estimate
*   pt      - Point characters
*   full_name  -  Full site name with occ and pt added
*   chocc   - Character version of occupancy (version 1.0)
 
      character*8 code, stype, estu, full_name
      character*4 pt, chocc

* MOD TAH 210505: Added support for multi-frequency satellite
*     antenna offsets.
      integer*4 svant_off  ! Functon to return index in qparn_svs
                  ! array for antenna offset (L1/LC are in 
                  ! max_svs_elem - (2:0); Other are put in IC and
                  ! radiation parameter slots.
      integer*4 ao   ! Offset value returned from svant_off based
                  ! on pt string. 

*   found  - Set true if parameter correctly idenified
      logical found
 
****  Start reading the block
      end_block = .false.
      ierr = 0
      nr = 0
      do i = 1, max_glb_parn
        itoo(i) = 0
        qsig(i) = 0.d0
      end do
      
      do while ( ierr.eq.0 .and. .not.end_block)
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ') then
              
              if( cglb_vers.eq.5 ) then
                  read(line,210,iostat=jerr) 
     .                 jp, stype, code,pt,occ, cyr, cdoy, csec,
     .                 estu, cons_type, est, est_sig
 210              format(i6,1x,a4,1x,a4,1x,a2,1x,i4,1x,i2,1x,i3,1x,i5,
     .                  1x,a4,i2,d21.14,1x,e14.8)
              else
                  read(line,215,iostat=jerr) 
     .                 jp, stype, code,pt,chocc, cyr, cdoy, csec,
     .                 estu, cons_type, est, est_sig
 215              format(i6,1x,a6,1x,a4,1x,a2,1x,a4,1x,i2,1x,i3,1x,i5,
     .                  1x,a4,i2,d22.14,1x,e11.5)
     
*                 Check the occupancy value
                  if( chocc.eq.'----' ) chocc = '   1'
                  read(chocc,*,iostat=kerr) occ
                  if( kerr.ne.0 ) occ = 1                  
              end if
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                          'dec_soln_apr')
 
              call yds_to_jd( cyr, cdoy, csec, ref_jd)

* MOD TAH 970716: Check for blank pt names, replace with ' A'
              if( pt(2:2).eq.' ' .or. pt(2:2).eq.'-' ) pt = ' A'

****          OK, construct the full name for the site and
*             see if we find it.

              write(full_name, 240) code, occ, pt(2:2)
 240          format(a4,i3.3,a1)
              call casefold(full_name)
              indx = 1
              call get_cmd(full_name, qsite_names, qnum_sites,
     .                         sn, indx )

              found = .false. 
              if( stype.eq.'STAX' .and. sn.le.0 ) then
                  write(*,250) full_name(1:4), full_name(5:7),
     .                         full_name(8:8)
 250              format('**WARNING** No site/id entry for code ',
     .                   a4, ' Occ ',a3,' Point ',a1)
              end if

****          See if satellite axis offset
              if( stype(1:4).eq.'SATA' ) then
                  read(code,'(1x,i3)',iostat=kerr) svn
                  sn = 0
                  do i = 1, qnum_svs
                     if( qsvi_svn(i).eq.svn ) sn = i
                  end do
                  if( kerr.ne.0 .or.sn.eq.0 ) then
                     write(*,260) line(1:80)
 260                 format('***WARNING*** Unable to find SV in: ',a)
                  endif
* MOD TAH 210505: Extract out frequency for antenna offsets to handle
*                 TUG L1, L2, L5, L7,L8 values.  Save as offset to 
*                 max_svs_elem (see SINEX has not IC or Rad parameters
*                 use these slots for non-L1, LC values.
                  ao = svant_off(pt)
             endif

* MOD TAH 200727: See if we will mark this parameter for elimination
*             due to being on wrong day
              elim = .false.
              if( qowner(1:3).eq.'COD' .and. sngday .and.
     .            abs(ref_jd-qmid_epoch).gt.0.1 ) elim=.true.                                      
              if( stype.eq.'STAX' .and. sn.gt.0 ) then
                  nr = nr + 1
                  qparn_sites(1,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  rm_param(nr) = elim
* MOD TAH 200728: If not eleiminated, mark site as used.  Only
*                 check of X comp (YZ should be present anyway)/
                  if( .not.elim ) rm_site(sn) = .false.
                  site_pos(1,sn) = qsol(nr)
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'STAY' .and. sn.gt.0) then
                  nr = nr + 1
                  qparn_sites(2,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  rm_param(nr) = elim
                  site_pos(2,sn) = qsol(nr)
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'STAZ' .and. sn.gt.0) then
                  nr = nr + 1
                  qparn_sites(3,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  rm_param(nr) = elim
                  site_pos(3,sn) = qsol(nr)
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'VELX' .and. sn.gt.0) then
                  nr = nr + 1
                  qparn_vel(1,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  site_vel(1,sn) = qsol(nr)
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'VELY' .and. sn.gt.0) then
                  nr = nr + 1
                  qparn_vel(2,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  site_vel(2,sn) = qsol(nr)
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'VELZ' .and. sn.gt.0) then
                  nr = nr + 1
                  qparn_vel(3,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  site_vel(3,sn) = qsol(nr)
                  itoo(jp) = nr
                  found = .true.

* MOD TAH 050710: Addin satellite axis offsets
              else if( stype.eq.'SATA_X' .and. sn.gt.0
     .                 .and. est_sig.gt. 0 ) then
                  nr = nr + 1
C                 qparn_svs(max_svs_elem-2,sn) = nr
C MOD TAH 210505: Changed access for multi-frequency antenna
C                 offset (TUG SINEX)
                  qparn_svs(ao+1,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'SATA_Y' .and. sn.gt.0
     .                 .and. est_sig.gt. 0 ) then
                  nr = nr + 1
C MOD TAH 210505: Changed access for multi-frequency antenna
C                 offset (TUG SINEX)
                  qparn_svs(ao+2,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'SATA_Z' .and. sn.gt.0
     .                 .and. est_sig.gt. 0 ) then
                  nr = nr + 1
C MOD TAH 210505: Changed access for multi-frequency antenna
C                 offset (TUG SINEX)
                  qparn_svs(ao+3,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.

*****         Check for Earth orientation parameters
              else if( stype.eq.'LOD ' .or. stype.eq.'LODR' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.  Convert LOD to UT1 Rate
*                 (mas per day is the unit)
                  found = .true.
                  nr = nr + 1
                  qnum_mul_pmu(2,3) = qnum_mul_pmu(2,3) + 1
                  ne = qnum_mul_pmu(2,3)
                  qscale(nr) = -15.d0
*                 Special mod for JPL: MOD REMOVED TAH 030214.
c                  if( est.lt.0 ) then
c                      est = -est
c                      qscale(nr) = -15.d0
c                  end if
                  
*                 See if we need to restore short-period UT1 corrections.
*                 Need to CHANGE JPL Date when problem fixed.
                  if ( (stype.eq.'LODR' .and. qowner(1:3).ne.'NOA')
     .		      .or.  (qowner(1:3).eq. 'JPL' .and. 
     .                 ref_jd.lt.  2451544.50d0) )    then
                       call short_period_lod( ref_jd, dut1, dlod)
*                      Add correction; convert sec to ms.                       
                       est = est + dlod*1.d3
                  end if
                  
                  qparn_mul_pmu(2,3,ne) = nr
                  qsol(nr) = est
                  qpmu_mul_apr(2,3,ne) = -est*15.d0
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  rm_param(nr) = elim
                  itoo(jp) = nr

              else if( stype.eq.'UT  ' .or. stype.eq.'UT1 ' .or.
     .                 stype.eq.'UTR ' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  found = .true.
                 
                  nr = nr + 1
                  qnum_mul_pmu(1,3) = qnum_mul_pmu(1,3) + 1
                  ne = qnum_mul_pmu(1,3)
                  qscale(nr) = 15.d0

* MOD TAH 991207: Check GFZ values.  Sign seems to be wrong
* MOD TAH 140106: Removed magnitude test
C                 if( qowner(1:3).eq. 'GFZ' .and. est.gt. 1000.d0 ) then
                  if( qowner(1:3).eq. 'GFZ' ) then
                      write(*,*) 'Changing SIGN of GFZ UT values'
                      est = -est
                  end if
                  
*                 See if we need to restore short-period UT1 corrections.
*                 Need to CHANGE JPL Date when problem fixed.
                  if ( stype.eq.'UTR ' .or. 
     .                (qowner(1:3).eq. 'JPL' .and. 
     .                 ref_jd.lt.  2451544.50d0) )    then
                       call short_period_lod( ref_jd, dut1, dlod)
*                      Add correction; convert sec to ms.                       
                       est = est + dut1*1.d3
                  end if
                  
                  qparn_mul_pmu(1,3,ne) = nr
                  qsol(nr) = est
                  qpmu_mul_apr(1,3,ne) = est*15.d0
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  rm_param(nr) = elim
                  itoo(jp) = nr
              else if( stype.eq.'XPO ' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  found = .true.
                  nr = nr + 1
                  qnum_mul_pmu(1,1) = qnum_mul_pmu(1,1) + 1
                  ne = qnum_mul_pmu(1,1)
                  qscale(nr) = 1.d0
                  qparn_mul_pmu(1,1,ne) = nr
                  qsol(nr) = est
                  qpmu_mul_apr(1,1,ne) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  rm_param(nr) = elim
                  itoo(jp) = nr
              else if( stype.eq.'XPOR' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.  Sinex Rate is
                  found = .true.
                  nr = nr + 1
                  qnum_mul_pmu(2,1) = qnum_mul_pmu(2,1) + 1
                  ne = qnum_mul_pmu(2,1)
                  qscale(nr) = 1.d0     
                  qpmu_mul_apr(2,1,ne) = est
                  if( estu.eq.'ma/s' ) then
                       qscale(nr) = 86400.d0
                       qpmu_mul_apr(2,1,ne) = est*86400.d0
                  end if
                  qparn_mul_pmu(2,1,ne) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  rm_param(nr) = elim
                  itoo(jp) = nr
              else if( stype.eq.'YPO ' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  found = .true.
                  nr = nr + 1
                  qnum_mul_pmu(1,2) = qnum_mul_pmu(1,2) + 1
                  ne = qnum_mul_pmu(1,2)
                  qscale(nr) = 1.d0
                  qparn_mul_pmu(1,2,ne) = nr
                  qsol(nr) = est
                  qpmu_mul_apr(1,2,ne) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  rm_param(nr) = elim
                  itoo(jp) = nr

              else if( stype.eq.'YPOR' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  found = .true.
                  nr = nr + 1
                  qnum_mul_pmu(2,2) = qnum_mul_pmu(2,2) + 1
                  ne = qnum_mul_pmu(2,2)
                  qscale(nr) = 1.d0
                  qpmu_mul_apr(2,2,ne) = est
                  if( estu.eq.'ma/s' ) then
                       qscale(nr) = 86400.d0
                       qpmu_mul_apr(2,2,ne) = est*86400.d0
                  end if
                  qparn_mul_pmu(2,2,ne) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  rm_param(nr) = elim
                  itoo(jp) = nr
              else if((stype.eq.'TRANX ' .or.
     .                 stype.eq.'XGC   ').and. chocc.eq.'   1' ) then
                  nr = nr + 1
                  qparn_tran(1,1) = nr
                  qscale(nr) = 1.d0
                  if( stype.eq.'XGC   ' ) qscale(nr) = -1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if((stype.eq.'TRANY ' .or.
     .                 stype.eq.'YGC   ').and. chocc.eq.'   1' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_tran(2,1) = nr
                  qscale(nr) = 1.d0
                  if( stype.eq.'YGC   ' ) qscale(nr) = -1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if((stype.eq.'TRANZ ' .or.
     .                 stype.eq.'ZGC   ').and. chocc.eq.'   1'  ) then
                  nr = nr + 1 
                  qscale(nr) = 1.d0
                  qparn_tran(3,1) = nr
                  qscale(nr) = 1.d0
                  if( stype.eq.'ZGC   ' ) qscale(nr) = -1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANRX' .or.
     .                 stype.eq.'XGCR  ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  if( stype.eq.'XGCR  ' ) qscale(nr) = -1.d0
                  qparn_tran(1,2) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANRY' .or.
     .                 stype.eq.'YGCR  ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  if( stype.eq.'YGCR  ' ) qscale(nr) = -1.d0
                  qparn_tran(2,2) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANRZ' .or.
     .                 stype.eq.'ZGCR  ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  if( stype.eq.'ZGCR  ' ) qscale(nr) = -1.d0
                  qparn_tran(3,2) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'SCALE ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_scale(1) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
               else if( stype.eq.'SCALER' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_scale(2) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
       
*                 Thats all we know about for the moment
              end if
          end if
          if( line(1:1).eq.'-' ) then
              end_block = .true.
          end if
      end do

****  Write out status
      write(*,400) nr, np
 400  format(' Found ',i4,' usable parameters, expected ',i4,
     .       ' parameters in dec_soln_est')
 
****  Thats all
      return
      end
 
 
CTITLE DEC_SOLN_APR
 
      subroutine dec_soln_apr(unit,line, np, cov_parm, sol_parm)
 
      implicit none

*     Routine to decode the solution apriori part of the sinex
*     files.  These values are written directly into the site_pos
*     vector.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
*   line        - Line read from file file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   sn      - Site number from list of sites
*   jp      - Parameter number read from the input list
*   cyr, cdoy, csec - Reference yr, doy and sec
*   indx        - Position in string
*   nr          - "real" parameter number in final covariance matrix
*   occ         - Occupation number
 
      integer*4 ierr, sn, jp, cyr, cdoy, csec, indx, nr, occ, i, jerr
      integer*4 svn, kerr
 
*   est     - Estimate from file
*   est_sig - Sigma of estimate.  Needed when correlations are given
*   ref_jd  - JD of reference epoch for value
 
      real*8 est, ref_jd, est_sig

* dut1, dlod -- Short period corrections to UT1 amd LOD

      real*8 dut1, dlod 
      
 
*   end_block   - Indicates end of block found
 
      logical end_block
 
*   code    - Site code read from Sinex file
*   stype   - Parameter type read from sinex file
*   estu        - Units for the estimate
*   pt      - Point characters
*   full_name  -  Full site name with occ and pt added
*   chocc   - Character occupancy
 
      character*8 code, stype, estu, full_name
      character*4 pt, chocc

* MOD TAH 210505: Added support for multi-frequency satellite
*     antenna offsets.
      integer*4 svant_off  ! Functon to return index in qparn_svs
                  ! array for antenna offset (L1/LC are in 
                  ! max_svs_elem - (2:0); Other are put in IC and
                  ! radiation parameter slots.
      integer*4 ao   ! Offset value returned from svant_off based
                  ! on pt string. 

 
****  Start reading the block
      end_block = .false.
      ierr = 0
      do i = 1, max_glb_parn
         atos(i) = 0
         atoo(i) = 0
         qapr_sig(i) = 0.d0
      end do
 
      do while ( ierr.eq.0 .and. .not.end_block)
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ' ) then

              if( cglb_vers.eq.5 ) then 
                  read(line,210,iostat=jerr) jp, stype, code,pt,occ, 
     .                 cyr, cdoy, csec, estu, est, est_sig
 210              format(i6,1x,a4,1x,a4,1x,a2,1x,i4,1x,i2,1x,i3,1x,i5,
     .                   1x,a4,2x,d21.14,1x,e14.8)
              else
                  read(line,215,iostat=jerr) jp, stype, code,pt,chocc, 
     .                 cyr, cdoy, csec, estu, est, est_sig
 215              format(i6,1x,a6,1x,a4,1x,a2,1x,a4,1x,i2,1x,i3,1x,i5,
     .                   1x,a4,2x,d22.14,1x,e11.5)
     
*                 Check the occupancy value
                  if( chocc.eq.'----' ) chocc = '   1'
                  read(chocc,*,iostat=jerr) occ 
                                   
              end if
              call report_error('IOSTAT',ierr,'decod',line,0,
     .                          'dec_soln_apr')
 
              call yds_to_jd( cyr, cdoy, csec, ref_jd)

* MOD TAH 970716: Check for blank pt names, replace with ' A'
              if( pt(2:2).eq.' ' .or. pt(2:2).eq.'-' ) pt = ' A'

****          OK, construct the full name for the site and
*             see if we find it.
              write(full_name, 240) code, occ, pt(2:2)
 240          format(a4,i3.3,a1)
              call casefold(full_name)
              indx = 1
              call get_cmd(full_name, qsite_names, qnum_sites,
     .                        sn, indx )

****          See if satellite axis offset
              if( stype(1:4).eq.'SATA' ) then
                  read(code,'(1x,i3)',iostat=kerr) svn
                  sn = 0
                  do i = 1, qnum_svs
                     if( qsvi_svn(i).eq.svn ) sn = i
                  end do
                  if( kerr.ne.0 .or.sn.eq.0 ) then
                     write(*,260) line(1:80)
 260                 format('***WARNING*** Unable to find SV in: ',a)
                  endif
* MOD TAH 210505: Extract out frequency for antenna offsets to handle
*                 TUG L1, L2, L5, L7,L8 values.  Save as offset to 
*                 max_svs_elem (see SINEX has not IC or Rad parameters
*                 use these slots for non-L1, LC values.
                  ao = svant_off(pt)
              endif
              if( stype.eq.'STAX' .and. sn.gt.0 ) then
                  site_pos(1,sn) = est
                  nr = qparn_sites(1,sn)
                  qapr(nr) = est
                  qapr_sig(jp) = est_sig
                  qapr_cov(1,1,sn) = est_sig**2
                  atos(jp) = sn
              else if( stype.eq.'STAY' .and. sn.gt.0 ) then
                  site_pos(2,sn) = est
                  nr = qparn_sites(2,sn)
                  qapr(nr) = est
                  qapr_sig(jp) = est_sig
                  qapr_cov(2,2,sn) = est_sig**2
                  atos(jp) = sn
              else if( stype.eq.'STAZ' .and. sn.gt.0 ) then
                  site_pos(3,sn) = est
                  nr = qparn_sites(3,sn)
                  qapr(nr) = est
                  qapr_sig(jp) = est_sig
                  qapr_cov(3,3,sn) = est_sig**2
                  atos(jp) = sn
              else if( stype.eq.'VELX' .and. sn.gt.0 ) then
                  site_vel(1,sn) = est
                  nr = qparn_vel(1,sn)
                  if( nr.gt.0 ) then
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atos(jp) = -sn 
                  end if
              else if( stype.eq.'VELY' .and. sn.gt.0 ) then
                  site_vel(2,sn) = est
                  nr = qparn_vel(2,sn)
                  if( nr.gt.0 ) then
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atos(jp) = -sn 
                  end if
              else if( stype.eq.'VELZ' .and. sn.gt.0 ) then
                  site_vel(3,sn) = est
                  nr = qparn_vel(3,sn)
                  if( nr.gt.0 ) then
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atos(jp) = -sn 
                  end if

* MOD TAH 050710: Addin satellite axis offsets
              else if( stype.eq.'SATA_X' .and. sn.gt.0 
     .                 .and. est_sig.gt. 0 ) then
C MOD TAH 210505:    Changed access for multi-frequency antenna
C                    offset (TUG SINEX)
C                 nr = qparn_svs(max_svs_elem-2,sn)
                  nr = qparn_svs(ao+1,sn)
                  if( nr.gt.0 ) then
                     qapr(nr) = est
C MOD TAH 210505:    Changed access for multi-frequency antenna
C                    offset (TUG SINEX)
C                    svs_pos(max_svs_elem-2,sn) = est
                     svs_pos(ao+1,sn) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  end if
              else if( stype.eq.'SATA_Y' .and. sn.gt.0
     .                 .and. est_sig.gt. 0 ) then
C MOD TAH 210505:    Changed access for multi-frequency antenna
C                    offset (TUG SINEX)
                  nr = qparn_svs(ao+2,sn)
                  if( nr.gt.0 ) then
                     qapr(nr) = est
C MOD TAH 210505: Changed access for multi-frequency antenna
C                 offset (TUG SINEX)
                     svs_pos(ao+2,sn) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  end if
              else if( stype.eq.'SATA_Z' .and. sn.gt.0
     .                 .and. est_sig.gt. 0 ) then
C MOD TAH 210505:    Changed access for multi-frequency antenna
C                    offset (TUG SINEX)
                  nr = qparn_svs(ao+3,sn)
                  if( nr.gt.0 ) then 
                     qapr(nr) = est
C MOD TAH 210505: Changed access for multi-frequency antenna
C                 offset (TUG SINEX)
                     svs_pos(ao+3,sn) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  end if                     

*             See if Earth orientation results
              else if ( stype.eq.'LOD ' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.  Convert LOD to UT1 Rate
*                 (mas per day is the unit)
                  nr = qparn_mul_pmu(2,3,occ)
                  if( nr.gt.0 ) then
                     qapr(nr) = -est*15.d0
*                    Special JPL Mod
                     if( est.lt.0 ) then
                         qapr(nr) = -qapr(nr)
                     end if
                     qapr_sig(jp) = est_sig*15.d0
                     qpmu_mul_apr(2,3,occ) = qapr(nr)
                  end if
              else if( stype.eq.'UT  ' .or. stype.eq.'UTR  ') then

* MOD TAH 991207: Check GFZ values.  Sign seems to be wrong
                  if( qowner(1:3).eq. 'GFZ' .and. est.gt. 1000.d0 ) then
                      write(*,*) 'Changing SIGN of GFZ UT values'
                      est = -est
                  end if
                  
*                 See if we need to restore short-period UT1 corrections.
*                 Need to CHANGE JPL Date when problem fixed.
                  if ( stype.eq.'UTR ' .or. 
     .                (qowner(1:3).eq. 'JPL' .and. 
     .                 ref_jd.lt.  2451544.50d0) )    then
                       call short_period_lod( ref_jd, dut1, dlod)
*                      Add correction; convert sec to ms.                       
                       est = est + dut1*1.d3
                  end if

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  nr = qparn_mul_pmu(1,3,occ)
                  if( nr.gt.0 ) then
                     qapr(nr) = est*15.d0
                     qpmu_mul_apr(1,3,occ) = qapr(nr)
                     qapr_sig(jp) = est_sig*15.d0
                  end if
              else if( stype.eq.'XPO ' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  nr = qparn_mul_pmu(1,1,occ)
                  if( nr.gt.0 ) then
                     qapr(nr) = est
                     qpmu_mul_apr(1,1,occ) = qapr(nr)
                     qapr_sig(jp) = est_sig
                  end if
              else if( stype.eq.'XPOR' ) then
                  nr = qparn_mul_pmu(2,1,occ)
                  if( nr.gt.0 ) then
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     if( estu.eq.'ma/s' ) then
                         qapr_sig(jp) = est_sig*86400.d0
                         qapr(nr) = est*86400.d0
                     endif
                     qpmu_mul_apr(2,1,occ) = qapr(nr)
                  end if

              else if( stype.eq.'YPO ' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  nr = qparn_mul_pmu(1,2,occ)
                  if( nr.gt.0 ) then
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                  end if
                  qpmu_mul_apr(1,2,occ) = qapr(nr)
              else if( stype.eq.'YPOR' ) then
                  nr = qparn_mul_pmu(2,2,occ)
                  if( nr.gt.0 ) then
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     if( estu.eq.'ma/s' ) then
                         qapr(nr) = est*86400.d0
                         qapr_sig(jp) = est_sig*86400.d0
                     endif
                     qpmu_mul_apr(2,2,occ) = qapr(nr)
                  end if

               else if((stype.eq.'TRANX ' .or.
     .                 stype.eq.'XGC   ').and. chocc.eq.'   1' ) then
                  nr = qparn_tran(1,1)
                  if( nr.gt.0 ) then 
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  endif
              else if((stype.eq.'TRANY ' .or.
     .                 stype.eq.'YGC   ').and. chocc.eq.'   1' ) then
                  nr = qparn_tran(2,1)
                  if( nr.gt.0 ) then 
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  endif
              else if((stype.eq.'TRANZ ' .or.
     .                 stype.eq.'ZGC   ').and. chocc.eq.'   1'  ) then
                  nr = qparn_tran(3,1)
                  if( nr.gt.0 ) then 
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  endif
              else if( stype.eq.'TRANRX' .or.
     .                 stype.eq.'XGCR  ' ) then
                  nr = qparn_tran(1,2) 
                  if( nr.gt.0 ) then 
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  endif
              else if( stype.eq.'TRANRY' .or.
     .                 stype.eq.'YGCR  ' ) then
                  nr = qparn_tran(2,2) 
                  if( nr.gt.0 ) then 
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  endif
              else if( stype.eq.'TRANRZ' .or.
     .                 stype.eq.'ZGCR  ' ) then
                  nr = qparn_tran(3,2) 
                  if( nr.gt.0 ) then 
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  endif
              else if( stype.eq.'SCALE ' ) then
                  nr = qparn_scale(1) 
                  if( nr.gt.0 ) then 
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  endif
               else if( stype.eq.'SCALER' ) then
                  nr = qparn_scale(2) 
                  if( nr.gt.0 ) then 
                     qapr(nr) = est
                     qapr_sig(jp) = est_sig
                     atoo(jp) = nr
                  endif
               else
                  print *,'APR: Drop through ',trim(line)                               
*             Thats all we know about for the moment
              end if
          end if
          if( line(1:1).eq.'-' ) then
              end_block = .true.
          end if
      end do
 
****  Thats all
      return
      end
 
CTITLE DEC_SOLN_MAT
 
      subroutine dec_soln_mat(unit, line, np, cov_parm, sol_parm )
 
      implicit none

*     Routine to decode the covariance matrix of the estimates from
*     the sinex files.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
*   line        - Line read from file file
 
      character*(*) line
 
* LOCAL VARAIBLES
 
 
*   ierr        - IOSTAT error
*   np1, np2    - Parameter numbers
*   nr1, nr2    - Parameter numbers in the output solution
*   i,j     - Loop counters
*   lim     - Number of elements to copy with each read
*   max_np  - Largest parameter number expected when reading
*             upper triangular forms.  (Found from itoo array).
 
      integer*4 ierr, np1, np2, i,j, k, lim, nr1, nr2, max_np, jerr
 
*   cov(3)      - Up to thress values read from file.
 
      real*8 cov(3)
 
*   end_block   - Indicates end of block found
*   corr_mat    - Set true if this is a correlation matrix rather than
*                 covariance matrix
*   lower       - Set true to indicate lower triangular form
*   noreport    - Set false once we have reported sigmas on diagonal.
 
 
      logical end_block, corr_mat, lower, noreport
 
****  Start reading the block
      end_block = .false.
      ierr = 0
      max_np = 0

*     See if correlation or covariance matrix
      corr_mat = .false.
      if( index(line,'CORR').gt.0 ) corr_mat = .true.
      lower    = .true. 
      if( index(line,'U C').gt.0  ) lower = .false.

*     If this is a correlation matrix; set the diagonal first
      if( corr_mat ) then
          do i = 1, np
             cov_parm(i,i) = qsig(i)**2
          end do
      end if

*     IF not a lower form then see maxparam expected
      if( .not.lower ) then
          max_np = 0
          i = max_glb_parn
          do while ( i.gt.1 .and. max_np.eq.0 )
              i = i - 1
              if( itoo(i).ne.0 ) max_np = i
          end do
      end if
 
      noreport = .true.

      do while ( ierr.eq.0 .and. .not.end_block)
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 ) then
 
              if( line(1:1).eq.' ' ) then
* MOD TAH 960903: Cleared the cov array incase there are
*                 blanks in line which will retain there
*                 old value.             
                  do j = 1,3
                      cov(j) = 0.d0
                  end do
                  read(line, * ,iostat=jerr) np1, np2, (cov(j), j=1,3)
cc 210              format(2i6,3d22.13)

*****             Allow for -1 errors on read with free format (to support
*                 GSZ files)
                  if( ierr.eq.-1 ) ierr = 0
                  if( ierr.ne.0 .and. ierr.ne. -1 ) then
                      write(*,240) ierr, line      
 240                  format('IOSTAT error ',i5,' Reading COVA block',/,
     .                       'Line:',a)
                      stop 'Error reading COVA'
                  end if
                  
* MOD TAH 141114: Changed to read correlation matrix as if it were a covariance
*                 matrix and the convert at the end.  Needed for U CORR matrices
*                 since we don't have full variance available as file is read.
*                 Now assign the elements to the matrix
C                 if( corr_mat ) then
C                     lim = np1 - np2
C                     
* MOD TAH 960614: Check added for Ver 1.0 sinex.
*                     See if sigma is on the diagonal
C                     if( lim.lt.3 .and. cov(lim+1).gt.0 .and.
C    .                    cglb_vers.gt.5  ) then
C                         nr1 = itoo(np1)
*                         Make sure this parameter still being
*                         used.                          
C                         if( nr1.gt.0 ) then
* MOD TAH 071011: Re-normalize the covariances which were originally
*                 compute with lower accuracy qsig
C                            do j = 1,nr1-1
C                                cov_parm(nr1,j) = cov_parm(nr1,j)*
C    .                                             cov(lim+1)/qsig(nr1)
C                            end do
C                            qsig(nr1) = cov(lim+1)
C                            cov_parm(nr1,nr1) = qsig(nr1)**2
C                         end if
C                      end if 
C                      print *,'Corr ', lim, np1, np2, nr1, 
C    .                                cov_parm(nr1,1:nr1-1), ' | ', 
C    .                                cov(lim+1),qsig(nr1)
C                   
C                 else
                      if ( lower ) then
                          lim = np1 - np2 + 1
                      else

* MOD TAH 960822: Set limit on max parameter number based on max
*   in input covaraince matrix.
                          lim = max_np - np2 + 1
                      end if
C                 end if
                  if( lim.gt.3 ) lim = 3

                  do j = 1, lim
*                    See if we need to convert to covariance
                     nr1 = itoo(np1)
                     nr2 = itoo(np2+j-1)
c                     write(*,710) np1, np2, nr1,nr2, cov(j)
 710                 format('N ',4i5,3E16.6)
* MOD TAH 1411414: Removed special treatement pof correlation matrix
C                   if( corr_mat .and. nr1.gt.0 .and. nr2.gt.0 ) then

* MOD TAH 981006:  See if diagonal element.
C                         if( nr1.eq.nr2 ) then
C                             if( cov(j).ne.1.d0 ) then 
*                                 see if variance or sigma
C                                 if( abs(cov(j)/qsig(nr1)-1.d0).
C    .                                             lt.1.d-6 ) then
C                                     qsig(nr1) = cov(j)
C                                     cov(j) = 1.d0
C                                     if ( noreport ) then
C                                         write(*,*) 
C    .                                    '***Sigmas found on diagonal'
C                                         noreport = .false.
C                                     endif
C                                 end if
C                             end if
C                         end if
C                         cov_parm(nr1, nr2) = cov(j)*qsig(nr1)*
C    .                                                qsig(nr2)
C                    end if
C                    if( .not.corr_mat .and. 
C    .                   nr1.gt.0 .and. nr2.gt.0 ) then
                     if( nr1.gt.0 .and. nr2.gt.0 ) then
                          cov_parm(nr1, nr2) = cov(j)
                          cov_parm(nr2, nr1) = cov(j)
                     end if
                  end do
              else
                  if( line(1:1).eq.'-' ) then
                      end_block = .true.
                  end if
              end if
          end if
      end do

* MOD TAH 141114: If this is correlation matrix ; convert to covariance
      if( corr_mat ) then
          write(*,'(a)') ' Converting correlations to covariances'
          if( cov_parm(1,1).eq.1.d0 ) then
             write(*,'(a)') ' Diagonal not variance using qsig only'
          else   ! Variances or sigmas on diagonaal
             noreport = .true.
             do i = 1, np 
                 if( abs(cov_parm(i,i)/qsig(i)-1.d0).
     .                    lt.1.d-5 ) then
                     qsig(i) = cov_parm(i,i)
                     cov_parm(i,i) = qsig(i)**2
                     if ( noreport ) then
                         write(*,*) 
     .                   '***Sigmas found on diagonal NP ',i
                         noreport = .false.
                     endif
                 else     ! Save sigma
                     if( abs(cov_parm(i,i)/qsig(i)**2-1.d0)
     .                   .gt.1.d-5 ) then
                          write(*,*) 
     .                   '*WARNING: qsig and cov_parm not consistent',
     .                        i, sqrt(cov_parm(i,i)), qsig(i)
                     end if 
                     qsig(i) = sqrt(cov_parm(i,i))   
                 end if
               end do
          end if  
*         Now covert correlations to covariances
          do i = 1, np
             do j = 1, np
                if( i.ne.j ) cov_parm(i,j) = cov_parm(i,j)*
     .                                       qsig(i)*qsig(j) 
             end do
          end do
      end if
 
****  Now copy the matrix into the upper triangular section
c      do i = 1, np
c         write(*,999) i, qsig(i)
c 999     format(i4,e22.10) 
c      end do
      
      do i = 1, np - 1
          do j = i+1, np
              cov_parm(i,j) = cov_parm(j,i)
          end do
      end do

      if( 1.eq.2 ) then 
         write(*,997) (i,itoo(i),i=1,np)
 997     format(550(i3,1x,i3,' | '))
* MOD TAH 141114: Changed correlation matrix to be imaged.
         do i = 1, np
            write(*,998 ) i,(int(100*(cov_parm(i,j)/
     .                    sqrt(cov_parm(i,i)*cov_parm(j,j)))),j=1,np)
 998        format('C ',i5,3000(1x,I4))
         end do
      end if

****  Now check the block diagonal in the covariance to make sure
*     no negative elements on diagonal and not correclations greater
*     than +-1.
      call check_cov_htog(cov_parm, np )


****  Thats all
      return
      end

CTITLE CHECK_COV_HTOG

      subroutine check_cov_htog( cov_parm, np )

      implicit none

*     Routine to check block diagonal on cov_parm.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   np      - number of parameters expected
 
      integer*4 np
 
* cov_parm(np,np) - Covarince matrix
 
 
      real*8 cov_parm(np,np)

 
* LOCAL VARAIBLES
      integer*4 i,j,k  ! Loop counters
     .,         npx    ! Parameter number for x-coordinate
      character*1  err(3,3) ! Set to 1 for each bad element

      logical bad   ! Set true if covariance is bad.

      real*8 corr   ! Correlation coefficient
      
      character*256 mess  ! Error message

****  Loop over the sites
      do i = 1, qnum_sites
         npx = qparn_sites(1,i)
         do j = 1,3
            do k = 1,3
               err(j,k) = '-'
            enddo
         enddo
         bad = .false.
         if( npx.gt.0 ) then 
            do j = 0,2
               if( cov_parm(npx+j,npx+j).le.0 ) then
                   bad = .true.
                   err(j+1,j+1) = 'N'
               end if
*              if not bad check off diagonal
               if( .not.bad ) then
                   do k = 0,2
                      if( k.ne.j ) then
                         corr = cov_parm(npx+j,npx+k)/
     .                          sqrt(cov_parm(npx+j,npx+j)*
     .                              cov_parm(npx+k,npx+k))
                         if( abs(corr).gt.1.d0 ) then
                             err(k+1,j+1) = 'C'
                         endif
                      endif
                   enddo
               endif
            enddo
            do j = 1,3
               do k = 1,3
                  if( err(j,k).eq.'C' .or.  err(j,k).eq.'N' ) 
     .                bad = .true.
               end do
            end do  
         endif
         

*        See if element is bad
         if ( bad ) then
            write(mess,220) qsite_names(i)(1:4), npx, err
 220        format('BAD_COVAR: Site ',a4,' Parameter ',i5,
     .             ' Neg. Diag. or correlation >1. Elements ',
     .               3(3a1,1x,':'))
            call report_stat('WARNING','HTOGLB','check_cov_htog',
     .             hfile, mess,0)
            write(*,420) npx, ((cov_parm(npx+k,npx+j),k=0,2),j=0,2)
 420        format('B ',i5,1x,9(e11.5,1x))
            do j = 0,2
               cov_parm(npx+j,npx+j) = 100.d0
            end do
         end if
      end do

***** Thats all
      return
      end
    
 
 
CTITLE DEC_SOLN_CON
 
      subroutine dec_soln_con(unit, line, np, cov_parm, sol_parm )
 
      implicit none

*     Routine to decode the apriori contraint matrix.  If the 
*     contraints are tight then they are removed from the solution.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
*   line        - Line read from file file
 
      character*(*) line
 
* LOCAL VARAIBLES
 
 
*   ierr        - IOSTAT error
*   np1, np2    - Parameter numbers
*   nr1, nr2    - Parameter numbers in the output solution
*   i,j     - Loop counters
*   lim     - Number of elements to copy with each read
*   l1, l2  - Elements in 3*3 matrix for site constraints
*   sn      - Site number
*   fsp     - First parameter for a station (used to get appropriate
*             coordinate).
 
      integer*4 ierr, np1, np2, i,j,k, lim, nr1, nr2,
     .          l1,l2, sn, fsp, jerr, jp
 
*   cov(3)      - Up to thress values read from file.
 
      real*8 cov(3)
 
*   end_block   - Indicates end of block found
*   corr_mat    - Set true if this is a correlation matrix rather than
*                 covariance matrix
*   first_done  - Set true once the diagonal has been done.  Constraints
*                 are assumed to be either diagonal or 3*3 values.
 
 
      logical end_block, corr_mat, first_done
 
****  Start reading the block
      end_block = .false.
      first_done = .false.
      ierr = 0

****  Start by inverting the covariance matrix to generate "normal equations"
*     use the difference between the apriori and estimates in sol_parm

      do i = 1, qnum_sites
         do j = 1,3
            do k = 1, 3
               qapr_cov(j,k,i) = 0.d0
               qvel_cov(j,k,i) = 0.d0
            end do
         end do
      end do

*     See if correlation or covariance matrix
      corr_mat = .false.
      if( index(line,'CORR').gt.0 ) corr_mat = .true.

*     If correlation matrix, copy over the diagonal elements
      if( corr_mat ) then
         do i = 1, np , 3
            sn = atos(i)    
            if( sn.gt.0 ) then
                qnum_apr_cov = qnum_apr_cov + 1
                do j = 1,3
                   qapr_cov(j,j,sn) = qapr_sig(i+j-1)**2
                end do
            else if ( sn.lt.0 ) then
                qnum_vel_cov = qnum_vel_cov + 1 
                do j = 1,3
                   qvel_cov(j,j,-sn) = qapr_sig(i+j-1)**2 
                end do
            end if
         end do
      end if
      
****  Get first site coordinate.
      fsp = -1
      i   = 0
      do while ( fsp.lt.0 .and. i.lt.qnum_parn )
          i = i + 1
          if( atos(i).gt.0 ) fsp = i - 1
      end do

****  Now read the entries to see what is in the full matrix      
      do while ( ierr.eq.0 .and. .not.end_block)
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 ) then
 
              if( line(1:1).eq.' ' ) then
                  read(line, * ,iostat=jerr) np1, np2, (cov(j), j=1,3)
 210              format(2i6,3d22.13)
*                 Allow for free format errors
                  if( jerr.eq.-1 ) jerr = 0
                  if( jerr.ne.0 ) then
                      write(*,240) jerr, line      
 240                  format('IOSTAT error ',i5,' Reading COVA block',/,
     .                       'Line:',a)
                      stop 'Error reading COVA'
                  end if

*                 Now assign the elements to the matrix.  See if parameter
*                 corresponds to site number.
                  sn = atos(np1)
                  lim = np1 - np2 + 1
                  if( lim.gt.3 ) lim = 3
*                 If this is not for a site position, then set lim to zero
*                 so that do loop will not iterate.
                  if( sn.eq.0 ) lim = 0
                  do j = 1, lim
*                    See if we need to convert to covariance
                     nr1 = np1
                     nr2 = np2+j-1
*                    The following assumes that the apriori's were given
*                    in groups of three
                     l1  = mod(nr1-fsp,3)    
                     if( l1.eq.0 ) l1 = 3
                     l2  = mod(nr2-fsp,3)
                     if( l2.eq.0 ) l2 = 3

*                    Assign values if valid parameter numbers and correlation
*                    matrix.
                     if( corr_mat .and. nr1.gt.0 .and. nr2.gt.0 ) then
                          if( sn.gt.0 ) then
                              if( nr1.ne.nr2 ) 
     .                            qapr_cov(l1, l2 ,sn) = cov(j)*
     .                                            qapr_sig(nr1)*
     .                                            qapr_sig(nr2)
* MOD TAH 960614: Diagonal check for Ver 1.0 CORR files.
*                             Check if sigma on the diagonal.
                              if( nr1.eq.nr2 .and. cglb_vers.gt.5 )
     .                            qapr_cov(l1, l2 ,sn) = cov(j)**2
                             
                          else
                              qvel_cov(l1, l2 ,-sn) = cov(j)*
     .                                        qapr_sig(nr1)*
     .                                        qapr_sig(nr2)
                          end if
                     end if
                     
*                    Now do covariance matrix, if parameter
*                    numbers are not zero.                     
                     if( .not.corr_mat .and. 
     .                   nr1.gt.0 .and. nr2.gt.0 ) then

*                         Count number only for first element
                          if( l1.eq.1 .and. sn.gt.0 ) 
     .                             qnum_apr_cov = qnum_apr_cov+1
                          if( l1.eq.1 .and. sn.lt.0 ) 
     .                             qnum_vel_cov = qnum_vel_cov+1
                          if( sn.gt.0 ) then
                              qapr_cov(l1, l2, sn) = cov(j)
                              qapr_cov(l2, l1, sn) = cov(j)
                          else
                              qvel_cov(l1, l2,-sn) = cov(j)
                              qvel_cov(l2, l1,-sn) = cov(j)
                          endif
                     end if
                  end do
**** MOD TAH 070111: See if this is translation, scale or satellite antenna
*                 offset
                  if( lim.eq.0 ) then
                      lim = np1 - np2 + 1
*                     See if we can find this parameter.  Too far off
                      if ( lim.ge.1 .and. lim.le.3 ) then
                         jp = np1
*                        Only save diagonal
cd                            write(*,510) jp, atoo(jp), qapr_sig(jp),
cd   .                             cov(lim)  
cd 510                     format('Replacing apr param ',i5,
cd     .                    ' Soln parn ',i5,' from ',e23.12,
cd     .                    ' to ',e23.12,e23.12)
                         if( corr_mat ) then
                             qapr_sig(jp) = cov(lim)
                         else
cd                            write(*,510) jp, atoo(jp), qapr_sig(jp),
cd    .                             sqrt(cov(lim))
                             qapr_sig(jp) = sqrt(cov(lim))
                         endif
                      endif 
                  end if

              else
                  if( line(1:1).eq.'-' ) then
                      end_block = .true.
                  end if
              end if
          end if
      end do
 
 
****  Thats all
      return
      end

CTITLE DEC_REMV_CON
 
      subroutine dec_remv_con(unit, line, np, cov_parm, sol_parm,  
     .           cov_copy, eigvec, eigval, bwork, zwork)
 
      implicit none

*     Routine to decode the apriori contraint matrix.  If the 
*     contraints are tight then they are removed from the solution.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
* cov_copy(np,np) - Copy of covariance matrx
* eigvec(np,np)  - MAtrix of eigenvectors 
* eigval(np)      - Eigen values; checked for negative values
* bwork(np), zwork(np) - Work arrays used by jacobi
 
 
      real*8 cov_parm(np,np), sol_parm(np),
     .       cov_copy(np,np+1), eigvec(np,np+1), eigval(np),
     .       bwork(np), zwork(np)

      integer*4 eigen_bins(15), bin 
 
 
*   line        - Line read from file file
 
      character*(*) line
 
* LOCAL VARAIBLES
 
*   nr1, nr2    - Parameter numbers in the output solution
*   i,j     - Loop counters
*   ipivot(max_glb_parn)  - Pivot elements for matrix inverse
*   nrot    - Number of rotations
*   eigen_iter  - Number of iterations through eigenroutine.
 
      integer*4  i,j,k,  nr1, nr2, ipivot(max_glb_parn), nrot,
     .           eigen_iter, l, m, p1, p2
 
*   cov(3)      - Up to thress values read from file.
*   scale(max_glb_parn) - Scale factors for matrix inverse
*   cov_apr(3,3)  - 3*3 matrix to invert for to remove from normal
*                   equations
*   snv         - snx_constraint variance
 
      real*8 cov(3), scale(max_glb_parn), cov_apr(3,3), snv,
     .       rot_var, tran_var

*   min_eigval, max_eigval - Min and max Eigenvalues

      real*8 min_eigval, max_eigval, pmu_rbrot(3,3,max_glb_sites)

*   num_nege  - Number of negative eigenvalues
      integer*4 num_nege

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      data I8 / 1 /

* MOD TAH 180202: Moved to -d=TR constrainted its owwn subroutine so that
*     cam also be used on GAMIT hfiles

      call decon_strapp(decon_str, np, cov_parm, sol_parm)

* MOD TAH 991206: See if user wants to compute eigenvalues
      if( comp_eigen ) then
* MOD TAH 180202: Moved copy of covariannce matrix here when decon_strapp
*     subroutine introced.  (To allow application to GAMIT BASELINE h-files).
         call dwmov8(cov_parm, 1, cov_copy, 1, (I8*np)*np )
         call jacobi(cov_copy,np,np,eigval,eigvec,nrot,
     .               bwork,zwork)

         write(*,*) 'Eigen values took ',nrot,' rotations'
     
         min_eigval = 1.d20
         max_eigval = -1.d20
         num_nege   = 0 
         do i = 1,15
            eigen_bins(i) = 0
         end do

         do i = 1, np 

*           Get the bin
            bin = int(log10(eigval(i)))+11
            if( bin.le.0 ) bin = 1
            if( bin.ge.15 ) bin = 15
            eigen_bins(bin) = eigen_bins(bin) + 1 

            if( eigval(i).lt.min_eigval ) min_eigval = eigval(i)
            if( eigval(i).gt.max_eigval ) max_eigval = eigval(i)
            if( eigval(i).lt.0.d0 ) num_nege = num_nege + 1
c           write(*,450) i, eigval(i), sqrt(cov_parm(i,i))
c 450       format('NE ',I4,' Eigen ',d20.8,' Sigma ',F20.8)
         end do
         write(*,460) min_eigval, max_eigval, num_nege
  460    format(' Min and Max Eigenvalues ',2d16.4, ' Negs ',i4)
         write(*,470) (i,i=-10,4)
  470    format('Eigenvalue Distribution: powers of 10',/, 
     .          'Eigenvalue Distribution ', 15I4) 
         write(*,475) (eigen_bins(i),i=1,15)
  475    format('Eigenvalue Distribution ', 15i4)   
      else
         write(*,*) 'Eigenvalues not checked'
      end if

****  Now copy the matrix into the upper triangular section
      if( snx_constraint.ne.0 ) then
         do i = 1, np

*           Make sure we have an apriori.  If we don't use the the estimate
*           as the aproiri
            if( qapr(i).eq.0.d0 ) then
                qapr(i) = qsol(i)
                write(*,110) i, qsol(i)
 110            format('*** WARNING *** No apriori for parameter ',i4,
     .                 ' Using estimate ',f20.4)
            end if

            sol_parm(i) = qsol(i) - qapr(i)
         end do
         call invert_vis(cov_parm, sol_parm, scale, ipivot, np, np, 1)
C     end if

C     if( snx_constraint.ne.0 ) then
         write(*,500) snx_constraint
 500     format(' Replacing User constraints with ',F10.4,' m')

         snv = snx_constraint**2
         qnum_apr_diag = 0

         do i = 1, qnum_sites

****        See if constraints are significant
            do j = 1,3
               cov(j) = qapr_cov(j,j,i)
               do k = 1,3
                  cov_apr(j,k) = qapr_cov(j,k,i)
               end do
            end do

            if((cov(1).lt.snv.and.cov(1).gt.0.d0) .or.
     .         (cov(2).lt.snv.and.cov(2).gt.0.d0) .or.
     .         (cov(3).lt.snv.and.cov(3).gt.0.d0)    ) then
                write(*,510) qsite_names(i), cov
 510            format('Site ',a8,' Constraints ',3E15.4)

*               Invert the apriori covariance matrix
                call invert_vis(cov_apr,cov, scale, ipivot,3,3,0)

*               Now remove
                do j = 1,3
                   nr1 = qparn_sites(j,i)
                   qnum_apr_diag = qnum_apr_diag + 1
                   do k = 1,3
                      nr2 = qparn_sites(k,i)
                      cov_parm(nr1,nr2) = cov_parm(nr1,nr2) - 
     .                                    cov_apr(j,k)

*                     Check the diagonal and add new constraint
                      if ( k.eq.j ) then
                          if( cov_parm(nr1,nr2).lt.0 ) then
                              write(*,550) i, j, cov_parm(nr1,nr2),
     .                              cov_apr(j,k)                  
 550                          format('**ERROR** Site ',i4,' CMP ',i2,
     .                               ' Diagonal < 0 ',f21.12,
     .                               1x,' Apr ',f21.12)
                          end if
                          cov_parm(nr1,nr2) = cov_parm(nr1,nr2) + 1/snv
                          qapr_cov(j,k,i) = snv  
                      else
                          qapr_cov(j,k,i) = 0.d0
                      end if
                   end do
                end do
            end if
         end do

*****    Now see if we velocity constraints that need to be removed.

         do i = 1, qnum_sites

****        See if constraints are significant
            do j = 1,3
               cov(j) = qvel_cov(j,j,i)
               do k = 1,3
                  cov_apr(j,k) = qvel_cov(j,k,i)
               end do
            end do
            if((cov(1).lt.snv.and.cov(1).gt.0.d0) .or.
     .         (cov(2).lt.snv.and.cov(2).gt.0.d0) .or.
     .         (cov(3).lt.snv.and.cov(3).gt.0.d0)    ) then
                write(*,610) qsite_names(i), cov
 610            format('Site ',a8,' Velconst    ',3F15.8)

*               Invert the apriori covariance matrix
                call invert_vis(cov_apr,cov, scale, ipivot,3,3,0)

*               Now remove
                do j = 1,3
                   nr1 = qparn_vel(j,i)
                   qnum_apr_diag = qnum_apr_diag + 1
                   do k = 1,3
                      nr2 = qparn_vel(k,i)
                      cov_parm(nr1,nr2) = cov_parm(nr1,nr2) - 
     .                                    cov_apr(j,k)

*                     Check the diagonal and add new constraint
                      if ( k.eq.j ) then
                          if( cov_parm(nr1,nr2).lt.0 ) then
                              write(*,550) i, j, cov_parm(nr1,nr2)
                          end if
                          cov_parm(nr1,nr2) = cov_parm(nr1,nr2) + 1/snv
                          qvel_cov(j,k,i) = snv  
                      else
                          qvel_cov(j,k,i) = 0.d0
                      end if
                   end do
                end do
            end if
         end do

****     Now convert the "normal equations" back to covariance matrix.
         call invert_vis(cov_parm, sol_parm, scale, ipivot, np, np, 1)
         do i = 1, np
            qsol(i) = sol_parm(i) + qapr(i)
            if( cov_parm(i,i).lt.0 ) write(*,580) i, cov_parm(i,i), 
     .            sol_parm(i)
 580        format(' Parameter ',i4,' Neg. Var ',2d20.10)
         end do
         
*        Now check the matrix for negative eignvals
         eigen_iter = 0
*        Set number of negative eigenvalues >1 so that we get into
*        test loop.
         num_nege = 1
         do while ( eigen_iter.lt.5 .and. num_nege.gt.0 )

             eigen_iter = eigen_iter + 1

             call dwmov8(cov_parm, 1, cov_copy, 1, (I8*np)*np )

             call jacobi(cov_copy,np,np,eigval,eigvec,nrot,
     .                   bwork,zwork)
c             call eigen( cov_copy,np, np, 1, eigvec, nrot, eigval, 
c    .            1.d-6, zwork)

*            Now check the eigenvalues
             min_eigval = 1.d20
             max_eigval = -1.d20
             num_nege   = 0
             do i = 1, np
                if( eigval(i).lt.min_eigval ) min_eigval = eigval(i)
                if( eigval(i).gt.max_eigval ) max_eigval = eigval(i)
                if( eigval(i).le.0.d0 ) num_nege = num_nege + 1
             end do
             write(*,590) eigen_iter, min_eigval, max_eigval
 590         format('At Iter ',i2,': Smallest and largest ',
     .              'Eigenvalues: ',2d14.5)
             if( num_nege.gt.0 ) then
                 write(*,595) eigen_iter, num_nege
 595             format('At Iter ',i2,': There are ',i4,
     .                  ' Negative Eigenvalues, Rescaling diagonal')
                 do i = 1,np
                    cov_parm(i,i) = cov_parm(i,i) + 1.d-6*cov_parm(i,i)
                 end do
             end if


          end do

*         Make sure we ended up with no negative eigenvalues
          if( num_nege.gt.0 ) then 
              write(*,597) num_nege
 597          format('**WARNING** ',i4,' Negative Eigenvalues Remain')
          end if


*     ELSE

*        At this point count the number of site constraints saved at this
*        time
         qnum_apr_diag = 0
         qnum_apr_cov  = 0
         qnum_vel_cov  = 0
         do i = 1, qnum_sites

****        See if constraints are significant
            if( qapr_cov(1,1,i).gt.0 .or. qapr_cov(2,2,i).gt.0 .or.
     .          qapr_cov(3,3,i).gt.0 ) qnum_apr_cov = qnum_apr_cov + 1
            if( qvel_cov(1,1,i).gt.0 .or. qvel_cov(2,2,i).gt.0 .or.
     .          qvel_cov(3,3,i).gt.0 ) qnum_vel_cov = qnum_vel_cov + 1
         end do

      ENDIF

****  See if antenna offsets and CoM constraints to be replaced
      if( index(decon_str,'C').gt.0 .or.index(decon_str,'A').gt.0 ) then

*        Remove apriori 
         do i = 1, np

*           Make sure we have an apriori.  If we don't use the the estimate
*           as the aproiri
            if( qapr(i).eq.0.d0 ) then
                qapr(i) = qsol(i)
cd               write(*,710) i, qsol(i)
cd 710            format('*** WARNING *** No apriori for parameter ',i4,
cd     .                 ' Using estimate ',f20.4)
            end if

            sol_parm(i) = qsol(i) - qapr(i)
         end do
         call invert_vis(cov_parm, sol_parm, scale, ipivot, np, np, 1)

****     Now see which constraints to remove.
         if( index(decon_str,'C').gt. 0 ) then
            do i = 1,3
               p1 = qparn_tran(i,1)
               if( p1.gt.0 ) then
                  do j = 1, np
                     if( atoo(j).eq.p1 ) p2 = j
                  end do
               end if
               if ( p1.gt.0 .and. p2.gt.0 ) then
                   write(*,720) i, p1,p2, cov_parm(p1,p1), 
     .                          1.d0/qapr_sig(p2)**2
 720               format('Removing CoM ',i1,' NP ',2i4,
     .                 'Val and con ',2e23.12)
                   if( cov_parm(p1,p1).gt.1.d0/qapr_sig(p2)**2 )
     .             cov_parm(p1,p1)=cov_parm(p1,p1)-1.d0/qapr_sig(p2)**2
* MOD TAH 200805: Don't apply any additional constraint (not needed)
*    .                    + 0.01d0   ! 10 meter constraint 
               end if
            end do
         end if

****     See for axis offsets
         if( index(decon_str,'A').gt. 0 ) then
            do i = 1,qnum_svs
            do k = 1,3
               p1 = qparn_svs(max_svs_elem-3+k,i)
               if( p1.gt.0 ) then
                  do j = 1, np
                     if( atoo(j).eq.p1 ) p2 = j
                  end do
               end if
               if ( p1.gt.0 .and. p2.gt.0 ) then
                   write(*,740) i, k, p1,p2, cov_parm(p1,p1), 
     .                          1.d0/qapr_sig(p2)**2
 740               format('Removing Sat PC Pos ',i3,i4,' NP ',2i5,
     .                 ' Val and con ',2e15.6)
                   if( cov_parm(p1,p1).gt.1.d0/qapr_sig(p2)**2 )
     .             cov_parm(p1,p1)=cov_parm(p1,p1)-1.d0/qapr_sig(p2)**2
* MOD TAH 200805: Don't apply additonal constraint (biases resolts)
*    .                 + 0.01d0   ! 10 meter constraint 
               end if
            end do
            end do
         end if
****     Now convert the "normal equations" back to covariance matrix.
         call invert_vis(cov_parm, sol_parm, scale, ipivot, np, np, 1)

         do i = 1, np
            qsol(i) = sol_parm(i) + qapr(i)
            if( cov_parm(i,i).lt.0 ) write(*,780) i, cov_parm(i,i), 
     .            sol_parm(i)
 780        format('C+A decon ',i4,' Neg. Var ',2d20.10)
         end do
   
      end if
 
****  Thats all
      return
      end

CTITLE DEC_SOLN_VEC
 
      subroutine dec_soln_vec(unit,line, np, cov_parm, sol_parm)
 
      implicit none

*     Routine to decode the solution vector from a normal equations
*     system.   The values are read into sol_parm and later inverted
*     when all information is available.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
*   line        - Line read from file file
 
      character*(*) line
 
* LOCAL VARIABLES
 
*   ierr        - IOSTAT error
*   sn      - Site number from list of sites
*   jp      - Parameter number read from the input list
*   cyr, cdoy, csec - Reference yr, doy and sec
*   indx        - Position in string
*   nr          - "real" parameter number in final covariance matrix
*   occ         - Occupation number
*   cons_type   - Constraint type on station (0,1 or 2)
 
      integer*4 ierr, sn, jp, cyr, cdoy, csec, indx, nr, occ, i,
     .          cons_type, jerr, ne, kerr, svn
 
*   est     - Estimate from file
*   est_sig - Sigma of estimate.  Needed when correlations are given
*   ref_jd  - JD of reference epoch for value
 
      real*8 est, ref_jd
      
* dut1, dlod -- Short period corrections to UT1 amd LOD

      real*8 dut1, dlod 
      
*   end_block   - Indicates end of block found
 
      logical end_block
 
*   code    - Site code read from Sinex file
*   stype   - Parameter type read from sinex file
*   estu        - Units for the estimate
*   pt      - Point characters
*   full_name  -  Full site name with occ and pt added
*   chocc   - Character version of occupancy (version 1.0)
 
      character*8 code, stype, estu, full_name
      character*4 pt, chocc

* MOD TAH 210505: Added support for multi-frequency satellite
*     antenna offsets.
      integer*4 svant_off  ! Functon to return index in qparn_svs
                  ! array for antenna offset (L1/LC are in 
                  ! max_svs_elem - (2:0); Other are put in IC and
                  ! radiation parameter slots.
      integer*4 ao   ! Offset value returned from svant_off based
                  ! on pt string. 

*   found  - Set true if parameter correctly idenified
      logical found
 
****  Start reading the block
      end_block = .false.
      ierr = 0
      nr = 0
      do i = 1, max_glb_parn
        itoo(i) = 0
        qsig(i) = 0.d0
      end do
      
      do while ( ierr.eq.0 .and. .not.end_block)
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 .and. line(1:1).eq.' ') then
              
              if( cglb_vers.eq.5 ) then
                  read(line,210,iostat=jerr) 
     .                 jp, stype, code,pt,occ, cyr, cdoy, csec,
     .                 estu, cons_type, est
 210              format(i6,1x,a4,1x,a4,1x,a2,1x,i4,1x,i2,1x,i3,1x,i5,
     .                  1x,a4,i2,d21.14,1x,e14.8)
              else
                  read(line,215,iostat=jerr) 
     .                 jp, stype, code,pt,chocc, cyr, cdoy, csec,
     .                 estu, cons_type, est
 215              format(i6,1x,a6,1x,a4,1x,a2,1x,a4,1x,i2,1x,i3,1x,i5,
     .                  1x,a4,i2,d22.14,1x,e11.5)
     
*                 Check the occupancy value
                  if( chocc.eq.'----' ) chocc = '   1'
                  read(chocc,*,iostat=kerr) occ
                  if( kerr.ne.0 ) occ = 1                  
              end if
              call report_error('IOSTAT',jerr,'decod',line,0,
     .                          'dec_soln_apr')
 
              call yds_to_jd( cyr, cdoy, csec, ref_jd)

* MOD TAH 970716: Check for blank pt names, replace with ' A'
              if( pt(2:2).eq.' ' .or. pt(2:2).eq.'-' ) pt = ' A'

****          OK, construct the full name for the site and
*             see if we find it.

              write(full_name, 240) code, occ, pt(2:2)
 240          format(a4,i3.3,a1)
              call casefold(full_name)
              indx = 1
              call get_cmd(full_name, qsite_names, qnum_sites,
     .                         sn, indx )

              found = .false. 
              if( stype.eq.'STAX' .and. sn.le.0 ) then
                  write(*,250) full_name(1:4), full_name(5:7),
     .                         full_name(8:8)
 250              format('**WARNING** No site/id entry for code ',
     .                   a4, ' Occ ',a3,' Point ',a1)
              end if

****          See if satellite axis offset
              if( stype(1:4).eq.'SATA' ) then
                  read(code,'(1x,i3)',iostat=kerr) svn
                  sn = 0
                  do i = 1, qnum_svs
                     if( qsvi_svn(i).eq.svn ) sn = i
                  end do
                  if( kerr.ne.0 .or.sn.eq.0 ) then
                     write(*,260) line(1:80)
 260                 format('***WARNING*** Unable to find SV in: ',a)
                  endif

* MOD TAH 210505: Extract out frequency for antenna offsets to handle
*                 TUG L1, L2, L5, L7,L8 values.  Save as offset to 
*                 max_svs_elem (see SINEX has not IC or Rad parameters
*                 use these slots for non-L1, LC values.
                  ao = svant_off(pt)
              endif
                                      
              if( stype.eq.'STAX' .and. sn.gt.0 ) then
                  nr = nr + 1
                  qparn_sites(1,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'STAY' .and. sn.gt.0) then
                  nr = nr + 1
                  qparn_sites(2,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'STAZ' .and. sn.gt.0) then
                  nr = nr + 1
                  qparn_sites(3,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'VELX' .and. sn.gt.0) then
                  nr = nr + 1
                  qparn_vel(1,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'VELY' .and. sn.gt.0) then
                  nr = nr + 1
                  qparn_vel(2,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'VELZ' .and. sn.gt.0) then
                  nr = nr + 1
                  qparn_vel(3,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.

* MOD TAH 050710: Addin satellite axis offsets
              else if( stype.eq.'SATA_X' .and. sn.gt.0) then
                  nr = nr + 1
C MOD TAH 210505: Changed access for multi-frequency antenna
C                 offset (TUG SINEX)
C                 qparn_svs(max_svs_elem-2,sn) = nr
                  qparn_svs(ao+1,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'SATA_Y' .and. sn.gt.0 ) then
                  nr = nr + 1
C MOD TAH 210505: Changed access for multi-frequency antenna
C                 offset (TUG SINEX)
                  qparn_svs(ao+2,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.
              else if( stype.eq.'SATA_Z' .and. sn.gt.0 ) then
                  nr = nr + 1
C MOD TAH 210505: Changed access for multi-frequency antenna
C                 offset (TUG SINEX)
                  qparn_svs(ao+3,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
                  found = .true.

*****         Check for Earth orientation parameters
              else if( stype.eq.'LOD ' .or. stype.eq.'LODR' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.  Convert LOD to UT1 Rate
*                 (mas per day is the unit)
                  found = .true.
                  nr = nr + 1
c                  qnum_mul_pmu(2,3) = qnum_mul_pmu(2,3) + 1
c                  ne = qnum_mul_pmu(2,3)
                  qscale(nr) = -15.d0
*                 Special mod for JPL: MOD REMOVED TAH 030214.
c                  if( est.lt.0 ) then
c                      est = -est
c                      qscale(nr) = -15.d0
c                  end if
                                    
c                  qparn_mul_pmu(2,3,ne) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr

              else if( stype.eq.'UT  ' .or. stype.eq.'UT1 ' .or.
     .                 stype.eq.'UTR ' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  found = .true.
                 
                  nr = nr + 1
c                  qnum_mul_pmu(1,3) = qnum_mul_pmu(1,3) + 1
c                  ne = qnum_mul_pmu(1,3)
                  qscale(nr) = 15.d0

* MOD TAH 991207: Check GFZ values.  Sign seems to be wrong
                  if( qowner(1:3).eq. 'GFZ' .and. est.gt. 1000.d0 ) then
                      write(*,*) 'Changing SIGN of GFZ UT values'
                      est = -est
                  end if
                                    
c                  qparn_mul_pmu(1,3,ne) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'XPO ' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  found = .true.
                  nr = nr + 1
c                  qnum_mul_pmu(1,1) = qnum_mul_pmu(1,1) + 1
c                  ne = qnum_mul_pmu(1,1)
                  qscale(nr) = 1.d0
c                  qparn_mul_pmu(1,1,ne) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'XPOR' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.  Sinex Rate is
                  found = .true.
                  nr = nr + 1
c                  qnum_mul_pmu(2,1) = qnum_mul_pmu(2,1) + 1
c                  ne = qnum_mul_pmu(2,1)
                  qscale(nr) = 1.d0     
c                  qparn_mul_pmu(2,1,ne) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'YPO ' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  found = .true.
                  nr = nr + 1
c                  qnum_mul_pmu(1,2) = qnum_mul_pmu(1,2) + 1
c                  ne = qnum_mul_pmu(1,2)
                  qscale(nr) = 1.d0
c                  qparn_mul_pmu(1,2,ne) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr

              else if( stype.eq.'YPOR' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  found = .true.
                  nr = nr + 1
c                  qnum_mul_pmu(2,2) = qnum_mul_pmu(2,2) + 1
c                  ne = qnum_mul_pmu(2,2)
                  qscale(nr) = 1.d0
c                  qparn_mul_pmu(2,2,ne) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if((stype.eq.'TRANX ' .or.
     .                 stype.eq.'XGC   ').and. chocc.eq.'   1' ) then
                  nr = nr + 1
                  qparn_tran(1,1) = nr
                  qscale(nr) = 1.d0
                  if( stype.eq.'XGC   ' ) qscale(nr) = -1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if((stype.eq.'TRANY ' .or.
     .                 stype.eq.'YGC   ').and. chocc.eq.'   1' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_tran(2,1) = nr
                  qscale(nr) = 1.d0
                  if( stype.eq.'YGC   ' ) qscale(nr) = -1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if((stype.eq.'TRANZ ' .or.
     .                 stype.eq.'ZGC   ').and. chocc.eq.'   1'  ) then
                  nr = nr + 1 
                  qscale(nr) = 1.d0
                  qparn_tran(3,1) = nr
                  qscale(nr) = 1.d0
                  if( stype.eq.'ZGC   ' ) qscale(nr) = -1.d0
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANRX' .or.
     .                 stype.eq.'XGCR  ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  if( stype.eq.'XGCR  ' ) qscale(nr) = -1.d0
                  qparn_tran(1,2) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANRY' .or.
     .                 stype.eq.'YGCR  ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  if( stype.eq.'YGCR  ' ) qscale(nr) = -1.d0
                  qparn_tran(2,2) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANRZ' .or.
     .                 stype.eq.'ZGCR  ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  if( stype.eq.'ZGCR  ' ) qscale(nr) = -1.d0
                  qparn_tran(3,2) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'SCALE ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_scale(1) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
               else if( stype.eq.'SCALER' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_scale(2) = nr
                  qsol(nr) = est
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
       
*                 Thats all we know about for the moment
              end if
          end if
          if( line(1:1).eq.'-' ) then
              end_block = .true.
          end if
      end do

****  Write out status
      write(*,400) nr, np
 400  format(' Found ',i4,' usable parameters, expected ',i4,
     .       ' parameters in dec_soln_vec')
c      do i = 1,nr
c         write(*,420) qsol(i)
c 420     format('VEC ',E20.12)
c      enddo
 
****  Thats all
      return
      end 


CTITLE DEC_SOLN_NEQ
 
      subroutine dec_soln_neq(unit, line, np, cov_parm, sol_parm )
 
      implicit none

*     Routine to decode the normal equations of the estimates from
*     the sinex files.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
* PASSED VARIABLES
 
*   unit        - Unit number
*   np      - number of parameters expected
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix/Normal equations
* sol_parm(np)    - Solution vector
 
 
      real*8 cov_parm(np,np), sol_parm(np)
 
*   line        - Line read from file file
 
      character*(*) line
 
* LOCAL VARAIBLES
 
 
*   ierr        - IOSTAT error
*   np1, np2    - Parameter numbers
*   nr1, nr2    - Parameter numbers in the output solution
*   i,j     - Loop counters
*   lim     - Number of elements to copy with each read
*   max_np  - Largest parameter number expected when reading
*             upper triangular forms.  (Found from itoo array).
 
      integer*4 ierr, np1, np2, i,j, k, lim, nr1, nr2, max_np, jerr
 
*   cov(3)      - Up to thress values read from file.
 
      real*8 cov(3)
 
*   end_block   - Indicates end of block found
*   corr_mat    - Set true if this is a correlation matrix rather than
*                 covariance matrix
*   lower       - Set true to indicate lower triangular form
*   noreport    - Set false once we have reported sigmas on diagonal.
 
 
      logical end_block, corr_mat, lower, noreport
 
****  Start reading the block
      end_block = .false.
      ierr = 0
      max_np = 0

*     See if correlation or covariance matrix
      lower    = .true. 
* MOD TAH 210503: Changed test for Upper from 'U I' to 'X U'
*+SOLUTION/NORMAL_EQUATION_MATRIX U
C     if( index(line,'U I').gt.0  ) lower = .false.
      if( index(line,'X U').gt.0  ) lower = .false.

*     IF not a lower form then see maxparam expected
      if( .not.lower ) then
          max_np = 0
          i = max_glb_parn
          do while ( i.gt.1 .and. max_np.eq.0 )
              i = i - 1
              if( itoo(i).ne.0 ) max_np = i
          end do
      end if
 
      noreport = .true.

      do while ( ierr.eq.0 .and. .not.end_block)
 
          read(unit,'(a)', iostat=ierr) line
          if( ierr.eq.0 ) then
 
              if( line(1:1).eq.' ' ) then
* MOD TAH 960903: Cleared the cov array incase there are
*                 blanks in line which will retain there
*                 old value.             
                  do j = 1,3
                      cov(j) = 0.d0
                  end do
                  read(line, * ,iostat=jerr) np1, np2, (cov(j), j=1,3)
cc 210              format(2i6,3d22.13)

*****             Allow for -1 errors on read with free format (to support
*                 GSZ files)
                  if( ierr.eq.-1 ) ierr = 0
                  if( ierr.ne.0 .and. ierr.ne. -1 ) then
                      write(*,240) ierr, line      
 240                  format('IOSTAT error ',i5,' Reading COVA block',/,
     .                       'Line:',a)
                      stop 'Error reading COVA'
                  end if
                  
                  if ( lower ) then
                      lim = np1 - np2 + 1
                  else

* MOD TAH 960822: limit on max parameter number based on max
*   in input covare matrix.
                      lim = max_np - np2 + 1
                  end if
                  if( lim.gt.3 ) lim = 3

                  do j = 1, lim
*                    See if we need to convert to covariance
                     nr1 = itoo(np1)
                     nr2 = itoo(np2+j-1)
c                     write(*,710) np1, np2, nr1,nr2, cov(j)
 710                 format('N ',4i5,3E16.6)
                     if( nr1.gt.0 .and. nr2.gt.0 ) then
                          cov_parm(nr1, nr2) = cov(j)
                          cov_parm(nr2, nr1) = cov(j)
                     end if
                  end do
              else
                  if( line(1:1).eq.'-' ) then
                      end_block = .true.
                  end if
              end if
          end if
      end do

c      do i = 1,np
c         write(*,720) (cov_parm(i,j),j=1,np)
c 720     format('NEQ ',400E20.12)
c      enddo
 
****  Now copy the matrix into the upper triangular section
c      do i = 1, np
c         write(*,999) i, qsig(i)
c 999     format(i4,e22.10)
c      end do
      
      do i = 1, np - 1
          do j = i+1, np
              cov_parm(i,j) = cov_parm(j,i)
          end do
      end do

      if( 1.eq.2 ) then 
         write(*,997) (i,itoo(i),i=1,np)
 997     format(550(i3,1x,i3,' | '))
         do i = 1, np
            write(*,998 ) i,(cov_parm(i,j),j=1,np)
 998        format('C ',i3,500(1x,G13.5))
         end do
      end if

****  Now check the block diagonal in the covariance to make sure
*     no negative elements on diagonal and not correclations greater
*     than +-1.
      call check_cov_htog(cov_parm, np )


****  Thats all
      return
      end

CTITLE ADD_NEQ_CON

      subroutine add_neq_con( cov_parm, np) 

      implicit none

*     Routine to add constraints to normal eqautions

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
  
* PASSED VARIABLES
 
*   np      - number of parameters expected
 
      integer*4 np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector

 
      real*8 cov_parm(np,np)

* MOD TAH 141119: Add scale effect for COD (NEQ entries are too large
*     so apply appropriate scale factor

      real*8 S   ! Scale factor for COD (1.d4), 1.00 for other
* MOD TAH 210503: Added AC dependent "loose" constraint due to different
*     sizes of elelements in NEQ systems.  Coded for COD, ESA and TUG.
      real*8 neq_diag  ! Value to be added to NEQ diag as weak constraint.
      real*8 MeanD(2)  ! Mean NEQ diagonal for Postion and EOP
      integer*4 MeanN(2)   ! Number of values
 
* LOCAL VARIABLES

      integer*4 i,j, nq, k, l, p1, p2

****  Add station position constraints
      S = 1.d0
      MeanD = 0.0
      neq_diag = 1.d-2   ! Value for COD and EDA
      if( qowner(1:3).eq.'COD' ) then
          S = 1.0d-4
      elseif ( qowner(1:3).eq.'TUG' ) then
          neq_diag = 1.d0 
      endif 
      do i = 1, qnum_sites
         do j = 1, 3
            nq = qparn_sites(j,i)
            if( nq.gt.0 ) then
                MeanN(1) = MeanN(1) + 1
                MeanD(1) = MeanD(1) + cov_parm(nq,nq)
                cov_parm(nq,nq) = cov_parm(nq,nq) + S*neq_diag
            end if
         end do
      end do

****  Now do EOP: Just constrain UT1 J = 3; Constraint allows for
*     over all rotation.
C     j = 3
C     do k = 1, qnum_mul_pmu(1,j)   ! Number of UT1 values
C        p1 = qparn_mul_pmu(1,j,k)
C        do l = 1, qnum_mul_pmu(1,j)
C           p2 = qparn_mul_pmu(1,j,l)
C           if( p1.gt.0 .and. p2.gt. 0 ) then
C               MeanN(2) = MeanN(2) + 1
C               MeanD(2) = MeanD(2) + cov_parm(nq,nq)
C               cov_parm(p1,p2) = cov_parm(p1,p2) + S*1.d-2
C           end if 
C        end do
C     end do
* MOD TAH 200728: Just add loose constraint to diagonal for EOP
      do i = 1,2   ! Offset and rate terms
         do j = 1,3   ! X,Y, UT1
           do k = 1, qnum_mul_pmu(i,j)  ! Number of terms
               nq = qparn_mul_pmu(i,j,k)
               MeanN(2) = MeanN(2) + 1
               MeanD(2) = MeanD(2) + cov_parm(nq,nq)
               cov_parm(nq,nq) = cov_parm(nq,nq) + S*neq_diag
           enddo
         enddo
      enddo

      write(*,220) qowner(1:3), S, MeanD/MeanN
 220  format('NEQ Constraint ',a,' Scale ',G8.1,' Mean Diag ',
     .       2E12.3)

****  Now add constraints to antenna offsets
C     if( index(decon_str,'A').eq.0 ) then ! No option to 
*             deconstrain antenna offsets, so constrain
*             to 0.1 mm
         do i = 1, qnum_svs
* MOD TAH 210503: Scan accross all elements to handle multiple
*           frequency savtellite antenna offsets in TUG solutions
C           do j = 0,2
            do j = 1, max_svs_elem 
C              nq = qparn_svs(max_svs_elem-j,i)
               nq = qparn_svs(j,i)
               if( nq.gt.0 ) then
* MOD TAH 210503: See if A option in decon_str. Use 'M' for 
*                 multi-frequency PCO values
                  if( index(decon_str,'A').gt.0 .or.
     .                index(decon_str,'M').gt.0  ) then 
*                     Place weak constraint. (or no constraint)
                     cov_parm(nq,nq) = cov_parm(nq,nq) + 0*neq_diag
                  else
*                    Satellite PCO should be constrained
                     cov_parm(nq,nq) = cov_parm(nq,nq) + 1.d10*neq_diag
                  endif 
C                 if( index(decon_str,'A').eq.0 ) then
*                    Satellite PCO should be constrained
C                    cov_parm(nq,nq) = cov_parm(nq,nq) + 1.d10*neq_diag
C                 else   ! Place weak constraint.
C                    cov_parm(nq,nq) = cov_parm(nq,nq) + neq_diag
C                 endif
               endif
            enddo
         enddo
C     endif

****  Thats all
      return
      end

CTTITE SVANT_OFF

      integer*4 function svant_off( pt )
 
      implicit none

*     Function 210505 to extract out frequency for antenna offsets to handle
*     TUG L1, L2, L5, L7,L8 values.  Save as offset to 
*     max_svs_elem (see SINEX has not IC or Rad parameters
*     use these slots for non-L1, LC values.
*     Slot specifically to fit in radiation parameter:
*     Allow L2, L5, L7 and L8 onl7

      include '../includes/kalman_param.h'

* VARIABLES
      character*(*) pt   ! PT code from SATA line (LC, L1, L2, L5-8)
C     integer*4 svant_off ! Offset for antenna offset parameters at
                          ! different frequencies.  Only LC anad L1
                          ! values are passed to binary h-file.

*     Test of pt code to see what is presemt.

      if ( pt(1:2).eq.'LC' ) then
         svant_off = max_svs_elem-3  ! Save with +1,+2,+3
      elseif ( pt(1:2).eq.'L1' ) then ! Same as LC
         svant_off = max_svs_elem-3  ! Save with +1,+2,+3
      elseif ( pt(1:2).eq.'L2' ) then
         svant_off = 8  ! Map to 1,2,3
      elseif(  pt(1:2).eq.'L5' ) then
         svant_off = 11 
      elseif(  pt(1:2).eq.'L7' ) then
         svant_off = 14 
      elseif(  pt(1:2).eq.'L8' ) then
         svant_off = 17
      else
*        Generate fatal error if non-numeric
         call report_error('IOSTAT',-1,'unknow band',pt,1,'SVANT_OFF')
         svant_off = 0
      endif

      RETURN
      end


