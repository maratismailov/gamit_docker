CTITLE DECSNX_SOLN
 
      subroutine decsnx_soln(unit, line, np, cov_parm, sol_parm)
 
*     Routine to decode the solution  blocks from SINEX:
 
      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'
 
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

      logical estread, aprread, epread

* estread, aprread -- Set true once the matrix est and apriori
*     have been read.  We then can remove the constraints.
      data  estread / .false. /, aprread / .false. /, epread / .false. /
      
* LOCAL VARIABLES
 
*   indx        - Position in string
*   ierr        - IOSTAT error
*   trimlen     - Length of string
 
      integer*4 indx, trimlen, i
     
 
*       block_found - True if we know how to decode
*           - the block
 
 
      logical block_found
 
****  Start looking for blocks
      block_found = .false.
      estread = .false.
      aprread = .false.
      epread  = .false.
 
      indx = index(line,'SOLUTION/EPOCHS')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_soln_epochs(unit,line)
          write(*,'(a)') 'Done doing SOLUTION/EPOCHS'
      end if
      
      indx = index(line,'SOLUTION/STATISTICS')
      if( indx.gt.0 ) then
          block_found = .true.
          call dec_soln_stats(unit,line)
          write(*,'(a)') 'Done doing SOLUTION/STATISTICS'
      end if
 
      indx = index(line,'SOLUTION/ESTIMATE')
      if( indx.gt.0 ) then 
          block_found = .true.
          if ( .not. epread ) then
             call dec_soln_est(unit,line, np, cov_parm, sol_parm)
             do i = 1, np
                sol_parm(i) = qsol(i)
             end do
             epread = .true.
          else
             call snx_finbl(unit)
          endif
          write(*,'(a)') 'Done doing SOLUTION/ESTIMATE'
  
      end if
 
      indx = index(line,'SOLUTION/APRIORI')
      if( indx.gt.0 ) then 
          block_found = .true.
          if( epread ) then
             call dec_soln_apr(unit,line, np, cov_parm, sol_parm)
             apr_missed = .false.
          else
              call snx_finbl(unit)
          endif
          write(*,'(a)') 'Done doing SOLUTION/APRIORI'
          
      end if
 
      indx = index(line,'SOLUTION/MATRIX_ESTIMATE')
      if( indx.gt.0 ) then
          block_found = .true.
          estread = .true.
          cov_parm(1,1) = 1
          call dec_soln_mat(unit,line, np, cov_parm, sol_parm)
          write(*,'(a)') 'Done doing SOLUTION/MATRIX_ESTIMATE'
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
          write(*,'(a)') 'Done doing SOLUTION/MATRIX_APRIORI'
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
 
*     Routine to decode the solution epochs to get the start and
*     stop times for each station
 
      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'
 
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
              call get_command(full_name, qsite_names, qnum_sites,
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
 
*     Routine to decode the solution stats block to get number of
*     data and chiqsq 
 
      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'
 
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
 210              format(32x,i22)
                  cnum_obs = snum_obs
              else if( index(line,'VARIANCE FACTOR').gt.0 ) then
                  read(line, 220, iostat=jerr) cchisq
 220               format(32x,f22.15)
              endif
              
          end if
          if( line(1:1).eq.'-' ) then
              end_block = .true.
              write(*,240) cnum_obs, cchisq
 240          format(' Analysis with ',i8,' Obs, Chi**2 ',f22.6)
          end if
      end do
 
****  Thats all
      return
      end

CTITLE DEC_SOLN_EST
 
      subroutine dec_soln_est(unit,line, np, cov_parm, sol_parm)
 
*     Routine to decode the solution estimates part of the sinex
*     files.  These values are written directly into the sol_parm
*     vector.
 
      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'
 
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
     .          cons_type, jerr, ne
 
*   est     - Estimate from file
*   est_sig - Sigma of estimate.  Needed when correlations are given
*   ref_jd  - JD of reference epoch for value
 
      real*8 est, est_sig, ref_jd
      
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
     .                  1x,a4,i2,d21.14,1x,d14.8)
              else
                  read(line,215,iostat=jerr) 
     .                 jp, stype, code,pt,chocc, cyr, cdoy, csec,
     .                 estu, cons_type, est, est_sig
 215              format(i6,1x,a6,1x,a4,1x,a2,1x,a4,1x,i2,1x,i3,1x,i5,
     .                  1x,a4,i2,d22.14,1x,d11.5)
     
*                 Check the occupancy value
                  if( chocc.eq.'----' ) chocc = '   1'
                  read(chocc,*,iostat=ierr) occ                  
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
              call get_command(full_name, qsite_names, qnum_sites,
     .                         sn, indx )

              found = .false. 
              if( stype.eq.'STAX' .and. sn.le.0 ) then
                  write(*,250) full_name(1:4), full_name(5:7),
     .                         full_name(8:8)
 250              format('**WARNING** No site/id entry for code ',
     .                   a4, ' Occ ',a3,' Point ',a1)
              end if                 
              if( stype.eq.'STAX' .and. sn.gt.0 ) then
                  nr = nr + 1
                  qparn_sites(1,sn) = nr
                  qscale(nr) = 1.d0
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
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
*                 Special mod for JPL:
                  if( est.lt.0 ) then
                      est = -est
                      qscale(nr) = -15.d0
                  end if
                  
*                 See if we need to restore short-period UT1 corrections.
*                 Need to CHANGE JPL Date when problem fixed.
                  if ( stype.eq.'LODR' .or. 
     .                (qowner(1:3).eq. 'JPL' .and. 
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
                  itoo(jp) = nr

              else if( stype.eq.'UT  ' .or.stype.eq.'UTR ' ) then

*                 Check the epoch.  If it matches the solution epoch then
*                 save, otherwise ignore it.
                  found = .true.
                 
                  nr = nr + 1
                  qnum_mul_pmu(1,3) = qnum_mul_pmu(1,3) + 1
                  ne = qnum_mul_pmu(1,3)
                  qscale(nr) = 15.d0

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
                  
                  qparn_mul_pmu(1,3,ne) = nr
                  qsol(nr) = est
                  qpmu_mul_apr(1,3,ne) = est*15.d0
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
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
                  itoo(jp) = nr
              else if( stype.eq.'TRANX ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_tran(1,2) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANY ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_tran(2,2) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANZ ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_tran(3,1) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANRX' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_tran(1,2) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANRY' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
                  qparn_tran(2,2) = nr
                  qsol(nr) = est
                  qsig(nr) = est_sig
                  qref_ep(nr) = ref_jd
                  itoo(jp) = nr
              else if( stype.eq.'TRANRZ' ) then
                  nr = nr + 1
                  qscale(nr) = 1.d0
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
     .       ' parameters')
 
****  Thats all
      return
      end
 
 
CTITLE DEC_SOLN_APR
 
      subroutine dec_soln_apr(unit,line, np, cov_parm, sol_parm)
 
*     Routine to decode the solution apriori part of the sinex
*     files.  These values are written directly into the site_pos
*     vector.
 
      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'
 
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
     .                   1x,a4,2x,d21.14,1x,d14.8)
              else
                  read(line,215,iostat=jerr) jp, stype, code,pt,chocc, 
     .                 cyr, cdoy, csec, estu, est, est_sig
 215              format(i6,1x,a6,1x,a4,1x,a2,1x,a4,1x,i2,1x,i3,1x,i5,
     .                   1x,a4,2x,d22.14,1x,d11.5)
     
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
              call get_command(full_name, qsite_names, qnum_sites,
     .                        sn, indx )
  
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
 
*     Routine to decode the covariance matrix of the estimates from
*     the sinex files.
 
      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'
 
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
 
      integer*4 ierr, np1, np2, i,j, lim, nr1, nr2, max_np, jerr
 
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
 210              format(2i6,3d21.13)

*****             Allow for -1 errors on read with free format (to support
*                 GSZ files)
                  if( ierr.eq.-1 ) ierr = 0
                  if( ierr.ne.0 .and. ierr.ne. -1 ) then
                      write(*,240) ierr, line      
 240                  format('IOSTAT error ',i5,' Reading COVA block',/,
     .                       'Line:',a)
                      stop 'Error reading COVA'
                  end if

*                 Now assign the elements to the matrix
                  if( corr_mat ) then
                      lim = np1 - np2
                      
* MOD TAH 960614: Check added for Ver 1.0 sinex.
*                     See if sigma is on the diagonal
                      if( lim.lt.3 .and. cov(lim+1).gt.0 .and.
     .                    cglb_vers.gt.5  ) then
                          nr1 = itoo(np1)
*                         Make sure this parameter still being
*                         used.                          
                          if( nr1.gt.0 ) then
                             qsig(nr1) = cov(lim+1)
                             cov_parm(nr1,nr1) = qsig(nr1)**2
                          end if
                       end if                     
                  else
                      if ( lower ) then
                          lim = np1 - np2 + 1
                      else

* MOD TAH 960822: Set limit on max parameter number based on max
*   in input covaraince matrix.
                          lim = max_np - np2 + 1
                      end if
                  end if
                  if( lim.gt.3 ) lim = 3

                  do j = 1, lim
*                    See if we need to convert to covariance
                     nr1 = itoo(np1)
                     nr2 = itoo(np2+j-1)
                     if( corr_mat .and. nr1.gt.0 .and. nr2.gt.0 ) then

* MOD TAH 981006:  See if diagonal element.
                          if( nr1.eq.nr2 ) then
                              if( cov(j).ne.1.d0 ) then 
*                                 see if variance or sigma
                                  if( abs(cov(j)/qsig(nr1)-1.d0).
     .                                             lt.1.d-6 ) then
                                      qsig(nr1) = cov(j)
                                      cov(j) = 1.d0
                                      if ( noreport ) then
                                          write(*,*) 
     .                                    '***Sigmas found on diagonal'
                                          noreport = .false.
                                      endif
                                  end if
                              end if
                          end if
                          cov_parm(nr1, nr2) = cov(j)*qsig(nr1)*
     .                                                qsig(nr2)
                     end if
                     if( .not.corr_mat .and. 
     .                   nr1.gt.0 .and. nr2.gt.0 ) then
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
             
****  Thats all
      return
      end
 
CTITLE DEC_SOLN_CON
 
      subroutine dec_soln_con(unit, line, np, cov_parm, sol_parm )
 
*     Routine to decode the apriori contraint matrix.  If the 
*     contraints are tight then they are removed from the solution.
 
      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'
 
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
     .          l1,l2, sn, fsp, jerr
 
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
 210              format(2i6,3d21.13)
*                 Allow for free format errors
                  if( ierr.eq.-1 ) ierr = 0
                  if( ierr.ne.0 ) then
                      write(*,240) ierr, line      
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
 
*     Routine to decode the apriori contraint matrix.  If the 
*     contraints are tight then they are removed from the solution.
 
      include '../includes/kalman_param.h'
      include '../includes/htoglb_comm.h'
 
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
     .       rot_var

*   min_eigval, max_eigval - Min and max Eigenvalues

      real*8 min_eigval, max_eigval, pmu_parts(3,3,max_glb_sites)

*   num_nege  - Number of negative eigenvalues
      integer*4 num_nege

*
*     Check the eiegn values value of the original matrix
      if( (index(hfile,'emr').gt.0 .or.  index(hfile,'gfz').gt.0 .or.
     .     index(hfile,'ngs').gt.0 .or.  index(hfile,'cod').gt.0 .or.
     .     index(hfile,'six').gt.0 .or.  index(hfile,'esa').gt.0 .or.
     .     index(hfile,'jpx').gt.0)  .and.
     .    snx_constraint.eq.0               ) then 
         rot_var = 10.0d0**2
         write(*,*) 'Applying 10 mas Rotational Deconstraint '
         do i = 1, qnum_sites
*                     ! XYZ
             do j = 1,3
                 call pmu_main_part(i, pmu_parts(1,j,i), site_pos(1,i),
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
     .                       pmu_parts(m,j,i)*pmu_parts(m,l,k)*
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
         do j = 1, 2                 ! Loop over XY pole
            do k = 1, qnum_mul_pmu(1,j)
               p1 = qparn_mul_pmu(1,j,k)
               do i = 1, qnum_sites  ! Loop over Sites
                  do l = 1,3         ! Loop over XYZ at site
                     p2 = qparn_sites(l,i)
                     cov_parm(p1,p2) = cov_parm(p1,p2) +
     .                       pmu_parts(j,l,i)*(rot_var) 
                     cov_parm(p2,p1) = cov_parm(p1,p2)
                  end do
               end do
            end do
         end do

      end if

* MOD TAH 990513: Add a little noise to the UT1 values in the esa
*    sinex files
      if( index(hfile,'esa').gt.0 ) then
          write(*,*) 'Adding noise to ESA pole and UT1'
          do j = 1, 3
             do i = 1, qnum_mul_pmu(1,j)
                p1 = qparn_mul_pmu(1,j,i)
                if( p1.gt.0 ) then
                   cov_parm(p1,p1) = cov_parm(p1,p1) + 1.d-6
                end if
             end do
          end do
      end if

      call dwmov(cov_parm, 1, cov_copy, 1, np*np )

* MOD TAH 991206: See if user wants to compute eigenvalues
      if( comp_eigen ) then
         call jacobi(cov_copy,np,np,eigval,eigvec,nrot,
     .               bwork,zwork)
c         call eigen( cov_copy,np, np, 1, eigvec, nrot, eigval, 
c    .               1.d-1, zwork)

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
      end if

      if( snx_constraint.ne.0 ) then
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
 510            format('Site ',a8,' Constraints ',3F15.8)

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

             call dwmov(cov_parm, 1, cov_copy, 1, np*np )

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
 
****  Thats all
      return
      end
 
