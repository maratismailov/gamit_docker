 
CTITLE MAKE_SNX
 
      subroutine make_snx( unit, line, np, cov_parm, sol_parm,
     .                     cov_copy, eigvec, eigval, bwork, zwork )

      implicit none
 
*     Routine to read and make a binary hfile for SINEX files:
*     Current version for Sinex 1.00
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
 
* unit   - Unit for reading the file
* np     - Number of parameters in this solution
 
 
      integer*4 unit, np
 
* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector

* cov_copy(np,np) - Copy of covariance matrix for making eigvlaues
* eigvec(np,np)   - Eignvectors
* eigval(np)      - Eigenvalues
* bwork(np), zwork(np) - Working arrays for eigenvectors.
 

* MOD TAH 200728: Add column to cov_parm because of remapping
*     of CODE 3-day SINEX files. cov_parm(:,np+1) is same
*     memory location as sol_parm(:) 
      real*8 cov_parm(np,np+1), sol_parm(np),
     .       cov_copy(np,np), eigvec(np,np), eigval(np), bwork(np),
     .       zwork(np)
 
* line   - Line read from file

      real*8 snx_meaneq  ! Mean epoch of CWU solutions (to be GIPSY bug)
 
 
      character*(*) line
 
* LOCAL VARIABLES
 
* i,j      - Loop counters
* num_in_file = number of solutions in file (1 only for snx solns)
* ierr     - IOSTAT error
* trimlen  - Length of string
 
 
      integer*4 i, num_in_file, ierr, trimlen, j, na, k

* name_unique - Set true is non _GPS name is the only one
*     and site will be named back

      logical name_unique   
* MOD TAH 200727.
      logical out  ! True to output multiple site names. Set
                   ! false for COD and sngday.
 
*  cyr, cdoy, csec - SNX creation start year, doy and ssecs
*  syr, sdoy, ssec - SNX data start year, doy and ssecs
*  eyr, edoy, esec - SNX data end year, doy and secs
 
 
      integer*4 syr, sdoy, ssec, eyr, edoy, esec, indx_q,
     .          cyr, cdoy, csec

      real*8 run_time, sectag, snx_vers

*  hname - String to hold the name of file for getting the first
*     part of name from

      character*64 hname
      character*4  sys_type
 
*
****  This is snx file:  Read the number of parameters
*     Read the start and stop times from the SNX line
* MOD TAH 960430: Add reading version.
      read(line,120, iostat=ierr) snx_vers, ccreator, cyr, cdoy, csec, 
     .               cowner, syr, sdoy, ssec,
     .               eyr, edoy, esec, sys_type, 
     .               ccons_type, canal_type
 120  format(5x,F5.2,1x, a3,1x,i2,1x,i3,1x,i5,1x,
     .       a3,1x,i2,1x,i3,1x,i5,1x,i2,1x,i3,1x,i5,1x,
     .       a1,1x,5x,1x,i1,1x,a5)
      
      call report_error('IOSTAT',ierr,'decod',line,0,'make_snx')

* MOD TAH 960530: Save the snx version in the globk style (ie. integer
*     multiplied by 100).
      cglb_vers = snx_vers*100
          
*     Set the data type for this solution.
* MOD TAH 970922: Set the system type for the first one.  Normally
*     the history entries will overwrite this.
      if( sys_type(1:1).eq.'P') qssys_type(1) = 2
      if( sys_type(1:1).eq.'R') qssys_type(1) = 1
      if( sys_type(1:1).eq.'L') qssys_type(1) = 4
      if( sys_type(1:1).eq.'C') qssys_type(1) = 7
      if( sys_type(1:1).eq.'M') qssys_type(1) = 8
*     This is combined and really don't know so we set all main
*     bits
      if( sys_type(1:1).eq.'C') csys_type = 7

      call yds_to_jd(cyr, cdoy, csec, run_time )
      call jd_to_ymdhms( run_time, qrun_time, sectag)
      qrun_time(6) = nint(sectag)
      qrun_time(7) = 0
 
*     Now convert values to JD's
*     Now convert values to JD's
* MOD TAH 161008: See if syr == 0 -- Bad JPL sinex files
      if( syr.gt.0 .and. sdoy.gt.0 ) then
          call yds_to_jd( syr, sdoy, ssec, qstart_epoch)
      else
          write(*,220) trim(hfile), sepoch
 220      format('*** WARNING *** Sinex file ',a,
     .           ' has no start date. Using JD ',F12.2)
          qstart_epoch = sepoch
      endif
 
* MOD TAH 000308: Added check on all components of date since
*      for 2000 the year will be zero.
      if( eyr.gt.0 .or. edoy.gt.0 .or.esec.gt.0 ) then
          call yds_to_jd( eyr, edoy, esec, qend_epoch)
      else
* MOD TAH 161008: Assume 1-day file
          qend_epoch = qstart_epoch + 1
      end if

* MOD TAH 140223: Fix the error in GIPSY start and stop times from CWU
*     PBO solutions.
      if ( index(hfile,'CWU').gt.0 .and. 
     .     (qend_epoch-qstart_epoch).gt.1.d0 ) then
*         Adjust start and stop to be 1-day.
          snx_meaneq = (qend_epoch+qstart_epoch)/2
          qend_epoch   = snx_meaneq + 0.5d0
          qstart_epoch = snx_meaneq - 0.5d0
*         Re-set sepoch (read in preread_snx).
          sepoch = snx_meaneq
      endif

* MOD TAH 200727: Get mean epoch COD solutions with 3-days
      if( ccreator(1:3).eq.'COD' .and. sngday ) then
*         Adjust start and stop to be 1-day.
          snx_meaneq = (qend_epoch+qstart_epoch)/2
          qend_epoch   = snx_meaneq + 0.5d0
          qstart_epoch = snx_meaneq - 0.5d0
          qmid_epoch   = snx_meaneq
*         Re-set sepoch (read in preread_snx).
          sepoch = snx_meaneq
      endif

*     NOTE: The solution epoch is found in preread_snx

      hfile_type = 'SNX'

*     Save the number of parameters . (This was done in preread_snx in htoglb)
      qnum_parn = np
* MOD TAH: 050910: Use the maximum number of satellites since we don't know
*     how many yet.
      call init_qparn( qnum_sites,  max_glb_svs )
      call init_site(site_pos, site_vel, qnum_sites)

      qsnum_save = 0
      estread = .false.
      aprread = .false.
      epread  = .false.
      neq_sys = .false.
      neq_read = .false.
      vec_read = .false.
      neq_inv  = .false.

      sgamit_mod = 0   ! For SINEX we don't know what these are.
      sload_mod = 0 
 
****  Now start looping over the blocks in the Sinex file
      apr_missed = .true.
      
      do while ( ierr.eq.0 )
 
          read(unit,'(a)', iostat=ierr ) line
          if( ierr.eq.0 ) then
 
*             File read OK, see what we have.  Look
*             for start of basic blocks
              if( line(1:1).eq.'+' ) then
 
*                 This start of block. Call the main block
*                 decoder
                  call decode_snxb(unit, line, np, cov_parm, sol_parm)
              else if( line(1:4).eq.'%END' ) then
 
*                 End of sinex file
                  ierr = -1
              end if
          end if
      end do

*     See if we need to go back and get the apriori blocks.
      if( apr_missed ) then
          write(*,*) 'Rewinding to get apriori blocks'
          na = 0
          rewind(unit)
          ierr = 0
          do while ( ierr.eq.0 )
 
              read(unit,'(a)', iostat=ierr ) line
              if( ierr.eq.0 ) then
 
*                 File read OK, see what we have.  Look
*                 for start of basic blocks
                  if( line(1:17).eq.'+SOLUTION/APRIORI' .or.
     .                line(1:24).eq.'+SOLUTION/MATRIX_APRIORI' ) then
                      na = na + 1
 
*                     This start of block. Call the main block
*                     decoder
                      call decode_snxb(unit, line, np, 
     .                                 cov_parm, sol_parm)
*                     See if we have read both apriori blocks.  If so
*                     stop reading file.     
                      if( na.eq. 2 ) ierr = -1
                  else if( line(1:4).eq.'%END' ) then
 
*                     End of sinex file
                      ierr = -1
                  end if
              end if
          end do
      end if
                      

****  Now remove the constraints
* MOD TAH 971022: Pass arrays needed for eigenvalues.
      if( .not.neq_sys ) then 
         call dec_remv_con(unit,line, np, cov_parm, sol_parm,
     .           cov_copy, eigvec, eigval, bwork, zwork )

*        QSOL is computed in dec_remv_con and has changes to sol_parm
*        accounted for. 
         do i = 1, np
            sol_parm(i) = qsol(i)
         end do
      else

* MOD TAH 200723: See if we are removing parameter.  The remove_neq
*        implicitly solve the NEQ for parmaters to be removed.  They
*        then removed physically with the pack_eq routine.  This code
*        will also remove any PCO constraints if -d=A option used.

         call remove_neq( np, cov_parm, sol_parm)

* MOD TAH 200728: Now compact the normal equations and remap
*        the parmeter lists.  After the pack, soln_parm can't
*        be used; use cov_parm(1,np+1) which is where it was
*        packed to. (RETURN only; not fully working yet).
*       (Scaling code below is still a problem becuase of change 
*        in cov_parm dimensions).
         call pack_neq( np, cov_parm, sol_parm)

*        Now inverr norm eqn to get convariance matrix.
C        call inv_normeq( np, cov_parm, sol_parm)
         call inv_normeq( np, cov_parm, cov_parm(1,np+1))
      endif

* MOD TAH 970321: Check the names of the stes and if we have
*     unique solutions numbers rename back to _GPS
* MOD TAH 200727: Set out to be true of multiple sites unless
*     this is COD and sngday is true (default)
      out = .true.
      if( qowner(1:3).eq.'COD' .and. sngday ) out = .false.

      do i = 1, qnum_sites
         if( qsite_names(i)(5:7).ne.'001' ) then

*        See if the first part of the name is uniquue
             name_unique = .true.
             do j = 1, qnum_sites
                if( j.ne.i ) then
C                   if( qsite_names(i)(1:8).eq.
C    .                  qsite_names(j)(1:8)  ) then
                    if( qsite_names(i)(1:4).eq.qsite_names(j)(1:4)
     .                                     .and.
     .                  qsite_names(i)(8:8).eq.qsite_names(j)(8:8) 
     .                  ) then
                        name_unique = .false.
* MOD TAH 200727: Only report if we are not removeing extra days.
                        if( out )
     .                  write(*,250) j, qsite_names(j), i,
     .                               qsite_names(i)  
 250                    format(' Multiple names: Site ',i4,2x,a8,
     .                         ' and site ',i4,
     .                          2x,a8,' are both needed')  
                     end if
                end if
             end do
 
*            If the name is unqiue then sumply rename back to
*            the standard name.
             if( name_unique .and. .not.igs_ptname ) then
                 qsite_names(i)(5:7) = '001' 
             end if
          end if
      end do
 
*     Write out some information
      write(*,300) qnum_sites, hfile(1:trimlen(hfile))
 300  format(i5,' sites found in ',a)
 
****  Now rescale the matrix (meter everywhere except solar radiation
*     which will be ratio to direct effect)
      do i = 1, qnum_parn
          call dwsmy(qscale(i), cov_parm(1,i),1,
     .                          cov_parm(1,i),1,  qnum_parn)
          call dwsmy(qscale(i), cov_parm(i,1),qnum_parn,
     .                          cov_parm(i,1),qnum_parn, qnum_parn)
 
* MOD TAH 200728: Due to possible remapping apply scale to
*        np+1 column of cov_parm (noramlly same memory location
*        as sol_parm 
C        sol_parm(i) = sol_parm(i)*qscale(i)
         cov_parm(i,np+1) = cov_parm(i,np+1)*qscale(i)
c          write(*,320) 'SCALE ',i,qscale(i), sol_parm(i)
c 320      format(a,1x,i5,1x,F19.6,1x,F19.6)
      end do
 
*     Generate the Fake KalObs file name
      sKalobs_file = hfile(1:trimlen(hfile))
 
*     Get Mfile name (claim this is database)
      call get_indx_q( hfile, indx_q)
      if( index(hfile,'_X').gt.0 ) then
          hname = 'V' // hfile(indx_q+3:indx_q+4) //
     .            hfile(indx_q-1:indx_q-1) // 'V.DBAS'
      else
          hname = hfile(indx_q-8:indx_q-6) // 
     .            hfile(indx_q-1:indx_q-1) // 'W.' // 
     .            hfile(indx_q-4:indx_q-1)
      end if 
      sdata_base = hname
 
      ssvs_epoch   = sepoch
 
****  Now generate the apriori and parameter codes
 
      call qgen_codes
 
***** Now write out the global file in GLOBK format.
      hf_ext = 'gls'

      call gen_gname( hname, sepoch, num_in_file, glb_dir,
     .                glb_file, hf_ext, name_format,expt_code)
      write(*,600) glb_file(1:trimlen(glb_file))
 600  Format(' Creating Binary file ',a)
 
***** Write out the header (create file if need be)
*     Update any information which might not have been in
*     soln information records.
      call check_soln_rec
      
      call qw_glb_header
 
*     Write out the names
      call qw_glb_names
 
*     Write out the full names
      call qw_glb_full
 
*     Write out the solution description
      call qw_description
 
*     Write out the codes
      call qw_codes
      
*     Write out the multi-parameter epoch records
      call qw_mul_ep
 
*     Write out the apriori values
      call qw_aprioris
 
*     Write out the solution and covariance matrix
      call qw_soln( cov_parm )
 
***** Thats all for this solution.  Now see if there is another
*     one
      return
      end

CTITLE CHECK_SOLN_REC

      subroutine check_soln_rec
      
      implicit none

*     Routine to copy any %=SNX line information to the sln
*     information records.
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

      integer*4 cn, orig_cn
      real*8 sectag, run_time

****  Copy over the information read from the %=SNX line
      orig_cn = cnum_soln_recs
      if( cnum_soln_recs.eq.0 ) cnum_soln_recs = 1
      cn = cnum_soln_recs
      qscons_type(cnum_soln_recs) = ccons_type
      qssys_type(cnum_soln_recs)  = csys_type
      qsowner(cnum_soln_recs)     = cowner
      qscreator(cnum_soln_recs)   = ccreator
      qsanal_type(cnum_soln_recs) = canal_type
      if( orig_cn.eq.0 ) then
          qsepoch_start(cn) = qstart_epoch
          qsepoch_end(cn)   = qend_epoch
          sectag = qrun_time(6)
          call ymdhms_to_jd( qrun_time, sectag, run_time)
          qsrun_time(cn)    = run_time   
          qsglb_ver(cn)     = cglb_vers
      end if
      
****  Thats all
      return
      end

CTITLE INV_NORMEQ

      subroutine inv_normeq(np, cov_parm, sol_parm)
 
      implicit none

*     Routine to invert normal equations to make covariance
*     matrix and solution 

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* PASSED
  
      integer*4 np   ! Number of parameters

* cov_parm(np,np) - Covarince matrix
* sol_parm(np)    - Solution vector
 
      real*8 cov_parm(np,np), sol_parm(np)
 
* LOCAL
* ipivot(max_glb_parn)  - Pivot elements for matrix inverse

      integer*4 i, j, ipivot(max_glb_parn)

* scale(max_glb_parn) -- Scale for matrix inversion

      real*8 scale(max_glb_parn)


****  If this is normal eqaution system and all is found; convert to
*     solution and covariance matrix
      if( neq_sys .and. neq_read .and. vec_read 
     .   .and. .not.apr_missed ) then
          print *,'Inverting NEQ system'
*         First apply some constraints
          call add_neq_con( cov_parm, np) 

*         Now add the Pre-inversion 
          call invert_vis(cov_parm, sol_parm, scale, ipivot, np, np, 1)
*
* MOD TAH 141119: Scale COD covariance matrix to make more reasonable
          if( qowner(1:3).eq.'COD' ) then
              write(*,'(a)') 'Scaling COD COV by 1.d-4'
              do i = 1, np
                  do j = 1, np
                      cov_parm(i,j) = cov_parm(i,j)*1.d-4
                  end do
              end do
          endif
 
* MOD TAH 210505: Change 2 to 1 to see inverted normal equations. 
         if( 1.eq.2 ) then
             write(*,'(a,1x,I6)') 'INVERT NEQ NP ',np

             write(*,510) (i,qapr(i)/qscale(i),sol_parm(i), 
     .            sol_parm(i)+qapr(i)/qscale(i),
     .            sqrt(cov_parm(i,i)),qscale(i),i=1,np) 
 510         format(('PN# ',i6,1x,3(F15.5,1x),1x,F10.5,1x,E10.2))
         endif
*
*         Now add the apriori back to the solution vector
          do i = 1, np
             sol_parm(i) = qapr(i)/qscale(i) + sol_parm(i)
             qsol(i) = sol_parm(i)
          enddo
          neq_sys = .false.
      else
          write(*,520) neq_sys , neq_read , vec_read, apr_missed
  520     format('ERROR: NEQ SYSTEM Missing parts: neq_sys ',L1,
     .           ' NEQ_READ ',L1,' VEC_READ ',L1,' APR_MISSED ',L1)
      endif

****  That all
      return
      end

 
    
