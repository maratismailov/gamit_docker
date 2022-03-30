CTITLE READ_GLB

      subroutine read_glb(cov_parm, sol_parm)
      
      implicit none

*     This routine will read the global file and save the information
*     needed for the sinex files      
      
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'
      
* PASSED VARIABLES
     
*   cov_parm(cnum_parn, cnum_parn)    - Covariance matrix
*   sol_parm(cnum_parn)          - Solution vector
 
      real*8 cov_parm(cnum_parn, cnum_parn), sol_parm(cnum_parn)

*
****  Start be reading the names from the binary files
      call sr_names
      
***** Get the full names
      call sr_full

***** Get the solution description
      call sr_description
      
***** Get the codes for the estimated parameters
      call sr_codes

***** Get the aprori parameter codes
      call sr_aprioris( cov_parm )
     
*     Finally read the solution and the aprioris
      call sr_soln( cov_parm, sol_parm)
      call sr_apr_cov 

***** Thats all
      return
      end
      
ctitle sr_names

      subroutine sr_names
      
      implicit none

*     This routine reads the names of the sites/satellites and 
*     radio sources (for VLBI).

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'

* LOCAL VARIABLES
*   i,j         - Loop counters
*   ierr        - FmpRead error.
*   iname(128)  - Integer equivalence to names_buffer (used
*               - for FmpRead.
*   len_read    - Length of record read from file (Not used)
*   lnum        - Number of current site or source in global
*               - list
*   site_boundary   - Record at which we shift from site names
*               - source names
*   cn, cs, cv  - Dummy counters for numbers of sites, sources
*                 and satellites
 
      integer*4 i,j, ierr, iname(128), len_read, lnum, 
     .          cn, cs, cv
 
 
*   names_buffer(32)    - Buffer containing the names to be
*               - read or written to the file.
*   next_name   - Next name read from the names buffer
 
      character*8 names_buffer(64), next_name
 
      equivalence ( names_buffer, iname )

***** Do site names and source names together.
      cn = 0
      cs = 0
      cv = 0
      do i = 1,cnum_names
         call readd(cglb_dcb,ierr, iname, 128, len_read,
     .              crec_names+i-1)
 
*****     Now loop over the names in the names buffer.
          do j = 1, cnum_sites
              call get_name_from_buffer( names_buffer, j, i, 0,
     .                next_name)
              call add_name_to_list( next_name, lnum,
     .                qsite_names, cs )
          end do
 
*****     Get the source names
          do j = 1, cnum_sources
              call get_name_from_buffer( names_buffer, j, i,
     .                cnum_sites, next_name)
              call add_name_to_list( next_name, lnum,
     .                qsource_names, cs )
          end do

*****     Get the SVS names
          do j = 1, cnum_svs
              call get_name_from_buffer( names_buffer, j, i,
     .                cnum_sites+cnum_sources, next_name)
              call add_name_to_list( next_name, lnum,
     .                qsvs_names, cv )
          end do
*                     ! Looping over records
      end do

***** Thats all
      return 
      end
      
CTITLE SR_FULL 

      subroutine sr_full
      
      implicit none

*     Routine to read the full names of the sites.

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'

* LOCAL VARIABLES
*   i,j,k       - Loop counters
*   ierr        - File error flag
*   full_list(1)   - Integer alias for cnames_list so that we can
*               - write to a file.
*   trimlen     - Length of string

      integer*4 i, ierr, full_list(max_glb_sites*8), len_read

*   cfull_list(1)  - List of site and sources names (saved in
*               - scr_common temporarily

      character*32 cfull_list(max_glb_sites)

      equivalence ( cfull_list, full_list )


***** Make the list of site and source names

      do i = 1, cnum_full*128
         full_list(i) = 8224*65536 + 8224
      end do

***** Now read record(s)

      call readd(cglb_dcb,ierr, full_list,128*cnum_full,
     .           len_read, crec_full)
      call report_error('FmpRead',ierr,'read','full block',1,
     .                  'sr_full')

      do i = 1, cnum_sites
         qfull_names(i) = cfull_list(i) 
      end do

***** Thats all
      return 
      end
    
CTITLE SR_DESCRIPTION

      subroutine sr_description     
     
      implicit none

*     This routine will read the descriptions of the indivual
*     solutions that went into the binary h-file.

*     Routine to read the full names of the sites.

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'
      
* LOCAL VARIABLES     
*   i,j,k       - Loop counters
*   ierr        - File error flag
*   in_rec     - Record number for input
*   rec_copy(128) - Record for station information
*   len_read    - Length of record read
      integer*4 i,j, ierr, in_rec, ns, rec_copy(128), len_read
      
****  Start reading over the records
      in_rec = crec_solutions
      
* MOD TAH 020628: Check to make value is not too large
      if( cnum_soln_recs.gt.max_sln_save ) then
*         Too many solution records, tell user what to do
          write(*,140)  max_sln_save,  cnum_soln_recs
 140      format(/,'**DISASTER** Number of saved solutions exceeds ',
     .           ' limit of ',i6,/,
     .           '             Modify the parameter max_sln_save ',
     .           ' in gg/kf/htoglb/htoglb_comm.h to at least ',I7,/,
     .           '             Use make to re-make the glbtosnx.')
          call report_stat('FATAL','glbtosnx','sr_description',' ',
     .           'Too many solution records',0)
      end if
      	  	      
      do ns = 1, cnum_soln_recs
      
          call readd(cglb_dcb,ierr,sdelete_count,128,len_read, in_rec)
         
*         Now save the information

          qtai_utc = stai_utc 
      
          qut1_apr(1) = sut1_apr(1)
          qut1_apr(2) = sut1_apr(2)

          do i = 1,2
              do j = 1,2
                  qwob_apr(j,i) = swob_apr(j,i) 
                  qnut_ang_apr(j,i) = snut_ang_apr(j,i) 
              end do
          end do

          qsepoch_start(ns) = sepoch_start 
          qsepoch_end(ns)   = sepoch_end  
          qsrun_time(ns)    = srun_time    
          qsnum_parn(ns)    = snum_parn    
          qscons_type(ns)   = scons_type   
          qssys_type(ns)    = ssys_type 
          qsglb_ver(ns)     = sglb_vers 
          call sub_null( sowner )
          qsowner(ns)       = sowner
          call sub_null( screator )
          qscreator(ns)     = screator 
          call sub_null( sprog_gen ) 
          qsprog_gen(ns)    = sprog_gen 
          call sub_null ( sanal_type )
          qsanal_type(ns)   = sanal_type  
          qskalobs_file(ns) = sKalObs_file
          qsexpt_title(ns)  = sexpt_title 

          in_rec = in_rec + 1
      end do
      
***** Now read the station information
      in_rec = crec_sinf
      do i = 1, cnum_sites, 2
          call readd(cglb_dcb, ierr, rec_copy, 128, len_read, in_rec)
          call sr_sinf_rec(i,rec_copy(1))
          call sr_sinf_rec(i+1, rec_copy(65))
          in_rec = in_rec + 1
      end do

****  Now read in the satellite information (if available)
      in_rec = crec_svinf
      if( in_rec.gt.0 ) then
          do i = 1, cnum_svs, 4
             call readd(cglb_dcb, ierr, rec_copy, 128, len_read, in_rec)
             call sr_svinf_rec(i,rec_copy(1))
             call sr_svinf_rec(i+1, rec_copy(33))
             call sr_svinf_rec(i+2, rec_copy(65))
             call sr_svinf_rec(i+3, rec_copy(97))
             in_rec = in_rec + 1
          end do
      endif 

      
***** Thats all
      return
      end
      
CTITLE SR_SINF_REC

      subroutine sr_sinf_rec(i, record )
      
      implicit none

*     Routine to decode the station information records

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'
      
* PASSED VALUES

* i  -  Station number
* record(64)  - Record for the Station information

      integer*4 i, j, record(64)


*     Start saving the information
      call wmov( record, 1, ssdata_st, 1, 64 )

*     Save the values
*     Code already known, so don't save
*     qsite_names(i) = sscode 
      qdata_st(i) = ssdata_st 
      qdata_en(i) = ssdata_en  
      qrecv_st(i) = ssrecv_st 
      qrecv_en(i) = ssrecv_en 
      qante_st(i) = ssante_st 
      qante_en(i) = ssante_en 
      do j = 1,3
         qarp_ecc(j,i) = ssarp_ecc(j)
         qL1a_ecc(j,i) = ssL1a_ecc(j) 
         qL2a_ecc(j,i) = ssL2a_ecc(j) 
      end do
* MOD TAH 200205: Save the antenna azimth offset 
      qantdaz(i) = ssantdaz
      qelev_cut(i) = sselev_cut 
      qnum_zen(i)  = ssnum_zen 
      call sub_null(ssant_mod )
      if( ichar(sti_antmod(1:1)).ne.0 ) then  ! New full name
         qant_mod(i)  = sti_antmod
      else
         qant_mod(i)  = ssant_mod 
      end if 
      call sub_null(ssante_sn)
      qante_sn(i)  = ssante_sn
      call sub_null( ssante_ty)      
      qante_ty(i)  = ssante_ty
      call sub_null ( ssradome_ty )
      qradome_ty(i) = ssradome_ty 
      call sub_null(ssrecv_sn)       
      qrecv_sn(i)  = ssrecv_sn 
      call sub_null( ssrecv_ty)    
      qrecv_ty(i)  = ssrecv_ty
      if( ichar(ssrecv_ty_end(1:1)).ne.0 ) then
         qrecv_ty(i) = ssrecv_ty // ssrecv_ty_end
      end if
      call sub_null(ssrecv_fw)      
      qrecv_fw(i)  = ssrecv_fw  

* MOD TAH 190627: Save the local site number to see if used
      gtol_sites(i) = ss_lnum

      return
      end

CTITLE SR_SVINF_REC

      subroutine sr_svinf_rec(i, record )
      
      implicit none

*     Routine to decode the satellite information records

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/svinf_def.h'
      
* PASSED VALUES

* i  -  Station number
* record(32)  - Record for the Station information

      integer*4 i, j, record(32)


*     Start saving the information
      call wmov( record, 1, svi_prn, 1, 32 )

*     Save the values
      qsvi_prn(i)    = svi_prn
      qsvi_svn(i)    = svi_svn
      qsvi_block(i)  = svi_block

      qsvi_antmod(i) = svi_antmod
      qsvi_ocode(i)  = svi_ocode

      do j = 1,3
         qsvi_antpos(j,1,i) = svi_antpos(j,1)
         qsvi_antpos(j,2,i) = svi_antpos(j,2)
      end do
      qsvi_launch(i) = svi_launch

* MOD TAH 19627: Save to local satellite number to see if used
*     (-1 if not used)
      gtol_svs(i) = svi_lnum
                             
      return
      end
      
CTITLE SR_APRIORIS 
     
      subroutine sr_aprioris( apr_values )    
      
      implicit none

*     This routine decodes the codes for the apriori values
*     and saves their values.  Currently, only station positions,
*     velocities.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'

* Passed variables

*   apr_values(*)   - array for read the apriori's into

      real*8 apr_values(*)
      
*   apr_codes(128)  - One record worth of apriori codes
*   i,j,k           - Loop counters
*   ierr            - File manipluation error routines
*   indx            - Index for the parameter being decoded eg. could
*                   - site number, source number etc.
*   len_read        - Length from readd.
*   type            - Type of parameter being decoded. See
*                   - GLB_HEADER_DEF.FTNI for descriptions.

*   rec_soln        - Number of records before the first
*                     record with the solution vector in
*                     it.
*   len_soln        - Number of records that need to be
*                     read to get the solution vector
*   start_soln      - The start of the solution in the first
*                     record that has been read
*   nea(2,3)        - Countersr for the epoch number of the multi-pmu
*                     parameters.
*   ne              - Individual value for nea.
  
      integer*4 apr_codes(128), i,j, ierr, indx, len_read, type,
     .          orp, ep_ent, ne, nea(2,3), orb_el, sv_num
  
 
***** Read all of the apriori values into the ema area.
      call readd(cglb_dcb , ierr, apr_values, 128*cnum_apr_vals,
     .           len_read, crec_apr_vals)
      call report_error('VREAD',ierr,'read','apriori values',
     .                  0,'GET_APRIORIS')
      if( ierr.ne.0 ) RETURN
      
*     Clear the apriori multi-pmu parameter array
      do i = 1,2
         do j = 1,3
            nea(i,j) = 0 
         end do
      end do
      
*     See if we need to read the epoch records for multi-pmu
*     paramaters
      if( cent_apr_ep.gt.0 ) then
          call readd(cglb_dcb, ierr, qmul_apr_ep, 128*cnum_apr_ep,
     .               len_read, crec_apr_ep)
      end if

***** Loop over the apriori code records and get the code values.
*     These are then decoded and the corresponding values are
*     read from apriori values block
      do i = 1, cnum_apr_types
 
          call readd(cglb_dcb, ierr, apr_codes, 128, len_read,
     .               crec_apr_types+i-1) 
*****     Now start decoding and saving the corresponding values.
          do j = 1, 128
 
*             Only decode if a valid value
              if( j+(i-1)*128.le. cnum_apr_codes ) then
                  call decode_code ( apr_codes(j), type, indx )
                  
*                 See if this is a type that we want to decode
                  if( type.eq. 7 ) then  
                      site_pos(1,indx) = apr_values(j+(i-1)*128)
                  else if( type.eq. 8 ) then
                      site_pos(2,indx) = apr_values(j+(i-1)*128)
                  else if( type.eq. 9 ) then
                      site_pos(3,indx) = apr_values(j+(i-1)*128)
                  else if( type.eq. 42 ) then
                      site_vel(1,indx) = apr_values(j+(i-1)*128)
                  else if( type.eq. 43 ) then
                      site_vel(2,indx) = apr_values(j+(i-1)*128)
                  else if( type.eq. 44 ) then
                      site_vel(3,indx) = apr_values(j+(i-1)*128)
* MOD TAH 050622: Check the apriri satellite axis offset
                  else if ( type.eq.51 ) then
* MOD TAH 190610: No need to change for ECOMC because PCO not changed.
                      call decode_code( indx, orb_el, sv_num ) 
                      if( orb_el.ge.21 .and.orb_el.le.23 ) then
                         svs_pos(max_svs_elem+(orb_el-21)-2,sv_num) =
     .                       apr_values(j+(i-1)*128)
                      endif

                  else if( type.eq. 56 ) then
                      call decode_code( indx, orp, ep_ent )
                      qnum_mul_pmu(orp,1) = qnum_mul_pmu(orp,1)+1
                      nea(orp,1) = nea(orp,1) + 1
                      ne = nea(orp,1)
                      qpmu_mul_apr(orp,1,ne) = apr_values(j+(i-1)*128)
                   else if( type.eq. 57 ) then
                      call decode_code( indx, orp, ep_ent )
                      qnum_mul_pmu(orp,2) = qnum_mul_pmu(orp,2)+1
                      nea(orp,2) = nea(orp,2) + 1
                      ne = nea(orp,2)
                      qpmu_mul_apr(orp,2,ne) = apr_values(j+(i-1)*128)
                   else if( type.eq. 58 ) then
                      call decode_code( indx, orp, ep_ent )
                      qnum_mul_pmu(orp,3) = qnum_mul_pmu(orp,3)+1
                      nea(orp,3) = nea(orp,3) + 1
                      ne = nea(orp,3)
                      qpmu_mul_apr(orp,3,ne) = apr_values(j+(i-1)*128)/
     .                               15.d0
                      if( orp.eq.2 ) then
                          qpmu_mul_apr(orp,3,ne)= 
     .                                         -qpmu_mul_apr(orp,3,ne) 
                      endif
                     
                  end if
              end if
C              print *,'APR: ',i,j,j+(i-1)*128, type,indx, 
C     .                 apr_values(j+(i-1)*128)
          end do
      end do
 

****  Thats all
      return
      end
 
CTITLE SR_CODES 

      subroutine sr_codes 

      implicit none

*     Routine to save the codes for the estimated parameters
*     Currently only station positions/velocities and eop
*     parameters are saved.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'

* LOCAL PARAMETERS
      
*   i,j,k           - Loop counters
*   ierr            - File manipluation error routines
*   indx            - Index for the parameter being decoded eg. could
*                   - site number, source number etc.
*   len_read        - Length from readd.
*   type            - Type of parameter being decoded. See
*                   - GLB_HEADER_DEF.FTNI for descriptions.

*   rec_soln        - Number of records before the first
*                     record with the solution vector in
*                     it.
*   len_soln        - Number of records that need to be
*                     read to get the solution vector
*   start_soln      - The start of the solution in the first
*                     record that has been read
*   ep_ent          - Epoch Entry number for multi-pmu values
*   orp             - Offset or rate index
*   orb_el          - Orbit element number
*   sv_num          - Satellite number
 
      integer*4 i,j, ierr, indx, len_read, type,
     .          no, np, ep_ent, orp, ne, k, orb_el, sv_num

* MOD TAH 190610; ECOMC model update
      integer*4 new_code ! Remapped satellire parametert code/
      integer*4 svcode_to_code ! Function to return new code given old one.
  
  
***** Loop over the parameter code records and get the code values.
*     These are then decoded and the corresponding values are
*     read from apriori values block
      do i = 1, cnum_parn
         itoo(i) = 0
         atos(i) = 0
         qscale(i) = 1.d0
      end do
      
      do i = 1,2
         do j = 1,3
            qnum_mul_pmu(i,j) = 0 
         end do
      end do
      
*     See if we need to read the epoch records for multi-pmu
*     paramaters
      if( cent_par_ep.gt.0 ) then
          call readd(cglb_dcb, ierr, qmul_par_ep, 128*cnum_par_ep,
     .               len_read, crec_par_ep)
      end if
      
      
*     Read in all the parameter codes and then save as individual
*     data types       
      call readd(cglb_dcb, ierr, qglb_codes, 128*cnum_par_types, 
     .           len_read, crec_par_types)
 
***** Now start decoding and saving the corresponding values.
      do j = 1, cnum_parn
*         
          np = j   
          call decode_code ( qglb_codes(j), type, indx )
 
*         See if this is a type that we want to decode
          if( type.eq. 7 ) then  
* MOD TAH 190627: Only save if site used in solution (indx is
*             site number). (0 test is for backwads compatabilty)
*             Only test on position since this arises from back solutions.
              if ( gtol_sites(indx).ge.0 ) then
                  qparn_sites(1,indx) = np
              endif 
          else if( type.eq. 8 ) then
              if ( gtol_sites(indx).ge.0 ) qparn_sites(2,indx) = np
          else if( type.eq. 9 ) then
              if ( gtol_sites(indx).ge.0 ) qparn_sites(3,indx) = np
          else if( type.eq. 42 ) then
              qparn_vel(1,indx) = np
          else if( type.eq. 43 ) then
              qparn_vel(2,indx) = np
          else if( type.eq. 44 ) then
              qparn_vel(3,indx) = np
          else if( type.eq. 13) then
              if( indx.eq.1 ) qparn_pmu(1,1) = np
              if( indx.eq.2 ) qparn_pmu(1,2) = np
              if( indx.eq.3 ) qparn_pmu(2,1) = np
              if( indx.eq.4 ) qparn_pmu(2,2) = np
* MDO TAH 070824: Units for rates are per-day so factor
*             1.d0/86400.d0 not needed.
              if( indx.eq.3 .or. indx.eq.4 ) 
     .                        qscale(np) = 1.d0
          else if( type.eq. 14) then
              if( indx.eq.1 ) qparn_pmu(1,3) = np
              if( indx.eq.1 ) qscale(np) = 1.d0/15.0d0
              if( indx.eq.2 ) qparn_pmu(2,3) = np
              if( indx.eq.2 ) qscale(np) = -1.d0/15.0d0
* MOD TAH 050622: Check if satellite antenna offsets are
*         estimated
* NOTE TAH 190606: Only Satellite PCO treated so don't need
*         to test from 107 cglb_vers change to radiation 
*         parameter mapping.
          else if ( type.eq.51 ) then
* MOD TAH 1910610: See if we need to reassign code
* MOD TAH 191201: Removed update here becuase causes problems in 
*          hfuopd with the aprioris (all apriori slotes are written
*          so change to 107 wouold beed aproiri vector and covariance
*          re-arrangement,
C             if( cglb_vers.lt.107 .and. type.eq.51 ) then
C                 new_code = svcode_to_code(qglb_codes(j))
C                 qglb_codes(j) = new_code
*                 Now re-code it.
C                 call decode_code ( qglb_codes(j), type, indx )
C             endif
              call decode_code( indx, orb_el, sv_num ) 
              if( orb_el.ge.21 .and.orb_el.le.23 ) then
* MOD TAH 190627: Only save the parameter number if the satellite 
*                 was used. (0 test is for backwads compatabilty)
                  if( gtol_svs(sv_num).ge.0 ) then
                     qparn_svs(max_svs_elem+(orb_el-21)-2,sv_num) = np
                  endif
              endif

* MOD TAH 050622: Change sign of translation and scale so that
*         they are consistent with IGS standard
          else if( type.eq.52 ) then
              qparn_tran(indx,1) = np
              qscale(np) = -1.d0
          else if( type.eq.53 ) then
              qparn_tran(indx,2) = np
              qscale(np) = -1.d0
          else if( type.eq.54 ) then
              qscale(np) = -1.d0
              qparn_scale(1)    = np
          else if( type.eq.55 ) then
              qparn_scale(2)    = np
              qscale(np) = -1.d0
* MOD TAH 981020: Added reading of multiple PMU parameters
          else if( type.eq.56 ) then
              call decode_code( indx, orp, ep_ent )
              qnum_mul_pmu(orp,1) = qnum_mul_pmu(orp,1)+1
              ne = qnum_mul_pmu(orp,1)
              qparn_mul_pmu(orp,1,ne) = np 
              qref_ep(np) = qmul_par_ep(ep_ent)
          else if( type.eq.57 ) then
              call decode_code( indx, orp, ep_ent )
              qnum_mul_pmu(orp,2) = qnum_mul_pmu(orp,2)+1
              ne = qnum_mul_pmu(orp,2)
              qparn_mul_pmu(orp,2,ne) = np 
              qref_ep(np) = qmul_par_ep(ep_ent)
          else if( type.eq.58 ) then
              call decode_code( indx, orp, ep_ent )
              qnum_mul_pmu(orp,3) = qnum_mul_pmu(orp,3)+1
              ne = qnum_mul_pmu(orp,3)
              qparn_mul_pmu(orp,3,ne) = np 
              qref_ep(np) = qmul_par_ep(ep_ent) 
              if( orp.eq.1 ) qscale(np) = 1.d0/15.0d0
              if( orp.eq.2 ) qscale(np) = -1.d0/15.0d0
          end if
      end do

****  Now get the mapping of the output parameters to the input ones.
      no = 0
      do i = 1, cnum_sites
         do j = 1, 3
            if( qparn_sites(j,i).ne.0 ) then
                no = no + 1
                otoi(no) = qparn_sites(j,i)
                itoo(qparn_sites(j,i)) = no
                atos(qparn_sites(j,i)) = i 
            end if
         end do
         do j = 1, 3
            if( qparn_vel(j,i).ne.0 ) then
                no = no + 1
                otoi(no) = qparn_vel(j,i)
                itoo(qparn_vel(j,i)) = no
                atos(qparn_vel(j,i)) = i
            end if
         end do
      end do

* MOD TAH 0506022: Check satelliet offsets
      do i = 1, cnum_svs
         do j = 1,3
            if( qparn_svs(max_svs_elem-3+j,i).ne.0 ) then
                no = no + 1
                otoi(no) = qparn_svs(max_svs_elem-3+j,i)
                itoo(qparn_svs(max_svs_elem-3+j,i)) = no
                atos(qparn_svs(max_svs_elem-3+j,i)) = i
            end if
          end do
      end do   
                
****  Do translation and rate
      do j = 1,2 
         do i = 1,3
            if( qparn_tran(i,j).ne.0 ) then
                no = no + 1
                otoi(no) = qparn_tran(i,j)
                itoo(qparn_tran(i,j)) = no
            end if
         end do
      end do

*     Now do the scale and rate of change    
      do i = 1,2
         if( qparn_scale(i).ne.0 ) then
             no = no + 1
             otoi(no) = qparn_scale(i)
             itoo(qparn_scale(i)) = no
         end if
      end do

****  Now do the PMU parameters
      do i = 1, 3
         do j = 1, 2
            if( qparn_pmu(j,i).ne.0 ) then
                no = no + 1
                otoi(no) = qparn_pmu(j,i)
                itoo(qparn_pmu(j,i)) = no
            end if
         end do
      end do
      
****  Now do the multi-PMU parameters
      do i = 1, 3
         do j = 1, 2
            do k = 1, qnum_mul_pmu(j,i) 
               if( qparn_mul_pmu(j,i,k).ne.0 ) then
                   no = no + 1
                   otoi(no) = qparn_mul_pmu(j,i,k)
                   itoo(qparn_mul_pmu(j,i,k)) = no
               end if
            end do
         end do
      end do

***** Save the number of parameters to be output to the SINEX file
      qnum_parn = no

****  Read the apriori codes.  The actual aprioris are read
*     in sr_aprioris
      call readd(cglb_dcb, ierr, qapr_codes, 128*cnum_apr_types, 
     .           len_read, crec_apr_types)
      
***** Thats all the apriori values we need
      return
      end

CTITLE SR_SOLN 

      subroutine sr_soln( cov_parm, sol_parm) 

      implicit none

*     Routine to read the covariance matrix for this solution

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'
      
* PASSED VARIABLES
     
*   cov_parm(cnum_parn, cnum_parn)    - Covariance matrix.  At one to
*                                     - dimension for soln
*   sol_parm(cnum_parn)          - Solution vector
 
      real*8 cov_parm(cnum_parn, cnum_parn+1), sol_parm(cnum_parn)

*
* LOCAL PARAMETERS

      integer*4 len_read, ierr, i, j

      integer*8 I8

      data I8 / 1 /
      
***** Read the full covariance matrix and then scale the results
* COMMENT TAH 190520: Only a margin of 16 in cnum_par_vals not
* overflowing I*4 when number of paramters is greater than 32767.
* MOD TAH: Change to writd8 to allow for >32767 parameters. 
      if( cnum_parn.gt.32767 ) then 
         call readd8(cglb_dcb,ierr,cov_parm, (I8*128)*cnum_par_vals,
     .              len_read, crec_par_vals)
      else
         call readd(cglb_dcb,ierr,cov_parm, 128*cnum_par_vals, 
     .              len_read, crec_par_vals)
      endif

      do i = 1, cnum_parn
          sol_parm(i) = sol_parm(i) * qscale(i)
          do j = 1, cnum_parn
             cov_parm(i,j) = cov_parm(i,j)* qscale(i)*qscale(j)
          end do
      end do
      
****  Thats all 
      return
      end
      
      
ctitle sr_apr_cov

      subroutine sr_apr_cov
      
      implicit none

*     This routine will read the apriori covariance matrix used
*     in the solution.
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'
 

      integer*4 i,j, len_read, record(128), np1, np2, no1, no2,
     .          pir, ierr, j1, j2, is, lastj, k, l,m 

      logical found

***** Loop over the records saving the values that need to be
*     saved

      do i = 1, cnum_acvc_recs
          call readd(cglb_dcb, ierr, record, 128, len_read, 
     .               crec_acvc+i-1)
     
****      Now loop over the triplets in the record and see if we
*         should save the values
          if( i.eq.cnum_acvc_recs) then
              lastj = (cnum_acvc - (i-1)*128-1)/4 + 1
          else
              lastj = 32
          end if

          do j = 1, lastj
              pir = (j-1)*4 + 1
              np1 = record(pir)
              np2 = record(pir+1)

*             Now find if these are in the output parameter list
              if( np1.le.max_glb_parn .and. 
     .            np2.le.max_glb_parn .and. 
     .            np1.gt.0 .and.np2.gt.0 ) then
              no1 = itoo(np1)
              no2 = itoo(np2)
              if( no1.ne.0 .and. no2.ne.0 ) then

*                 This parameter is in the output so save the value
*                 See what the parameter is.
                  found = .false.
*                 Now see if atos points to site or satellite.
                  is = atos(np1)
                  if( is.gt.0 ) then
                     j1 = 0
                     j2 = 0 
                     do k = 1, 3
                        if( np1.eq. qparn_sites(k,is) ) then
                           j1 = k
                        end if
                        if( np2.eq. qparn_sites(k,is) ) then
                            j2 = k
                        end if
                     end do
                     if( j1.gt.0 .and.j2.gt.0 ) then
                         call dwmov(record(pir+2),1,
     .                              qapr_cov(j1,j2,is),1,1)
                         call dwmov(record(pir+2),1,
     .                              qapr_cov(j2,j1,is),1,1)
                         found = .true.
                     end if
*                    If not found see if velocity
                     if( .not.found ) then
                        do k = 1, 3
                           if( np1.eq. qparn_vel(k,is) ) then
                              j1 = k
                           end if
                           if( np2.eq. qparn_vel(k,is) ) then
                               j2 = k
                           end if
                        end do
                        if( j1.gt.0 .and.j2.gt.0 ) then
                            call dwmov(record(pir+2),1,
     .                                 qapr_vel(j1,j2,is),1,1)
                            call dwmov(record(pir+2),1,
     .                                 qapr_vel(j2,j1,is),1,1)
                            found = .true.
                        end if
                     end if
* MOD TAH 050622:  If not found then see if satellite.  Constraints
*                    here are diagonal so check just the one index
                     if( .not.found ) then
                         do k = 1,3
                            if( np1.eq.
     .                          qparn_svs(max_svs_elem-3+k,is) ) then
                                j1 = k
                            end if
                            if( j1.gt.0  ) then
                                call dwmov(record(pir+2),1,
     .                            qapr_svs_ant(j1,is),1,1)
                               found = .true.
                            end if
                         end do
                     end if

                  else
*                    Must be polar motion and therefore diagonal
                     if( np1.ne.np2 ) then
                         write(*,300) np1, np2
 300                     format('**WARNING** Expecting diagonal',
     .                          ' apriori covariance P1 P2 ',2i5)
                     end if

*                    Do normal polar motion
                     do k = 1,3
                        do l = 1, 2
                           if( qparn_pmu(l,k).eq.np1 ) then
                               call dwmov(record(pir+2),1,
     .                                    qapr_pmu(l,k),1,1)
                               if( k.eq.3 ) then
                                   qapr_pmu(l,k) = qapr_pmu(l,k)/
     .                                             15.d0**2
                               end if
                           end if
                        end do
                     end do   

*                    Do multi-epoch pmu.  We need to consider how to
*                    do multiple offdiagonal elements.
                     do k = 1,3
                        do l = 1, 2
                           do m = 1, qnum_mul_pmu(l,k)
                              if( qparn_mul_pmu(l,k,m).eq.np1 ) then
                                  call dwmov(record(pir+2),1,
     .                                       qpmu_mul_asig(l,k,m),1,1)
                                  if( k.eq.3 ) then
                                      qpmu_mul_asig(l,k,m) = 
     .                                           qpmu_mul_asig(l,k,m)/
     .                                           15.d0**2
                                  endif
                               end if
                           end do
                        end do
                     end do  
*                    Do translation and rate of change                     
                     do l = 1,2
                        do k = 1,3
                           if( qparn_tran(k,l).eq.np1 ) then
                               call dwmov(record(pir+2),1, 
     .                                    qapr_tran(k,l),1,1)
                           end if
                        end do
                     end do
*                    Do scale and rate of change                     
                     do k = 1,2
                        if( qparn_scale(k).eq.np1 ) then
                            call dwmov(record(pir+2),1, 
     .                                 qapr_scale(k),1,1)
                        end if
                     end do
                 end if
              end if
              end if
          end do
      end do

***** Thats all
      return
      end 

