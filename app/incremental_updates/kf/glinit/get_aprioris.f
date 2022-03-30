CTITLE GET_APRIORIS
 
      subroutine get_aprioris( ierr, apr_values, crunjd, grunjd )
 
      implicit none
 
*     Routine to read the apriori codes and then read the values
*     associated with each of these and save the results in the
*     GLOBK common block.
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
      include '../includes/svinf_def.h'
      include '../glinit/qsvi_rec.h'
 
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
 
      integer*4 apr_codes(128), i,j, ierr, indx, len_read, type,
     .          len_soln, k 
      integer*4 rec_soln , start_soln
      integer*4 rec_copy(128)   ! Buffer for reading svs information
                                ! records
      integer*4 in_rec          ! SVS record number to read

*   apr_values(64554j)   - apriori values read from file (stored
*                    - temporaryily in ema). Use large size to
*                    - convince compiler to used I*4)
      real*8 apr_values(*)
 
      real*8 crunjd, grunjd  ! Run time of current hfile being
              ! read (crunjd), run time of latest solution in this
              ! epoch of data (grunjd) and sectag for conversion.
              ! Values are used to see if SV antenna offsets should
              ! by updated. 

* MOD TAH 190610: Updated the satellite radition parameters codes
*     accomiodate the ECOMC model
      integer*4 new_code  ! Remapped satellite code (51).
     .,         svcode_to_code ! Function to return new code given old one.

      integer*4 incode, inindx, intyp, insv
      integer*4 oorb_el, osv_num, norb_el, nsv_num
      integer*4 nind  ! Indices in apr_values array to get
                      ! remapped rad parameters
     .,         off   ! Offset to get correct aprioiri for 107 rad model  
* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

* MOD TAH 200130: Test to see if we have satellite PCO.
      logical svs_pco_set  ! Set true when PCO in apriori list
                           ! if not value from svs_inf is saved

      data I8 / 1 /

***** Read all of the apriori values into the ema area.
      svs_pco_set = .false.
 
      call readd(cglb_dcb , ierr, apr_values, 128*cnum_apr_vals,
     .           len_read, crec_apr_vals)
      call report_error('VREAD',ierr,'read','apriori values',
     .                  0,'GET_APRIORIS')
      if( ierr.ne.0 ) RETURN
 
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
                  incode = apr_codes(j)
                  inindx = indx
* MOD TAH 200123: Added decoding of satellite index to generate
*                 the orbit element number when lastesr cglb_vers
*                 files are used.
                  call decode_code( indx,   norb_el, nsv_num )

* MOD TAH 1910610: See if we need to reassign code
                  off = 0 
                  if( cglb_vers.lt.107 .and. type.eq.51 ) then
                      new_code = svcode_to_code(apr_codes(j))
                      apr_codes(j) = new_code
*                     Now re-code it.
                      call decode_code ( apr_codes(j), type, indx )
* MOD TAH 191201: Also ssap the apriori value as well.  All aprioris
*                 are in memory so should be straight forward.
                      call decode_code( indx,   norb_el, nsv_num )
                      call decode_code( inindx, oorb_el, osv_num )
                      nind = (i-1)*128+j 
                      if( oorb_el.eq.9 ) then  ! First re-mapped rad parameters
                                 ! *      9  - B AXIS BIAS       - Z AXIS BIAS (never used) 0
* MID TAH 191201: Remap the apriori values
                          apr_values(nind) = apr_values(nind+1)
                          do k = 1,6
                              apr_values(nind+k)=apr_values(nind+k+2)
                          enddo
*                         Clear remaining values
                          do k = 7,7+4
                              apr_values(nind+k)=0
                          enddo                        

                      endif
                      if( norb_el.eq.9 ) off = 1
                      if( norb_el.ge.10 .and. norb_el.lt.20 ) off = 2
                  endif

*                 Now save the apriori value if needed.
* MOD TAH 950830: Turned off checking of the aprioris so that apriori
*                 station coordinates can be passed in the global file.
*                 (Original code only fixed parameters were put in the
*                 apriori list, now stations and satellites are put
*                 there).
                  if( type.ne.51 ) then
                     call save_apr_val( type, indx, (i-1)*128+j,
     .                                  apr_values, 0, 0)
* MOD TAH 191201: Only save value if norb_el is > 0 (i.e., valid)
                  elseif( norb_el.gt.0 ) then
* MOD TAH 200429: Removed off from index becauase 191201 code remaps the
*                 apriories to the correct place.
                     off = 0
                     call save_apr_val( type, indx, (i-1)*128+j-off,
     .                                  apr_values, 0, 0)
                     if( norb_el.ge.21 ) svs_pco_set = .true.
                  endif 

              end if
          end do
      end do
 
***** Now read in the solution vector (and covariance matrix) and
*     get approximate apriori values based on estimates.

*     Compute the record number we need to start at to get the 
*     solutiuon vector
      rec_soln = 2*(I8*cnum_parn)*cnum_parn/128 + crec_par_vals
      len_soln = cnum_par_vals - (rec_soln - crec_par_vals)
      start_soln = (I8*cnum_parn)*cnum_parn - 
     .             (rec_soln-crec_par_vals)*64 + 1

      call readd(cglb_dcb , ierr, apr_values, 128*len_soln,
     .           len_read , rec_soln )
 
      call report_error('VREAD',ierr,'read','parameter values',
     .                  0,'GET_APRIORIS')
      if( ierr.ne.0 ) RETURN
 
***** Loop over the parameter code records and get the code values.
*     These are then decoded and the corresponding values are
*     read from apriori values block
      do i = 1, cnum_par_types
 
          call readd(cglb_dcb, ierr, apr_codes, 128, len_read,
     .               crec_par_types+i-1)
 
*****     Now start decoding and saving the corresponding values.
          do j = 1, 128
 
*             Only decode if valid code
              if( j+(i-1)*128.le. cnum_parn ) then
                  call decode_code ( apr_codes(j), type, indx )
 
* MOD TAH 1910610: See if we need to reassign code
                  if( cglb_vers.lt.107 .and. type.eq.51 ) then
                      new_code = svcode_to_code(apr_codes(j))
                      apr_codes(j) = new_code
*                     Now re-code it.
                      call decode_code ( apr_codes(j), type, indx )
                  endif

*                 Now save the apriori value if needed.
*                 ! Point to start of solution vector
                  call save_apr_val( type, indx, (i-1)*128+j,
     .                 apr_values(start_soln), 0, 1)
              end if
          end do
      end do
 
***** Now save some other apriori values based on values written in
*     header.  We only save these values for the first experiment
 
      if( num_glb_sol.eq.1 ) then
          apr_val_gamma = 1.d0
 
*                              ! These values are adjustments to the
          apr_val_wob(1) = 0
*                              ! apriori values from the the global file
          apr_val_wob(2) = 0
*                              ! (At the moment anyway)
          apr_val_ut1(1) = 0
          apr_val_ut1(2) = 0
 
      end if
 
***** Now save earth paramters. Later these could be extracted from
*     apr_code values saved in the apriori block or solution block
*     of the global matrix.  Not much need for this at the moment
*     since we always use the same values.
 
      do i = 1, gnum_sites
          do j = 1,3
              apr_val_tid(j,i) = cetd_apr(j)
          end do
      end do

* MOD TAH 131111: Now check the satellite antenna offset values.  Here we 
*     read the records and update if this hfile has a run-time greater than
*     any previously processed ones. 
*     Read over the satellite information records
      if( crec_svinf.gt.0 ) then   ! We have satellite records
          in_rec = crec_svinf
          do i = 1, cnum_svs, 4
             call readd(cglb_dcb, ierr, rec_copy, 128, len_read, in_rec)
             call gr_svinf_rec(i,rec_copy(1))
             call gr_svinf_rec(i+1, rec_copy(33))
             call gr_svinf_rec(i+2, rec_copy(65))
             call gr_svinf_rec(i+3, rec_copy(97))
             in_rec = in_rec + 1
          end do

*         See if we should update based on crunjd (current) and saved grunjd
* MOD TAH 191221: Removed test on run-time.  Original idea was to latest values
*         but that is not needed anymore.
C         if( crunjd.ge.grunjd ) then
*         copy the satellite antenna offsets into the apr_val_svs arrays.
          grunjd = crunjd
          do i = 1, cnum_svs 
             do j = 1,3  !  XYZ offset
*               Only copy (L1/LC entry into apriori; assume the 2 frequency
*               values are the same).
* MOD TAH 200130: Only reset value to svs_inf value if we have not saved it.
                if( .not. svs_pco_set )
     .          apr_val_svs(max_svs_elem-3+j,ltog_svs(i)) = 
     .                                    qsvi_antpos(j,1,i)
             enddo
          enddo
      endif 
  
 
***** Thats all the apriori values we need
      return
      end
 
