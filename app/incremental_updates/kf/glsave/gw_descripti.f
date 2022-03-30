CTITLE GW_DESCRIPTION

      subroutine gw_description( type )

      implicit none 
 
 
*     Routine to create the solution description and write the record
*     to the global file.
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/globk_common.h'

* MOD TAH 190621: Added type option for BACK solution writes.  During
*     back solution only need information on current glb_inp_file rather
*     re-reading all the input files in the srt_file.  (This file is already
*     open in glbak so we can't directly re-read it).

* PASSED VARIABLES
      character*(*) type   ! When set to BACK glb_inp_file used 
      
* LOCAL VARIABLES 
*   i,j,k       - Loop counters
*   ierr        - File error flag
*   isort_dcb(144)  - Sort file DCB buffer
*   iglb_dcb(16)    - DCB buffer of reading each of the global
*                   - files in this solution
*   len_read        - Length of record read from input file
*   out_rec         - Record number to write with the solution
*                   - description.
*   rec_copy(128)   - First record of the global files
*                   - read in and a place to read the solution
*                   - records into.
*                   - ***NOTE*** If this global files
*                   - header is changed then this routine may
*                   - need updating.
*   sol_rec         - Record for the solution description part
*                   - of the global file being copied
*   sol_num         - Number of solution records to copy
*   oglb_dcb(16)    - Original DCB buffer.  Copied here so that
*                     we can read the headers.
 
      integer*4 i,j,k, ierr, oglb_dcb(16),
     .    len_read, out_rec, rec_copy(128), sol_rec,
     .    sol_num
     

* MOD TAH 980519: Added reading and writing of forward and 
*   back chi**2/f to srt_file.

*   for_chi, bak_chi -- Forward and backwarsd chi**2/f

      real*4 for_chi, bak_chi
             
*   old_srec_sinf - Starting record in the individual globals
*         for the sation info records
*   sinf_rec      - Current counter for the solution information
*         records

      integer*4 old_srec_sinf, sinf_rec, trimlen
      logical kbit

*   var_scale       - Solution scaling varinace

      real*8 var_scale

* MOD TAH 190621: Added upper limit on loop for reading files
*     when type BACK is passed
      integer*4 loc_glb_sol  ! Local number of solutions to read

*   glb_name        - Name of the global file to be copied
 
      character*(sort_recl) glb_name
 
* MOD TAH 190621: See when type of solution
      if ( type.eq.'BACK' ) then
          ! back_type = .true.   ! Already set
          loc_glb_sol = 1
      else     ! Default stanard type run
          ! back_type = .false.  ! Already set
          loc_glb_sol = num_glb_sol
      endif

****  Copy the output global dcb buffere
      do i = 1,16
         oglb_dcb(i) = cglb_dcb(i)
      end do

*     Advance the unit number
* MOD TAH 190621: Increment by 2 to avoid clash with glb_sol_file unit
*     (only an issue when glsave run from back solution)
      cglb_dcb(1) = cglb_dcb(1) + 2
 
****  Open the sort file containing the names of all of the solutions
*     in this global solution
* MOD TAH 190621: Only read sort_file if not BACK type
      if( .not. back_type ) then 
         open(100, file = sort_file, status='old', iostat=ierr, 
     .             access='direct', recl=sort_recl+8+8)
         call report_error('IOSTAT',ierr,'open',sort_file,0,
     .                     'GW_DESCRIPTION')
 
         if( ierr.ne.0 ) then
*            No use trying to copy records so return
             RETURN
         end if
      endif 
 
***** Now loop over the file names in the sort file and copy
*     the solution descriptions from each into combined global file
 
*                               ! Initialize the output record number
      out_rec = crec_solutions
*                               ! This value will be incremented for
*                               ! each record written
*     Save the location of the solution information records initally
      sinf_rec = crec_sinf

* MOD TAH 980518: Initialize the array that says the ECC values are
*     known.

      do i = 1, max_glb_site_wrds
          gecc_known(i) = 0
      end do
      do i = 1, max_glb_svs 
           gsvi_prn(i) = 0
           gsvi_svn(i) = 0
      enddo


* MOD TAH 190621: Change loop to used loc_glb_sol instead of num_glb_sol 
      do i = 1, loc_glb_sol

* MOD TAH 190621: Only read file if not BACK type
          if( back_type ) then
             glb_name = glb_inp_file
             for_chi = for_chi_save
             bak_chi = bak_chi_save
             var_scale = 1.0         ! Not needed here
          else
             read(100,iostat=ierr,rec=i) glb_name, var_scale, 
     .                                   for_chi, bak_chi
          endif

*****     See if this file not used
          if( ichar(glb_name(1:1)).gt.128 ) then
             write(*,200) char(ichar(glb_name(1:1))-128),
     .             glb_name(2:trimlen(glb_name)), min(for_chi,999.99),
     .             min(bak_chi,999.99)
 200         format('Global ',a1,a,' Not used Chis ',2f8.2)
          else 
*            Open the global file
C            write(*,210) cglb_dcb(1), trim(glb_name)
C210         format('Opening on unit ',i5,1x,a)
             call FmpOpen( cglb_dcb,ierr,glb_name,'RO',0)
             call report_error('FmpOpen',ierr,'open',glb_name,0,
     .                      'GW_DESCRIPTION')
             if( ierr.eq.0 ) then
                 call rw_glb_header('R', ierr) 
*                MOD TAH 150825: Get the PRN to SVN information
                 call read_svinf_rec

                 call rw_names_block('R')
                 glb_inp_file = glb_name
                 call eq_name_change('NO')
             end if

*            Copy if no errors
             if( ierr.ge.0 ) then

*                Get and save the long names of the sites (the last long 
*                name will be saved).

                 call save_full_names 
 
*                Get solution records number and the number of solution
*                records
                 sol_rec = crec_solutions    !  rec_copy(12)
                 sol_num = cnum_soln_recs    !  rec_copy( 5)

*                Now read the solutions records.  If the verion number
*                is not zero, then copy the station information records 
*                as well.
 
                 do j = 1, sol_num
                     call readd( cglb_dcb,ierr,sdelete_count,128,
     .                                    len_read,  sol_rec+j-1)

*                    Save the old record number of the station information
*                    records
                     old_srec_sinf = srec_sinf
C                    srec_sinf = sinf_rec
*                    Set the station information record number to zero. Only
*                    last solution record has station information.
                     srec_sinf = 0
                     call writd( oglb_dcb,ierr,sdelete_count,128,
     .                                         out_rec)
                     call report_error('FmpWrite',ierr,'writ',
     .                           'Solution records',0,'GW_DESCRIPTION')
    
*                                         ! GET out if write error
                     if( ierr.lt.0 ) RETURN
                     out_rec = out_rec + 1

*                    Now check the version number.  Below use the valid
*                    version number range.  Also accept ver 0.05 SNX. 
*                    MOD TAH 030415: Accept any version >= 0
                     if( ((sglb_vers.ge.100 .and. sglb_vers.le.110) .or.
     .                    (sglb_vers.ge.0 ) ) .and.
     .                   j.eq.sol_num ) then
                         do k = 1, cnum_sinf
                            call readd( cglb_dcb,ierr,rec_copy,128,
     .                           len_read, crec_sinf+k-1)
   
****                        Saves the information about which sites have
*                           been used.

                            if( kbit(guse_site,ltog_sites(k*2-1)) )
     .                      call acc_sinf( k*2-1, rec_copy(1))
*                           The last part of the record may be blank at the
*                           end of the block.  Check to see.
* MOD TAH 980328: Corrected bug with ltog_site lookup for k*2-1 instead
*                           of k*2 as current shown.
                            if( 2*k.le.cnum_sites ) then
                                 if( kbit(guse_site,ltog_sites(k*2)) )
     .                           call acc_sinf( k*2, rec_copy(65))
                            end if

****                        Save the record number where we will write
*                           these records in the combined global file.
C                           call writd( cglb_dcb,ierr,rec_copy,128,
C    .                                  sinf_rec )       
C                           sinf_rec = sinf_rec + 1
                         end do

* MOD TAH 050622: Accumulate the information about satellites.
                         do k = 1,cnum_svinf
                            call readd( cglb_dcb,ierr,rec_copy,128,
     .                           len_read, crec_svinf+k-1)
*                           Accumulation code handles the end of record
*                           (0 for PRN and so we just check all records
*                           even if some do not exist)
                            call acc_svinf( (k-1)*4+1, rec_copy(1))
                            call acc_svinf( (k-1)*4+2, rec_copy(33))
                            call acc_svinf( (k-1)*4+3, rec_copy(65))
                            call acc_svinf( (k-1)*4+4, rec_copy(97))
                         end do                            

                        end if
                 end do
*                     ! No error opening input file
              end if
          end if       ! File actually used
*                     ! Looping over all of the global files in this
      end do
*                     ! solution

*     Now create and write out the solution record for the combined
*     global file itself.
*     Copy the DCB buffer back into place
      do i = 1,16
         cglb_dcb(i) = oglb_dcb(i)
      end do
      call rw_glb_header('R', ierr)
      call gw_soln_rec( sinf_rec )
 
***** Thats all
* MOD TAH 190621: Unit 100 only opened here is not BACK type
      if( .not. back_type ) close(100)
      return
      end

CTITLE GW_SOLN_REC 
    
      subroutine gw_soln_rec( sinf_rec ) 

      implicit none  
 
*     This rouitine creates the new soln record for the combined 
*     global file if there was more than sone input solution.
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/globk_common.h'

* PASSED VARIABLES
* sinf_rec - Current record number to write the station info
*                 rec. (Since previous station records were read
*                 value at the end of the writing these is correct)
* svinf_rec - Current record for satellite records

       integer*4 sinf_rec, svinf_rec

*   i,j,k       - Loop counters
*   rec_copy(128)  - Record for the station infomatation records
*                   (two stations per record)

      integer*4 i, j, k, l, rec_copy(128), ierr
      real*8 sectag
      logical kbit, finished
 

****  Loop over the variables solution record variables setting the
*     new values
      do i = 1, max_edit_types
         sdelete_count(i) = gdelete_count(i)
      end do

      sversion = -999
      stai_utc = gtai_utc

* MOD TAH 190628: If back type save number of obs from this global file
      if( back_type ) then  
         snum_obs = cnum_obs
      else
         snum_obs = gnum_obs
      endif 
      snum_sites = cnum_sites
 
      sepoch_end   = gepoch_end
      if( sort_direction.eq.1 ) then
          sepoch       = gepoch_end 
      else
          sepoch       = gepoch_start 
      end if
      sepoch_start = gepoch_start
      ssvs_epoch   = csvs_epoch

      sectag = grun_time(6)
      call ymdhms_to_jd( grun_time, sectag, srun_time)
 
      do i = 1,2
          snut_ang_apr(i,1) = gnut_ang_apr(i)
          snut_ang_apr(i,2) = gnut_ang_apr(i+2)
      end do
 
      sut1_apr(1) = gut1_apr(1)
*     Convert UT1 rate back to mas/day
      sut1_apr(2) = gut1_apr(2)
 
      do i = 1,2
          swob_apr(i,1) = gwob_apr(i)
          swob_apr(i,2) = gwob_apr(i+2)
      end do

****  Add the new values

      sglb_vers = glbf_version
      scons_type = gcons_type
      ssys_type  = gsys_type
      call get_institute( sowner, default_institute)
      call get_institute( screator, default_institute)
      sprog_gen = '=GLK'
      sanal_type = canal_type
      skalobs_file = glb_out_file
      sexpt_title = gdescription

***** Save infor about where the station information will be
      srec_sinf = crec_sinf
      snum_sinf = (cnum_sites-1)/2+1
 
***** Now write out the header records to the file
 
      call writd( cglb_dcb, ierr, sdelete_count, 128, crec_comb_soln)
      call report_error('FmpWrite',ierr,'writ','combined soln rec',1,
     .                  'GW_soln_rec')

***** Now create and write the station information records
*     (Do two at a time since there are two stations per physical
*      record).
      sinf_rec = crec_sinf
      i = 0
      k = 0
      finished = .false.
      do while ( .not.finished ) 
         i = i + 1
         if( i.eq.gnum_sites ) finished = .true.
         if( kbit(guse_site,i) ) then
             call cr_sinf_rec(i, rec_copy(k*64+1))
             k = k + 1
             if( k.eq.2 ) then
                 call writd(cglb_dcb, ierr, rec_copy, 128, sinf_rec)
                 sinf_rec = sinf_rec + 1
                 k = 0
             end if
         end if
      end do
      if( k.eq.1 ) then
          call writd(cglb_dcb, ierr, rec_copy, 128, sinf_rec) 
      end if
             
***** Now create and write the satellite information records
*     (Do four at a time since there are four satellites per physical
*      record).
      svinf_rec = crec_svinf
      i = 0
      k = 0
      finished = .false.
      if ( gnum_svs.eq.0 ) finished = .true.
      do while ( .not.finished ) 
         i = i + 1
         if( i.eq.gnum_svs ) finished = .true.
         l = i
         k = k + 1
         call cr_svinf_rec(l, rec_copy((k-1)*32+1))
         if( k.eq.4 ) then
              call writd(cglb_dcb, ierr, rec_copy, 128, svinf_rec)
              svinf_rec = svinf_rec + 1
              k = 0
              do j = 1,128     ! Clear the next record
                 rec_copy(j) = 0
              end do
         end if
      end do
      if( k.ne.0 ) then
          call writd(cglb_dcb, ierr, rec_copy, 128, svinf_rec) 
      end if


***** Thats all
      return
      end
 
CTITLE CR_SINF_REC

      subroutine cr_sinf_rec( i, record )

      implicit none 

*     Routine to create the solution description and write the record
*     to the global file.
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'
      include '../includes/globk_common.h'
 
* PASSED VALUES

* i  -  Station number
* record(64)  - Record for the Station information

      integer*4 i, j, record(64)
      integer*4 ne   ! Earthquake number

* LOCAL VARIABLES
* MOD TAH 190627: Added new function to get ss_lnum value
      integer*4 gtol_map   ! Function to map global site number
                     ! to local number number to see if it is used.

*
*     Start saving the information
      sscode = gsite_names(i)
      ssdata_st = gdata_st(i) 
      if( ssdata_st.eq.0.d0 ) ssdata_st = gepoch_start
      ssdata_en = gdata_en(i) 
      if( ssdata_en.eq.0.d0 ) ssdata_en = gepoch_end
      ssrecv_st = grecv_st(i) 
      if( ssrecv_st.eq.0.d0 ) ssrecv_st = gepoch_start
      ssrecv_en = grecv_en(i) 
      if( ssrecv_en.eq.0.d0 ) ssrecv_en = gepoch_end
      ssante_st = gante_st(i) 
      if( ssante_st.eq.0.d0 ) ssante_st = gepoch_start
      ssante_en = gante_en(i) 
      if( ssante_en.eq.0.d0 ) ssante_en = gepoch_end
      do j = 1,3
         ssarp_ecc(j) = garp_ecc(j,i)
         ssL1a_ecc(j) = gL1a_ecc(j,i)
         ssL2a_ecc(j) = gL2a_ecc(j,i)
         satmload(j) =  gatmload(j,i) 
         shydload(j) =  ghydload(j,i) 
      end do
* MOD TAH 200205: Save the antdaz antenna azimith value
      ssantdaz   = gantdaz(i)
      sselev_cut = gelev_cut(i)
      ssnum_zen  = gnum_zen(i)
      ssant_mod  = gant_mod(i)
      sti_antmod = gant_mod(i)
      ssante_sn  = gante_sn(i)
      ssante_ty  = gante_ty(i)
      ssradome_ty = gradome_ty(i)
      ssrecv_sn  = grecv_sn(i)
      ssrecv_ty  = grecv_ty(i)
      ssrecv_ty_end = grecv_ty(i)(17:20)
      ssrecv_fw  = grecv_fw(i)

* MOD TAH 030615: Save the JD and time constain for log fits at this
*     station.  Fid the qarthquake for this site
      ne = 0
      do j = 1, num_eq
         if( gsite_names(i)(7:8).eq.eq_codes(j)(1:2) ) ne = j
      end do

*     If ne is zero, there is no earthquake associated with this
*     site.
      if( ne.gt.0 ) then
          slog_ep = eq_epoch(ne)
          slog_tau = eq_log_tau(ne)
      else
          slog_ep = 0.d0
          slog_tau = 0.0
      end if

* MOD TAH 190627: See if this site is actually used when GLX file created
*     from glbak run
      if( back_type ) then   ! See if site is local
          ss_lnum = gtol_map( i, ltog_sites, gnum_sites )
      else
          ss_lnum = i
      endif
           

****  Now copy the values to the data record
      call wmov( ssdata_st, 1, record, 1, 64 )
      return
      end

CTITLE CR_SVINF_REC

      subroutine cr_svinf_rec( i, record )

      implicit none 

*     Routine to create the satellite information records for output
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/svinf_def.h'
      include '../includes/globk_common.h'
 
* PASSED VALUES

* i  -  Satellite number
* record(32)  - Record for the Station information

      integer*4 i, j, record(32)

* LOCAL VARIABLES
      integer*4 gtol_map   ! Function to return local satellite 
            ! number from global one.  Only used in BACK GLX files.

*
*     Start saving the information
      svi_prn = gsvi_prn(i)
      svi_svn = gsvi_svn(i)
      svi_block = gsvi_block(i)
* MOD TAH 190627: If this is back solution assign svi_lnum 
      if( back_type ) then   ! Creating GLX file in glbak
         svi_lnum = gtol_map( i, ltog_svs, gnum_svs )
      else      ! Just save infoex
         svi_lnum = i
      endif

      do j = 1,3
         svi_antpos(j,1) =  gsvi_antpos(j,1,i)
         svi_antpos(j,2) =  gsvi_antpos(j,2,i)
      end do
      svi_launch = gsvi_launch(i)
      svi_antmod = gsvi_antmod(i)
      svi_ocode = gsvi_ocode(i)
      do j = 1,9    ! Make sure we clear the rest of the record

         svi_fill(j) = 0
      end do

****  Now copy the values to the data record
      call wmov( svi_prn, 1, record, 1, 32 )

      return
      end
 
CTITLE ACC_SINF    

      subroutine acc_sinf( is, record )

      implicit none 

*     Routine to create the solution description and write the record
*     to the global file.
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'
      include '../includes/globk_common.h'
 
* PASSED VALUES

* is  -  Station number
* record(64)  - Record for the Station information

      integer*4 is, i, j,  record(64),  indx
      logical kbit

*
*     Find out which station this referrs to
      indx = 0
      call wmov( record, 1,  ssdata_st, 1, 64 )
      if( sscode(1:1).eq.char(0) ) RETURN
      i = ltog_sites(is)
C     call get_cmd( sscode, gsite_names, gnum_sites, i, indx)
      if(  i.gt.0 ) then

*         See if values have been set yet.  
* MOD TAH 190628: Also reset values if back_type is true so that
*         values for current experiment only are saved.
          if( gdata_st(i).eq.0.d0.or.back_type ) gdata_st(i) = ssdata_st 
          if( gdata_en(i).eq.0.d0.or.back_type ) gdata_en(i) = ssdata_en 
          if( grecv_st(i).eq.0.d0.or.back_type ) grecv_st(i) = ssrecv_st 
          if( grecv_en(i).eq.0.d0.or.back_type ) grecv_en(i) = ssrecv_en 
          if( gante_st(i).eq.0.d0.or.back_type ) gante_st(i) = ssante_st 
          if( gante_en(i).eq.0.d0.or.back_type ) gante_en(i) = ssante_en 

*         Now see limits shoulb be changed
          if( gdata_st(i).gt. ssdata_st ) then
              gdata_st(i) = ssdata_st 
          end if
          if( gdata_en(i).lt. ssdata_en ) then
              gdata_en(i) = ssdata_en 
          end if
          if( grecv_st(i).gt. ssrecv_st ) then
              grecv_st(i) = ssrecv_st 
          end if
          if( grecv_en(i).lt. ssrecv_en ) then
              grecv_en(i) = ssrecv_en 
          end if
          if( gante_st(i).gt. ssante_st ) then
              gante_st(i) = ssante_en 
          end if
          if( gante_en(i).lt. ssante_en ) then
              gante_en(i) = ssante_en 
          end if

*         If the site has been renamed, check to see if we
*         changed the eccentricity of the mark (if so update
*         the values)
          if( kbit(rn_name_changed,i) ) then
              call update_arp(i)
          end if
          
          if ( .not. kbit(gecc_known,i) ) then
             do j = 1,3
                garp_ecc(j,i) = ssarp_ecc(j) 
                gL1a_ecc(j,i) = ssL1a_ecc(j) 
                gL2a_ecc(j,i) = ssL2a_ecc(j)
                gatmload(j,i) = satmload(j)
                ghydload(j,i) = shydload(j)
             end do
* MOD TAH 200205: Save the antdaz antenna azimith value
             gantdaz(i)   = ssantdaz
             gelev_cut(i) = sselev_cut 
             gnum_zen(i)  = ssnum_zen  
             gant_mod(i)  = ssant_mod 
             if( ichar(sti_antmod(1:1)).ne.0 ) then
                gant_mod(i) = sti_antmod
             end if 
             gante_sn(i) = ssante_sn  
             gante_ty(i) = ssante_ty 
             gradome_ty(i) = ssradome_ty
             grecv_sn(i) = ssrecv_sn 
             grecv_ty(i) = ssrecv_ty
             if( ichar(ssrecv_ty_end(1:1)).ne.0 ) then
                 grecv_ty(i) = ssrecv_ty // ssrecv_ty_end
             end if
             grecv_fw(i) = ssrecv_fw  
*
             call sbit(gecc_known,i,1)
          else

* MOD TAH 980518: Check to see if values are OK and
*            update if they do not appear to be complete.
* MOD TAH 101030: Test direction.  If going forward in time
*            update with latest values; if going backwards
*            only update if we seem to mussing data (original
*            code).
             if( sort_direction.lt.0 ) then
                if( gant_mod(i)(1:4).eq.'----'.or. 
     .              ichar(gant_mod(i)(1:1)).eq.0 ) then
                    gant_mod(i) = ssant_mod
                    if( ichar(sti_antmod(1:1)).ne.0 ) 
     .                  gant_mod(i) = sti_antmod
                end if
                if( gante_sn(i)(1:4).eq.'----'.or. 
     .              ichar(gante_sn(i)(1:1)).eq.0 )
     .          gante_sn(i) = ssante_sn  
                if( gante_ty(i)(1:4).eq.'----'.or. 
     .              ichar(gante_ty(i)(1:1)).eq.0 )
     .          gante_ty(i) = ssante_ty
                if( gradome_ty(i)(1:4).eq.'----'.or. 
     .              ichar(gradome_ty(i)(1:1)).eq.0 )
     .          gradome_ty(i) = ssradome_ty
                if( grecv_sn(i)(1:4).eq.'----'.or. 
     .              ichar(grecv_sn(i)(1:1)).eq.0 )
     .          grecv_sn(i) = ssrecv_sn 
* MOD TAH 200803: Removed test on the extended part
*               of name (17-20) because it may naturally
*               be blank and not an indication to update
                if( grecv_ty(i)(1:4).eq.'----'.or. 
     .              ichar(grecv_ty(i)(1:1)).eq.0 ) then  ! .or.
!    .              grecv_ty(i)(17:17).eq.' ' ) then
                    grecv_ty(i) = ssrecv_ty
                    if( ichar(ssrecv_ty_end(1:1)).ne.0 )
     .                  grecv_ty(i) = ssrecv_ty //  ssrecv_ty_end
                endif
             else   !  Running forward so use the latest information
* MOD TAH 101030: Check the direction of data so that that latest 
*               information
                if( ssant_mod(1:4).ne.'----' .and.
     .              ichar(ssant_mod(1:1)).ne.0 ) 
     .              gant_mod(i) = ssant_mod
                if( ichar(sti_antmod(1:1)).ne.0 ) 
     .              gant_mod(i) = sti_antmod
                if( ssante_sn(1:4).ne.'----'.and. 
     .              ichar(ssante_sn(1:1)).ne.0 )
     .          gante_sn(i) = ssante_sn  
                if( ssante_ty(1:4).ne.'----'.and. 
     .              ichar(ssante_ty(1:1)).ne.0 )
     .          gante_ty(i) = ssante_ty
                if( ssradome_ty(1:4).ne.'----'.and. 
     .              ichar(ssradome_ty(1:1)).ne.0 )
     .          gradome_ty(i) = ssradome_ty
                if( ssrecv_sn(1:4).ne.'----'.and. 
     .              ichar(ssrecv_sn(1:1)).ne.0 )
     .          grecv_sn(i) = ssrecv_sn 
                if( ssrecv_ty(1:4).ne.'----'.and. 
     .              ichar(ssrecv_ty(1:1)).ne.0 ) then
                    grecv_ty(i) = ssrecv_ty
                    if( ichar(ssrecv_ty_end(1:1)).ne.0 )
     .                  grecv_ty(i) = ssrecv_ty //  ssrecv_ty_end
                endif
             end if
             if( grecv_fw(i)(1:4).eq.'----'.or. 
     .           ichar(grecv_fw(i)(1:1)).eq.0 )
     .       grecv_fw(i) = ssrecv_fw  
          endif
      else
          write(*,300) sscode
 300      format('**WARNING** Code ',a,' from SINF record not',
     .           ' found in global site list')
      end if

***** Thats all
      return
      end

CTITLE ACC_SVINF    

      subroutine acc_svinf( is, record )

      implicit none 

*     Routine to accumulate the satellite information records
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/svinf_def.h'
      include '../includes/globk_common.h'
 
* PASSED VALUES

* is  -  satellite number
* record(32)  - Record for the satellite information

      integer*4 is, i, j,  record(32),  indx
*
*     Find out which station this referrs to
      indx = 0
      i = ltog_svs(is)
      call wmov( record, 1,  svi_prn, 1, 32 )
      if( svi_prn.eq.0  ) RETURN

      if(  i.gt.0 ) then
          gsvi_prn(i)    = svi_prn
          gsvi_svn(i)    = svi_svn
          gsvi_block(i)  = svi_block
* MOD TAH 101030: Only update this information if running forward
*         in time and new informaion is OK; or when running backward
*         if current is valid
          if( sort_direction.gt.0 ) then  ! Forward; update if new is OK
             if( svi_antmod(1:4).ne.'UNKN' .and. 
     .           ichar(svi_antmod(1:1)).ne.0 ) then
                gsvi_antmod(i) = svi_antmod
                gsvi_ocode(i)  = svi_ocode
                do j = 1, 3
                   gsvi_antpos(j,1,i) = svi_antpos(j,1)
                   gsvi_antpos(j,2,i) = svi_antpos(j,2)
                end do
                gsvi_launch(i) = svi_launch
             end if
          else        ! Only update if current is not valid
             if( gsvi_antmod(i)(1:4).eq.'UNKN' .or. 
     .           ichar(gsvi_antmod(i)(1:1)).eq.0 ) then
                gsvi_antmod(i) = svi_antmod
                gsvi_ocode(i)  = svi_ocode
                do j = 1, 3
                   gsvi_antpos(j,1,i) = svi_antpos(j,1)
                   gsvi_antpos(j,2,i) = svi_antpos(j,2)
                end do
                gsvi_launch(i) = svi_launch
             end if
          end if
      end if

***** Thats all
      return
      end

CTITLE SAVE_FULL_NAMES

      subroutine save_full_names

      implicit none 

*     Routine to read the full names and save them

      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/globk_common.h'
 
      include '../includes/glb_hdr_def.h'
 
*   i,j,k       - Loop counters
*   ierr        - File error flag
*   full_list(1)   - Integer alias for cnames_list so that we can
*               - write to a file.
 
      integer*4 i, ierr, full_list(max_glb_sites*8), len_read,
     .          k, gtos(max_glb_sites)
      logical kbit
 
*   cnames_list(max_names)  - List of site and sources names
 
      character*32 cnames_list(max_glb_sites)
 
      equivalence ( cnames_list, full_list )
 
***** Make the list of site and source names

      call readd(cglb_dcb,ierr, full_list,128*cnum_full ,len_read,
     .           crec_full)
      call report_error('FmpRead',ierr,'read','Full  block',0,
     .                  'save_full_names')

****  Get the mapping from the global site numbers to saved
*     site numbers in case we have not used a site in the solution.
      k = 0 
      do i = 1, gnum_sites
         if( kbit(guse_site,i) ) then
             k = k + 1
             gtos(i) = k
         else
             gtos(i) = 0
         end if 
      end do

      do i = 1, cnum_sites
*         Only save the long site name if we are still using
*         this site.
          if( gtos(ltog_sites(i)).gt.0 ) then
              gsite_full(gtos(ltog_sites(i))) = cnames_list(i)
          end if
      end do
 
***** Thats all
      return
      end
 
 
CTITLE UPDATE_ARP  

      subroutine update_arp(is)          

      implicit none 

*     Routine to update the ARP eccentricity based on the positons
*     changes in the rename command. The site number passed (is) is
*     the global site number.
 
      include '../includes/kalman_param.h'
      include '../includes/glsave_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'
      include '../includes/globk_common.h'
 
* PASSED VALUES

*  is   - Global site number

      integer*4 is

* LOCAL VARIABLES

*   i,j  - loop counters

      integer*4 i, j

*   loc_coord(3)  - Local coordinates
*   rot_matrix(3,3) - Rotation from XYZ to NEU
*   dNEU(3)       - Chnage in NEU coordinates

      real*8 loc_coord(3) , rot_matrix(3,3), dNEU(3)

****  Find the correct rename entry.
      do i = 1, num_renames
         if( rn_codes(2,i).eq.gsite_names(is) .and.
     .       ssdata_st    .ge.rn_times(1,i)   .and.
     .       ssdata_en    .le.rn_times(2,i)  ) then

*            OK found the correct entry.  Now rotate
*            XYZ in rn_dpos to NEU and remove from
*            ssarp_ecc entry.
             if( rn_dpos(1,i)**2 + rn_dpos(2,i)**2 + 
     .           rn_dpos(3,i)**2 .gt.0 ) then
                 call rotate_geod(rn_dpos(1,i), dNEU,rn_types(i),'NEU',
     .                    apr_val_site(1,1,is), loc_coord, rot_matrix)

                 write(*,120) gsite_names(is), dNEU
 120             format('Changing ARP Eccentricity of ',a,' by ',
     .                  3F10.4,' m NEU')

                 do j = 1,3
                    ssarp_ecc(j) = ssarp_ecc(j) - dNEU(j)
*                   Do not need to do for the L1/L2 to apr point. 
C                   ssL1a_ecc(j) = ssL1a_ecc(j) - dNEU(j)
C                   ssL2a_ecc(j) = ssL2a_ecc(j) - dNEU(j)
                 end do
             end if
         end if
      end do

****  Thats all
      return
      end 

