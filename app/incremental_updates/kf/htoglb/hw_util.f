CTITLE QW_GLB_HEADER
 
      subroutine qw_glb_header
 
      implicit none 

 
*     Routine to create and write the global header record
*     for the global solution file.  This routine will also purge
*     and create the file to which the solution will be written.
*
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
*   AddressOf   - Gets the address of a variable (used to find
*               - out how large the global header records are)
*   i,j         - Loop counter
*   ierr        - File error flag
 
* MOD TAH 190511: Changed to I*8 for 64-bit memory
      integer*8 AddressOf
      integer*4 i, ierr
 
*   jsize       - Size of the file in blocks
 
      integer*4 jsize, trimlen
 
*   Full_file_name  - Full name of the global file (with size
*                   - Information)
 
      character*128 Full_file_name

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      data I8 / 1 /
 
***** First compute the record boundaries of all the blocks for the
*     global file
 
      cnum_header = ( AddressOf(last_glb_header_wrd) -
     .                AddressOf(glb_header) )/(4*128)   + 1
      cnum_names  = ( qnum_sites+qnum_svs-1 )*2/128   + 1
*                         ! One record per solution description
      cnum_full   = ( qnum_sites*8 - 1 )/ 128 + 1
      cnum_soln_recs = qnum_soln_recs
      
*     Compute number of records needed for the site information 
*     for each solution.  The site information is allocated 
*     256 bytes per station.  The solution records themselves
*     contain the record number of the station information.  
*     Initialize number of site info records with the number
*     needed for this combined solution.
*     We only save these for the final solution
      cnum_sinf = ((qnum_sites-1)/2+1)
* MOD TAH 050622: Compute number of records for satellites
      cnum_svinf = ((qnum_svs-1)/4+1)
      cnum_comb = qnum_comb
      
c     cnum_par_types = (qnum_par_codes  -1)/128 + 1
c     cnum_apr_types = (qnum_apr_codes  -1)/128 + 1
      cnum_par_types = (qnum_parn       -1)/128 + 1
      cnum_apr_types = (qnum_apr_codes  -1)/128 + 1
      cnum_apr_vals  = (2*qnum_apr_codes-1)/128 + 1
* MOD TAH 190603: Make calculation does not overflow. (Introduce
*     I8 to make integer*8
      cnum_par_vals  = (2*qnum_parn*(qnum_parn+I8)-1)/128 + 1
      
* MOD TAH 981020: Compute number of records needed for the parameter
*     epoch values
      cent_par_ep = qent_par_codes
      cent_apr_ep = qent_apr_codes
      cnum_par_ep = (cent_par_ep - 1)/64 + 1
      cnum_apr_ep = (cent_apr_ep - 1)/64 + 1
      

****  Save the number of entries on the apriori constraints.  This
*     is num_glb_parn + off diagonals for NEU constraints (3 per one)
      cnum_acvc      = 4*(qnum_apr_diag+6*qnum_apr_cov+
     .                    21*qnum_apr_svs)
      cnum_acvc_recs = (cnum_acvc - 1)/128 + 1

      crec_names     = 1              + cnum_header
      crec_full      = crec_names     + cnum_names
      crec_solutions = crec_full      + cnum_full 
      crec_sinf      = crec_solutions + cnum_soln_recs
*     Save the record for the solution record that referrs to the
*     combined global file.
      crec_comb_soln = crec_sinf      - 1

* MOD TAH 050622: Add in count for satellite records      
      crec_svinf     = crec_sinf      + cnum_sinf
      crec_par_types = crec_svinf     + cnum_svinf
      crec_apr_types = crec_par_types + cnum_par_types
      
* MOD TAH 981020: Add the records for the multi-epoch polar motion
      crec_par_ep    = crec_apr_types + cnum_apr_types
      crec_apr_ep    = crec_par_ep    + cnum_par_ep
      
      crec_apr_vals  = crec_apr_ep    + cnum_apr_ep
      crec_par_vals  = crec_apr_vals  + cnum_apr_vals
      
*     Compute record number for start of apriori constraint block.
*     The apriori constraint values and the number of records will
*     computed later as the records are written.
      crec_acvc      = crec_par_vals +  cnum_par_vals
 
      jsize  = crec_acvc      + cnum_acvc_recs - 1

***** Purge any file with the same name as we will create
 
      call clean_file( glb_file )
 
      call FullFileName( glb_file, 1, jsize, 128, full_file_name)
 
      call FmpOpen( cglb_dcb, ierr, Full_File_Name, 'wc', 1)
      call report_error('FmpOpen',ierr,'creat',Full_file_name,1,
     .                  'qw_glb_header')
 
***** Now set up the rest of the header blocks
 
      cfile_type   = 'GLOBAL'
      cexpt_title  = qexpt_title
      cnum_apr_codes = qnum_apr_codes
      cnum_parn    = qnum_parn
      cnum_sites   = qnum_sites
      cnum_sources = 0          
      cnum_svs     = qnum_svs

      do i = 1, 7
          crun_time(i) = qrun_time(i)
      end do
 
      ctai_utc     = qtai_utc
      cgpst_utc    = qgpst_utc
      
      do i = 1, max_edit_types
          cdelete_count(i) = 0               
      end do
 
      cnum_obs = snum_obs
 
      cepoch_end   = qend_epoch
      cepoch_expt  = sepoch
      cepoch_start = qstart_epoch
      csvs_epoch   = ssvs_epoch
 
      do i = 1,3
          cetd_apr(i) = 0          
      end do

      do i = 1,2
          cnut_ang_apr(i) = qnut_ang_apr(i,1)
          cnut_ang_apr(i+2) = qnut_ang_apr(i,2)*365.25
      end do

      cut1_apr(1) = qut1_apr(1)
      cut1_apr(2) = qut1_apr(2)
 
      do i = 1,2
          cwob_apr(i) = qwob_apr(i,1)
          cwob_apr(i+2) = qwob_apr(i,2)
      end do
 

****  Add the new values

      cglb_vers = glbf_version
      ccons_type = qscons_type(cnum_soln_recs)
      csys_type  = qssys_type(cnum_soln_recs)
      cowner     = qsowner(cnum_soln_recs)
      call get_institute(ccreator, default_institute )   
      cprog_gen  = '=GLK'
      canal_type = qsanal_type(cnum_soln_recs)
      if( trimlen(sgtime).eq.0 .or. ichar(sgtime(1:1)).eq.0 ) then
          sgtime = 'GPST'
      end if
      if( trimlen(sgframe).eq.0 .or. ichar(sgframe(1:1)).eq.0 ) then
          sgframe = 'B1950'
      end if
      if( trimlen(sgprec).eq.0 .or. ichar(sgprec(1:1)).eq.0 ) then
          sgprec = 'IAU68'
      end if
      if( trimlen(sgsrpmod).eq.0 .or. ichar(sgsrpmod(1:1)).eq.0 ) then
          sgsrpmod = 'SPHRC'
      end if
      cgtime   = sgtime
      cgframe  = sgframe
      cgprec   = sgprec
      cgsrpmod = sgsrpmod  
* MOD TAH 981112:  Save the gamit models used in the header definition
*     block.
      cgamit_mod = sgamit_mod
      cload_mod  = sload_mod
      
***** Now write out the header records to the file 
      call writd( cglb_dcb, ierr, glb_header, 128*cnum_header, 1)
      call report_error('FmpWrite',ierr,'writ','global header',1,
     .                  'qw_glb_header')
 
***** Thats all
      return
      end
 
CTITLE QW_GLB_NAMES
 
      subroutine qw_glb_names
 
      implicit none 

 
*     Routine to create the list of site and source names and save
*     these in the global file.
*
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
 
*   i,j,k       - Loop counters
*   ierr        - File error flag
*   names_list(1)   - Integer alias for cnames_list so that we can
*               - write to a file.
*   lname       - Length of site name
*   trimlem     - Length of string
 
      integer*4 i,j, ierr, names_list(1), trimlen, lname

      integer*4 jerr ! OCC decode error
     .,         occ  ! occupancy code

      character*1 cocc  ! One character occ code
     .,         pt      ! Pt code for site  

      character*4 code  ! 4-char site code
 
*   cnames_list(1)  - List of site and sources names (saved in
*               - scr_common temporarily
 
      character*8 cnames_list(max_glb_sites+max_glb_svs)
      character*8 qorg_name 
 
      equivalence ( cnames_list, names_list )
 
 
      common cnames_list
 
***** Make the list of site and source names

* MOD TAH 911119: Modded names to add a _GPS  to end of name.
*     This is to avoid confusion between VLBI site names and GPS
*     site_names
 
      write(*,110) cnum_sites, hfile(1:trimlen(hfile))
 110  format(/,'There are ',i5,' sites in ',a,/,
     .         'Number    Name  Code     | Long name ')

***   See if we need to add qfull names to sites where they were
*     added during decoding.
      do i = 1, cnum_sites
         if( trimlen(qfull_names(i)).eq.0 ) then
*           Scan sites to find original name
            do j = 1, i-1
               if( qsite_names(i)(1:4).eq.qsite_names(j)(1:4) .and.
     .             qsite_names(i)(8:8).eq.qsite_names(j)(8:8) ) then
                  qfull_names(i) = qfull_names(j)
                  exit
               end if
            end do
         end if
      end do

****  Now generate globk type names

      do i = 1, cnum_sites
          qorg_name = qsite_names(i) 
          lname = trimlen(qsite_names(i))
          if( lname.lt.5 .and. hfile_type(1:3).ne.'TER' .and.
     .        hfile_type(1:3).ne.'VLB' ) then
              qsite_names(i)(lname+1:) = '_' // hfile_type
          end if
* MOD TAH 100511: Use igs_ptname to see if orignal conversion of 
*         name should be used (igs_ptname false) or if the direct
*         mapping should be used.
          if( hfile_type(1:3).eq.'SNX' .and. .not. igs_ptname ) then
              if( qsite_names(i)(7:8).eq.'1A' ) then
                  if( qfull_names(i)(32:32).eq.'R' ) then
                     qsite_names(i)(5:) = '_VLB'
                  else if( qfull_names(i)(32:32).eq.'S' ) then
                     qsite_names(i)(5:) = '_SLR'
                  else if( qfull_names(i)(32:32).eq.'D' ) then
                     qsite_names(i)(5:) = '_DOR'
                  else if( qfull_names(i)(32:32).eq.'E' ) then
                     qsite_names(i)(5:) = '_EDM'
                  else if( qfull_names(i)(32:32).eq.'C' ) then
*                     See if DORIS is in station name 
                      if( index(qfull_names(i),'DORIS').gt.0 ) then
                          qsite_names(i)(5:) = '_DOR'
                      else
                          qsite_names(i)(5:) = '_COM'
                      end if                 
                  else
                     qsite_names(i)(5:) = '_GPS'
                  end if
*                 MOD TAH 090504: Added extra resolution to name for sites
*                 with OCC values greater than 9
              elseif( qsite_names(i)(6:6).eq.'0' ) then
                  if( qfull_names(i)(32:32).eq.'R' ) then
                     qsite_names(i)(5:6) = '_V'
                  else if( qfull_names(i)(32:32).eq.'S' ) then
                     qsite_names(i)(5:6) = '_S'
                  else if( qfull_names(i)(32:32).eq.'D' ) then
                     qsite_names(i)(5:6) = '_D'
                  else
                     qsite_names(i)(5:6) = '_G'
                  end if
              end if
          end if
* MOD TAH 210308: Patch to fix bad CWU sinex files where names SITE/ID 
*         block are not written correctly (wrong columns so non _GPS extents
*         can be added.
          if( qowner(1:3).eq.'CWU' ) qsite_names(i)(5:) = '_GPS'

          if( hfile_type(1:3).eq.'SNX' .and. igs_ptname ) then
* MOD TAH 100511: Use a direct mapping to names
*             Split the occ and pt code from qsite_name.
              code = qsite_names(i)(1:4)
              read(qsite_names(i)(5:7),*,iostat=jerr) occ  
              call report_error('IOSTAT',jerr,'Occ decod',
     .             qsite_names(i),0)
              if( jerr.ne.0 ) occ = 35
              pt   = qsite_names(i)(8:8)
*             Now form name
* MOD SCM 100513: There is a bug in this code if >= 17 renames occur in the SNX file since the "G" rename will get created twice!
              if( occ.eq.1 ) then
                 cocc = 'G'
              elseif( occ.le.9 ) then
                 write(cocc,'(I1)') occ
              elseif( occ.le.35 ) then
                 cocc = char(65+(occ-10))
              else
*                Report error:
                 write(*,140) i, qsite_names(i), qfull_names(i)
 140             format('ERROR: Site ',i4,' Names ',a,1x,a,
     .                  'has too large an OCC code cannot convert',/,
     .                  '       Decode without the -s option')
                 stop 'OCC code too large in SNX name'
              end if
*             Now do character 7
              if( pt.eq.'A' ) then
                  qsite_names(i)(5:8) = '_' // cocc(1:1) // 'PS'
              else
                  qsite_names(i)(5:8) = '_' // cocc(1:1) // pt(1:1) //
     .                                  'S'
              end if
          end if

*         For Terrestrial sites use the first 8 characters of the
*         full name
          if( hfile_type(1:3).eq.'TER' ) then
              qsite_names(i) = qfull_names(i)(1:8)
          end if
          cnames_list(i) = qsite_names(i)
          if( ichar(qfull_names(i)(1:1)).eq. 0 ) then
              qfull_names(i) = qsite_names(i) // '_NoLongName'
          endif 
          write(*,160) i, qsite_names(i), qorg_name, qfull_names(i)
 160      format(i5,1x,a8,2x,a8,' | ',a32) 
      end do

* MOD TAH 200728: Added use_site line to be used to remove sites
*     on days that are not center day.
      if( qowner(1:3).eq.'COD' ) then
         write(*,'(a,$)') 'CODE use_site clear'
         j = 0
         do i = 1, cnum_sites
            if( .not. rm_site(i) ) then
               j = j + 1
               if( mod(j-1,5).eq.0  ) then
                  write(*,210) qsite_names(i)
 210              format(/,'CODE use_site ',a,1x,$)
               else
                  write(*,215) qsite_names(i)
 215              format(1x,a,1x,$)
               endif
            endif
         enddo
         write(*,'(1x)')
      endif

 
      do i = 1, cnum_svs
          cnames_list(i+cnum_sites) = qsvs_names(i)
      end do
 
***** Now write out record(s)
 
      call writd(cglb_dcb,ierr, names_list,128*cnum_names,crec_names)
      call report_error('FmpWrite',ierr,'writ','Names block',1,
     .                  'qw_glb_Names')
 
***** Thats all
      return
      end

CTITLE QW_GLB_FULL  

      subroutine qw_glb_full

      implicit none 


*     Routine to create the full names for sites and save
*     these in the global file.
*

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'

*   i,j,k       - Loop counters
*   ierr        - File error flag
*   full_list(1)   - Integer alias for cnames_list so that we can
*               - write to a file.
*   trimlen     - Length of string

      integer*4 i, ierr, full_list(max_glb_sites*8)

*   cfull_list(1)  - List of site and sources names (saved in
*               - scr_common temporarily

      character*32 cfull_list(max_glb_sites)

      equivalence ( cfull_list, full_list )


      common cfull_list

***** Make the list of site and source names

      do i = 1, cnum_full*128
         full_list(i) = 8224*65536 + 8224
      end do

      do i = 1, cnum_sites
          cfull_list(i) = qfull_names(i)
      end do

***** Now write out record(s)

      call writd(cglb_dcb,ierr, full_list,128*cnum_full,crec_full)
      call report_error('FmpWrite',ierr,'writ','full block',1,
     .                  'qw_glb_full')

***** Thats all
      return
      end

CTITLE QW_DESCRIPTION
 
      subroutine qw_description

      implicit none 
 
 
*     Routine to create the solution description and write the record
*     to the global file.
*
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
 
*   i,j,k       - Loop counters
*   ierr        - File error flag
*   out_rec     - Record number for output
*   rec_copy(128) - Record for station information
 
      integer*4 i,j, ierr, out_rec, ns, rec_copy(128)
 
***** Copy the information we need.  Most of this information was obtained
*     while the hfile was read
*                               ! Initialize the output record number
      out_rec = crec_solutions

      do ns = 1, cnum_soln_recs 
          do i = 1, max_edit_types
              sdelete_count(i) = 0
          end do
 
          stai_utc = qtai_utc
      
          sut1_apr(1) = qut1_apr(1)
          sut1_apr(2) = qut1_apr(2)

          do i = 1,2
              do j = 1,2
                  swob_apr(j,i)     = qwob_apr(j,i)
                  snut_ang_apr(j,i) = qnut_ang_apr(j,i)
              end do
          end do

          sepoch_start = qsepoch_start(ns)
          sepoch_end   = qsepoch_end(ns)
          srun_time    = qsrun_time(ns)
          snum_parn    = qsnum_parn(ns)
          scons_type   = qscons_type(ns)
          ssys_type    = qssys_type(ns)
          sglb_vers    = qsglb_ver(ns)
          sowner       = qsowner(ns)
          screator     = qscreator(ns)
          sprog_gen    = qsprog_gen(ns)
          sanal_type   = qsanal_type(ns)
          sKalObs_file = qskalobs_file(ns)
          sexpt_title  = qsexpt_title(ns)

*         There is station information only for the last 
*         combined solution file so set values zero except for
*         the last one
          if( ns.lt. cnum_soln_recs ) then
              srec_sinf = 0
              snum_sinf = 0 
              snum_sites = 0
          else
               srec_sinf = crec_sinf
               snum_sinf = (qnum_sites-1)/2+1
               snum_sites = qnum_sites
          end if
                  
*****     Now write out record
 
          call writd(cglb_dcb,ierr,sdelete_count,128,out_rec)
          call report_error('FmpWrite',ierr,'writ','sln block',1,
     .                  'qw_description')
          out_rec = out_rec + 1
      end do
      
****  Now write out the station information records
      out_rec = crec_sinf
      do i = 1, qnum_sites, 2
         call qr_sinf_rec(i, rec_copy(1))
         call qr_sinf_rec(i+1, rec_copy(65))
         call writd(cglb_dcb, ierr, rec_copy, 128, out_rec)
         out_rec = out_rec + 1
      end do

****  Now write out the satellite records
      out_rec = crec_svinf
      do i = 1, qnum_svs, 4
         call qr_svinf_rec(i,   rec_copy(1))
         call qr_svinf_rec(i+1, rec_copy(33))
         call qr_svinf_rec(i+2, rec_copy(65))
         call qr_svinf_rec(i+3, rec_copy(97))

         call writd(cglb_dcb, ierr, rec_copy, 128, out_rec)
         out_rec = out_rec + 1
      end do
         
         
***** Thats all
      return
      end
      
CTITLE QR_SINF_REC 
  
      subroutine qr_sinf_rec( i, record)

      implicit none 

*     Routine to create the solution description and write the record
*     to the global file.
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'
 
* PASSED VALUES

* i  -  Station number
* record(64)  - Record for the Station information

      integer*4 i, j, record(64)
      

* LOCAL VARIABLES

*
*     Start saving the information
      sscode = qsite_names(i)
      ssdata_st = qdata_st(i) 
      if( ssdata_st.eq.0.d0 ) ssdata_st = sepoch_start
      ssdata_en = qdata_en(i) 
      if( ssdata_en.eq.0.d0 ) ssdata_en = sepoch_end
      ssrecv_st = qrecv_st(i) 
      if( ssrecv_st.eq.0.d0 ) ssrecv_st = sepoch_start
      ssrecv_en = qrecv_en(i) 
      if( ssrecv_en.eq.0.d0 ) ssrecv_en = sepoch_end
      ssante_st = qante_st(i) 
      if( ssante_st.eq.0.d0 ) ssante_st = sepoch_start
      ssante_en = qante_en(i) 
      if( ssante_en.eq.0.d0 ) ssante_en = sepoch_end
      do j = 1,3
         ssarp_ecc(j) = qarp_ecc(j,i)
         ssL1a_ecc(j) = qL1a_ecc(j,i)
         ssL2a_ecc(j) = qL2a_ecc(j,i)
         satmload(j)  = qatmload(j,i)
         shydload(j)  = qhydload(j,i)
      end do
* MOD TAH 200205: Save the antenna azimuth
      ssantdaz    = qantdaz(i)
      sselev_cut  = qelev_cut(i)
      ssnum_zen   = qnum_zen(i)
      ssant_mod   = qant_mod(i)
      sti_antmod  = qant_mod(i)  ! Saved twice with full name in sti_antmod
      ssante_sn   = qante_sn(i)
      ssante_ty   = qante_ty(i)
      ssradome_ty = qradome_ty(i)
      ssrecv_sn   = qrecv_sn(i)
      ssrecv_ty   = qrecv_ty(i)(1:16)
      ssrecv_fw   = qrecv_fw(i)
      ssrecv_ty_end = qrecv_ty(i)(17:20)
                 
      sscons_size = -1.d0
             
****  Now copy the values to the data record
      call wmov( ssdata_st, 1, record, 1, 64 )
      return
      end
 
CTITLE QR_SVINF_REC 
  
      subroutine qr_svinf_rec( i, record)

      implicit none 

*     Routine to create the Satellite information record
 
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/svinf_def.h'
 
* PASSED VALUES

* i  -  Station number
* record(32)  - Record for the satellite information

      integer*4 i, j, record(32)
      
*
*     Start saving the information
      svi_prn = qsvi_prn(i)
      svi_svn = qsvi_svn(i)
      svi_block = qsvi_block(i)
* MOD TAH 190627: Save index i tp svi_lnum to show satellite used.
      svi_lnum = i
      do j = 1,3
         svi_antpos(j,1) =  qsvi_antpos(j,1,i)
         svi_antpos(j,2) =  qsvi_antpos(j,2,i)
      end do
      svi_launch = qsvi_launch(i)
      svi_antmod = qsvi_antmod(i)
      svi_ocode = qsvi_ocode(i)
      do j = 1,9    ! Make sure we clear the rest of the record

         svi_fill(j) = 0
      end do
             
****  Now copy the values to the data record
      call wmov( svi_prn, 1, record, 1, 32 )
      return
      end
 

CTITLE QW_CODES
 
      subroutine qw_codes

      implicit none 
 
 
*     Routine to write the list of parameter and aproiri codes
*     to the global file.
*
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
 
*   i,j,k       - Loop counters
*   ierr        - File error flag
 
 
      integer*4 ierr
 
***** Now write out record(s)

      call writd(cglb_dcb,ierr,qglb_codes,128*cnum_par_types,
     .           crec_par_types)
      call report_error('FmpWrite',ierr,'writ','Global codes',1,
     .                  'qw_CODES')
 
      call writd(cglb_dcb,ierr,qapr_codes,128*cnum_apr_types,
     .           crec_apr_types)
      call report_error('FmpWrite',ierr,'writ','Aproiri codes',1,
     .                  'qw_CODES')
 
***** Thats all
      return
      end
 
CTITLE QW_MUL_EP
 
      subroutine qw_mul_ep
 
      implicit none 

 
*     Routine to write the epochs of the multi-epoch parameter
*     and apriori values.
*
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
 
*   i,j,k       - Loop counters
*   ierr        - File error flag
 
 
      integer*4 ierr
 
***** Now write out record(s)
      if( cnum_par_ep.gt.0 ) then
         call writd(cglb_dcb,ierr,qmul_par_ep,128*cnum_par_ep,
     .               crec_par_ep)
          call report_error('FmpWrite',ierr,'writ',
     .                      'Parameter multi-epochs',1,
     .                      'qw_CODES')
      end if
      if( cnum_apr_ep.gt.0 ) then
          call writd(cglb_dcb,ierr,qmul_apr_ep,128*cnum_apr_ep,
     .               crec_apr_ep)
          call report_error('FmpWrite',ierr,'writ',
     .                      'Apriori multi-epochs',1,
     .                      'qw_CODES')
      end if
 
 
***** Thats all
      return
      end
 
CTITLE QW_APRIORIS
 
      subroutine qw_aprioris
 
      implicit none 

 
*     Routine to make the list of apriori values and write these
*     to the global file.
*
      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
 
*   ierr        - File error flag

      integer*4 ierr

***** Now write out record(s)
 
      call writd(cglb_dcb,ierr,qapr_list,128*cnum_apr_vals,
     .           crec_apr_vals)
      call report_error('FmpWrite',ierr,'writ','aproiri vals',1,
     .                  'qw_aproiris')
 
***** Thats all
      return
      end
 
CTITLE QW_SOLN
 
      subroutine qw_soln ( cov_parm )
 
      implicit none 

 
*     Routine to write the covariance matrix and solution vector
*     to the disk file.  To do this the VWRIT routine is used
*     to write directly from ema to disk (via the scratch common area)
 
      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include 'htoglb_comm.h'
 
*   ierr            - Fmp Error
*   i,j,k           - Loop counters
*   FmpClose        - File closing emualtio routine
 
      integer*4 ierr, FmpClose
 
*   cov_parm(1)     - Covariance matrix with solution at the
*                   - end
 
      real*8 cov_parm(*)

* MOD TAH 190520: Mod to allow more than 32767x32767 matrices
      integer*8 I8   ! Needed for large numbers of parameters

      data I8 / 1 /
 
 
*     Now write out the matrices
* MOD TAH 190603: Change to writd8 to allow for >32767 parameters. 
      if( qnum_parn.gt.32767 ) then 
         call writd8(cglb_dcb, ierr, cov_parm, (I8*128)*cnum_par_vals,
     .              crec_par_vals)
         call report_error('WRITD8',ierr,'writ','COV_PARM',0,'qw_SOLN')
      else
         call writd(cglb_dcb, ierr, cov_parm, 128*cnum_par_vals,
     .              crec_par_vals)
         call report_error('WRITD',ierr,'writ','COV_PARM',0,'qw_SOLN')
      endif 
      
*     Write out the apriori
      if ( cnum_acvc.gt.0 ) then
          call qw_cons
      end if
 
***** Now close the output file
 
      ierr = FmpClose(cglb_dcb)
      call report_error('FmpClose',ierr,'clos', glb_file,0,'qw_SOLN')
 
* DEBUG
 
C     write(1,50) (crun_time(i), i=1,5), glb_file
C 50  format(/' Run ',5i3,' Global ',a20,/,
C    .        ' Parameter values ')
 
***** Thats all
      return
      end
  
CTITLE QW_CONS

      subroutine qw_cons
 
      implicit none 

 
*     Routine to write the apriori covariance matrix used in the
*     analysis into the global file.
 
      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
 
      include 'htoglb_comm.h'
 

*   record(128) -  Record to be output (paired as p1,p2, covariance)
*                  Therefore 32 triplets per record.
*   rec_acvc    -  Current record number of the apriori covariance
*                  matrix
*   num_in_rec  -  Number of triplets in current record.  When this
*                  equals 32 it is written to the file.
*   pir          - Position in record of first parameter number

      integer*4 record(128),rec_acvc, num_in_rec, pir, i,j,k,
     .          ierr, np
 
      
****  Initialize the record counters
      rec_acvc =  crec_acvc 
      num_in_rec = 0
      cnum_acvc  = 0   
      
***** Now add the NEU station constrains.  Also count the number of
*     3x3 constraint matrices put on the solution.  This is used when
*     write the apriori constraints to combined global files.

      do i = 1, qnum_sites

*        Loop the sites
         np = qparn_sites(1,i)
         if( np.ne.0.and. qapr_cov(1,1,i).gt.0 ) then

*            Write out the values            
             do j = 1,3
                 do k = 1,j
                    num_in_rec = num_in_rec + 1
                    pir = (num_in_rec-1)*4 + 1
                    record(pir) = np+j-1
                    record(pir+1) = np+k-1
                    
*                   Move two values for real*8 variable                      
                    call dwmov(qapr_cov(j,k,i), 1, record(pir+2),1,1)
                      

*****               Now see if we need to write the record
                    if( num_in_rec.eq.32 ) then
                       call writd(cglb_dcb, ierr, record, 128,
     .                            rec_acvc)
                       rec_acvc = rec_acvc + 1
                       cnum_acvc = cnum_acvc + num_in_rec
                       num_in_rec = 0
                    end if
                 end do
             end do          
*                       ! Parameter estimated
          end if
*                       ! Looping over the sites.
      end do
      
****  Now do the satellite elements
      do i = 1, qnum_apr_svs
         np = qparn_svs(1,i)
         if( np.gt.0 ) then
             do j = 1,6
                do k = 1,j
                    num_in_rec = num_in_rec + 1
                    pir = (num_in_rec-1)*4 + 1
                    record(pir) = np+j-1
                    record(pir+1) = np+k-1
                    
*                   Move two values for real*8 variable                      
                    call dwmov(qapr_svs(j,k,i), 1, record(pir+2),1,1)
                      

*****               Now see if we need to write the record
                    if( num_in_rec.eq.32 ) then
                       call writd(cglb_dcb, ierr, record, 128,
     .                            rec_acvc)
                       rec_acvc = rec_acvc + 1
                       cnum_acvc = cnum_acvc + num_in_rec
                       num_in_rec = 0
                    end if
                 end do
             end do          
          end if
      end do
      
****  Now do the diagonal elements
      if ( qnum_apr_diag.gt.0 ) then
          do i = 1, qnum_parn
             if( qapr_diag(i).gt.0 ) then
                 num_in_rec = num_in_rec + 1
                 pir = (num_in_rec-1)*4 + 1
                 record(pir) = i        
                 record(pir+1) = i       
                    
*                Move two values for real*8 variable                      
                 call dwmov(qapr_diag(i), 1, record(pir+2),1,1)
                      

*****            Now see if we need to write the record
                 if( num_in_rec.eq.32 ) then
                    call writd(cglb_dcb, ierr, record, 128,
     .                         rec_acvc)
                    rec_acvc = rec_acvc + 1
                    cnum_acvc = cnum_acvc + num_in_rec
                    num_in_rec = 0
                 end if
             end if
          end do
      end if
                   
****  See if we have a residual amount of record to write to
*     file
      if( num_in_rec.gt.0 ) then

*         Clear the trailing part of the record      
          do i = num_in_rec*4+1, 128
              record(i) = 0
          end do         
          cnum_acvc = cnum_acvc + num_in_rec
          call writd(cglb_dcb, ierr, record, 128,
     .               rec_acvc)
          call report_error('VWRIT',ierr,'writ','COV_PARM',
     .                      0,'GW_CONS')
      end if
 
***** Thats all
      return
      end

CTITLE mk_sln_gamit

      subroutine mk_sln_gamit

      implicit none 

*     This routine makes the sln_def record for a gamit solution
*     for those entries that are not set during the run.

      include '../includes/kalman_param.h'
      include 'htoglb_comm.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

      integer*4 ns
      real*8 sectag

      ns = 1 
      qsepoch_start(ns) = qstart_epoch 
      qsepoch_end(ns)   = qend_epoch  
      sectag = qrun_time(6)
      call ymdhms_to_jd(qrun_time, sectag, qsrun_time(ns))
      qsnum_parn(ns)    = qnum_parn 
      qscons_type(ns)   = 2           
      qssys_type(ns)    = 2
      qsglb_ver(ns)     = glbf_version
      qsowner(ns)       = cowner       
      call get_institute(qscreator(ns), default_institute)
      qsprog_gen(ns)    = '+GLK'
      qsanal_type(ns)   = canal_type   
      qskalobs_file(ns) = sKalobs_file
      qsexpt_title(ns)  = qexpt_title   

*     Thats all
      return
      end
