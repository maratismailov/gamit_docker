      program swaph

      implicit none 

*     Program to read a globk binary hfile and swap the byte order if necessary.
*     This is to allow hfiles generated on the linux and Dec systems to used on
*     the HP and Sun sytems and visa versa.  The version number of the hfile
*     is used to indicate that the bytes have been swapped.

* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'

* VARIABLES

* len_run  -- Length of runstring
* nf       -- File number being processed
* rcpar    -- Routine to read the runstring
* ierr, kerr -- Error returns.
* FmpClose   -- Rotuine to close binary files.

      integer*4 len_run, nf, rcpar, ierr, kerr,  FmpClosE

* glb_inp_file -- Name of the global file being processed.

      character*128 glb_inp_file 


***** Start, loop over the names of the hfiles  from the runstring.
      len_run = 1
      nf = 0
      do while ( len_run.gt.0 ) 
          nf = nf + 1
          len_run = rcpar(nf, glb_inp_file)
          if( len_run.le.0 ) then
              if( nf.eq.1 )  
     .        call proper_runstring('swaph.hlp','swaph',1)
          else

*             Open the hfile and read the header
              call fmpopen(cglb_dcb,kerr,glb_inp_file,'rwo',0 )
              call report_error('FmpOpen',kerr,'open', glb_inp_file,
     .                          0,'SWAPH')

*             Read the header and see what we have
              if( kerr.eq.0 ) then 
                  call rw_glb_header('R',kerr)

*                 Check the file type
                  if( cfile_type.ne.'GLOBAL' ) then
                      write(*,110)  nf,glb_inp_file(1:len_run)
 110                  format(' File ',i3,1x,a,': NOT a GLOBAL FILE')
                      kerr = -1
                  end if 
               end if

****           See if file still OK to process.
               if( kerr.eq.0 .or. kerr.eq.2001 ) then 

*                 Check the version number
                  if( cglb_vers.eq.0 ) then
                      write(*,120) nf,glb_inp_file(1:len_run)
 120                  format(' File ',i3,1x,a,
     .                       ' is version 0.0: Cannot convert')  
                      ierr = FmpClose(cglb_dcb)
                  else if( cglb_vers.lt.32767) then
                      write(*,140) nf,glb_inp_file(1:len_run),
     .                       cglb_vers/100.0
 140                  format(' File ',i3,1x,a,' is version ',
     .                       F4.2,': Does not need conversion')
                      ierr = FmpClose(cglb_dcb)
                  else

*                     Start the converion:
                      call swap_head
                      write(*,160) nf,glb_inp_file(1:len_run),
     .                       cglb_vers/100.0
 160                  format(' File ',i3,1x,a,' is version ',
     .                       F4.2,': Converting')

*                     Write the header out
                      call rw_glb_header('W',kerr)
 
*                     Now run over the blocks in the file converting
*                     each type
                      call swap_write_blocks

                      ierr = FmpClose(cglb_dcb)
                 end if
             end if
         end if
      end do

****  Thats all
      end


CTITLE SWAP_HEAD

      subroutine swap_head

      implicit none 

*     This routine swaps the header information for the global binary
*     hfiles.  (This is also used to read the rest of the file).

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
 

****  Start, just run through each variable swapping as need be.

      call swap_bytes(4, cnum_header, 1)
      call swap_bytes(4, cnum_names, 1) 
      call swap_bytes(4, cnum_full, 1) 
      call swap_bytes(4, cnum_soln_recs, 1) 
      call swap_bytes(4, cnum_par_types, 1) 
      call swap_bytes(4, cnum_apr_types, 1)
      call swap_bytes(4, cnum_apr_vals, 1) 
      call swap_bytes(4, cnum_par_vals, 1) 
      call swap_bytes(4, crec_names, 1) 
      call swap_bytes(4, crec_full, 1) 
      call swap_bytes(4, crec_solutions, 1)
      call swap_bytes(4, crec_par_types, 1) 
      call swap_bytes(4, crec_apr_types, 1) 
      call swap_bytes(4, crec_apr_vals, 1) 
      call swap_bytes(4, crec_par_vals, 1)
c     cfile_type  
c     cexpt_title  
      call swap_bytes(4, cnum_apr_codes, 1) 
      call swap_bytes(4, cnum_parn, 1)
      call swap_bytes(4, cnum_sites, 1) 
      call swap_bytes(4, cnum_sources, 1) 
      call swap_bytes(4, cchisq, 1) 
      call swap_bytes(4, cnum_svs, 1) 
      call swap_bytes(4, crun_time, 7)
      call swap_bytes(4, ctai_utc, 1) 
      call swap_bytes(4, cdelete_count, max_edit_types)
      call swap_bytes(4, cnum_obs, 1) 

*     Real*8 values
      call swap_bytes(8, cepoch_end, 1) 
      call swap_bytes(8, cepoch_expt, 1) 
      call swap_bytes(8, cepoch_start, 1) 
      call swap_bytes(8, cetd_apr, 3)
      call swap_bytes(8, cnut_ang_apr, 8)
      call swap_bytes(8, cut1_apr, 6) 
      call swap_bytes(8, cwob_apr, 8)
      call swap_bytes(8, csvs_epoch, 1)

      call swap_bytes(4, cglb_vers, 1) 
      call swap_bytes(4, cnum_sinf, 1) 
      call swap_bytes(4, crec_sinf, 1) 
      call swap_bytes(4, cnum_acvc, 1) 
      call swap_bytes(4, cnum_acvc_recs, 1) 
      call swap_bytes(4, crec_acvc, 1) 
      call swap_bytes(4, cgpst_utc, 1) 
      call swap_bytes(4, crec_comb_soln, 1) 
      call swap_bytes(4, cnum_comb, 1)
      call swap_bytes(4, ccons_type, 1) 
      call swap_bytes(4, csys_type, 1)

* MOD TAH 981020: Add new variables for multi-pmu parameters
      call swap_bytes(4, cnum_par_ep, 1) 
      call swap_bytes(4, cent_par_ep, 1)
      call swap_bytes(4, crec_par_ep, 1)
      call swap_bytes(4, cnum_apr_ep, 1)
      call swap_bytes(4, cent_apr_ep, 1)
      call swap_bytes(4, crec_apr_ep, 1)
    
* MOD TAH 991130: Swap the gamit_mod variable.  (Make sure that it has
*     not been forgotten to be swappped before.  If value > 2^24 then
*     not previously swapped, so leave untouched.
* MOD TAH 030517: Still a problem with original code.  If the bytes
*     need swapping, then value should be > 2**30 since the EOP bits 
*     will be in the top byte.  Therefore if greater than 2**30, 
*     swap the bytes, otherwise leave untouched.  Same change made to
*     sgamit_mod below.
C     if( cgamit_mod.lt. 2**30 )  call swap_bytes(4, cgamit_mod, 1)  
      if( cgamit_mod.gt. 2**30 )  call swap_bytes(4, cgamit_mod, 1) 
* MOD TAH 050622: New load model bit map.
      call swap_bytes(4, cload_mod, 1)

* MOD TAH 050622: New variables for satellite information records
      call swap_bytes(4,cnum_svinf,1)
      call swap_bytes(4,crec_svinf,1)

****  Thats all
      return
      end
 
      
CTITLE SWAP_WRITE_BLOCKS

      subroutine Swap_write_blocks

      implicit none 

*     Routine to swap the remaining bytes in a hfiles and write out
*     records with the swapped bytes in them


* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'


****  Loop over the blocks in the global hfile

      call sw_glb_names

      call sw_description

      call sw_glb_full

      call sw_codes

      call sw_aprioris 

      call sw_mul_pmu

      call sw_soln

      call sw_cons

***** Thats all
      return 
      end

CTITLE SW_GLB_NAMES

      subroutine sw_glb_names

      implicit none 

*     Routine to swap the names block

* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

****  Since all the records here are character strings their is no need 
*     to do anything

      return
      end

CTITLE SW_DESCRIPTION

      subroutine sw_description 

      implicit none 

*     Routine to write out the description arrays.  Here we need to account
*     for the mix of variables in the description blocks.


* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* LOCAL VARIABLES

* i  -- Loop counter
* rec -- record number to read and wrte
* ierr -- IOSTAT error
* len_read -- Length of record record

      integer*4 i, rec, ierr, len_read

* array(128)  -- Array to read the pair of station information records into 

      integer*4 array(128) 

****  Start by reading/swapping and writing the solution information records
*     There is should be one of these per input global file. 

      do i = 1, cnum_soln_recs
         rec = crec_solutions + i - 1
         call readd(cglb_dcb,ierr,sdelete_count,128,len_read,rec) 

*        Now convert the values in this record
          
         call swap_bytes(4, sdelete_count, max_edit_types)
         call swap_bytes(4, snum_obs, 1)
         call swap_bytes(4, stai_utc, 1) 
         call swap_bytes(4, sversion, 1)
         
         call swap_bytes(4, snum_parn, 1)
         call swap_bytes(4, secc_change, 3*max_sln_sites )   
         call swap_bytes(4, sgpst_utc, 1)  
         
C        call swap_bytes(4, sgframe  
C        call swap_bytes(4, sgprec  
C        call swap_bytes(4, sgsrpmod  
C        call swap_bytes(4, sgtime 
         
C        call swap_bytes(4, sprog_gen  
         call swap_bytes(8, srun_time, 1) 
         
         call swap_bytes(8, sepoch, 1) 
         call swap_bytes(8, ssvs_epoch, 1) 
         
         call swap_bytes(8, sut1_apr, 2) 
         call swap_bytes(8, swob_apr, 4) 
         call swap_bytes(8, snut_ang_apr, 4) 
C        call swap_bytes(4, sdata_base 
          
C        call swap_bytes(4, sKalObs_file  
C        call swap_bytes(4, sanal_type  
          
         call swap_bytes(8, sepoch_end, 1)
         call swap_bytes(8, sepoch_start, 1)
          
         call swap_bytes(4, sglb_vers, 1)
         call swap_bytes(4, srec_sinf, 1)
         call swap_bytes(4, snum_sinf, 1)
         call swap_bytes(4, scons_type, 1)
         call swap_bytes(4, ssys_type, 1)
         
         call swap_bytes(4, snum_sites, 1)
C        call swap_bytes(4, sowner  
C        call swap_bytes(4, screator 
          
C        call swap_bytes(4, sexpt_title 

* MOD TAH 991130: Swap the gamit_mod variable.  (Make sure that it has
*     not been forgotten to be swappped before.  If value > 2^24 then
*     not previously swapped, so leave untouched.
*     not previously swapped, so leave untouched.
* MOD TAH 030517: Still a problem with original code.  If the bytes
*     need swapping, then value should be > 2**30 since the EOP bits 
*     will be in the top byte.  Therefore if greater than 2**30, 
*     swap the bytes, otherwise leave untouched.  Same change made to
*     cgamit_mod above.
C        if( sgamit_mod.lt. 2**30 )  call swap_bytes(4, sgamit_mod, 1)  
         if( sgamit_mod.gt. 2**30 )  call swap_bytes(4, sgamit_mod, 1)
* MOD TAH 050622: New load model bit map
         call swap_bytes(4, sload_mod, 1) 


*        Now write out record
         call writd(cglb_dcb,ierr,sdelete_count,128,rec)
      end do

****  Now do the station information records
      do i = 1, cnum_sites, 2
          rec = crec_sinf + (i-1)/2 
          call readd(cglb_dcb, ierr, array, 128, len_read, rec)
          call sw_sinf(array(1))
          call sw_sinf(array(65))
          call writd(cglb_dcb, ierr, array, 128, rec)
      end do

***** Now swap the satellite information records if there are present
      if( crec_svinf.gt.0 ) then
          do i = 1, cnum_svs, 4
             rec = crec_svinf + (i-1)/4
             call readd(cglb_dcb, ierr, array, 128, len_read, rec)
             call sw_svinf(array(1))
             call sw_svinf(array(33))
             call sw_svinf(array(65))
             call sw_svinf(array(97))
             call writd(cglb_dcb, ierr, array, 128, rec)
          end do
      endif
 
 
****  Thats all
      return
      end

CTITLE SW_SINF

      subroutine sw_sinf ( array )

      implicit none 

*     Routine to swap the bytes in the station information array.

* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/sinf_def.h'

* PASSED VARIABLES

* array(64)  -- Array containing the record

      integer*4 array(64)

***** Copy the record copy into to the allocated common so that we know
*     what things are

      call wmov( array, 1, ssdata_st, 1, 64 )

      call swap_bytes(8, ssdata_st, 1)
      call swap_bytes(8, ssdata_en, 1)
      call swap_bytes(8, ssrecv_st, 1)
      call swap_bytes(8, ssrecv_en, 1)
      
      call swap_bytes(8, ssante_st, 1)
      call swap_bytes(8, ssante_en, 1)
      
      call swap_bytes(4, ssarp_ecc, 3)
      call swap_bytes(4, ssL1a_ecc, 3)
      call swap_bytes(4, ssL2a_ecc, 3)
      call swap_bytes(4, sselev_cut, 1)
      call swap_bytes(4, sscons_size, 1)
      
      call swap_bytes(4, ssnum_zen, 1)

* MOD TAH 050622: New variables
      call swap_bytes(8, slog_ep, 1)
      call swap_bytes(4, slog_tau, 1)
      call swap_bytes(4, satmload, 3)
      call swap_bytes(4, shydload, 3)

C     call swap_bytes(4, ssant_mod  
C     call swap_bytes(4, sscode  
C     call swap_bytes(4, ssante_sn  
C     call swap_bytes(4, ssrecv_sn 
      
C     call swap_bytes(4, ssrecv_ty  
C     call swap_bytes(4, ssante_ty  
C     call swap_bytes(4, ssrecv_fw  
C     call swap_bytes(4, ssres

*     Now move the values back into the array 
      call wmov( ssdata_st, 1, array, 1, 64 )

****  Thats all 
      return 
      end

CTITLE SW_SVINF

      subroutine sw_svinf ( array )

      implicit none 

*     Routine to swap the bytes in the station information array.

* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
      include '../includes/svinf_def.h'

* PASSED VARIABLES

* array(32)  -- Array containing the record

      integer*4 array(32)

***** Copy the record copy into to the allocated common so that we know
*     what things are

      call wmov( array, 1, svi_prn, 1, 32 )

      call swap_bytes(4, svi_prn,   1)
      call swap_bytes(4, svi_svn,   1)
      call swap_bytes(4, svi_block, 1)

      call swap_bytes(8, svi_antmod, 6)
      call swap_bytes(8, svi_launch, 1)

*     Now move the values back into the array 
      call wmov( svi_prn, 1, array, 1, 32 )

****  Thats all 
      return 
      end

CTITLE SW_GLB_FULL 

      subroutine sw_glb_full 

      implicit none 

*     Routine to swap the full names.  Since everything here is character
*     strings there is no need to do anthing.


* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

****  Since all the records here are character strings their is no need 
*     to do anything

      return
      end

CTITLE SW_CODES 

      subroutine sw_codes 

      implicit none 

*     Routine to swap the bytes in the codes block.  All the variables
*     here are integer*4 so we simply need to read the records, swap the
*     bytes and write out the records


* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* LOCAL VARIABLES

* rec  -- Record number to read/write file
* i    -- Loop counter over the records in the parameter and apriori
*         codes block.
* ierr -- IOSTAT error on read.
* len_read -- Length actually read from file (not used).

* array(128) -- Record read from file

      integer*4  rec, i, ierr, len_read, array(128)

***** Loop over the records, reading, swapping and writting as we go.

      do i = 1, cnum_par_types
         rec = crec_par_types + i - 1
         call readd(cglb_dcb,ierr,array, 128, len_read, rec )

         call swap_bytes(4, array, 128 )

         call writd(cglb_dcb,ierr,array, 128, rec )

      end do

      do i = 1, cnum_apr_types
         rec = crec_apr_types + i - 1
         call readd(cglb_dcb,ierr,array, 128, len_read, rec )

         call swap_bytes(4, array, 128 )

         call writd(cglb_dcb,ierr,array, 128, rec )

      end do

****  Thats all
      return
      end


CTITLE SW_APRIORIS 

      subroutine sw_aprioris 

      implicit none 

*     Routine to read the apriori values block, swap and write

* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* LOCAL VARIABLES

* rec  -- Record number to read/write file
* i    -- Loop counter over the records in the parameter and apriori
*         codes block.
* ierr -- IOSTAT error on read.
* len_read -- Length actually read from file (not used).

* array(64) -- Record read from file (real*8 array)

      integer*4  rec, i, ierr, len_read
      real*8     array( 64)

***** Loop over the records, reading, swapping and writting as we go.

      do i = 1, cnum_apr_vals 
         rec = crec_apr_vals + i - 1
         call readd(cglb_dcb,ierr,array, 128, len_read, rec )

         call swap_bytes(8, array,  64 )

         call writd(cglb_dcb,ierr,array, 128, rec )

      end do

****  Thats all
      return
      end

CTITLE SW_MUL_PMU 

      subroutine sw_mul_pmu 

      implicit none 

*     Routine to read the multi-pmu epochs and swap the bytes in
*     the record.

* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* LOCAL VARIABLES

* rec  -- Record number to read/write file
* i    -- Loop counter over the records in the parameter and apriori
*         codes block.
* ierr -- IOSTAT error on read.
* len_read -- Length actually read from file (not used).

* array(64) -- Record read from file (real*8 array)

      integer*4  rec, i, ierr, len_read
      real*8     array( 64)

***** Loop over the records, reading, swapping and writting as we go.
*     Swap the parameter estimate epochs
      do i = 1, cnum_par_ep 
         rec = crec_par_ep + i - 1
         call readd(cglb_dcb,ierr,array, 128, len_read, rec )

         call swap_bytes(8, array,  64 )

         call writd(cglb_dcb,ierr,array, 128, rec )

      end do

*     Swap the apriori epochs
      do i = 1, cnum_apr_ep 
         rec = crec_apr_ep + i - 1
         call readd(cglb_dcb,ierr,array, 128, len_read, rec )

         call swap_bytes(8, array,  64 )

         call writd(cglb_dcb,ierr,array, 128, rec )

      end do

****  Thats all
      return
      end



CTITLE SW_SOLN     

      subroutine sw_soln 

      implicit none 

*     Rouitne to read/swap and write the solution (plus covariance matrix)
*     of the binary hfiles.

* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* LOCAL VARIABLES

* rec  -- Record number to read/write file
* i    -- Loop counter over the records in the parameter and apriori
*         codes block.
* ierr -- IOSTAT error on read.
* len_read -- Length actually read from file (not used).

* array( 64) -- Record read from file

      integer*4  rec, i, ierr, len_read
      real*8     array( 64)

***** Loop over the records, reading, swapping and writting as we go.

      do i = 1, cnum_par_vals
         rec = crec_par_vals + i - 1
         call readd(cglb_dcb,ierr,array, 128, len_read, rec )

         call swap_bytes(8, array,  64 )

         call writd(cglb_dcb,ierr,array, 128, rec )

      end do

****  Thats all
      return 
      end 

CTITLE SW_CONS  

      subroutine sw_cons

      implicit none 

*     Routine to read/swap and write the constraints block in the 
*     binary hfiles.  For this block, the parameters are arranged
*     as a pair of I*4 parameter numbers, then a real*8 value.


* INCLUDE FILES

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'

* LOCAL VARIABLES

* rec  -- Record number to read/write file,
* num_recs -- Number so records to be read.
* i    -- Loop counter over the records in the parameter and apriori
*         codes block.
* ierr -- IOSTAT error on read.
* len_read -- Length actually read from file (not used).

* array(128) -- Record read from file.  These are 2 I*4 values and then
*     a real*8

      integer*4  rec, num_recs, i, j, ierr, len_read
      integer*4  array(128)

***** Loop over the records, reading, swapping and writting as we go.

      num_recs = cnum_acvc/32 + 1
      if( num_recs.eq.0 ) RETURN 

      do i = 1, num_recs 
         rec = crec_acvc + i - 1
         call readd(cglb_dcb,ierr,array, 128, len_read, rec )

*        Now work through the array
         do j = 1, 32
            call swap_bytes(4, array((j-1)*4+1),  2 )
            call swap_bytes(8, array((j-1)*4+3),  1 )
         end do 

         call writd(cglb_dcb,ierr,array, 128, rec )

      end do

****  Thats all
      return 
      end 
