CTITLE GET_PARAM_CODES
 
      subroutine get_param_codes ( codes )

      implicit none  
 
*     Routine to read and save the parameter codes from the global
*     file.
 
      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/globk_common.h'
 
*   buffer(128)     - buffer for reading codes into before
*                   - transferring to CODES
*   codes(1)        - Array to hold the codes
*   i,j             - Loop counter
*   ierr            - Error flag
*   len             - Length read from file (not used)
*   nbuf            - Number of values to be copied from
*                   - buffers we have to read
 
 
      integer*4 buffer(128), codes(1), i,j, ierr, len, nbuf

*   in_rec -- Record number to read from h-file
*   curr_max_ent -- Number of entries currently avalaible
*   type  -- Type of parmeter (lower 8-bits)
*   indx  -- sub-type of parameter 

      integer*4 in_rec, curr_max_ent, type, indx 

* MOD TAH 190610: Updated the satellite radition parameters codes
*     accomiodate the ECOMC model
      integer*4 svcode_to_code ! Function to return new code given old one.
 
****  Read the records from the file and transferr into CODES
*     Note that here we read the codes one record at a time.  This
*     if an artifact of the HP version not being able to directly
*     read into VMA.  Code has been keep this way incase the max
*     number of global parameters (and thus the dimension of code)
*     is not a multiple of 128)
      do i = 1,cnum_parn
        codes(i) = 0
      end do
      do i = 1, cnum_par_types
 
        call readd(cglb_dcb,ierr,buffer,128,len,crec_par_types+i-1)
 
*                                           ! Only copy part of buffer
        if( i.eq.cnum_par_types ) then
              nbuf = cnum_parn - (i-1)*128  
*                                             ! Copy all of buffer
          else
              nbuf = 128
          end if
 
*         Copy values to codes
          do j = 1, nbuf
* MOD TAH 190610: Update any satellite radition codes if <107 GLX file
!             codes((i-1)*128+j) = buffer(j)
              codes((i-1)*128+j) = svcode_to_code( buffer(j))
          end do
      end do 

****  Now see if any of the codes referr to multiple echo parameters
*     (So far these are only types 56, 57 and 58.
*     Clear the number point array and initialize the record number
* MOD TAH 080422: Based on Bei Jinzhong email, size of ij loop reversed
*     to be consistent with dimensions.
      do i = 1, 2
         do j = 1, 3
            cnum_mul_pmu(i,j) = 0
         end do
      end do

      if( cnum_par_ep.gt.0 ) then 
          in_rec = crec_par_ep - 1
          curr_max_ent = 0
    
          do i = 1, cnum_parn
             call decode_code( codes(i), type, indx)
             if( type.ge.56 .and. type.le.58 ) then

*                These are multi-epoch parameters types.  We now need to 
*                read and save the epochs that correspond to these
*                parameters,
                 call get_mul_ep( type, indx, in_rec, curr_max_ent )
             end if
          end do 
      end if

***** Thats all
      return
      end
 
CTITLE GET_MUL_EP

      subroutine get_mul_ep ( type, indx, in_rec, curr_max_ent )

      implicit none 

*     Routine to read the multiple epoch parameter epoch values from
*     the input hfile.

      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/globk_common.h'

* PASSED Variables
* ----------------
*  type  -- Type of parameter (currently 56-58)
*  indx  -- Contains pointer to the sub-type and the entry number in the
*           parameter epoch block corresponding to this entry
*  in_rec -- Current record to read for the epoch block
*  curr_max_ent -- Current maximum entry number

      integer*4 type, indx, in_rec, curr_max_ent, len_read

* LOCAL Variables
*----------------
* m1, m2  -- Split of the indx in sub-type and epoch numeber
* values(128) -- Values read from epoch block

      integer*4 m1, m2, ierr
      real*8 values(128)


****  Split the sub-type (offset or rate) from the indx and get the entry
*     number for the epoch number
      call decode_code( indx, m1, m2)
      cnum_mul_pmu(m1,type-55) = cnum_mul_pmu(m1,type-55)+1

*     Now try to the epoch.  See if entry in current list
      if( m2.gt. curr_max_ent ) then

*         Need to read the next block of values from the binary file
          in_rec = in_rec + 1
          curr_max_ent = curr_max_ent + 64
          call readd(cglb_dcb, ierr,values,128,len_read, in_rec)
      end if

*     Save  the value
      cmul_pmu_ep(m2) = values(mod(m2-1,64)+1)

****  Thats all
      return
      end 

