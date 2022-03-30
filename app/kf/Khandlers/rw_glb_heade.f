CTITLE RW_GLB_HEADER
 
      subroutine rw_glb_header( option, ierr )
 

      implicit none
 
*     Routine to read or write the header records from a global
*     file.  The file is assumed to be already open.

* MOD TAH 930611: added reading of the first solution description.
*     (if their is more than solution, then the sdata_base name is
*      set to CombinedG (this is for output).

* MOD TAH 970609: Checked if files appears byte swapped.  Ierr returns
*     as 2001 if it appears to be. 
 
 
      include '../includes/kalman_param.h'
      include '../includes/glb_hdr_def.h'
      include '../includes/sln_def.h'
* MOD TAH 140223: Needed to add globk_markov.h to fix error with
*     times in CWU GIPSY sinex files.
      include '../includes/globk_markov.h'
*   ierr        - FmpRead error
*   len_read    - Number of words read with readd
 
      integer*4 ierr, len_read

*   num_head -- Number of header records (after possible swap).
*   rec_soln -- Record number for solution after possible swap
*   num_soln_rec - Number of solution recs after possible swap.
*   glb_vers -- File version after possible swap 
*   swap     -- Logical to indicate that we need to swap

      integer*4 num_head, rec_soln, num_soln_rec, glb_vers
      logical swap 
 
*   option      - One character option passed by user
*   opt         - Upper case version of option.
 
      character*1 option, opt
 
***** Convert option to upper case
 
      opt = option
      call casefold( opt )
 
***** Read or write the file based on opt
 
*                                 ! Read.
      if( opt.eq.'R' ) then
 
*         read first block to get the size of the header
          call readd(cglb_dcb, ierr, glb_header, 128, len_read,1)

*         Now read the full amount based on size.
          if( ierr.ge.0 ) then

* MOD TAH 970609: Check that this is global file and see if byte swapped
              if( cfile_type.eq.'GLOBAL' ) then
                  num_head = cnum_header
                  swap = .false. 
                  if( num_head.gt.32767) then
                      call swap_bytes(4, num_head, 1)
                      swap = .true. 
                  end if 
                  call readd(cglb_dcb, ierr, glb_header, 128*num_head,
     .                   len_read, 1 )
              else
                  ierr = -1
              end if 
          end if
 
          call report_error('FmpRead',ierr,'read','glb_header',0,
     .                      'RW_GLB_HEADER')
*                 ! Read option.

*         If there is only one experiment in this solution get the
*         data base name (of mfile for GPS) otherwize set to combined
*         global

*         Check for possible byte swapping (NOTE: the swapped version is
*         not saved.  This is done latter in routine to swap the file.
          glb_vers = cglb_vers
          rec_soln = crec_solutions
          num_soln_rec = cnum_soln_recs
          if( swap ) then
             call swap_bytes(4, rec_soln, 1)
             call swap_bytes(4, num_soln_rec, 1) 
             call swap_bytes(4, glb_vers, 1)
          end if 

          if( glb_vers.lt.100 .and. num_soln_rec.gt.1 ) then
             sdata_base = 'CombinedG'
          else
             call readd(cglb_dcb,ierr, sdelete_count,128,len_read,
     .                  rec_soln+num_soln_rec-1)
          end if

          if( swap ) ierr = 2001 

* MOD TAH 140223: Added correction to start and stop times for CWU and 
*         COM (pbo) solutions (based on hfile name).  
*         Check for 7-day error exactly (since is only case
*         that needs fixing).

          if( (index(glb_inp_file,'cwu').gt.0 .or. 
     .        index(glb_inp_file,'com').gt.0) .and.
     .        abs((cepoch_end-cepoch_start)-7.d0).lt.1.d-3 )
     .                                                 then
              cepoch_end   = cepoch_end   - 3.0d0 
              cepoch_start = cepoch_start + 3.0d0
           endif

      end if
 
      if( opt.eq.'W' ) then
          call writd(cglb_dcb,ierr, glb_header, 128*cnum_header, 1)
          call report_error('FmpWrite',ierr,'writ','glb_header',0,
     .                      'RW_GLB_HEADER')
 
*                 ! Write option.
      end if
 
***** Thats all
      return
      end
 
