CTITLE CREATE_SOLVK
 
      subroutine create_solvk( CO_name )
 

      implicit none
c
c     This routine creates 'SOLVK common' with the name in CO_name.
c
c     Here we set 'control_READ' true so that the RW_SOLVK_Block routine
c     knows that the record lengh information is available
 
c
c     Include files
c     -------------
      include '../includes/kalman_param.h'
      include '../includes/solvk_common.h'
c
      integer*4 total_blocks, err
 
c
      integer*4 num_con_words, num_mar_words
 
* MOD TAH 190511: Changed to I*8 for 64-bit memory
      integer*8 AddressOf
 
c
*   full_len    - Length of the full Solvk file name.
*   root_len    - Length of the root part of the Solvk file name
*   trimlen     - HP utility to return length of string
 
      integer*4 full_len, trimlen
 
*   CO_name     - Name of the SOLVK common file
 
      character*(*) CO_name
 
*   full_solvk_filename     - the full name of the Solvk file
*                           - including the file type and size.
 
 
      character*128 full_solvk_filename
 
c.... Get info about the size of the common blocks (NOTE: we remove
c     one extra word because we do not need to save last_@_word for
c     each block)
 
      num_con_words    = (AddressOf(last_control_word)
     .                 - AddressOf(control) - 1)/4 + 1
 
      num_mar_words    = (AddressOf(last_markov_word)
     .                 - AddressOf(markov) - 1)/4 + 1
 
c
c.... Number of blocks for the values common (NOTE: we do not add extra
c     block becuase there is dummy space for odd amount of block at the
c     end of each block.)
      num_control_blocks= (num_con_words - 1) / 128 + 1
 
*     Get the start of the markov block and the number of blocks in it
      start_markov      = num_control_blocks + 1
      num_markov_blocks = (num_mar_words-1) / 128 + 1
 
c.... Total number of blocks in the file
      total_blocks = start_markov + num_markov_blocks - 1
c
c.... Create the full file name
 
      call FullFileName(CO_name, 1,total_blocks, 128,
     .                  full_solvk_filename )
 
      full_len = trimlen(Full_solvk_filename)
 
      call FmpOpen(isoldc, err, full_solvk_filename, 'RWOC', 1)
 
c.... Error?
      if( err .lt. 0 ) then
c
c....   Say not OK
        write(icrt  ,200) err,full_solvk_filename(1:full_len)
        if( ilog.ne.icrt )
     .  write(ilog  ,200) err,full_solvk_filename(1:full_len)
  200   format(/,' FMGR',I4.3,' on create of ',A,/ )
 
        stop ' SOLVK_CREATE: Error creating SOLVK common '
c
*             ! Indicate control read (This ensures that the block size
      else
*             ! information is known)
         control_read = .true.
 
      end if
c
      end
 
