CTITLE CREATE_KALOBS
 
      subroutine create_KalObs

      implicit none
c
c     This routine creates 'KalObs' with the name in KalObs_name
c
c MOD TAH 8701028 Modified to use FmpOpen rather than OPEN directly
C     This allows both CI and FMGR cartridges to be accessed and
c     avoid converting the file name to integer.
c     NOTE: The full file name should be entered in the database list file.
c
c MOD TAH 891117 Modified alrgorithm such that all integers are
c     integer*4
c
c     Here we set 'VALUES_READ' true so that the RW_KalObs_Block routine
c     knows that the record lengh information is available
 
c
c     Include files
c     -------------
      include '../includes/kalman_param.h'
      include '../includes/readin_user.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
c
      integer*4 total_blocks, err, i
 
c
      integer*4 num_obs_words, num_val_words
 
      integer*4 num_apr_words, num_names_words

* MOD TAH 190511: Changed to I*8 for 64-bit memory 
      integer*8 AddressOf
 
c
*   full_len    - Length of the full KalObs file name.
*   trimlen     - HP utility to return length of string
 
      integer*4 full_len, trimlen
 
*   full_kalobs_filename    - the full name of the KalObs file
*                           - including the file type and size.
      character*128 full_kalobs_filename
 
c*******cek modifications 88/6/14 to zero DATA block::
c   /DATA_BLOCK/ common block in include file OBS_DATA.FTNI
c
*            save_filename - local variable to save
      character*128 save_filename
 
c                              full_kalobs_filename during zeroing
C                              just in case
      integer*4 zero_data(1)
 
c
      equivalence (zero_data(1),site)
c
c******end of modifications for now**********************
c.... Get info about the size of the common blocks (NOTE: we remove
c     one extra word because we do not need to save last_@_word for
c     each block)
c     AddressOf returns the byte address of the start of a word.
c     Simce our "arrays" are now integer*4, we need to divide these
c     values by 4.  num_?_words will be just enough words to cover
c     all bytes up to the one before last_?_word
      num_val_words    = (AddressOf(last_values_word)
     .                 - AddressOf(values) - 1)/4 + 1
 
      num_names_words  = (AddressOf(last_names_word)
     .                 - AddressOf(names) - 1)/4 + 1
 
      num_apr_words    = (AddressOf(last_apr_word)
     .                 - AddressOf(aprioris) - 1)/4 + 1
 
      num_obs_words    = (AddressOf(last_data_word)
     .                 - AddressOf(site) - 1)/4 + 1
c
c  cek modification 88/6/14 zero out data block
c
      save_filename=full_kalobs_filename
c
c
      do i=1,num_obs_words
        zero_data(i)=0
      end do
c
      full_kalobs_filename=save_filename
c
c  end of modifications cek 88/6/14
c
c
c.... Number of blocks for the values common (NOTE: we do not add extra
c     block becuase there is dummy space for odd amount of block at the
c     end of each block.)
      num_values_blocks= (num_val_words - 1) / 128 + 1
 
*     Get the start of the names block and the number of blocks in it
      start_names      = num_values_blocks + 1
      num_names_blocks = (num_names_words-1) / 128 + 1
 
*     Get the start of the aprioris block and the number of blocks in it
      start_apr        = start_names + num_names_blocks
      num_apr_blocks   = (num_apr_words - 1) / 128 + 1
 
*     Get the start of the data blocks and the number of blocks
      start_data       = start_apr   + num_apr_blocks
      obs_rec_len      = (num_obs_words - 1) / 128 + 1
c
c.... Total number of blocks in the file
      total_blocks = start_data + num_obs * obs_rec_len - 1
c
c.... Create the full file name
 
      call FullFileName(KalObs_name, 1,total_blocks, 128,
     .                  full_KalObs_filename )
 
 
      full_len = trimlen(Full_KalObs_filename)
 
      call FmpOpen(ko_dcb, err, full_kalobs_filename, 'CW', 1)
*                                                     ! Create/write access
 
c.... Error?
      if (err .ge. 0) then
c
c....   Say OK
        write( *   ,150) full_kalobs_filename(1:full_len)
        if( log_lu.ne. 6 .and. log_lu.ne.0  )
     .  write(log_lu,150) full_kalobs_filename(1:full_len)
  150   format(/,1X,A,' successfully created.')
c
c....   Set created flag
        created =     .true.
        values_read = .true.
c
      else
c
c....   Say not OK
        write( *   ,200) err,full_KalObs_filename(1:full_len)
        if( log_lu.ne. 6 .and. log_lu.ne.0 )
     .  write(log_lu,200) err,full_KalObs_filename(1:full_len)
  200   format(/,' FMGR',I4.3,' on create of ',A,
     .         /,' Continue to next experiment.')
c
c....   Set created flag to false
        created = .false.
c
      end if
c
      end
 
