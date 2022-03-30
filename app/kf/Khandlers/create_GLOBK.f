CTITLE CREATE_GLOBK
 
      subroutine create_GLOBK( CO_name )

      implicit none
c
c     This routine creates 'GLOBK common' with the name in CO_name.
c
c     Here we set 'glb_con_read' true so that the RW_GLOBK_Block routine
c     knows that the record lengh information is available
 
c
c     Include files
c     -------------
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
c
      integer*4 total_blocks, err
 
c
      integer*4 num_con_words, num_mar_words, num_ema_words
 
* MOD TAH 190511: Changed to I*8 for 64-bit memory
      integer*8 AddressOf
 
c
*   full_len    - Length of the full Globk file name.
*   root_len    - Length of the root part of the Globk file name
*   trimlen     - HP utility to return length of string
 
      integer*4 full_len, trimlen
 
*   CO_name     - Name of the GLOBK common file
 
      character*(*) CO_name
 
*   full_globk_filename     - the full name of the Globk file
*                           - including the file type and size.
 
 
      character*128 full_globk_filename
 
c.... Get info about the size of the common blocks (NOTE: we remove
c     one extra word because we do not need to save last_@_word for
c     each block)
 
      num_con_words    = (AddressOf(last_glb_control)
     .                 - AddressOf(glb_control) - 1)/4 + 1
 
      num_mar_words    = (AddressOf(last_glb_markov)
     .                 - AddressOf(glb_markov) - 1)/4 + 1
 
c     num_ema_words    = EmaSize(glb_ema, last_glb_ema)
      num_ema_words    = (AddressOf(last_glb_ema) -
     .                    AddressOf(glb_ema) - 1)/4 + 1
 
c
c.... Number of blocks for the values common (NOTE: we do not add extra
c     block becuase there is dummy space for odd amount of block at the
c     end of each block.)
      gnum_control_sec= (num_con_words - 1) / 128 + 1
 
*     Get the start of the markov block and the number of blocks in it
      grec_markov_sec = gnum_control_sec + 1
      gnum_markov_sec = (num_mar_words-1) / 128 + 1
 
*     Get the start of the ema blovk and the number of blocks in it
      grec_ema_sec    = gnum_markov_sec + grec_markov_sec
      gnum_ema_sec    = (num_ema_words-1) / 128 + 1
 
c.... Total number of blocks in the file
      total_blocks = grec_ema_sec + gnum_ema_sec - 1
c
c.... Create the full file name
 
*                                         ! use the default name
      if( trimlen(CO_name).eq.0 ) then
          CO_name = glb_com_default
      end if
 
      call FullFileName(CO_name, 1,total_blocks, 128,
     .                  full_globk_filename )
 
      full_len = trimlen(Full_globk_filename)
 
      call FmpOpen(glb_com_dcb, err, full_globk_filename, 'RWOC', 1)
 
c.... Error?
      if( err .lt. 0 ) then
c
c....   Say not OK
        call report_error('FmpOpen', err,'creat',full_globk_filename,
*                                             ! Kill if error
     .                    1,'GLOBK_CREATE')
c
*             ! Indicate control read (This ensures that the block size
      else
*             ! information is known)
         glb_con_read = .true.
 
      end if
c
      end
 
