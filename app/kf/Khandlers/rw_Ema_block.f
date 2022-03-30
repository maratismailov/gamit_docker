CTITLE RW_EMA_BLOCK
 
      subroutine rw_Ema_block( option, block_name, block_address,
     .                         ierr)

      implicit none
 
 
 
*     Routine to read or write a logical ema block of the GLOBK common file.
*     There is only one such block:
*     EMA       -- EMA Part contains REAL*4 and REAL*8 variables
*
*     The First record of CONTROL must be available when this routine
*     is called so that the lengths of the blocks and the start record
*     numbers of the blocks can be found.
*
*                                  08:50 PM TUE., 28 Jul., 1987
*
      include '../includes/kalman_param.h'
      include '../includes/globk_cntl.h'
 
*   block_address(1)    - Start word for reading/writing file
*               - (assumed to correspond to the correct variable
*               - name)
*   ierr        - FMGR  error flag
*   len         - Length of record returned by VREAD
*   num_words   - number of words to read
*   start_rec   - Starting record number to read
 
      integer*4 block_address(1), ierr, len, num_words, start_rec
 
*   option      - Option R or W.
*   opt         - Option casefolded
 
      character*1 option, opt
 
*   block_name  - Name of block to be read/written
 
      character*(*) block_name
 
*   block_found - Check to see if block name found
 
      logical block_found
 
 
***** Firstly see if values block has been record (this is set true either
*     when GLOBK common is opened or created.  This ensures the record length
*     and start records are available
 
      if( .not.glb_con_read ) then
          call readd( glb_com_dcb, ierr, glb_control, 128, len, 1)
          glb_con_read = .true.
      end if
 
***** Now for the block name get the start record and length
 
      start_rec = grec_ema_sec
      num_words = gnum_ema_sec*128
*                                 ! There is only one block
      block_found = .true.
 
      opt = option
      call casefold( opt )
 
***** Now if the block found, then read or write file
 
      IF( block_found ) THEN
 
*                                     ! READ block
          if( opt.eq.'R' ) then
              call readd( glb_com_dcb, ierr, block_address, num_words,
     .                    len, start_rec )
          endif
 
*                                     ! Write block
          if( opt.eq.'W' ) then
              call writd( glb_com_dcb, ierr, block_address, num_words,
     .                    start_rec )
          endif
 
          if( opt.ne.'W' .and. opt.ne.'R' ) then
              call bad_option(option, 'RW_EMA_BLOCK')
              ierr = -1007
          end if
 
      END IF
 
***** Thats all
      return
      end
 
