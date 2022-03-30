CTITLE RW_SOLVK_BLOCK
 
      subroutine rw_solvk_block( option, block_name, block_address,
     .                            ierr)

      implicit none
 
*     Routine to read or write a logical block of the SOLVK common file.
*     These blocks are:
*     CONTROL   -- CONTROL part used for the filter programs
*     MARKOV    -- Contains the markov parameters and parameter number
*
*     The First record of CONTROL must be available when this routine
*     is called so that the lengths of the blocks and the start record
*     numbers of the blocks can be found.
*
*                                  9:15 AM  MON., 16  FEB., 1987
*
      include '../includes/kalman_param.h'
      include '../includes/solvk_cntl.h'
 
*   block_address(1)    - Start word for reading/writing file
*               - (assumed to correspond to the correct variable
*               - name)
*   ierr        - FMGR  error flag
*   len         - Length of record returned by readd
*   num_words   - number of words to read
*   start_rec   - Starting record number to read
 
      integer*4 block_address(1), ierr, len, num_words, start_rec
 
*   option      - Option R or W.
*   opt         - Option casefolded
 
      character*1 option, opt
 
*   bname       - Shortened and casefolded block name
 
      character*2 bname
 
*   block_name  - Name of block to be read/written
 
      character*(*) block_name
 
*   com_copy    - Copy of the common file name so that we will
*               - will not overwrite it when the file is read.
 
      character*128 com_copy
 
*   block_found - Check to see if block name found
 
      logical block_found
 
***** Firstly see if values block has been record (this is set true either
*     when SOLVK common is opened or created.  This ensures the record length
*     and start records are available
 
      if( .not.control_read ) then
          call readd( isoldc, ierr, control, 128, len, 1)
          control_read = .true.
      end if
 
***** Now for the block name get the start record and length
 
      opt = option
      bname = block_name
 
      call casefold(opt)
      call casefold(bname)
 
      block_found = .false.
 
*                                 ! CONTROL block
      if( bname.eq.'CO' ) then
          start_rec = 1
          num_words = num_control_blocks*128
          block_found = .true.
      end if
 
*                                 ! MARKOV block
      if( bname.eq.'MA' ) then
          start_rec = start_markov
          num_words = num_markov_blocks*128
          block_found = .true.
      end if
 
***** Now if the block found, then read or write file
 
      IF( block_found ) THEN
 
*                                     ! READ block
          if( opt.eq.'R' ) then
              com_copy = com_file
              call readd( isoldc, ierr, block_address, num_words, len,
     .                    start_rec )
              com_file = com_copy
          endif
 
*                                     ! Write block
          if( opt.eq.'W' ) then
              call writd( isoldc, ierr, block_address, num_words,
     .                    start_rec )
          endif
 
          if( opt.ne.'W' .and. opt.ne.'R' ) then
              call bad_option(option, 'RW_SOLVK_BLOCK')
              ierr = -1007
          end if
 
*                         ! Unknown block name given
      ELSE
 
          call bad_option(block_name, 'RW_SOLVK_BLOCK')
          ierr = -1006
      END IF
 
***** Thats all
      return
      end
 
