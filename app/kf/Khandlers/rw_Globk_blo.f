CTITLE RW_GLOBK_BLOCK
 
      subroutine rw_Globk_block( option, block_name, block_address,
     .                            ierr)
 

      implicit none
 
 
*     Routine to read or write a logical block of the GLOBK common file.
*     These blocks are:
*     CONTROL   -- CONTROL part used for the filter programs
*     MARKOV    -- Contains the markov parameters and parameter number
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
*   len         - Length of record returned by readd
*   num_words   - number of words to read
*   start_rec   - Starting record number to read
*   trimlen     - Length of string
 
      integer*4 block_address(1), ierr, len, num_words, start_rec,
     .          trimlen
 
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
*     when GLOBK common is opened or created.  This ensures the record length
*     and start records are available
 
      if( .not.glb_con_read ) then
          com_copy = glb_com_file
          call readd( glb_com_dcb, ierr, glb_control, 128, len, 1)
          glb_com_file = com_copy
          glb_con_read = .true.
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
          num_words = gnum_control_sec*128
          block_found = .true.
      end if
 
*                                 ! MARKOV block
      if( bname.eq.'MA' ) then
          start_rec = grec_markov_sec
          num_words = gnum_markov_sec*128
          block_found = .true.
      end if
 
***** Now if the block found, then read or write file
 
      IF( block_found ) THEN
 
*                                     ! READ block
          if( opt.eq.'R' ) then
              call readd( glb_com_dcb, ierr, block_address, num_words,
     .                    len, start_rec )

* MOD TAH 961230: if CO read, check the version
              if( bname.eq.'CO' ) then
                  if( curr_glb_com_ver.ne.glb_com_ver ) then
* MOD TAH 190610: Allow 515 as well 516 (difference only affects the
*                     radition parameters so if no orbits all is OK
                      if( glb_com_ver.eq.515 .and. 
     .                    curr_glb_com_ver.eq.516 ) then
                          write(*,115) trim(com_copy), glb_com_ver/100.
 115                      format('**WARNING** Globk common file ',a,
     .                       ' previous version ',F5.2,'. OK except',
     .                       ' satellite rad parameters')
                      else

                          write(*,120) com_copy(1:trimlen(com_copy)), 
     .                             glb_com_ver/100., 
     .                             curr_glb_com_ver/100.
 120                      format('**ERROR** Globk common file ',a,
     .                         ' is version ',F5.2,/,
     .                         '          Current version is ',F5.2) 
                          call report_stat('FATAL','GLOBK',
     .                        'rw_globk_block', com_copy,
     .                        'Obsolete common file version', 
     .                        glb_com_ver)
                      endif
                  end if
              end if
          endif
 
*                                     ! Write block
          if( opt.eq.'W' ) then

* MOD TAH 961230: Save the current version for these files.
              if( bname.eq.'CO' ) glb_com_ver = curr_glb_com_ver
              call writd( glb_com_dcb, ierr, block_address, num_words,
     .                    start_rec )
          endif
 
          if( opt.ne.'W' .and. opt.ne.'R' ) then
              call bad_option(option, 'RW_GLOBK_BLOCK')
              ierr = -1007
          end if
 
*                         ! Unknown block name given
      ELSE
 
          call bad_option(block_name, 'RW_GLOBK_BLOCK')
          ierr = -1006
      END IF
 
***** Thats all
      return
      end
 
