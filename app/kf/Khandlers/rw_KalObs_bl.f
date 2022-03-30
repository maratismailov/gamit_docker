CTITLE RW_KALOBS_BLOCK
 
      subroutine rw_KalObs_block( option, block_name, block_address,
     .                            ierr, ob )
 
 
      implicit none

 
*     Routine to read or write a logical block of the KalObs file.
*     These blocks are:
*     VALUES -- Basic control information
*     NAMES  -- Names of various things in data file
*     APR    -- apriories values for things in data file
*     DATA   -- the data records. For this block 'ob' is used for
*               the data (logical) record number
*
*     The First record of VALUES must be available when this routine
*     is called so that the lengths of the blocks and the start record
*     numbers of the blocks can be found.
*
*                                  9:15 AM  MON., 16  FEB., 1987
*
      include '../includes/kalman_param.h'
      include '../includes/obs_values.h'
 
*   block_address(1)    - Start word for reading/writing file
*               - (assumed to correspond to the correct variable
*               - name)
*   ierr        - FMGR  error flag
*   len         - Length of record returned by readd
*   num_words   - number of words to read
*   ob          - DATA logical record number (only used for 'DATA'
*               - block
*   start_rec   - Starting record number to read
*   version_date(5) - The date for the current version of the KalObs
*               - files.
 
      integer*4 block_address(1), ierr, len, num_words, ob,
     .    start_rec, version_date(5)
 
*   seconds     - Dummy seconds tag for version date
*   version_epoch   - JD of current version of files.
 
      real*8 seconds, version_epoch
 
*   option      - Option R or W.
*   opt         - Option casefolded
 
      character*1 option, opt
 
*   bname       - Shortened and casefolded block name
 
      character*2 bname
 
*   block_name  - Name of block to be read/written
 
      character*(*) block_name
 
*   block_found - Check to see if block name found
 
      logical block_found
 
***** Firstly see if values block has been record (this is set true either
*     when KalObs is opened or created.  This ensures the record length
*     and start records are available
 
      if( .not.values_read ) then
          call readd( ko_dcb, ierr, values, 128, len, 1)
          values_read = .true.
 
*         Now check that the version of this file is OK (i.e. later than
*         version date
          version_date(1) = version_year
          version_date(2) = version_month
          version_date(3) = version_day
          version_date(4) = version_hour
          version_date(5) = 0
          seconds         = 0.d0
 
          call ymdhms_to_jd( version_date, seconds, version_epoch)
*                                                 ! OBSOLETE file
          if( version_epoch.gt.read_epoch) then
*                       ! Either user's terminal or system console,
              write( *  ,100) data_base, version
  100         format(/' ***DISASTER*** This is an obselete version',
     .                ' of the KalObs Files.',/,
     .                ' Re-read data base ',a,' version ',i3,
     .                ' with the latest version of READIN',/)
 
              stop ' OBSELETE Version of KalObs file. Terminating'
          end if
      end if
 
***** Now for the block name get the start record and length
 
      opt = option
      bname = block_name
 
      call casefold(opt)
      call casefold(bname)
 
      block_found = .false.
 
*                                 ! VALUES block
      if( bname.eq.'VA' ) then
          start_rec = 1
          num_words = num_values_blocks*128
          block_found = .true.
      end if
 
*                                 ! NAMES  block
      if( bname.eq.'NA' ) then
          start_rec = start_names
          num_words = num_names_blocks*128
          block_found = .true.
      end if
 
*                                 ! APR block
      if( bname.eq.'AP' ) then
          start_rec = start_apr
          num_words = num_apr_blocks*128
          block_found = .true.
      end if
 
*                                 ! DATA block
      if( bname.eq.'DA' ) then
          start_rec = start_data + (ob-1)*obs_rec_len
          num_words = obs_rec_len*128
          block_found = .true.
      end if
 
***** Now if the block found, then read or write file
 
      IF( block_found ) THEN
 
*                                     ! READ block
          if( opt.eq.'R' ) then
              call readd( ko_dcb, ierr, block_address, num_words, len,
     .                    start_rec )
          endif
 
*                                     ! Write block
          if( opt.eq.'W' ) then
              call writd( ko_dcb, ierr, block_address, num_words,
     .                    start_rec )
          endif
 
          if( opt.ne.'W' .and. opt.ne.'R' ) then
              call bad_option(option, 'RW_KALOBS_BLOCK')
              ierr = -1007
          end if
 
*                         ! Unknown block name given
      ELSE
 
          call bad_option(block_name, 'RW_KALOBS_BLOCK')
          ierr = -1006
      END IF
 
***** Thats all
      return
      end
 
