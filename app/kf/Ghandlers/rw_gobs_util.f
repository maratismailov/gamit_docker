CTITLE CLOSE_GOBS
 
      subroutine close_Gobs(ierr)

      implicit none
 
c     Close the Gobs file.  Here we check to make sure that the
*     current ephead block is current.  If it is not then we write
*     it out before closing the file.
 
      include '../includes/kalman_param.h'
      include '../includes/gobs_def.h'
 
*   FmpClose            - HP function to close file, return error code
 
      integer*4 ierr, FmpClose
 
****  See if we need to update the ephead block (write it to disk)
*     before we close the file
 
      if( .not.ephead_current .and. curr_epoch.ne.0 ) then
 
*         Write out the current ephead
          call rw_Gobs_blk( go_dcb, 'WRITE', 'EPHEAD', ephead,
     .            num_header_recs, num_epoch_recs, num_ephead_recs,
     .            num_eph_only_recs, num_data_recs, epochs, 
     .            curr_epoch, ierr )
          call report_error('Fmp Error',ierr,'writ','Gobs ephead',
     .                    0,'RW_GOBS_EPHEAD')
      end if
 
*     Now close the file.
 
      ierr = FmpClose(go_dcb)
 
      call report_error('FmpClose',ierr,'clos','Gobs file',0,
     .                  'CLOSE_GOBS')
 
      header_read = .false.
      curr_epoch = 0
      curr_obs   = 0
 
      end
 
CTITLE CREATE_GOBS
 
      subroutine create_Gobs( Gobs_name, num_ep, num_rcv, num_obs,
     .                        ierr)

      implicit none
c
c     This routine creates a SOLVG data file.  It computes the sizes
c     of the various blocks in the file and computes the record numbers
c     needed to address each block.
c     When this routine is called the total number of epochs must be
c     known (The value is in the Gob_header num_gepochs)
c
c     Here we set 'VALUES_READ' true so that the RW_Gobs_Blk routine
c     knows that the record lengh information is available
 
c
c     Include files
c     -------------
      include '../includes/kalman_param.h'
      include '../includes/gobs_def.h'
c
* PASSED VARIABLES
 
* Gobs_name  - Name of the Gobs file
* num_ep     - Number of epochs (must be correct)
* num_rcv    - Number of receivers (must be correct)
* num_obs    - Number of one-way data records.  (Need not be
*              correct, used just to get an estimate of the file
*              size.
* ierr       - Error return from open.  IOSTAT Error number.
 
      character*(*) Gobs_name
 
      integer*4 num_ep, num_rcv, num_obs, ierr
 
* LOCAL VARIABLES
 
* total_blocks - Total number of GO_RECL words in the files
*
*   total_blocks   - Total number of GO_RECL words in the files
*   total_mbytes   - Total size of file in Mbytes.
*   num_header_words    - Number of i*4 words in the header
*   num_epoch_words     - Number of i*4 words in the epoch block
*   num_ephead_words    - Number of I*4 words in the epoch
*                         header block
*   num_eph_buff_recs   - Number of GO_RECL length records in 
*                         in max_grcv dependent records
*   start_epochs        - Start record number for the epochs block
*   start_ephead        - Start record number for the ephead block
*   num_data_words      - Number of I*4 words in the one-way
*                         data records
 
*   trimlen             - Length of string routine

      real*8 total_mbytes 
      integer*4 total_blocks, num_header_words, num_epoch_words,
     .    num_ephead_words, num_eph_buff_recs, 
     .    num_data_words, trimlen,
     .    start_epochs, start_ephead

      integer*8 AddressOf   ! Function for address of variable
 
*   full_len    - Length of the full Gobs file name.
 
 
      integer*4 full_len
 
*   full_gobs_name    - the full name of the Gobs file
*                             including the file type and size.
 
      character*256 full_gobs_name
 
c.... Get info about the size of the common blocks (NOTE: we remove
c     one extra word because we do not need to save last_@_word for
c     each block).  The only exception to this is the num_epoch_words
c     which is computed from the number of epochs
c     AddressOf returns the byte address of the start of a word.
c     Simce our "arrays" are now integer*4, we need to divide these
c     values by 4.  num_?_words will be just enough words to cover
c     all bytes up to the one before last_?_word
      num_header_words    = (AddressOf(last_header_word)
     .                    -  AddressOf(header) - 1)/4 + 1
 
      num_epoch_words  = num_ep + 1
 
      num_ephead_words    = (AddressOf(last_ephead_word)
     .                    -  AddressOf(ephead) - 1)/4 + 1

      num_eph_buff_words  = (AddressOf(last_eph_rcv) 
     .                    -  AddressOf(medium_flag)-1)/
     .                       (4*max_grcv)
 
      num_data_words    = (AddressOf(last_data_word)
     .                  -  AddressOf(site) - 1)/4 + 1
c
c
c.... Number of records needed for each block type in the file.
c     (NOTE: we do not add extra block becuase there is dummy space
c     for odd amount of block at the end of each block.)
 
*     Header records
      num_header_recs = (num_header_words - 1) / GO_RECL + 1
 
*     epochs records
      num_epoch_recs = (num_epoch_words-1) / GO_RECL + 1
 
*     epoch header records
      num_eph_only_recs   = (num_ephead_words - 1) / GO_RECL + 1
      num_eph_buff_recs   = (num_eph_buff_words -1 )/ GO_RECL + 1 
      num_ephead_recs     = num_eph_only_recs +
     .                      num_rcv * num_eph_buff_recs 
 
*     one-way data records
      num_data_recs      = (num_data_words - 1) / GO_RECL + 1
 
*     Now get the start record numbers in the file based on the sizes
      start_epochs      = num_header_recs + 1
      start_ephead      = start_epochs + num_epoch_recs

* DEBUG:
      write(*, 90) num_rcv, num_ep, num_obs, go_recl, 
     .             num_header_recs, num_epoch_recs,
     .             num_eph_only_recs, num_eph_buff_recs,
     .             num_ephead_recs, num_ephead_recs*num_ep,
     .             num_data_recs , num_data_recs * num_obs
  90  format(' GOBS Size statistics:',/,
     .       ' Num receivers ',i4,' Num epochs ',i4,
     .       ' Num one-ways  ',i6,' GO_RECL ',I4,/,
     .       ' Records in header   ',i4,/,
     .       ' Records in epochs   ',i4,/,
     .       ' Records in eph only ',i4,/,
     .       ' Records in eph buff ',i4,/
     .       ' Records in ephead   ',i4,' Total ',i8,/,
     .       ' Records in data     ',i4,' Total ',i8   )

*     Note: we do not have a start_data record because this will need to
*           be computed based on number of ephead's and data records
*           which occurr before any data record
 
c
c.... Approximate Total number of records in the file.  This just used
*     get some idea of the total size (it should be an upper bound).
      total_blocks = start_ephead + num_ep*num_ephead_recs +
     .                              num_obs*num_data_recs
      total_mbytes = total_blocks*GO_RECL*4/(1024.d0*1024.d0)
c
c.... Create the full file name
 
      call FullFileName(Gobs_name, 2,total_blocks, GO_RECL,
     .                  full_gobs_name )
 
      full_len = trimlen(full_gobs_name)
      write(*,100) full_gobs_name(1:full_len), total_mbytes
 100  format(' Creating ',a,' at ',F7.2,' Mbytes' )
 
*     An error has occurred here if the return is less than zero.
      if( trimlen(Gobs_name).gt.0 ) then 
          call FmpOpen(go_dcb, ierr, full_gobs_name, 'CW', 1)
*                                                 ! Create/write access
 
c....     Error?
          if (ierr .ge. 0) then
c
c....         Say OK
              write( *   ,150) full_gobs_name(1:full_len)
  150         format(/,1X,A,' successfully created.')
c
c....         Set the header read flag so that the RW_ routines will not
*             try to automatically read it.
              header_read = .true.
              curr_epoch = 0
*             Set the errorto zero to show all is OK
              ierr = 0
c
          else
c
c....         Say not OK
              call report_error('FmpOpen',ierr,'creat',full_gobs_name,
     .            0,'Create_Gobs')
c
          end if
      end if

****  Thats all
      return	
c
      end
 
CTITLE OBS_TO_EPNUM
 
      subroutine obs_to_epnum( obs, epochs, num_ep, ep, ierr)
 
      implicit none

*     Routine to scan the list of epochs and return the epoch number
*     to which one-way observation obs belongs to.  If the epoch can
*     not be found then the epoch number is returned as -1.
 
* PASSED VARIABLES
 
*   obs     - Observation number
*   num_ep      - Number of epochs (used to make sure that we
*               - scan past the end of the epochs array
*   epochs(num_ep)  - list of observations numbers at the start
*                     of each epoch
 
*   ep      - Epoch number corresponding to obs.
*   ierr    - Error flag from trying to find epoch.
*           - Returns -1002 if epoch number can not be found.
 
      integer*4 obs, num_ep, epochs(num_ep), ep, ierr
 
* LOCAL VARIABLES
 
*       epoch_found     - Indicates that the epoch number has
*           - been found
 
 
 
      logical epoch_found
 
****  loop over the epochs array to find which epoch number our
*     data record is in.
      epoch_found = .false.
      ep = 0
      do while ( .not.epoch_found .and. ep.le.num_ep)
          ep = ep + 1
          if ( obs.lt.epochs(ep+1) ) then
*             We have have found the epoch number
              epoch_found = .true.
          end if
      end do
 
*     See if found epoch number
****  Make sure that we found the epoch number.  If we did not
*     then report error and exit from this routine.
      if( .not.epoch_found ) then
 
*         Tell user we have a problem
          write(*,100) obs
 100      format('**ERROR** Finding the epoch number corresponding',
     .           ' to data record ',i5,/,
     .           '          Not reading this data record')
          ierr = -1002
      end if
 
 
***** Thats all
      return
      end
 
 
 
CTITLE OPEN_GOBS
 
      subroutine open_Gobs(Gobs_name,ierr)

      implicit none
 
c     Routine to open the Gobs files using the File Management Protocol
c     FMP routines.  On open the error on open is only reported here
c     it is the higher level subroutines responsibility to take action
c     in response to the error.
 
c     The header_read variable is set .false. here to force the header
c     block to be read latter.
 
      include '../includes/kalman_param.h'
      include '../includes/gobs_def.h'
 
 
 
      character*(*) Gobs_name
 
*   ierr            - FMP error flag
 
      integer*4 ierr

*    full_gobs_name - Full name of file containing the record
*                     length and file type

      character*256 full_gobs_name
 
***** Open file using the Fmp routine (ensures that the data control
*     block (dcb) is set up correctly.
*     Add the file type and record length to the file name
      call FullFileName(Gobs_name, 2, 0, GO_RECL, full_gobs_name ) 
      call FmpOpen(go_dcb, ierr, full_gobs_name, 'RWO', 1)
 
      call report_error('FmpOpen', ierr,'open',full_gobs_name,
     .                   0,'OPEN_GOBS')
 
***** Set header_read false
 
*                             ! Force the first record of file to be
      header_read = .false.
*                             ! read, so that we will get the record
*                             ! lengths and start record numbers in
*                             ! RW_Gobs_Blk
      curr_epoch = 0
      curr_obs   = 0
      end
 
CTITLE RW_GOBS_BLK
 
      subroutine rw_Gobs_blk( go_dcb, option, block_name,
     .       block_address, num_header_recs, num_epoch_recs,
     .       num_ephead_recs, num_eph_only_recs, 
     .       num_data_recs, epochs, ent, ierr )
 

      implicit none
 
*     Routine to read or write a logical block of the Gobs file.  This
*     routine is set up so that the records to be written or the
*     destination of the records to be read do not need to be in the
*     standard common block locations.  (Most of the time they will be
*     but these routines will allow multiple gobs files to be read at
*     same time.)
*
*     These blocks are:
*     HEADER - The header block of the Gobs file
*     EPOCHS - The epochs block (contains a list the data record numbers
*              at the start of each block)
*     EPHEAD - Epoch header block
*     DATA   - One-way data record block
*
*     NOTE: This is routine does not ensure that all necessary parts of
*           the gobs files are read e.g. when a data record is read,
*           corresponding ephead records are not read
* INCLUDES

      include '../includes/kalman_param.h'
*
* PASSED VARIABLES
 
*   go_dcb(16)  - Data control block for the fmp file routines
*   block_address(*) - Variable which is located that the start of
*               - memory to be used to read the records into
*   num_header_recs - Number of records in the header
*   num_epoch_recs - Number of records in the epochs list
*   num_ephead_recs - Number of records in the epoch header
*   num_eph_only_recs - Number of records in just the fixed part
*                     of the ephoch headers.
*   num_data_recs   - number of records in the one-way data
*                     records
*   epochs(num_epoch_recs*GO_RECL) - List of the first observations
*                     numbers in each epoch.  Can be dummy array
*                     if the header or epochs blocks are read.
*                     Must be valid when ephead or data blocks
*                     are read.
*   ent             - Either epoch number when ephead block read
*                     or observation number when the data records
*                     are read.
*   ierr            - Fmp read routine error flag return.
 
      integer*4 go_dcb(16), block_address(*), num_header_recs,
     .    num_epoch_recs, num_ephead_recs, num_eph_only_recs, 
     .    num_data_recs, epochs(num_epoch_recs*GO_RECL), ent, ierr
 
*   option          - Option to either READ or WRITE block (only
*                     first character is checked)
*   block_name      - Name of block to operated on.  Upto 4
*                     characters of this option are used.
 
      character*(*) option, block_name
 
* LOCAL VARIABLES
 
*   len         - Length of record returned by readd
*   num_words   - number of words to read
*   start_rec   - Starting record number to read
*   ep          - Epoch number corresponding to observation number
*                 passed
 
 
      integer*4 len, num_words,  start_rec, ep
 
*   opt         - Option casefolded
 
 
      character*4 opt
 
*   bname       - Shortened and casefolded block name
 
 
      character*4 bname
 
*   block_found - Check to see if block name found
 
      logical block_found
 
***** Now for the block name get the start record and length
 
      opt = option
      bname = block_name
 
      call casefold(opt)
      call casefold(bname)
 
      block_found = .false.
 
*     Find the block name and compute the starting record number in
*     the file
 
*                                 ! HEADER block
      if( bname(1:2).eq.'HE' ) then
          start_rec = 1
          num_words = num_header_recs*GO_RECL
          block_found = .true.
      end if
 
*                                     ! EPOCHS block
      if( bname(1:4).eq.'EPOC' ) then
          start_rec = 1 + num_header_recs
          num_words = num_epoch_recs*GO_RECL
          block_found = .true.
      end if
 
*                                     ! EPHEAD block
      if( bname(1:4).eq.'EPHE' ) then
 
*         The start record for this type depends on which epoch number
*         we want to read.  Need to account for the start of the file
*         and all the ephead and data records upto the point we want
          start_rec = 1 + num_header_recs + num_epoch_recs +
     .                (ent-1)*num_ephead_recs +
     .                (epochs(ent)-1)*num_data_recs

*         Since the ephead block is made up of a constant length
*         block and then one buffer record per reciever we only
*         set the words here for the fix length portion.
          num_words = num_eph_only_recs*GO_RECL
          block_found = .true.
      end if
 
*                                     ! One-way DATA record
      if( bname(1:2).eq.'DA' ) then
 
*         Again here we need to accound for all of the ephead and
*         data records to the point we want.  In this case we also need
*         to find which epoch number the one-way data record is
*         associated with.
          call obs_to_epnum( ent, epochs, num_epoch_recs*GO_RECL,
     .                       ep, ierr)
          if( ierr.ne.0 ) RETURN
 
*         Now that we have found epoch number, compute the record
*         number of the observation
          start_rec = 1 + num_header_recs + num_epoch_recs +
     .                ep*num_ephead_recs +
     .                (epochs(ep)-1)*num_data_recs +
     .                (ent-epochs(ep))*num_data_recs
          num_words = num_data_recs*GO_RECL
          block_found = .true.
      end if
 
***** If we found the block name then read/write the records
      IF( block_found ) THEN
 
*                                     ! READ block
          if( opt(1:1).eq.'R' ) then
              call readd( go_dcb, ierr, block_address, num_words, len,
     .                    start_rec )
*                                           ! WRITE the block
          else if( opt(1:1).eq.'W' ) then
              call writd( go_dcb, ierr, block_address, num_words,
     .                    start_rec )
          else
*             Bad operation passed to this routine
              write(*,200) option
 200          format('**ERROR** Bad operation passed to RW_GOBS_BLK',
     .                ' Valid operations are READ or WRITE. ',a4,
     .                ' Passed ',/,
     .                '          No operation performed')
              ierr = -1003
          end if
*                         ! Unknown block name given
      ELSE
          write(*,250) option
 250      format('**ERROR** Bad operation passed to RW_GOBS_BLK',
     .            ' Valid operations are READ or WRITE. ',a4,
     .            ' Passed ',/,
     .            '          No operation performed')
          ierr = -1004
      END IF
 
***** Thats all
      return
      end
 
CTITLE RW_GOBS_DATA
 
      subroutine rw_Gobs_data(option, ob, ierr)

      implicit none
c
c     Routine to read or write the Gobs data records.  OPTION
c     is a character variable.  If OPTION(1:1) = 'W', the record
c     will be written, if 'R', read.  Assumes Gobs is open.  This
*     will read into the standard common locations for the data.
c
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/gobs_def.h'
 
* PASSED VARIABLES
 
*   ob          - observation number to read or write
*   ep          - Epoch number for this observation.
*   ierr        - FMP error routine returns
 
      integer*4 ob, ierr, ep
 
*   option      - Option: either READ or WRITE (only first
*               - character checked.
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   opt         - Local version of opt so that we can casefold
 
      character*4 opt
 
c
c.... Get first character of option
      opt = option
      call CaseFold(opt)
c
*     See if the first record of the header block has been read.  If it
*     has not then read then we can write any block so report error and
*     return.
 
      if( .not.header_read .and. opt(1:1).eq.'W' ) then
          write(*,100)
 100      format('**ERROR** Attempt to write a record of a Gobs file',
     .            ' with header not read',/,
     .            '         File NOT written')
          ierr = -1005
          RETURN
      end if
 
****  See if we should read the header and the epochs blocks (normally
*     this would be done reading data but it will work this way.
 
      if( .not.header_read ) then
 
*         Warning user in case not aware of the problem
          write(*,150)
 150      format('**WARNING** Attempt to read DATA before read Gobs',
     .            ' file HEADER',/,
     .            '            Reading header and epochs now')
          call rw_gobs_blk( go_dcb, 'READ', 'HEADER', header,
     .        num_header_recs, num_epoch_recs, num_ephead_recs, 
     .        num_eph_only_recs, num_data_recs, epochs, 0, ierr )
          if( ierr.eq.0 ) then
              header_read = .true.
          else
              call report_error('Fmp Error',ierr,option,'Gobs Header',
     .                            0,'RW_GOBS_HEADER')
          end if
          call rw_gobs_blk( go_dcb, 'READ', 'EPOCHS', epochs,
     .        num_header_recs, num_epoch_recs, num_ephead_recs, 
     .        num_eph_only_recs, num_data_recs, epochs, 0, ierr )
          call report_error('Fmp Error',ierr,option,'Gobs epochs',
     .                0,'RW_GOBS_HEADER')
      end if
 
****  See if we are in same epoch as the last data record read.
 
      call obs_to_epnum( ob, epochs, num_gepochs, ep, ierr)
      if( ep.ne.curr_epoch ) then
 
*         The epoch will change with this read.  First see if we need
*         the current ephead before reading the next,
          if( .not.ephead_current .and. curr_epoch.gt.0 ) then
              call rw_gobs_ephead('WRITE', curr_epoch, ierr )
          end if
 
*         Now get the new epoch header
          call rw_gobs_ephead('READ', ep, ierr )
      end if
 
 
      call rw_Gobs_blk( go_dcb, opt, 'DATA', site,
     .        num_header_recs, num_epoch_recs,
     .        num_ephead_recs, num_eph_only_recs, 
     .        num_data_recs, epochs, ob, ierr )
      call report_error('Fmp Error',ierr,option,'Gobs ephead',
     .                0,'RW_GOBS_EPHEAD')
 
*     Save the current epoch number
      curr_obs = ob
 
*     If we just wrote the data record then indicate that the ephead
*     block should be updated (This is not necessarily so, but this way
*     we know it will always be updated),
      if( opt(1:1).eq.'W' ) then
          ephead_current = .false.
      end if
 
***** Thats all
      return
      end
 
CTITLE RW_GOBS_EPHEAD
 
      subroutine rw_Gobs_ephead(option, ep, ierr)

      implicit none
c
c     Routine to read or write the Gobs ephead records.  OPTION
c     is a character variable.  If OPTION(1:1) = 'W', the record
c     will be written, if 'R', read.  Assumes Gobs is open.  This
*     will read into the standard common locations for the data.
c
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/gobs_def.h'
 
* PASSED VARIABLES
 
*   ep          - Epoch number to read
*   ierr        - FMP error routine returns
 
      integer*4 ep, ierr
 
*   option      - Option: either READ or WRITE (only first
*               - character checked.
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   opt         - Local version of opt so that we can casefold
 
      character*4 opt
 
c
c.... Get first character of option
      opt = option
      call CaseFold(opt)
c
*     See if the first record of the header block has been read.  If it
*     has not then read then we can write any block so report error and
*     return.
 
      if( .not.header_read .and. opt(1:1).eq.'W' ) then
          write(*,100)
 100      format('**ERROR** Attempt to write a record of a Gobs file',
     .            ' with header not read',/,
     .            '         File NOT written')
          ierr = -1005
          RETURN
      end if
 
****  See if we should read the header and the epochs blocks (normally
*     this would be done reading ephead but it will work this way.
 
      if( .not.header_read ) then
 
*         Warning user in case not aware of the problem
          write(*,150)
 150      format('**WARNING** Attempt to read EPHEAD before read Gobs',
     .            ' file HEADER',/,
     .            '            Reading header and epochs now')
          call rw_gobs_blk( go_dcb, 'READ', 'HEADER', header,
     .        num_header_recs, num_epoch_recs, num_ephead_recs, 
     .        num_eph_only_recs, num_data_recs, epochs, 0, ierr )
          if( ierr.eq.0 ) then
              header_read = .true.
          else
              call report_error('Fmp Error',ierr,option,'Gobs Header',
     .                            0,'RW_GOBS_HEADER')
          end if
          call rw_gobs_blk( go_dcb, 'READ', 'EPOCHS', epochs,
     .        num_header_recs, num_epoch_recs, num_ephead_recs, 
     .        num_eph_only_recs, num_data_recs, epochs, 0, ierr )
          call report_error('Fmp Error',ierr,option,'Gobs epochs',
     .                0,'RW_GOBS_HEADER')
      end if
 
****  If we are going to read a new ephead block, make sure the current
*     one is current.  If it is not write it out before reading the
*     next one
 
      if( .not.ephead_current .and. curr_epoch.ne.0 ) then
 
*         Write out the current ephead
          call rw_gobs_blk( go_dcb, 'WRITE', 'EPHEAD', ephead,
     .            num_header_recs, num_epoch_recs, num_ephead_recs,
     .            num_eph_only_recs, num_data_recs, epochs, 
     .            curr_epoch, ierr )
          call report_error('Fmp Error',ierr,'writ','Gobs ephead',
     .                    0,'RW_GOBS_EPHEAD')
          call rw_rem_eph( 'WRITE', ierr)
          call report_error('Fmp Error',ierr,option,'Gobs remain eph',
     .                0,'RW_GOBS_EPHEAD')
      end if
 
 
****  Read or write the ephead block and the epochs block
 
      call rw_gobs_blk( go_dcb, opt, 'EPHEAD', ephead,
     .        num_header_recs, num_epoch_recs, num_ephead_recs, 
     .        num_eph_only_recs, num_data_recs, epochs, ep, ierr )
      call report_error('Fmp Error',ierr,option,'Gobs ephead',
     .                0,'RW_GOBS_EPHEAD')
      call rw_rem_eph( opt, ierr)
      call report_error('Fmp Error',ierr,option,'Gobs remain eph',
     .                0,'RW_GOBS_EPHEAD')
 
*     Save the current epoch number and set the status of the epoch
*     ephead block as current(i.e., does not need to be updated)
      curr_epoch = ep
      ephead_current = .true.
 
***** Thats all
      return
      end
 
CTITLE RW_REM_EPH     
 
      subroutine rw_rem_eph( option, ierr)

      implicit none
c
*     Routine to read or write the remaining the ephemeris 
*     header buffer entries with the recwiver dependent 
*     information in them.  (Each reciever gets one record)
c
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/gobs_def.h'
 
* PASSED VARIABLES
 
*   ierr        - FMP error routine returns
 
      integer*4 ierr
 
*   option      - Option: either READ or WRITE (only first
*               - character checked.
 
      character*(*) option
 
* LOCAL VARIABLES

*   i           - Loop counter
*   len         - Length of record read.

      integer*4 i, len
 
*   opt         - Local version of opt so that we can casefold
 
      character*4 opt
 
c
c.... Get first character of option
      opt = option
      call CaseFold(opt)
c
      if( opt(1:1).eq.'R' ) then

*         Loop over the records reading the buffer and saving
*         the values.
          do i = 1, num_grcv 
              call readd( go_dcb, ierr, eph_rcv_buff, 
     .                    num_eph_buff_words, len, 0)

*             now copy the values into the appropriate arrays.
**WARNING**   The following copy assumes that only integer*4
*             real*4's are in the gobs_ephead_rcv common  
*             block.
              call GO_buff_to_comm( eph_rcv_buff, medium_flag(i),
     .             num_eph_buff_words, max_grcv )
          end do
      else if( opt(1:1).eq.'W' ) then

*         Loop over the records writing the buffer and saving
*         the values.
          do i = 1, num_grcv 

*             Copy common to buffer (see warning above).
              call GO_comm_to_buff( eph_rcv_buff, medium_flag(i),
     .             num_eph_buff_words, max_grcv )
              call writd( go_dcb, ierr, eph_rcv_buff, 
     .                    num_eph_buff_words, 0)

          end do
      end if

****  Thats all
      return
      end

CTITLE GO_COMM_TO_BUFF

      subroutine GO_comm_to_buff( buffer, start, num_words, 
     .                            max_grcv )

      implicit none

*     Routine to copy the specific site dependent values from
*     common to the buffer for writing out to the Gobs files.
*
* INCLUDES
*  None

* PASSED VARIABLES

*  num_words - Number of words to copied.
*  max_grcv  - Size of the first dimension in start (must be
*              max_grv (max number of recievers)
*  buffer(num_words) - buffer for putting the specific site
*              values into.
*  start(max_grcv, num_words) - First word in the common block
*              for this sites.
*
      integer*4 num_words, max_grcv , buffer(num_words),
     .          start(max_grcv, num_words)

* LOCAL VARIABLES

*  i   - Loop counter over the values

      integer*4 i

***** Loop over the number of words copying from common area 
*    (passed in this call) to the buffer

      do i = 1, num_words
         buffer(i) = start(1,i)
      end do

***** Thats all
      return
      end

CTITLE GO_BUFF_TO_COMM

      subroutine GO_buff_to_comm( buffer, start, num_words, 
     .                            max_grcv )

      implicit none

*     Routine to copy the specific site dependent values from
*     from the buffer to the common area when the Gobs file 
*     is read.
*
* INCLUDES
*  None

* PASSED VARIABLES

*  num_words - Number of words to copied.
*  max_grcv  - Size of the first dimension in start (must be
*              max_grv (max number of recievers)
*  buffer(num_words) - buffer for putting the specific site
*              values into.
*  start(max_grcv, num_words) - First word in the common block
*              for this sites.
*
      integer*4 num_words, max_grcv , buffer(num_words),
     .          start(max_grcv, num_words)

* LOCAL VARIABLES

*  i   - Loop counter over the values

      integer*4 i

***** Loop over the number of words copying from buffer area 
*     to the common area (passed in this call).

      do i = 1, num_words
         start(1,i) = buffer(i)
      end do

***** Thats all
      return
      end

CTITLE RW_GOBS_HEADER
 
      subroutine rw_Gobs_header(option,ierr)

      implicit none
c
c     Routine to read or write the Gobs header records.  OPTION
c     is a character variable.  If OPTION(1:1) = 'W', the record
c     will be written, if 'R', read.  Assumes Gobs is open.  This
*     will read into the standard common locations for the data.
c
* INCLUDES
 
      include '../includes/kalman_param.h'
      include '../includes/gobs_def.h'
 
* PASSED VARIABLES
 
*   ierr        - FMP error routine returns
 
      integer*4 ierr
 
*   option      - Option: either READ or WRITE (only first
*               - character checked.
 
      character*(*) option
 
* LOCAL VARIABLES
 
*   opt         - Local version of opt so that we can casefold
 
      character*4 opt
 
c
c.... Get first character of option
      opt = option
      call CaseFold(opt)
c
*     See if the first record of the header block has been read.  If it
*     has not then read then we can write any block so report error and
*     return.
 
      if( .not.header_read .and. opt(1:1).eq.'W' ) then
          write(*,100)
 100      format('**ERROR** Attempt to write a record of a Gobs file',
     .            ' with header not read',/,
     .            '         File NOT written')
          ierr = -1005
          RETURN
      end if
 
****  Read or write the header block and the epochs block
      call rw_gobs_blk( go_dcb, opt, 'HEADER', header,
     .        num_header_recs, num_epoch_recs, num_ephead_recs, 
     .        num_eph_only_recs, num_data_recs, epochs, 0, ierr )
      if( ierr.eq.0 ) then
          header_read = .true.
          curr_epoch  = 0
          curr_obs    = 0
      else
          call report_error('Fmp Error',ierr,option,'Gobs Header',
     .                        0,'RW_GOBS_HEADER')
      end if
      call rw_gobs_blk( go_dcb, opt, 'EPOCHS', epochs,
     .        num_header_recs, num_epoch_recs, num_ephead_recs, 
     .        num_eph_only_recs, num_data_recs, epochs, 0, ierr )
      call report_error('Fmp Error',ierr,option,'Gobs epochs',
     .                0,'RW_GOBS_HEADER')
 
 
***** Thats all
      return
      end
 
 
