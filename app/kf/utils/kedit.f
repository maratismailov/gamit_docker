      program kedit
 
 
*     Program to set the interactive SOLVK edit bit by reading an
*     input file of data to be edited and matching it to the data in
*     the KalObs file.
*
*     Runstring
*     CI> kedit KalObs_file input_list
*
*     where KalObs_file is the name of the KalObs file to be
*         updated
*           Input_list is the name of a file or Lu which contains the
*         the edited data in the form
*          Site 1 Site 2 edit epoch (yy mm day hr min)
*
*     with either site name being ALL will match all baselines
*     and <edit> is 1 or 0 depending on whether the bit is to be set or
*                reset.
*
*          max_interedits   - Maxium number of interactive edits
 
      integer*4 max_interedits
 
      parameter ( max_interedits = 2000 )
 
      include '../includes/kalman_param.h'
      include '../includes/obs_header.h'
      include '../includes/obs_data.h'
 
 
*   date(5)     - Dates and time read from input
*   i,j,k       - Loop counters
*   ierr        - IOSTAT error
*   iel,jel     - Site numbers from Get command
*   indx        - counter used in decoding lines
*   len_run     - Length of the runstring
*   num_interedits  - Number of edits entered by user
*   num_edits_in    - Number of data restored
*   num_edits_out   - Number of data removed
*   base_edits(2,max_interedits) - User sites for edits (value of
*               - 999999 for site number will match all sites)
*   type_edits(max_interedits)  - Type of edits (1 will set bit,
*               - 0 will reset bit)
*   trimlen     - HO function for length of string
*   rcpar       - Runstring reading function
*   lui         - Input file unit number
 
      integer*4 date(5), i,j,k, ierr, iel,jel, indx, len_run,
     .    num_interedits, num_edits_in, num_edits_out,
     .    base_edits(2,max_interedits), type_edits(max_interedits),
     .    trimlen, rcpar, lui
 
*   sec_tag     - Seconds tag
*   epoch_edits(max_interedits)  - Epochs of the user edits (these
*               - are excepted to be with +- one minute of value)
 
      real*8 sec_tag, epoch_edits(max_interedits)
 
*   cleared     - Indicates all edits are to be cleared
*   interactive - Indicates interactive input
*   kbit        - Checks bit status
*   write_rec   - Indicates that record should be written
 
      logical cleared, interactive, kbit, write_rec
 
*   input_file  - Name of the input file
*   KO_name     - Name of KalObs file
*   runstring   - Runstring read from command line
 
      character*64 input_file, KO_name, runstring
 
*   buffer      - Line read from user file
 
      character*80 buffer
 
****  Decode the runstring
 
      cleared = .false.
 
      len_run = rcpar(1,runstring)
      if( len_run.gt.0 ) then
          KO_name = runstring
      else
          call proper_runstring('kedit.hlp',6,1)
      end if
 
      len_run = rcpar(2,runstring)
      if( len_run.gt.0 ) then
          input_file = runstring
          interactive = .false.
      else
          input_file = '6'
          interactive = .true.
      end if
 
***** Try to open the KalObs file
      call open_KalObs( KO_name, ierr)
      if( ierr.ne.0 ) stop ' kedit: Error opening KalObs file'
 
      call rw_KalObs_header('R', ierr)
      call report_error('FmpRead',ierr,'read','KalObs header',1,
     .                  'kedit')
 
*     Try to open interactivee
      lui = 100
      call open_lu(lui, input_file, ierr, 'old')
      call report_error('IOSTAT',ierr,'open',input_file,1,
     .    'kedit')
 
***** Now get the epochs of the user supplied breaks
      if( interactive ) then
          write(*,100)
 100      format(/' KEDIT: Program to enter clock breaks',/,
     .            ' Give Baseline, type (0/1) and epoch ',
     .            ' (yy mm dd hr min) [End with <cntl D>]')
      end if
 
      do while ( ierr.eq.0 )
 
          if( interactive ) then
              write(*,120)
  120         format('?',$)
          end if
 
          read(lui,'(a)',iostat=ierr) buffer
          call casefold( buffer )
 
          if( ierr.eq.0 .and. buffer(1:1).eq.' ' .and.
*                                                            ! Process
     .        trimlen(buffer).gt.0             ) then
 
*             See if we can find name
              indx = 1
              call get_cmd(buffer, site_names, num_sites, iel, indx)
              call get_cmd(buffer, site_names, num_sites, jel, indx)
 
*                                                     ! We found the sites
              if( iel.gt.0 .and. jel.gt.0  ) then
*                                     ! Increment number
                  k = k + 1
*                                                 ! Too many
                  if( k.eq. max_interedits ) then
                      write(*,200) max_interedits
  200                 format(' kedit: Too many edits entered.',
     .                       ' (Maxiumim is ',i3,')')
                      ierr = -1
                  end if
 
                  base_edits(1,k) = iel
                  base_edits(2,k) = jel
 
*                 Get the type of edit
                  call read_line(buffer,indx,'I4', ierr, 
     .                           type_edits(k),' ')
 
*                 Now get date
                  call multiread(buffer,indx,'I4', ierr, date,' ',5)
                  call report_error('IOSTAT',ierr,'decod',buffer,0,
     .                              'kedit')
*                                         ! Save epoch
                  if( ierr.eq.0 ) then
                      sec_tag = 0.d0
                      call ymdhms_to_jd( date, sec_tag, epoch_edits(k))
*                                         ! Ignore this line
                  else
                      k = k - 1
                  end if
*                                         ! Tell use name not found
              else
 
*                 Check to see if CLEAR given
*                                                       ! Yes
                  if( index(buffer,'CLEAR').gt.0 ) then
                      cleared = .true.
                  else
                      write(*,210) buffer(1:max(1,trimlen(buffer)))
  210                format(' Sites not found in :',a)
                  end if
              end if
*                         ! No error readinf input
          end if
*                         ! Loop over the input
      end do
 
***** Now see if we can find observations close to the breaks
 
      if( k.eq.0 .and. .not.cleared ) then
          stop ' KEDIT: No valid edits breaks entered'
      end if
 
      num_interedits = k
      num_edits_in  = 0
      num_edits_out = 0
 
*                         ! Loop over KalObs file
      do i = 1, num_obs
 
          call rw_KalObs_block('R','DATA',site,ierr, i)
          call report_error('FmpRead',ierr,'read','KalObs file',1,
     .                      'kedit')
 
*         Check to see if this obs needs editing
          write_rec = .false.
 
*         See if all interactives are to be cleared
          if( cleared .and. Kbit(data_flag,7) ) then
              call sbit(data_flag,7,0 )
              write_rec = .true.
          end if
 
          do j = 1, num_interedits
 
              if( abs( epoch-epoch_edits(j)).lt. 60.d0/86400. ) then
*                                            ! Times match, check baseline
 
                  if( (site(1).eq. base_edits(1,j)  .or.
     .                 base_edits(1,j).eq.999999 )     .or.
     .                (site(2).eq. base_edits(1,j))      ) then
*                                            ! Site 1 matches
                      if( (site(1).eq.base_edits(2,j) .or.
     .                     base_edits(2,j).eq.999999)    .or.
     .                    (site(2).eq. base_edits(2,j))  ) then
*                                            ! Site 2 matches
*                                                          ! See if we need
*                                                          ! to change bit
                          if( kbit(data_flag,7) .and. 
     .                        type_edits(j).eq.0 ) then
                              call sbit(data_flag,7,0)
                              write_rec = .true.
                          else if( .not.kbit(data_flag,7) .and.
     .                             type_edits(j).ne.0 ) then
                              call sbit(data_flag,7,1)
                              write_rec = .true.
                          end if
                      end if
                  end if
              end if
          end do
 
****      See if we need to write record
          if( write_rec ) then
              call rw_KalObs_block('W','DATA',site,ierr, i)
              if( kbit(data_flag,7) ) then
                  num_edits_out = num_edits_out + 1
              else
                  num_edits_in = num_edits_in + 1
              end if
          end if
*                 ! Looping over KalObs file
      end do
 
***** Now report to users
      write(*,300) KO_Name(1:trimlen(KO_name)), num_edits_in,
     .                                          num_edits_out
 300  format(' For KalObs file ',a,/,i4,' data restored; '
     .       i4,' data removed')

      if( num_edits_in.ne.0 .or. num_edits_out.ne.0 ) then
          call increment_KalVer
          call rw_KalObs_header('W',ierr)
      end if
 
      call close_KalObs(ierr)
 
      end
 
 
