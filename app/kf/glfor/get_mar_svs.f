
CTITLE get_mar_svs

      subroutine get_mar_svs( epoch )

      implicit none 

*     This routine will open and read the svs_mar_file and get
*     the markov statistics on this experiment.  These are added to
*     the standard markov elements being used for the orbits (i.e.,
*     these are in addition to ones being used for rest of experiments).
 
      include '../includes/kalman_param.h'
      include '../includes/globk_common.h'
      include '../includes/glb_hdr_def.h'

* PASSED VARIABLES

*   epoch       - Time at which the elements are required

      real*8 epoch
 
*   ierr        - Fmpread error
*   trimlen     - Length of line
*   date(5)     - Date read from apriori file
*   indx,jndx   - Pointer for position in string
*   is          - sv number in list of svs
*   nt          - Parameter number for satellite parameter.
 
      integer*4  ierr, date(5), indx, jndx, is, trimlen, nt, i, j, k,
     .           jerr

*   found       - Logical which stays set until we have found some
*                 entries and then it is used for exit

      logical found, done

*   sectag      - second tag of date
*   tepoch      - Test julian date to see if with in 0.5 days of the
*                 center of our experiment.
*   duration    - Time from epoch in file to the end of this markov
*                 statistic days.
*   vals(max_svs_elem)   - Temporary storage of markov values.

      real*8 sectag, tepoch, duration, vals(max_svs_elem)

*   cprn       - PRN name read from file.
*   cval       - dummy for multiread

      character*8 cprn, cval

*   line        - Line read from file

      character*256 line
      character*64 word

***** Try to open the apriori satellite position file

*     Clear the number of markov parameters for this epoch
      num_mar_svs = 0
      if( ichar(svs_mar_file(1:1)).eq.0 ) svs_mar_file = ' '

      if( trimlen(svs_mar_file).gt.0 .and. gnum_svs.gt.0 ) then

          open(105, file=svs_mar_file, status='old', iostat=ierr )
          call report_error('IOSTAT',ierr,'open', svs_mar_file, 0,
     .                      'get_mar_svs')

*         if there was no open error, then try to read the values from
*         the file.

          found = .false.

          do while ( ierr.eq.0 )

              read(105,'(a)', iostat=ierr ) line
              if( line(1:1).eq.' ' .and. ierr.eq.0 .and. 
     .            trimlen(line).gt.0 )   then

*                 See if the date on this line matchs experimentr
                  indx = 1
                  call multiread(line, indx, 'I4', ierr, date, cval, 4)
                  call read_line(line, indx, 'R8', ierr, duration, cval) 
                  date(5) = 0
                  sectag = 0.d0
                  call ymdhms_to_jd( date, sectag, tepoch)
                  if( epoch-tepoch.lt. duration .and. 
     .                epoch-tepoch.ge.0               ) then

                      found = .true.

*                     we have found a time within the duration of the
*                     additional markov statistics
*                     Find the SV Name
                      call getword(line, cprn, indx )
                      jndx = 1
                      call get_cmd(cprn, gsvs_names, gnum_svs, is, 
     .                                 jndx)
                      if( is.gt.0 ) then

* MOD TAH 981215:        Added more flexiable options for reading file
*                        Work along a line see what we have 
                         do j = 1, max_svs_elem
                            vals(j) = 0.d0
                         end do
                         j = 0
                         done = .false.
                         do while ( j.lt. max_svs_elem .and. .not.done )    
                             call GetWord( line, word, indx )
                             call casefold(word)

*                            See if we are end of line
                             if ( trimlen(word).eq.0 ) then
                                done = .true.
                             else
*                               See what we have:
                                call check_num( word, jerr )
                                jndx = 1 

                                if( jerr.eq.0 ) then   ! Value is numeric
                                   j = j + 1
                                   call read_line(word,jndx,'R8',
     .                                            jerr,vals(j), cval)
                                else if ( word(1:1).eq.'!' .or. 
     .                                    word(1:1).eq.'#' ) then
                                   done = .true.
                                else if ( index(word,'R').gt.0 ) then
*                                  Remaining radiation parameters given,
*                                  get the value
                                   call sub_char(word,'R',' ')
                                   j = j + 1 
                                   call read_line(word,jndx,'R8',jerr,
     .                                            vals(j), cval)
                                   do k = j, max_svs_elem-3
                                       vals(k) = vals(j)
                                   end do
                                   j = max_svs_elem-3
                                else if ( index(word,'A').gt.0 ) then
*                                  Remaining radiation parameters given, 
*                                  get the value
                                   call sub_char(word,'A',' ')
                                   j = j + 1 
                                   call read_line(word,jndx,'R8',jerr,
     .                                            vals(j), cval)
                                   do k = j,max_svs_elem 
                                       vals(k) = vals(j)
                                   end do
                                   done = .true.
                                else 
                                   write(*,4610) word(1:max(1,
     .                                     trimlen(word))), 'MAR_SVSF'
 4610                              format('**WARNING** Error decoding ',
     .                                   a, ' in ',a,' command')
                                end if
                             endif
                         end do

*                        Now save the values
                         do i = 1, max_svs_elem
                             nt = parn_svs(i,is)
                             if( nt.gt.0 ) then
                                 num_mar_svs = num_mar_svs + 1
                                 cov_mar_svs(num_mar_svs)  = vals(i)
                                 parn_mar_svs(num_mar_svs) = nt
                             end if 
                         end do
                         write(*,220) gsvs_names(is)
 220                     format('Adding Process noise ',a8)
                      else
                           write(*,'('' PRN '',a,'' Found at correct '',
     .                               '' time but not user PRN list'')')
     .                               cprn
                      end if
                  else
*                     If we have already found svs and now time is
*                     wrong set ierr so that we get out
                      if( found ) ierr = -1
                  end if
              end if
         end do

         close(105)

      end if

***** Thats all
      return
      end 

     
 
