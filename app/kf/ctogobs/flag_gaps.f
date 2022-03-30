
CTITLE FLAG_GAPS
 
      subroutine flag_gaps(pass, ctol_cse, data_flag_cse, bf_type_cse )

      implicit none
 
*     This routine will scan the one-way data and put bias flags
*     on all the gaps in the data.  Once this has been the done the
*     program should be able to more reliable fix data when there is
*     gap on all satellites at the same time.  (Once one double has
*     been fixed all the rest should be be referred back to the
*     one with bias flag removed.)
 
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   bf_type_cse(num_chan, num_cfiles, num_ep) - Bias flag type, records
*                   - why a bias flag was set.  (Never reset even when
*                     the bias flag is removed)

*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number

*   pass            - pass number through possible iteration.  Results
*                     only reported on first pass. 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    bf_type_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep), pass
 
* LOCAL VARIABLES
 
*   i,j,k   - Counters of stations, satellites and epochs
*   ch      - Channel number at current epoch corresponding to
*           - list number lv.
*   ltoc    - function to convert sv list number to channel
*           - number.  Returns -1 if not observed.
*   first_ep    - Epoch number of first data point after the first
*           - bias flag
*   step_ep - Number of epochs between each data point for this
*           - site.
*   num_added - Number of bias flags added due to gaps
*   num_in_gap - Keeps track of number of data in gaps and this
*               value is tested against the number of epochs allowed
*               in the gap
 
      integer*4 i,j,k, ch, ltoc, first_ep, step_ep, num_added,
     .          num_in_gap
 
*   kbit    - Checks if bit is set.
*   data_OK - Logical function to test if data is OK or not
*   first_good  - Indicates that first good data point has been
*           - found
*   in_gap  - Indicates that we are currently in a gap and the
*             next good data point should be flagged with a bias
*             flag.
*   good_bf - Logical function to indicate that we have a good
*             data point with a bias flag.
 
      logical kbit, data_OK, first_good, in_gap
      
*     Say whatwe are doing
      if( pass.le.1 ) write(*,100)
 100  format(/,' ADDING BIAS FLAGS TO ALL GAPS')
 
*     Now loop over the one-way data
      do i = 1, num_cfiles
          do j = 1, num_sat
              first_good = .false.
 
****          Scan up to the first good data point and
*             step at the data sampling interval for this site
              k = 0
              do while( k.lt.num_ep .and. .not.first_good )
                  k = k + 1
                  ch = ltoc(ctol_cse(1,i,k), j, actual_max_chan)
*                                         ! SV observed at this time
                  if( ch.gt.0 ) then
 
C                     if( good_bf(data_flag_cse(ch,i,k),
C    .                                    0,phs_mask) ) then
                      if( data_ok(data_flag_cse(ch,i,k),
     .                                    0,phs_mask) ) then
                          first_good = .true.
                      end if
                  end if
              end do
 
*             Based on the sampling interval at this site compute
*             number of epochs bewteen data points.
              if( orig_sampling(i).gt.sampling ) then
                  step_ep = orig_sampling(i)/sampling
              else
                  step_ep = 1
              end if
 
*             Start at the first point after the first bias flag
*             after saving the first point

              first_ep = k
              in_gap = .false.
              num_in_gap = 0
              num_added = 0
              do k = first_ep, num_ep, step_ep

                  ch = ltoc(ctol_cse(1,i,k), j, actual_max_chan)
*                                         ! SV observed at this time
                  if( ch.gt.0 ) then

*                     See if data is any good 
                      if( data_OK(data_flag_cse(ch,i,k),0,
     .                            phs_mask)  ) then

*                         See if we don't have a bias flag and we
*                         in a gap
                          if( .not.(kbit(data_flag_cse(ch,i,k),31).or.
     .                          kbit(data_flag_cse(ch,i,k),32)) .and.
     .                        in_gap  ) then

*                             Set a bias flag
                              call sbit(data_flag_cse(ch,i,k),31,1)
                              call sbit(bf_type_cse(ch,i,k), 4,1)
                              num_added = num_added + 1
                          end if
 
*                         Set in gap flase since we have found good
*                         data
                          in_gap = .false.
                          num_in_gap = 0
                      else
*                         Mark in gap due to bad data if more than limit
*                         on isze of gaps
                          num_in_gap = num_in_gap + step_ep
                          if( num_in_gap.ge.gap_size(i) ) 
     .                                            in_gap = .true.
                      end if
                  else
*                     Mark in gap due to missing data if more than limit
*                     on isze of gaps
                      num_in_gap = num_in_gap + step_ep
                      if( num_in_gap.ge.gap_size(i) ) in_gap = .true.
                  end if
              end do
              if( pass.le.1 .and. num_added.gt.4 )
     .                     write(*,120) cf_codes(i), prn_list(j), 
     .                     first_ep, step_ep, num_added, gap_size(i)
 120          format(' Gaps flagged on ',a4,' PRN ',i2.2,' Start ',i4,
     .                    ' Step ',i3,' epochs, ',i4,' Bflags added',
     .                    ' Separated by ',i3,' epochs')
 
          end do
*                         ! Looping oversites
      end do
 
***** Thats all
      return
      end
