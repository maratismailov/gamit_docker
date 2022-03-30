
CTITLE REMOVE_FB
 
      subroutine remove_fb(ctol_cse, data_flag_cse )

      implicit none
 
*     This routine will scan the one-way data and remove the first
*     bias flag  on each one-way sequence of data.
 
* INCLUDES
 
      include '../includes/kalman_param.h'
      include 'ctogobs_com.h'
 
* PASSED VARIABLES
 
*   data_flag_cse(num_chan, num_cfiles, num_ep) - Data flag for each
*                   - measurement (same definition as in Gobs File)
*   ctol_cse(num_chan, num_cfiles, num_ep)  - Conversion from
*                   - channel number to satellite number
 
 
      integer*4 data_flag_cse(num_chan, num_cfiles, num_ep),
     .    ctol_cse(num_chan, num_cfiles, num_ep)
 
* LOCAL VARIABLES
 
*   i,j,k   - Counters of stations, satellites and epochs
*   ch      - Channel number at current epoch corresponding to
*           - list number lv.
*   ltoc    - function to convert sv list number to channel
*           - number.  Returns -1 if not observed.
 
      integer*4 i,j,k, ch, ltoc
 
*   kbit    - Checks if bit is set.
*   data_OK - Logical function to test if data is OK or not
*   first_bias  - Indicates that first bias flag has been
*           - found
*   good_bf - Function which returns true on a good data point
*             with a bias flag.

      logical first_bias, good_bf

*     Say whatwe are doing
      write(*,100)
 100  format(/,' REMOVING BIAS FLAGS FROM START OF ONE-WAYS')
 
*     Now loop over the one-way data
      do i = 1, num_cfiles
          do j = 1, num_sat
              first_bias = .false.
 
****          Scan up to the first bias flag (always marked) and
*             step at the data sampling interval for this site
              k = 0
              do while( k.lt.num_ep .and. .not.first_bias )
                  k = k + 1
                  ch = ltoc(ctol_cse(1,i,k), j, actual_max_chan)
*                                         ! SV observed at this time
                  if( ch.gt.0 ) then
 
                      if( good_bf(data_flag_cse(ch,i,k), 
     .                                           0,phs_mask) ) then 
                          first_bias = .true.

*                         Now remove the bias flag
                          call sbit(data_flag_cse(ch,i,k),31,0)
                          call sbit(data_flag_cse(ch,i,k),32,0)
                          write(*,200) cf_codes(i), prn_list(j), k
 200                      format(' Removing first bias flag at ',a,
     .                           ' PRN ',I2.2,' Epoch ',i5)
                      end if
                  end if
              end do
 
          end do
*                         ! Looping oversites
      end do
 
***** Thats all
      return
      end
